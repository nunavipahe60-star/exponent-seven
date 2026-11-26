"""Benchmark or run an exhaustive split-prime trace sieve on an NVIDIA GPU.

Requires Python 3, NumPy and an NVIDIA driver (Turing or later). The default
kernel uses PTX. Optional --nvrtc compiles the equivalent CUDA C++ kernel.
All final bounds and all local conditions use integer arithmetic.
"""
import argparse
import ctypes as C
import math
import os
import time
from pathlib import Path

import numpy as np


def primes(n):
    return [p for p in range(2, n + 1)
            if all(p % d for d in range(2, math.isqrt(p) + 1))]


WHEEL = primes(23)[1:]
M = math.prod(WHEEL)
TEST = [p for p in primes(251) if p > 23]
LANES = 32
LOW = 31 * 10**19


def local(p, ell):
    a = 4 * pow(p, 7, ell)
    return [(pow((r*r-a) % ell, (ell-1)//2, ell) != 1) for r in range(ell)]


def square_set_test(p, ell):
    nonzero_squares = {j*j % ell for j in range(1, (ell+1)//2)}
    return [(r*r - 4*p**7) % ell not in nonzero_squares for r in range(ell)]


def wheel_crt(p):
    # Direct CRT reconstruction from each pair of admissible residues.
    residues = np.array([0], dtype=np.int64)
    modulus = 1
    for ell in WHEEL:
        allowed = np.flatnonzero(square_set_test(p, ell))
        inv = pow(modulus, -1, ell)
        residues = np.concatenate([residues + modulus * (((a-residues)*inv) % ell)
                                   for a in allowed])
        modulus *= ell
    return np.asarray(residues, dtype=np.uint32)


def wheel(p):
    residues = np.array([0], dtype=np.uint32)
    modulus = 1
    for ell in WHEEL:
        lifts = (residues[None, :] +
                 np.arange(ell, dtype=np.uint32)[:, None] * modulus).ravel()
        residues = lifts[np.asarray(local(p, ell))[lifts % ell]]
        modulus *= ell
    return residues


def masks(p):
    out = []
    for ell in TEST:
        ok = local(p, ell)
        out.extend(sum(1 << j for j in range(LANES) if ok[(r+j*M) % ell])
                   for r in range(ell))
    return np.asarray(out, dtype=np.uint32)


def source():
    body = []
    offset = 0
    for i, ell in enumerate(TEST):
        body.append(f"    v &= __ldg(mask+{offset}+(r%{ell}u+"
                    f"(tile%{ell}u)*{(LANES*M)%ell}u)%{ell}u);")
        if i in (5, 10, 15, 23, 31, 40):
            body.append("    if (!v) return;")
        offset += ell
    return ('''extern "C" __global__ void sieve(const unsigned *res,
        unsigned count, const unsigned *mask, unsigned first,
        unsigned long long high, unsigned long long *hits, unsigned long long *out,
        unsigned capacity) {
    unsigned idx=blockIdx.x*blockDim.x+threadIdx.x;
    if (idx>=count) return;
    unsigned r=res[idx], tile=first+blockIdx.y;
    unsigned long long base=(unsigned long long)tile*32*MODULUS+r;
    if (base>high) return;
    unsigned v=~0u;
''' + '\n'.join(body) + '''
    while(v) {
        unsigned bit=__ffs(v)-1;
        unsigned long long x=base+(unsigned long long)bit*MODULUS;
        if (x<=high && x>0) {
            unsigned long long pos=atomicAdd(hits,1ull);
            if (pos<capacity) out[pos]=x;
        }
        v&=v-1;
    }
}''').replace("MODULUS", str(M)+"ull").encode()


U = C.c_uint
U64 = C.c_uint64
PTR = C.c_void_p


def check(code):
    if code:
        raise RuntimeError("CUDA call failed with code "+str(code))


def portable_ptx():
    # Use the GPU driver's built-in PTX assembler, including its constant
    # modulus optimizations. The same register names are reused per prime.
    checks = []
    offset = 0
    for i, ell in enumerate(TEST):
        checks.append(f"""
        rem.u32 a,r,{ell};
        rem.u32 b,tile,{ell};
        mad.lo.u32 a,b,{(LANES*M)%ell},a;
        rem.u32 a,a,{ell};
        add.u32 a,a,{offset};
        mul.wide.u32 address,a,4;
        add.u64 address,masks,address;
        ld.global.nc.u32 a,[address];
        and.b32 v,v,a;""")
        if i in (5, 10, 15, 23, 31, 40):
            checks.append("setp.eq.u32 cond,v,0; @cond bra done;")
        offset += ell
    return (f""".version 6.4
.target sm_75
.address_size 64
.visible .entry sieve(.param .u64 R, .param .u32 N, .param .u64 T,
 .param .u32 F, .param .u64 H, .param .u64 C, .param .u64 O, .param .u32 K) {{
 .reg .u32 idx,n,r,tile,a,b,v,cap;
 .reg .u64 address,base,hi,masks,outp,cnt,x;
 .reg .pred cond;
 mov.u32 a,%ctaid.x; mov.u32 b,%ntid.x; mov.u32 idx,%tid.x;
 mad.lo.u32 idx,a,b,idx;
 ld.param.u32 n,[N]; setp.ge.u32 cond,idx,n; @cond bra done;
 ld.param.u64 address,[R]; mul.wide.u32 base,idx,4;
 add.u64 address,address,base; ld.global.u32 r,[address];
 mov.u32 tile,%ctaid.y; ld.param.u32 a,[F]; add.u32 tile,tile,a;
 mul.wide.u32 base,tile,{M}; shl.b64 base,base,5; cvt.u64.u32 x,r; add.u64 base,base,x;
 ld.param.u64 hi,[H]; setp.gt.u64 cond,base,hi; @cond bra done;
 ld.param.u64 masks,[T]; mov.u32 v,0xffffffff;
 """ + "\n".join(checks) + f"""
 ld.param.u64 cnt,[C]; ld.param.u64 outp,[O]; ld.param.u32 cap,[K];
again:
 setp.eq.u32 cond,v,0; @cond bra done;
 brev.b32 a,v; clz.b32 a,a; mul.wide.u32 x,a,{M}; add.u64 x,base,x;
 setp.gt.u64 cond,x,hi; @cond bra skip; setp.eq.u64 cond,x,0; @cond bra skip;
 atom.global.add.u64 address,[cnt],1; cvt.u64.u32 cnt,cap;
 setp.ge.u64 cond,address,cnt; ld.param.u64 cnt,[C]; @cond bra skip;
 shl.b64 address,address,3; add.u64 address,outp,address; st.global.u64 [address],x;
skip:
 sub.u32 a,v,1; and.b32 v,v,a; bra again;
done: ret;
}}
 """).encode()


class Device:
    def __init__(self, path):
        if os.name == "nt":
            driver = C.WinDLL("nvcuda.dll")
        else:
            driver = C.CDLL("libcuda.so.1")
        self.d = driver
        check(driver.cuInit(0))
        unit = C.c_int()
        check(driver.cuDeviceGet(C.byref(unit), 0))
        self.context = PTR()
        check(driver.cuCtxCreate_v2(C.byref(self.context), 0, unit))
        ptx = self.compile_cuda(path) if path else C.create_string_buffer(portable_ptx())
        mod = PTR()
        check(driver.cuModuleLoadData(C.byref(mod), ptx))
        self.fn = PTR()
        check(driver.cuModuleGetFunction(C.byref(self.fn), mod, b"sieve"))

    def compile_cuda(self, path):
        if os.name == "nt":
            directory = Path(path).resolve().parent
            self.dll_path = os.add_dll_directory(str(directory))
            self.builtins = [C.WinDLL(str(dll)) for dll in directory.glob("nvrtc-builtins64_*.dll")]
            runtime = C.WinDLL(str(Path(path).resolve()))
        else:
            runtime = C.CDLL(path)
        prog = PTR()
        check(runtime.nvrtcCreateProgram(C.byref(prog), source(), b"sieve.cu", 0, None, None))
        opts = (C.c_char_p * 1)(b"--gpu-architecture=compute_75")
        err = runtime.nvrtcCompileProgram(prog, 1, opts)
        if err:
            n = C.c_size_t()
            runtime.nvrtcGetProgramLogSize(prog, C.byref(n))
            log = C.create_string_buffer(n.value)
            runtime.nvrtcGetProgramLog(prog, log)
            raise RuntimeError(log.value)
        n = C.c_size_t()
        check(runtime.nvrtcGetPTXSize(prog, C.byref(n)))
        ptx = C.create_string_buffer(n.value)
        check(runtime.nvrtcGetPTX(prog, ptx))
        check(runtime.nvrtcDestroyProgram(C.byref(prog)))
        return ptx

    def upload(self, array):
        ptr = U64()
        check(self.d.cuMemAlloc_v2(C.byref(ptr), C.c_size_t(array.nbytes)))
        check(self.d.cuMemcpyHtoD_v2(ptr, PTR(array.ctypes.data), C.c_size_t(array.nbytes)))
        return ptr

    def run(self, residues, table, bound, spans):
        capacity = 65536
        outputs = np.zeros(capacity, dtype=np.uint64)
        nhit = np.zeros(1, dtype=np.uint64)
        res, mask, hits, out = map(self.upload, (residues, table, nhit, outputs))
        start = time.perf_counter()
        for first, length in spans:
            args = (res, U(len(residues)), mask, U(first), U64(bound), hits, out, U(capacity))
            ptrs = (PTR * len(args))(*(C.cast(C.byref(a), PTR) for a in args))
            check(self.d.cuLaunchKernel(self.fn, (len(residues)+127)//128, length, 1,
                                        128, 1, 1, 0, None, ptrs, None))
            check(self.d.cuCtxSynchronize())
        seconds = time.perf_counter()-start
        check(self.d.cuMemcpyDtoH_v2(PTR(nhit.ctypes.data), hits, C.c_size_t(nhit.nbytes)))
        if nhit[0] > capacity:
            raise OverflowError("Hit buffer filled; subdivide this run")
        check(self.d.cuMemcpyDtoH_v2(PTR(outputs.ctypes.data), out, C.c_size_t(outputs.nbytes)))
        for ptr in (out, hits, mask, res):
            check(self.d.cuMemFree_v2(ptr))
        return sorted(map(int, outputs[:int(nhit[0])])), seconds


def main():
    global WHEEL, TEST, M, local
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("p", nargs="*", type=int, help="Individual least split primes to test")
    ap.add_argument("--nvrtc")
    ap.add_argument("--full", action="store_true", help="Exhaust every prime from 701 to 4597")
    ap.add_argument("--short", action="store_true", help="Exhaust every prime from 701 to 3169")
    ap.add_argument("--crosscheck", action="store_true",
                    help="Use a different wheel, square sets, CRT construction and mask order")
    ap.add_argument("--tiles", type=int, default=8,
                    help="Tile count per sample location; zero scans every tile")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--low", type=int, default=LOW)
    args = ap.parse_args()
    if args.full or args.short:
        if args.p or args.low != LOW or (args.full and args.short):
            ap.error("Specify a single fixed proof range with the default lower endpoint")
        args.p = [p for p in primes(3169 if args.short else 4597) if p>=701]
        args.tiles = 0
    if not args.p or args.tiles<0 or args.low<1:
        ap.error("Use --full or supply primes with a positive lower endpoint")
    make_wheel = wheel
    if args.crosscheck:
        WHEEL = [3, 5, 7, 11, 13, 17, 19, 29]
        M = math.prod(WHEEL)
        TEST = [p for p in reversed(primes(251)) if p>2 and p not in WHEEL]
        local = square_set_test
        make_wheel = wheel_crt
    gpu = Device(args.nvrtc)
    if args.selftest:
        # Project several prescribed positive test traces (including x=1, the
        # last lane, and upper-end clipping) into every device mask table.
        res = np.asarray([1, 4, 200, M-1], dtype=np.uint32)
        cap = 35*M+4
        targets = [1, M-1, 31*M+4, 32*M+200, cap, cap+1]
        table = []
        oksets = []
        for ell in TEST:
            ok = {x % ell for x in targets}
            oksets.append(ok)
            table.extend(sum(1 << j for j in range(32) if (r+j*M)%ell in ok)
                         for r in range(ell))
        expected = sorted(x for r in map(int,res) for j in range(64)
                          if (x:=r+j*M)<=cap and all(x%e in s for e,s in zip(TEST,oksets)))
        actual,_ = gpu.run(res, np.asarray(table,dtype=np.uint32), cap, [(0,2)])
        assert expected==actual and set(targets[:-1]).issubset(actual), (actual,expected)
        # A second oracle directly evaluates the genuine quadratic conditions.
        p=701
        top=50000
        r=make_wheel(p)
        for e in WHEEL+TEST:
            assert local(p,e)==[pow((x*x-4*p**7)%e,(e-1)//2,e)!=1 for x in range(e)]
        assert sorted(map(int,r[r<=top]))==[x for x in range(top+1) if
            all(pow((x*x-4*p**7)%e,(e-1)//2,e)!=1 for e in WHEEL)]
        expected=[x for x in range(1,top+1)
                  if all(pow((x*x-4*p**7)%e,(e-1)//2,e)!=1 for e in WHEEL+TEST)]
        actual,_=gpu.run(r,masks(p),top,[(0,1)])
        assert actual==expected
        print("Self-tests passed.",flush=True)
    print("p,upper_x,wheel_residues,total_tiles,measured_tiles,gpu_seconds,setup_seconds,hits",flush=True)
    for p in args.p:
        if p <= max(WHEEL+TEST) or p not in primes(p):
            ap.error("Each p must be a prime greater than the last test prime")
        delta = 4*p**7 - args.low
        if delta <= 0:
            continue
        bound = math.isqrt(delta)
        nt = bound//(32*M)+1
        start = time.perf_counter()
        residues,table = make_wheel(p),masks(p)
        setup = time.perf_counter()-start
        if args.tiles == 0 or nt <= 3*args.tiles:
            spans = [(i,min(8,nt-i)) for i in range(0,nt,8)]
        else:
            n = args.tiles
            spans = [(0,n),((nt-n)//2,n),(nt-n,n)]
        result,sec = gpu.run(residues,table,bound,spans)
        print(f"{p},{bound},{len(residues)},{nt},{sum(n for _,n in spans)},"
              f"{sec:.6f},{setup:.6f},{len(result)}", flush=True)
        for x in result:
            print(f"trace p={p} x={x}",flush=True)


if __name__ == "__main__":
    main()
