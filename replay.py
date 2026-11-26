"""Replay the analytic certificates and audit the three saved GPU logs.

Use --gpu to regenerate the complete direct trace sieve, or --cpu to run
the mathematically distinct two-prime computation in all three rank bands.
"""
import argparse
import ast
from datetime import datetime
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
CPU_COUNTS = {11: [542464, 493357, 13115, 0], 12: [44026, 39339, 5802, 0],
              13: [860, 743, 655, 0]}


def call(name, *args, **options):
    return subprocess.run([sys.executable, "-B", str(ROOT / name), *map(str,args)],
                          check=True, **options)


def verify_cpu_totals(output):
    summary = re.search(r"^TOTAL (\{[^\n]*\}) Wall seconds", output, re.M)
    if summary is None or ast.literal_eval(summary.group(1)) != CPU_COUNTS:
        raise ValueError("Missing or inconsistent CPU totals")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gpu", action="store_true")
    ap.add_argument("--cpu", action="store_true")
    opts = ap.parse_args()
    if not __debug__:
        ap.error("Run with Python assertions enabled")
    certificate = call("kernel/verify.py", stdout=subprocess.PIPE, text=True).stdout
    print(certificate, end="", flush=True)
    rows = certificate.splitlines()
    published = (ROOT/"kernel/certificate.txt").read_text(encoding="utf-8").splitlines()
    if len(rows) != 57 or rows != published:
        raise ValueError("The 57 recomputed Fourier certificates differ from the release")
    call("cpu/rank_certificate.py")
    call("cpu/check_table.py")
    copies = [ROOT/"logs/trace-primary.log", ROOT/"logs/trace-crosscheck.log"]
    if opts.gpu:
        path = Path.cwd() / ("trace-replay-" + datetime.now().strftime("%Y%m%d-%H%M%S-%f") + ".log")
        with path.open("x", encoding="utf-8") as out:
            call("gpu/trace_sieve.py", "--full", "--selftest", stdout=out)
        copies.append(path)
        print("New trace log", path, flush=True)
    call("gpu/audit_log.py", *copies)
    call("gpu/audit_log.py", "--upper-prime", 3169, ROOT/"logs/trace-short.log")
    verify_cpu_totals((ROOT/"logs/cpu-divisors.log").read_text(encoding="utf-8-sig"))
    print("Archived CPU totals agree in all three rank bands.", flush=True)
    if opts.cpu:
        call("cpu/positive_control.py")
        result = call("cpu/replay_pairs.py", capture_output=True, text=True)
        print(result.stdout, end="")
        verify_cpu_totals(result.stdout)
    print("Checks passed.", flush=True)


if __name__ == "__main__":
    main()
