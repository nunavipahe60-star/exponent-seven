from concurrent.futures import ProcessPoolExecutor,as_completed
from pair_check import run
from time import perf_counter
if __name__=='__main__':
    start=perf_counter()
    with ProcessPoolExecutor(max_workers=6) as ex:
        fut={ex.submit(run,r,None,i):(r,i) for r in (11,12,13) for i in range(1,5)}
        totals={r:[0,0,0,0] for r in (11,12,13)}
        for f in as_completed(fut):
            (r,i)=fut[f];data=f.result()
            for j in range(4):totals[r][j]+=data[j]
            print('Rank, p mod 5, stats:',r,i,data,flush=True)
    print('TOTAL',totals,'Wall seconds',perf_counter()-start,flush=True)
