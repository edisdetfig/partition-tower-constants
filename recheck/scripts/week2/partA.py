# A: exact p(N) mod l^D for N = l^b n + delta_b, saved to a JSON file. Multiprocessing over (l,b,n).
import sys, json, time
from multiprocessing import Pool
from flint import fmpz
D=60
def job(t):
    l,b,n=t; m=l**b; d=pow(24,-1,m); N=m*n+d
    t0=time.time(); v=int(fmpz(N).partitions_p()%(l**D))
    return (l,b,n,N,v,time.time()-t0)
jobs=[]
for b in range(1,16,2):
    for n in range(8): jobs.append((7,b,n))
for b in range(1,12,2):
    for n in range(8): jobs.append((11,b,n))
for n in range(3): jobs.append((11,13,n))
jobs.sort(key=lambda t:-(t[0]**t[1]*(t[2]+1)))     # largest first
if __name__=='__main__':
    out=[]
    with Pool(int(sys.argv[1])) as P:
        for r in P.imap_unordered(job,jobs):
            out.append(r); print(r[:4],'%.0fs'%r[5],flush=True)
            json.dump(out,open('pvals.json','w'))
    print('done',len(out))
