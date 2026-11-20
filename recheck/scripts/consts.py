# Independent recheck: constants of 5..31 from exact partition numbers (FLINT), pairs of levels two apart, n = 0.
import sys, time, json
from flint import fmpz
P=lambda n:int(fmpz(n).partitions_p())
CAP=float(sys.argv[1]) if len(sys.argv)>1 else 1.5e14
REF={5:[4,2,3,1,0,2,0,3,0,1,1,0,1,4,0,1,3,2,3,0,3,1,0,4,3,0,4,1,2,4,2,3,0],
7:[5,6,1,1,2,1,3,1,0,6,5,3,1,1,5,3,3,6,4,2],
11:[9,3,2,0,1,7,9,6,8,4,5,0],
13:[6,3,1,2,5,8,0,9,5,9,9,6,8,2,8,5,3,0,8,10,9,11,0,1,0,10,8,8,7,2],
17:[6,12,12,7,7,5,16,3,10,14,0,0,15,0,5,15,16,14,1,13,1,12,4,1,0,15,1,1,8,13],
19:[10,2,10,8,10,1,16,5,7,10,6,10,8,2,18,0,8,12,12,11,14,14,0,1,3,14,10,17,10,18],
23:[5,20,8,4,19,3,6,6,20,20,3,8,14,11,21,0,10,7,6,22,15,17,1,4,18,7,14,6,16,22],
29:[3,14,17,22,28,6,14,19,1,14,0,17,2,0,26,24,3,0,21,17,16,21,12,4,11,9,22,20,8,15],
31:[30,6,2,28,1,8,3,15,15,13,16,5,3,6,27,18,2,9,28,2,17,23,10,26,8,3,1,4,0,19]}
def guar(l,b): return {5:b,7:b//2+1,11:b}.get(l,0)
def val(x,l):
    k=0
    while x%l==0: x//=l;k+=1
    return k
jobs=[]
for l in REF:
    b=1
    while True:
        m=l**b; d=pow(24,-1,m)
        if d>CAP: break
        jobs.append((d,l,b)); b+=1
        if b>40: break
jobs.sort()
print('jobs',len(jobs),'largest',jobs[-5:],flush=True)
H={}; t0=time.time()
for d,l,b in jobs:
    t=time.time(); v=P(d); H[(l,b)]=v
    if time.time()-t>5: print('  p(%d) l=%d b=%d  %.0f s (total %.0f)'%(d,l,b,time.time()-t,time.time()-t0),flush=True)
res={}
for l in REF:
    K=len(REF[l]); M=l**(K+2); out=[]
    bs=sorted(b for (ll,b) in H if ll==l)
    for b in bs:
        if (l,b+2) not in H: continue
        x,y=H[(l,b)],H[(l,b+2)]
        vx,vy=val(x,l),val(y,l)
        gx,gy=guar(l,b),guar(l,b+2)
        assert vx>=gx and vy>=gy,(l,b,vx,vy)
        ux,uy=x//l**vx,y//l**vy
        r=uy*pow(ux,-1,M)%M
        shift=(vy-gy)-(vx-gx)       # extra powers beyond the guaranteed ones
        dig=[]; z=r
        for i in range(K): dig.append(z%l); z//=l
        agree=0
        if shift==0:
            while agree<K and dig[agree]==REF[l][agree]: agree+=1
        out.append((b,b+2,vx-gx,vy-gy,agree,dig[:min(K,agree+2)]))
    res[l]=out
    print('l=%d (reference digits known: %d)'%(l,K))
    for o in out: print('   levels %d,%d  extra valuations %d,%d  digits agreeing with the project: %d   mine: %s'%o)
    sys.stdout.flush()
json.dump({str(k):v for k,v in res.items()},open('consts.json','w'))
print('total seconds',round(time.time()-t0))
