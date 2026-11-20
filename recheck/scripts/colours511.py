# Independent recheck: constants of r-coloured partition towers from exact values (power series mod l^K).
import sys, time
from flint import nmod_poly
def euler(N,M):
    c=[0]*N; k=0
    while True:
        a=k*(3*k-1)//2; b=k*(3*k+1)//2
        if a>=N: break
        s=1 if k%2==0 else M-1
        c[a]=(c[a]+s)%M
        if b<N and k>0: c[b]=(c[b]+s)%M
        k+=1
    return nmod_poly(c,M)
def pr(r,N,M):
    return euler(N,M).inverse_series_trunc(N).pow_trunc(r,N)
def val(x,l,cap):
    if x==0: return cap
    k=0
    while x%l==0: x//=l;k+=1
    return k
# (r, l, top level B, reference digits (project's files), slope s, guaranteed valuation function)
CASES=[(5,11,6,[2,8,5,2,10,5,2,0,1,1,3],1)]
t0=time.time()
for r,l,B,ref,s in CASES:
    K=len(ref)+B+2
    while l**K>=2**63: K-=1
    M=l**K
    N=l**B+1
    P=pr(r,N,M); c=lambda i:int(P[i])
    print('r=%d l=%d: series to %d terms mod l^%d  (%.0f s)'%(r,l,N,K,time.time()-t0),flush=True)
    for n in (0,1):
        row=[]
        for b in range(1,B-1):
            if l**(b+2)*n+(r*pow(24,-1,l**(b+2)))%l**(b+2)>=N: continue
            x=c(l**b*n+(r*pow(24,-1,l**b))%l**b); y=c(l**(b+2)*n+(r*pow(24,-1,l**(b+2)))%l**(b+2))
            vx,vy=val(x,l,K),val(y,l,K)
            if vx>=K-2 or vy>=K-2: row.append((b,'zero')); continue
            if vy-vx!=s: row.append((b,'valuations %d,%d'%(vx,vy))); continue
            Mm=l**(K-vy)
            q=(y//l**vy)*pow(x//l**vx,-1,Mm)%Mm
            dig=[]; z=q
            for i in range(len(ref)): dig.append(z%l); z//=l
            ag=0
            while ag<len(ref) and ag<K-vy and dig[ag]==ref[ag]: ag+=1
            row.append((b,vx,vy,ag,dig[:ag+1]))
        print('   n=%d: (level b, v(H_b), v(H_b+2), digits of H_(b+2)/(l^s H_b) agreeing with the project, mine)'%n)
        for x in row: print('      ',x)
    sys.stdout.flush()
print('total',round(time.time()-t0))
