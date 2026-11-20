# Independent spot-check of the share test, from the README's definition only (exact partition numbers, no forms).
import sys
from flint import fmpz
from sympy import jacobi_symbol as js, factorint
def val(x,l,N):
    x%=l**N
    if x==0: return N
    v=0
    while x%l==0: x//=l; v+=1
    return v
def snf_vals(M,l,N):
    m=l**N; M=[[x%m for x in r] for r in M]; n=len(M); out=[]
    for s in range(n):
        best=min(((val(M[i][j],l,N),i,j) for i in range(s,n) for j in range(s,n)))
        v,i,j=best; out.append(v)
        if v>=N: out+= [N]*(n-s-1); break
        M[s],M[i]=M[i],M[s]
        for r in M: r[s],r[j]=r[j],r[s]
        u=pow(M[s][s]//l**v,-1,m)
        for i in range(s+1,n):
            f=(M[i][s]//l**v)*u%m
            M[i]=[(a-f*b)%m for a,b in zip(M[i],M[s])]
        for j in range(s+1,n):
            f=(M[s][j]//l**v)*u%m
            for r in M[s:]: r[j]=(r[j]-f*r[s])%m
    return sorted(out)
c=lambda M: int(fmpz((M+1)//24).partitions_p())
for l in [int(a) for a in sys.argv[1:]]:
    R=(l+12)//24; K=l//24; S=R+K+2; N=5; m=l**N
    Qs=[q for q in (5,7,11,13,17,19,23,29,31,37) if q!=l][:S]
    for b in (1,3):
        ts=[t for t in range(1,4000) if (l**b*t+1)%24==0 and t%l and all(e==1 for e in factorint(t).values())][:S]
        Mx=[[ (c(Q*Q*l**b*t)+int(js(l**b*t%Q,Q))*int(js(-3%Q,Q))*pow(Q*Q,-1,m)*c(l**b*t))%m for Q in Qs] for t in ts]
        print('l=%d level %d: grid %dx%d, R=%d, K=%d; valuations of the invariant factors: %s'%(l,b,S,S,R,K,snf_vals(Mx,l,N)),flush=True)
