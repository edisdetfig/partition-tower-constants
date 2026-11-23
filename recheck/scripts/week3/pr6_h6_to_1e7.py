# a_p(h6) for every prime p < 10^7, exactly, to test p | a_p.  h6 = the weight-6 newform on Gamma0(6).
# Series are computed modulo P = 2^61-1; |a_p| <= 2 p^(5/2) < 6.4e17 < P/2, so the centred residue is a_p itself.
import time, sys, numpy as np, cypari2
from flint import nmod_poly
t0=time.time(); N=10**7+1; P=2**61-1
def euler(d):
    c=[0]*N; k=0
    while True:
        e1=d*k*(3*k-1)//2; e2=d*k*(3*k+1)//2
        if e1>=N: break
        s=1 if k%2==0 else P-1
        c[e1]=s
        if e2<N and k>0: c[e2]=s
        k+=1
    return nmod_poly(c,P)
E={d:euler(d) for d in (1,2,3,6)}; print('euler products',round(time.time()-t0),flush=True)
def etaq(r):                       # prod eta(d z)^r_d as a q-series (leading power of q returned separately)
    sh=sum(d*e for d,e in r.items()); assert sh%24==0
    out=None
    for d,e in r.items():
        if e==0: continue
        f=E[d].pow_trunc(abs(e),N)
        if e<0: f=f.inverse_series_trunc(N)
        out=f if out is None else out.mul_low(f,N)
    return out.left_shift(sh//24).truncate(N)
# cusp forms of weight 6 on Gamma0(6) among eta quotients: all four cusp orders positive
def orders(a,b,c,d): return (6*a+3*b+2*c+d, 3*a+6*b+c+2*d, 2*a+b+6*c+3*d, a+2*b+3*c+6*d)
cands=[]
R=range(-12,25)
for a in R:
    for b in R:
        for c in R:
            d=12-a-b-c
            o=orders(a,b,c,d)
            if all(x>0 and x%24==0 for x in o) and (b+d)%2==0 and (c+d)%2==0: cands.append((a,b,c,d))
print("candidate eta quotients",cands,flush=True); cands=[(6,0,6,0),(0,6,0,6),(1,1,5,5)]; print("using",cands,flush=True)
pari=cypari2.Pari(); pari.allocatemem(10**9)
h6=[int(x) for x in pari.mfcoefs(pari.mfeigenbasis(pari.mfinit([6,6],0))[0],3000)]
forms=[etaq({1:a,2:b,3:c,6:d}) for (a,b,c,d) in cands]; print('forms',len(forms),round(time.time()-t0),flush=True)
# solve h6 = sum x_i forms_i on the first coefficients, modulo P
from flint import nmod_mat
m=len(forms); rows=[[int(f[n]) for f in forms] for n in range(1,m+1)]
A=nmod_mat(m,m,[x for r in rows for x in r],P); b=nmod_mat(m,1,[h6[n]%P for n in range(1,m+1)],P)
x=A.solve(b); xs=[int(x[i,0]) for i in range(m)]
H=None
for xi,f in zip(xs,forms):
    t=f*xi; H=t if H is None else H+t
co=H.coeffs(); co+= [0]*(N-len(co)); print('combination built',round(time.time()-t0),flush=True)
cen=lambda v:(int(v)-P if int(v)>P//2 else int(v))
bad=[n for n in range(1,3001) if cen(co[n])!=h6[n]]
print('agreement with PARI on a_1..a_3000:',len(bad)==0, '; coefficients x (centred):',[cen(v) for v in xs],flush=True)
sieve=np.ones(N,dtype=bool); sieve[:2]=False
for i in range(2,int(N**0.5)+1):
    if sieve[i]: sieve[i*i::i]=False
primes=np.nonzero(sieve)[0]; print('primes below 10^7:',len(primes),flush=True)
hits=[]; viol=0
for p in primes:
    p=int(p); a=cen(co[p])
    if abs(a)>2*p**2.5+1: viol+=1
    if p>=5 and a%p==0: hits.append(p)
print('Deligne bound violations:',viol)
print('checks: a5,a7,a11,a13 =',[cen(co[p]) for p in (5,7,11,13)])
print('primes 5 <= p < 10^7 with p | a_p(h6):',hits)
print('total seconds',round(time.time()-t0))
