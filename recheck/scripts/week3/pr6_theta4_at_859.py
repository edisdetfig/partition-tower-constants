# Independent check: T_l mod l on the sign class of M_{l-3}(Gamma0(6)), from q-expansions of eta quotients.
import sys, time, numpy as np
from fractions import Fraction
from flint import nmod_mat, nmod_poly
L=int(sys.argv[1]); p=L; k=L-3; n=k+1
PREC=p*(n+2)+16
FN=1<<int(np.ceil(np.log2(2*PREC)))
def mul(a,b):
    a2=a.astype(np.float64); a2[a2>p//2]-=p
    b2=b.astype(np.float64); b2[b2>p//2]-=p
    c=np.fft.irfft(np.fft.rfft(a2,FN)*np.fft.rfft(b2,FN),FN)[:PREC]
    return np.mod(np.rint(c).astype(np.int64),p)
def inv(a):
    x=np.zeros(PREC,dtype=np.int64); x[0]=pow(int(a[0]),-1,p); m=1
    while m<PREC:
        m*=2
        t=mul(a,x); t=(-t)%p; t[0]=(t[0]+2)%p
        x=mul(x,t); x[m:]=0
    return x
def power(a,e):
    if e<0: a=inv(a); e=-e
    r=np.zeros(PREC,dtype=np.int64); r[0]=1
    while e:
        if e&1: r=mul(r,a)
        e>>=1
        if e: a=mul(a,a)
    return r
E=np.zeros(PREC,dtype=np.int64); j=0
while True:
    done=True
    for m in ((j*(3*j-1))//2,(j*(3*j+1))//2):
        if m<PREC:
            E[m]=(E[m]+(1 if j%2==0 else -1))%p if not (j==0 and m==0 and E[0]) else E[m]; done=False
    if done: break
    j+=1
def Ed(d):
    r=np.zeros(PREC,dtype=np.int64); r[::d]=E[:len(r[::d])]; return r
# solve eta exponents: orders at cusps d=1,2,3,6 (in units 1/24)
import itertools
def orders(a,b,c,d): return (6*a+3*b+2*c+d, 3*a+6*b+c+2*d, 2*a+b+6*c+3*d, a+2*b+3*c+6*d)
def find(target,wt):
    R=range(-14,15)
    for a in R:
        for b in R:
            for c in R:
                d=2*wt-a-b-c
                if orders(a,b,c,d)==target: return (a,b,c,d)
F0e=find((48,0,0,0),2); te=find((-24,0,0,24),0)
print('eta exponents F0',F0e,'t',te); sys.stdout.flush()
def etaq(e,scale=1):
    r=np.zeros(PREC,dtype=np.int64); r[0]=1
    for d,x in zip((1,2,3,6),e):
        if x: r=mul(r,power(Ed(d),x*scale))
    return r
t0=time.time()
t=etaq(te); t=np.concatenate(([0],t[:-1]))   # q * product
G=etaq(F0e,k//2)
print('series ready',time.time()-t0); sys.stdout.flush()
idxL=np.arange(n)*p
low=np.zeros((n,5*n+5),dtype=np.int64); atL=np.zeros((n,n),dtype=np.int64)
b=G.copy()
for i in range(n):
    low[i]=b[:5*n+5]; atL[i]=b[idxL]
    if i<n-1: b=mul(b,t)
print('basis done',time.time()-t0); sys.stdout.flush()
assert all(low[i,i]==1 and not low[i,:i].any() for i in range(n))
B=nmod_mat(n,n,[int(x) for x in low[:,:n].flatten()],p); Binv=B.inv()
def op(V): return nmod_mat(n,n,[int(x) for x in V.flatten()],p)*Binv   # row i = coords of image of b_i
U2=op(low[:,0:2*n:2]); U3=op(low[:,0:3*n:3]); TL=op(atL)
V5=low[:,0:5*n:5].copy(); c5=pow(5,k-1,p)
for m in range(0,n,5): V5[:,m]=(V5[:,m]+c5*low[:,m//5])%p
T5=op(V5)
def leg(a): return 1 if pow(a%p,(p-1)//2,p)==1 else -1
def ident(r): return nmod_mat(r,r,[1 if i==j else 0 for i in range(r) for j in range(r)],p)
def genleft(A,c):           # rows spanning {x: x(A-c)^N=0}, in rref, with pivots
    r=A.nrows(); M=A-ident(r)*c
    for _ in range(11): M=M*M
    X,nul=M.transpose().nullspace()
    R=nmod_mat(nul,r,[int(X[i,j]) for j in range(nul) for i in range(r)],p)
    R,rk=R.rref(); piv=[next(j for j in range(r) if int(R[i,j])) for i in range(rk)]
    return R,piv
def restrict(A,R,piv):
    RA=R*A; r=R.nrows()
    return nmod_mat(r,r,[int(RA[i,j]) for i in range(r) for j in piv],p)
# the class h6 predicts: U2 eigenvalue +(2|l)2^((k-2)/2), U3 eigenvalue -(3|l)3^((k-2)/2)
c2=leg(2)*pow(2,(k-2)//2,p)%p; c3=(-leg(3)*pow(3,(k-2)//2,p))%p
R1,p1=genleft(U2,c2); U3r=restrict(U3,R1,p1); R2,p2=genleft(U3r,c3)
TLc=restrict(restrict(TL,R1,p1),R2,p2); d=TLc.nrows()
X,nul=TLc.transpose().nullspace()          # x with x*TLc = 0
print('class dim',d,'; kernel of T_l has dimension',nul)
x=nmod_mat(1,d,[int(X[i,0]) for i in range(d)],p)
v=x*R2*R1                                   # coordinates in the basis b_i of M_k(Gamma0(6))
vv=np.array([int(v[0,i]) for i in range(n)],dtype=object)
m=low.shape[1]
import cypari2; pari=cypari2.Pari()
M=4000 if m>4000 else m-1
qq=(np.array(vv,dtype=object)@low.astype(object))%p
h6=[int(c) for c in pari.mfcoefs(pari.mfeigenbasis(pari.mfinit([6,6],0))[0],M)]
b1=int(qq[1]); print('a_0, a_1 of the T_l-kernel form:',int(qq[0]),b1)
inv=pow(b1,-1,p)
ok=sum(1 for nn in range(1,M+1) if (pow(nn,4,p)*int(qq[nn])*inv-h6[nn])%p==0)
print('n^4 * a_n(f) = a_n(h6) mod %d for n = 1..%d: %d of %d'%(p,M,ok,M))
print('seconds',round(time.time()-t0))
