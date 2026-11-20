# Independent check of the hypothesis of E1: T_Q acts on the class S_l (mod l) as one number.
# Class = generalized eigenspace of U2 (value 1/4) and U3 (value 1/9) on M_{l-3}(Gamma0(6)) mod l, from eta quotients.
import sys, time, numpy as np
from flint import nmod_mat
L=int(sys.argv[1]); Qs=[int(x) for x in sys.argv[2].split(',')]; FULL=(len(sys.argv)>3 and sys.argv[3]=='full')
p=L; k=L-3; n=k+1
PREC=max(p*(n+2)+16, (max(Qs)*n+16) if FULL else (max(Qs)+16))
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
def orders(a,b,c,d): return (6*a+3*b+2*c+d, 3*a+6*b+c+2*d, 2*a+b+6*c+3*d, a+2*b+3*c+6*d)
def find(target,wt):
    R=range(-14,15)
    for a in R:
        for b in R:
            for c in R:
                d=2*wt-a-b-c
                if orders(a,b,c,d)==target: return (a,b,c,d)
F0e=find((48,0,0,0),2); te=find((-24,0,0,24),0)
def etaq(e,scale=1):
    r=np.zeros(PREC,dtype=np.int64); r[0]=1
    for d,x in zip((1,2,3,6),e):
        if x: r=mul(r,power(Ed(d),x*scale))
    return r
t0=time.time()
t=etaq(te); t=np.concatenate(([0],t[:-1]))
G=etaq(F0e,k//2)
low=np.zeros((n,3*n+5),dtype=np.int64); atQ={Q:np.zeros((n,n if FULL else 1),dtype=np.int64) for Q in Qs}
b=G.copy()
for i in range(n):
    low[i]=b[:3*n+5]
    for Q in Qs:
        atQ[Q][i]=b[np.arange(n)*Q] if FULL else b[Q:Q+1]
    if i<n-1: b=mul(b,t)
assert all(low[i,i]==1 and not low[i,:i].any() for i in range(n))
B=nmod_mat(n,n,[int(x) for x in low[:,:n].flatten()],p); Binv=B.inv()
def op(V): return nmod_mat(n,n,[int(x) for x in V.flatten()],p)*Binv
U2=op(low[:,0:2*n:2]); U3=op(low[:,0:3*n:3])
def ident(r): return nmod_mat(r,r,[1 if i==j else 0 for i in range(r) for j in range(r)],p)
def genleft(A,c):
    r=A.nrows(); M=A-ident(r)*c
    for _ in range(11): M=M*M
    X,nul=M.transpose().nullspace()
    R=nmod_mat(nul,r,[int(X[i,j]) for j in range(nul) for i in range(r)],p)
    R,rk=R.rref(); piv=[next(j for j in range(r) if int(R[i,j])) for i in range(rk)]
    return R,piv
def restrict(A,R,piv):
    RA=R*A; r=R.nrows()
    return nmod_mat(r,r,[int(RA[i,j]) for i in range(r) for j in piv],p)
def leg(a,q): return 1 if pow(a%q,(q-1)//2,q)==1 else -1
c2=pow(4,-1,p); c3=pow(9,-1,p)
R1,p1=genleft(U2,c2); U3r=restrict(U3,R1,p1); R2,p2=genleft(U3r,c3)
C=R2*R1; d=C.nrows()                      # class basis in coordinates of the eta-quotient basis
print('l=%d weight %d: class dimension %d (rank bound %d); series to %d terms (%.0f s)'%(L,k,d,(L+12)//24,PREC,time.time()-t0),flush=True)
a1=[sum(int(C[j,i])*int(low[i,1]) for i in range(n))%p for j in range(d)]
for Q in Qs:
    aQ=[sum(int(C[j,i])*int(atQ[Q][i][1 if FULL else 0]) for i in range(n))%p for j in range(d)]
    # is (a_Q(S_j)) = c * (a_1(S_j)) for one c ?
    cs=set()
    for x,y in zip(aQ,a1):
        if y: cs.add(x*pow(y,-1,p)%p)
        elif x: cs.add('none')
    Qm2=pow(Q,-2,p); m1=leg(-1,Q)
    line='   Q=%d: a_Q on the class basis %s, a_1 %s -> one number c: %s ; Q^-2=%d, (-1|Q)=%+d, so +(-1|Q)Q^-2=%d and -(-1|Q)Q^-2=%d'%(Q,aQ,a1,sorted(cs,key=str),Qm2,m1,(m1*Qm2)%p,(-m1*Qm2)%p)
    if FULL:
        V=atQ[Q].copy(); cq=pow(Q,k-1,p)
        for m in range(0,n,Q): V[:,m]=(V[:,m]+cq*low[:,m//Q])%p
        TQ=op(V); TQc=restrict(restrict(TQ,R1,p1),R2,p2)
        line+=' ; full matrix: charpoly %s'%[int(x) for x in TQc.charpoly().coeffs()]
        line+=' scalar' if all(int(TQc[i,j])==(int(TQc[0,0]) if i==j else 0) for i in range(d) for j in range(d)) else ' NOT scalar'
    print(line,flush=True)
print('seconds',round(time.time()-t0))
