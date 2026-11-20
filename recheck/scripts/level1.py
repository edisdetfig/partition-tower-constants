# Independent recheck of two "other forms" claims: 1/E4 at 7 and 13 (pole at rho), and j at 13.
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
def dil(P,d,N,M):   # P(q^d) truncated
    c=[0]*N
    co=P.coeffs()
    for i in range(0,(N+d-1)//d):
        if i<len(co): c[i*d]=int(co[i])
    return nmod_poly(c,M)
def val(x,l,cap):
    if x==0: return cap
    k=0
    while x%l==0: x//=l;k+=1
    return k
def digits(x,l,K): 
    out=[]
    for i in range(K): out.append(x%l); x//=l
    return out
def E4series(N,M):
    e1=euler(N,M); e2=dil(e1,2,N,M)
    A=e1.pow_trunc(16,N).mul_low(e2.pow_trunc(8,N).inverse_series_trunc(N),N)           # eta(z)^16/eta(2z)^8
    Bq=e2.pow_trunc(16,N).mul_low(e1.pow_trunc(8,N).inverse_series_trunc(N),N)          # eta(2z)^16/eta(z)^8 without its q
    Bq=nmod_poly([0]+[int(x) for x in Bq.coeffs()[:N-1]],M)
    E4=A+Bq*256
    # check against 240*sigma_3
    for n in range(1,40):
        s3=sum(d**3 for d in range(1,n+1) if n%d==0)
        assert int(E4[n])==(240*s3)%M,(n,int(E4[n]))
    return E4,e1
def unitroot_m(l,K,curve,m):
    # unit root u of Frobenius of y^2 = x^3 + A x + B at l, then u^m mod l^K
    A,Bc=curve; cnt=1
    sq={}
    for y in range(l): sq[y*y%l]=sq.get(y*y%l,0)+1
    for x in range(l): cnt+=sq.get((x**3+A*x+Bc)%l,0)
    a=l+1-cnt; M=l**K
    u=a%l
    for _ in range(K+2): u=(u-(u*u-a*u+l)*pow(2*u-a,-1,M))%M
    assert (u*u-a*u+l)%M==0 and u%l
    return a,pow(u,m,M)
t0=time.time()
for l,B,K in ((7,8,14),(13,6,10)):
    M=l**K; N=l**B+2
    E4,e1=E4series(N,M); F=E4.inverse_series_trunc(N)
    a,um6=unitroot_m(l,K,(0,1),-6)
    print('1/E4 at l=%d: a_l(y^2=x^3+1)=%d; u^-6 digits %s   (%.0f s)'%(l,a,digits(um6,l,10),time.time()-t0),flush=True)
    for n in (1,2,3):
        xs=[int(F[l**b*n]) for b in range(0,B+1) if l**b*n<N]
        G=[(xs[i+1]-xs[i])%M for i in range(len(xs)-1)]      # (x-1) applied: the Eisenstein constant removed
        row=[]
        for i in range(len(G)-1):
            v0,v1=val(G[i],l,K),val(G[i+1],l,K)
            if v0>=K-1 or v1>=K-1: row.append((i,'exhausted')); continue
            Mm=l**(K-max(v0,v1))
            if v1!=v0: row.append((i,'valuations',v0,v1)); continue
            q=(G[i+1]//l**v1)*pow(G[i]//l**v0,-1,Mm)%Mm
            dg=digits(q,l,K-v1); ref=digits(um6,l,K); ag=0
            while ag<len(dg) and dg[ag]==ref[ag]: ag+=1
            row.append((i,'v',v0,'agree',ag))
        print('   n=%d: G_(b+1)/G_b vs u^-6, by level b: %s'%(n,row),flush=True)
# j at 13: c(13^(b+1) n) = u * c(13^b n), n>=1
l,B,K=13,6,12; M=l**K; N=l**B+3
E4,e1=E4series(N,M)
J=E4.pow_trunc(3,N).mul_low(e1.pow_trunc(24,N).inverse_series_trunc(N),N)     # q*j
c=lambda n:int(J[n+1])
assert c(1)==196884%M and c(2)==21493760%M
ref=[8,10,11,9,11,0,9,7,4,7,7]
for n in (1,2,3):
    row=[]
    for b in range(0,B):
        if l**(b+1)*n+1>=N: break
        x,y=c(l**b*n),c(l**(b+1)*n); vx,vy=val(x,l,K),val(y,l,K)
        if vx!=vy: row.append((b,'valuations',vx,vy)); continue
        Mm=l**(K-vx); q=(y//l**vy)*pow(x//l**vx,-1,Mm)%Mm; dg=digits(q,l,K-vx); ag=0
        while ag<min(len(dg),len(ref)) and dg[ag]==ref[ag]: ag+=1
        row.append((b,'agree',ag,dg[:ag+1]))
    print('j at 13, n=%d: c(13^(b+1)n)/c(13^b n) vs the project\'s constant %s: %s'%(n,ref[:6],row),flush=True)
print('seconds',round(time.time()-t0))
