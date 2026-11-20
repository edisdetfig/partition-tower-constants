# Independent recheck: congruences at lucky primes on exact partition numbers (FLINT).
import sys, time
from flint import fmpz
from sympy import factorint, jacobi_symbol
P=lambda n:int(fmpz(n).partitions_p())
t0=time.time()
def classlaw(l,Q,eps,K=8):
    # 24N-1 = Q^2 * M, l | M, Q not| M, (M|Q)=eps  -> p(N) = 0 mod l ; the other sign is the control
    hits=[];ctrl=[]; m=1
    while len(hits)<K or len(ctrl)<K:
        M=l*m
        if (Q*Q*M+1)%24==0 and M%Q!=0:
            N=(Q*Q*M+1)//24; s=jacobi_symbol(M%Q,Q)
            if s==eps and len(hits)<K: hits.append((m,N,P(N)%l))
            if s==-eps and len(ctrl)<K: ctrl.append((m,N,P(N)%l))
        m+=1
    z=sum(1 for h in hits if h[2]==0); zc=sum(1 for c in ctrl if c[2]==0)
    print('l=%d Q=%d eps=%+d: members zero mod l: %d of %d (largest N=%d); controls (other sign) zero: %d of %d  [%.0f s]'%(l,Q,eps,z,len(hits),max(h[1] for h in hits),zc,len(ctrl),time.time()-t0),flush=True)
for l,Q,eps in ((37,6599,-1),(37,7541,1),(37,9547,-1),(43,4049,-1),(43,8941,-1)):
    classlaw(l,Q,eps)
# Ono type at 37, Q=2999: 24N-1 = 37*2999^3*m, 2999 not| m
l,Q=37,2999; res=[]
for m in (13,37,61):
    X=l*Q**3*m; assert (X+1)%24==0; N=(X+1)//24; res.append((m,N,P(N)%l))
print('l=37 Q=2999 (both multipliers 0): ',res,' [%.0f s]'%(time.time()-t0),flush=True)
# controls for 2999: 24N-1 = 37*2999^2*m (one power of Q fewer; no law expected)
res=[]
for m in (11,35,59,83):
    X=l*Q**2*m
    if (X+1)%24: continue
    N=(X+1)//24; res.append((m,N,P(N)%l))
print('   controls with 2999^2 only:',res,flush=True)
# the mod-67 example
mod=19728144816847346813; N0=6177904563614
f=factorint(mod); print('67 example: modulus factors',f,flush=True)
Q=[q for q in f if q!=67][0]
X=24*N0-1; print('   24*N0-1 =',factorint(X) if X<10**16 else X, ' divisible by Q^2:',X%(Q*Q)==0,' by 67:',X%67==0)
M=X//(Q*Q); print('   M=',M,' (M|Q)=',jacobi_symbol(M%Q,Q))
t=time.time(); v=P(N0)%67; print('   p(N0) mod 67 =',v,' [%.0f s]'%(time.time()-t),flush=True)
# a control in the same family with the other sign (next m)
eps=jacobi_symbol(M%Q,Q); m=M//67+1; found=0
while found<2:
    Mm=67*m
    if (Q*Q*Mm+1)%24==0 and Mm%Q and jacobi_symbol(Mm%Q,Q)==-eps and (Q*Q*Mm+1)//24<4e13:
        N=(Q*Q*Mm+1)//24; t=time.time(); print('   control m=%d N=%d p mod 67 = %d [%.0f s]'%(m,N,P(N)%67,time.time()-t),flush=True); found+=1
    m+=1
    if m>M//67+200: break
print('total',round(time.time()-t0))
