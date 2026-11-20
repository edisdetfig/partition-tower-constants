from flint import fmpz
from sympy import jacobi_symbol
import time
P=lambda n:int(fmpz(n).partitions_p())
l=67; Q=665279; N0=6177904563614; M0=(24*N0-1)//(Q*Q); eps=int(jacobi_symbol(M0%Q,Q))
print('M0 =',M0,'= 67 *',M0//67,' (M0|Q) =',eps)
found=0; m=1
while found<3:
    M=l*m
    if (Q*Q*M+1)%24==0 and M%Q:
        s=int(jacobi_symbol(M%Q,Q)); N=(Q*Q*M+1)//24
        if N<4e13 or True:
            t=time.time(); v=P(N)%l; print('m=%d N=%d (M|Q)=%+d  p(N) mod 67 = %d  [%.0f s]'%(m,N,s,v,time.time()-t),flush=True); found+=1
    m+=1
