# third-digit control at 37: law at levels (5,7,9), n=0, with the true P and with P changed in its third digit
from flint import fmpz
import time
P=lambda n:int(fmpz(n).partitions_p())
l=37; pol=[13715640031306785055086320854257,12718641364044840723509473763747,1]
def val(x,l):
    k=0
    while x%l==0: x//=l;k+=1
    return k
t=time.time()
hs=[P(pow(24,-1,l**b)) for b in (5,7,9)]
s=sum(c*h for c,h in zip(pol,hs))
s3=sum(c*h for c,h in zip([pol[0]+l*l]+pol[1:],hs))
s3b=sum(c*h for c,h in zip([pol[0],pol[1]+l*l,1],hs))
print('37, levels (5,7,9), n=0: true law l^%d ; third digit of c0 changed: l^%d ; third digit of c1 changed: l^%d ; v(H5)=%d v(H7)=%d  (%.0f s)'%(val(s,l),val(s3,l),val(s3b,l),val(hs[0],l),val(hs[1],l),time.time()-t))
