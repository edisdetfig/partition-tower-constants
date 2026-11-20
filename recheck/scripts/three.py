# the four-term law on lane l at 67, 71, 79, 83, with P_l mod l (project's U_l polynomials, sign-adjusted; equal to my own)
from flint import fmpz
import time
P=lambda n:int(fmpz(n).partitions_p())
POL={67:[15,18,1,1],71:[2,2,18,1],79:[30,55,52,1],83:[64,18,76,1]}
t0=time.time()
for l,pol in POL.items():
    res=[];ctrl=[]
    for n in (0,1,2):
        hs=[]
        for b in (1,3,5,7):
            m=l**b; N=m*n+pow(24,-1,m)
            if N>6e13: hs=None; break
            hs.append(P(N))
        if hs is None: continue
        res.append((n,sum(c*h for c,h in zip(pol,hs))%l))
        ctrl.append((n,sum(c*h for c,h in zip([pol[0]+1]+pol[1:],hs))%l))
    print('l=%d levels (1,3,5,7): law mod l by n: %s ; control (constant term +1): %s   [%.0f s]'%(l,res,ctrl,time.time()-t0),flush=True)
