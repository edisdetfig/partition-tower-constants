# Independent recheck: the laws for 37..73 on exact partition numbers, with the project's polynomials P_l (taken from its files).
import sys, time, json
from flint import fmpz
P=lambda n:int(fmpz(n).partitions_p())
POL={37:[13715640031306785055086320854257,12718641364044840723509473763747,1],
41:[102634213897242228709985159985173,33113711977509355733826871033800,1],
43:[250089566708442925599949223613440,456265882141526569232945444006852,1],
47:[2745333515965138293747935832325324,451629498878395392752190439239349,1],
53:[17579347254392346885913618476050148,17989064678683986962258330439435519,1],
59:[90287726400164071080001336355303676,82896254441341078936689312668002993,1],
61:[216620822251427698494782631324931280,214772260455710589155236120068070578,367333559700715629968018781117360526,1],
73:[640148276439656497,1391144417610593507,1274766274151810209,1]}
CAP1=float(sys.argv[1]) if len(sys.argv)>1 else 1.5e14   # largest single argument
CAPN=3e13                                               # arguments allowed for n>=1
def val(x,l,cap=30):
    if x==0: return cap
    k=0
    while x%l==0: x//=l;k+=1
    return k
cache={}
def H(l,b,n):
    m=l**b; d=pow(24,-1,m); N=m*n+d
    if N not in cache: cache[N]=P(N)
    return cache[N]
t0=time.time(); out={}
for l,pol in POL.items():
    R=len(pol)-1; rows=[]
    for b in range(1,8):
        top=b+2*R; m=l**top; d=pow(24,-1,m)
        if d>CAP1: continue
        ns=[n for n in range(0,6) if m*n+d<=(CAP1 if n==0 else CAPN)]
        vs=[]; vbad=[]; vshort=[]
        for n in ns:
            hs=[H(l,b+2*i,n) for i in range(R+1)]
            s=sum(c*h for c,h in zip(pol,hs)); vs.append(val(s,l))
            s2=sum(c*h for c,h in zip([pol[0]+l]+pol[1:],hs)); vbad.append(val(s2,l))   # control: second digit of the constant term changed
            s3=sum(c*h for c,h in zip(pol[1:],hs[:-1])); vshort.append(val(s3,l))        # control: a law one term shorter (drop top level, shift)
        rows.append((b,[b+2*i for i in range(R+1)],ns,vs,vbad,vshort))
        print('l=%d levels %s n=%s: law divisible by l^%s ; control (P changed in digit 2): l^%s ; (%.0f s)'%(l,[b+2*i for i in range(R+1)],ns,vs,vbad,time.time()-t0),flush=True)
    out[l]=rows
json.dump({str(k):v for k,v in out.items()},open('twins.json','w'))
print('total',round(time.time()-t0))
