# Dimension of the class at weights l-3, 2l-4, 3l-5, 4l-6: newforms on Gamma0(6) with U_p = sign*p^((k-2)/2),
# sign(p,k) = (p|l)^(1 + (k-(l-3))/(l-1)) as in fading_forms4.py, i.e. w_p = -sign. Counted with PARI's mf package.
import cypari2, sys
pari=cypari2.Pari(); pari.allocatemem(2*10**9, silent=True)
def leg(a,l): return 1 if pow(a,(l-1)//2,l)==1 else -1
def classdim(l,k):
    x=(k-(l-3))//(l-1)
    w2=-(leg(2,l)**(1+x)); w3=-(leg(3,l)**(1+x))
    mf=pari.mfinit([6,k],0)                       # newspace
    if pari.mfdim(mf)==0: return 0
    sp=pari.mfsplit(mf); dims=[len(pari.Vec(v)) for v in sp[0]] if False else None
    # count by Atkin-Lehner eigenvalues of each Galois orbit of newforms
    n=0
    orbits=pari.mfsplit(mf,0,1)          # [vectors, polys]
    polys=orbits[1]
    A2=pari.mfatkineigenvalues(mf,2); A3=pari.mfatkineigenvalues(mf,3)
    for i in range(len(polys)):
        d=int(pari.poldegree(polys[i]))
        e2=int(A2[i][0]); e3=int(A3[i][0])
        if e2==w2 and e3==w3: n+=d
    return n
out=[]
for l in [int(a) for a in sys.argv[1:]]:
    ds=[classdim(l,l-3+(l-1)*x) for x in range(2)]+[0,0]
    out.append((l,ds,[ds[i+1]-ds[i] for i in range(3)], round((l-1)/24,2)))
    print('l=%3d  class dims at l-3, 2l-4, 3l-5, 4l-6: %-18s gains per step: %-12s (l-1)/24 = %.2f  R_l=%d'%(l,ds,[ds[i+1]-ds[i] for i in range(3)],(l-1)/24,(l+12)//24),flush=True)
