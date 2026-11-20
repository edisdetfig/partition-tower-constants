# Traces of a2, a3, a5 for every newform orbit of level 6 at the weights l-3, l = 13 ... 53,
# to compare by eye with the LMFDB pages lmfdb.org/ModularForm/GL2/Q/holomorphic/6/<weight>/a/
# The orbit marked "class" has the signs w2 = -(2|l), w3 = -(3|l).
import cypari2
pari=cypari2.Pari(); pari.allocatemem(10**9)
for l in (13,17,19,23,29,31,37,41,43,47,53):
    k=l-3; w2=-int(pari.kronecker(2,l)); w3=-int(pari.kronecker(3,l))
    mf=pari.mfinit([6,k],0); eb=pari.mfeigenbasis(mf); pols=pari.mfsplit(mf)[1]
    for i,f in enumerate(eb):
        e2=int(pari.round(pari.real(pari.mfatkineigenvalues(mf,2)[i][0]))); e3=int(pari.round(pari.real(pari.mfatkineigenvalues(mf,3)[i][0])))
        P=pols[i]; d=int(pari.poldegree(P))
        tr=lambda n: int(pari.trace(pari.Mod(pari.lift(pari.mfcoef(f,n)),P))) if d>1 else int(pari.mfcoef(f,n))
        mark='  <-- class' if (e2,e3)==(w2,w3) else ''
        print('l=%d k=%d orbit dim %d signs (%+d,%+d) tr a2=%d tr a3=%d tr a5=%d  (tr a5 mod l = %d)%s'%(l,k,d,e2,e3,tr(2),tr(3),tr(5),tr(5)%l,mark))
