import cypari2,sys,time
pari=cypari2.Pari(); pari.allocatemem(3*10**9)
from sympy import primerange
lo,hi=int(sys.argv[1]),int(sys.argv[2])
for l in primerange(lo,hi+1):
    t=time.time(); k=l-3
    w2=-int(pari.kronecker(2,l)); w3=-int(pari.kronecker(3,l))
    mf=pari.mfinit([6,k],0)
    d=int(pari.mfdim(mf))
    A2=pari.mfatkininit(mf,2); A3=pari.mfatkininit(mf,3)
    M2=A2[1]; M3=A3[1]
    # normalise: W^2 = scalar; find scalar c with M^2 = c*Id, eigenvalues +-sqrt(c)
    T5=pari.mfheckemat(mf,5)
    I=pari.matid(d)
    c2=(M2*M2)[0,0]; c3=(M3*M3)[0,0]
    s2=pari.sqrtint(pari.numerator(c2))/pari.sqrtint(pari.denominator(c2))
    s3=pari.sqrtint(pari.numerator(c3))/pari.sqrtint(pari.denominator(c3))
    assert s2*s2==c2 and s3*s3==c3
    # eigenvalue of f is (matrix eigenvalue)/s up to A[2] scaling; kernel of (M - w*s*I)
    K=pari.matintersect(pari.matker(M2-w2*s2*I),pari.matker(M3-w3*s3*I))
    r=int(pari.matsize(K)[1])
    alt=None
    if r!=(l+12)//24:
        alt=r
        K=pari.matintersect(pari.matker(M2+w2*s2*I),pari.matker(M3+w3*s3*I)); r=int(pari.matsize(K)[1])
    # T5 restricted
    B=pari.matinverseimage(K,T5*K)
    cp=pari.charpoly(B)
    sq=pari.polresultant(pari.subst(cp,'x','t'),pari('x-t^2'),'t')
    m=pari.Mod(1,l); sqm=sq*m; sqm=sqm/pari.pollead(sqm)
    g=pari.gcd(sqm,pari.deriv(sqm))
    fa=pari.factormod(cp,l); fs=pari.factormod(sqm,l)
    fd=lambda M:sorted((int(pari.poldegree(M[j,0])),int(M[j,1])) for j in range(M.nrows()))
    print(l,'dim',d,'class',r,'R',(l+12)//24,'alt' if alt is not None else '','a5 degs',fd(fa),'a5^2 degs',fd(fs),'distinct',int(pari.poldegree(g))==0,'(%.0fs)'%(time.time()-t),flush=True)
