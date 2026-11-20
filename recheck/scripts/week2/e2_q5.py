import cypari2,sys,time
pari=cypari2.Pari(); pari.allocatemem(4*10**9)
def kron(a,p): return int(pari.kronecker(a,p))
primes=[13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97]
for l in primes:
    t=time.time(); k=l-3
    w2=-kron(2,l); w3=-kron(3,l)
    mf=pari.mfinit([6,k],0)
    sp=pari.mfsplit(mf)
    forms=sp[0]; pols=sp[1]
    tot=None; n=0; orb=[]
    eb=pari.mfeigenbasis(mf)
    for i,f in enumerate(eb):
        ev2=pari.mfatkineigenvalues(mf,2)[i]; ev3=pari.mfatkineigenvalues(mf,3)[i]
        if int(pari.round(pari.real(ev2[0])))!=w2 or int(pari.round(pari.real(ev3[0])))!=w3: continue
        P=pols[i]; d=int(pari.poldegree(P)); n+=d; orb.append(d)
        a5=pari.mfcoef(f,5)
        cp=pari.charpoly(pari.Mod(pari.lift(a5),P)) if d>1 else pari('x')-a5
        tot=cp if tot is None else tot*cp
    R=(l+12)//24
    # polynomial with roots a^2
    sq=pari.polresultant(pari.subst(tot,'x','t'),pari('x-t^2'),'t')
    if pari.pollead(sq)<0: sq=-sq
    det=pow(5,k-1,l); deti=pow(det,-1,l)
    # roots a^2/det: scale
    m=pari.Mod(1,l)
    sqm=pari.subst(sq*m,'x',pari('x')*pari.Mod(det,l))
    sqm=sqm/pari.pollead(sqm)
    g=pari.gcd(sqm,pari.deriv(sqm))
    fac=pari.factormod(tot,l)
    facs=pari.factormod(sqm,l)
    rts=[int(pari.lift(r)) for r in pari.polrootsmod(tot,l)]
    fd=lambda M:[(int(pari.poldegree(M[j,0])),int(M[j,1])) for j in range(M.nrows())]
    print(l,'k=',k,'signs',w2,w3,'forms',n,'R',R,'orbits',orb,'| a5 mod l: degs',fd(fac),'roots',rts,'| a5^2/5^(k-1): degs',fd(facs),'roots',[int(pari.lift(r)) for r in pari.polrootsmod(sqm,l)],'squarefree',int(pari.poldegree(g))==0,'(%.0fs)'%(time.time()-t),flush=True)
