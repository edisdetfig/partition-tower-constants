import cypari2, itertools
pari=cypari2.Pari(); pari.allocatemem(2*10**9)
kr=lambda a,p:int(pari.kronecker(a,p))
def orbits(N,k):
    mf=pari.mfinit([N,k],0); eb=pari.mfeigenbasis(mf); pols=pari.mfsplit(mf)[1]; out=[]
    for i,f in enumerate(eb):
        P=pols[i]; d=int(pari.poldegree(P))
        sg=tuple(int(pari.round(pari.real(pari.mfatkineigenvalues(mf,p)[i][0]))) for p in (2,3) if N%p==0)
        cp={q:(pari.charpoly(pari.Mod(pari.lift(pari.mfcoef(f,q)),P)) if d>1 else pari('x')-pari.mfcoef(f,q)) for q in (2,3,5,7,11,13)}
        out.append((N,d,sg,cp))
    return out
for l in (17,29,37):
    k=l-3; new=orbits(6,k); old=[o for N in (1,2,3) for o in orbits(N,k)]
    cls=(-kr(2,l),-kr(3,l))
    print('l=%d weight %d: level-6 orbits %s; lower-level orbits %s'%(l,k,[(o[1],o[2]) for o in new],[(o[0],o[1]) for o in old]))
    for o in new:
        hits=[]
        for q in old:
            # shared root mod l of the charpolys of T_5, T_7, T_11, T_13 (necessary for a congruence)
            sh=[int(pari.poldegree(pari.gcd(o[3][t]*pari.Mod(1,l),q[3][t]*pari.Mod(1,l)))) for t in (5,7,11,13)]
            if all(s>0 for s in sh):
                # level-raising condition at the primes new for o and absent from q's level
                cond=[]
                for p in (2,3):
                    if q[0]%p:
                        c=pari.polresultant(q[3][p], pari('x^2')-pari(p)**(k-2)*(p+1)**2)
                        cond.append((p,int(c%l)==0))
                hits.append(('level %d dim %d'%(q[0],q[1]),sh,cond))
        tag='  <-- the class S_l' if o[2]==cls else ''
        print('   new orbit dim %d signs %s%s: congruent lower-level orbits: %s'%(o[1],o[2],tag,hits if hits else 'none'))
