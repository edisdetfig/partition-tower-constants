import json, sys
files=sys.argv[1:]
P={}
for f in files:
    for l,b,n,N,v,t in json.load(open(f)): P[(l,b,n)]=v
D=60
def vl(x,l,cap):
    if x==0: return cap
    k=0
    while x%l==0 and k<cap: x//=l; k+=1
    return k
def digs(v,l,n): return [(v//l**i)%l for i in range(n)]
for l in (7,11):
    e=(lambda b:b//2+1) if l==7 else (lambda b:b)
    bs=sorted({b for (L,b,n) in P if L==l})
    ns=lambda b: sorted(n for (L,bb,n) in P if L==l and bb==b)
    H={}
    for (L,b,n),v in P.items():
        if L!=l: continue
        assert v%l**e(b)==0, (l,b,n,'not divisible')
        H[(b,n)]=(v//l**e(b))%l**(D-e(b))
    prec=D-e(bs[-1])      # digits of H known
    m=l**prec
    bt=bs[-1]; nn=[n for n in ns(bt) if H[(bt-2,n)]%l]
    print('l = %d: levels b = %s; H_b known mod %d^%d'%(l,bs,l,prec))
    for b in bs: print('     v(H_%d(n)), n=0..: %s'%(b,[vl(H[(b,n)],l,prec) for n in ns(b)]))
    # alpha from the deepest pair, for each n; and agreement between n
    al={n:H[(bt,n)]*pow(H[(bt-2,n)],-1,m)%m for n in nn}
    a0=al[nn[0]]
    agr=min(vl(al[n]-a0,l,prec) for n in nn)
    print('  alpha = H_%d/H_%d: the %d values of n agree on %d digits; alpha digits: %s'%(bt,bt-2,len(nn),agr,digs(a0,l,agr)))
    # convergence of the ratio by level, for reference
    for b in bs[1:]:
        u=[n for n in ns(b) if H[(b-2,n)]%l]
        r=[H[(b,n)]*pow(H[(b-2,n)],-1,m)%m for n in u]
        print('     H_%d/H_%d agrees with alpha on %s digits (n = %s)'%(b,b-2,[vl(x-a0,l,prec) for x in r],u))
    A=a0; Ap=agr    # alpha known to Ap digits
    E={}
    for b in bs[:-1]:
        for n in ns(b+2):
            E[(b,n)]=(H[(b+2,n)]-A*H[(b,n)])%m
    print('  E_b(n) = H_(b+2) - alpha H_b, valuations (alpha good to %d digits, so values >= %d are unresolved):'%(Ap,Ap))
    for b in bs[:-1]: print('     b=%2d: %s'%(b,[vl(E[(b,n)],l,prec) for n in ns(b+2)]))
    print('  E_(b+2)/E_b: valuation, and the unit-part digits shared by all n (only digits resolved by alpha\'s precision):')
    for b in bs[:-2]:
        rows=[]; vs=[]; ok=[]
        for n in ns(b+4):
            x,y=E[(b+2,n)],E[(b,n)]; v1,v0=vl(x,l,prec),vl(y,l,prec)
            vs.append(v1-v0)
            good=Ap-v1          # digits of E_(b+2) beyond its valuation that alpha resolves
            u=(x//l**v1)*pow(y//l**v0,-1,m)%m
            rows.append((u,good))
        g=min(r[1] for r in rows)
        sh=0
        if g>0:
            while sh<g and len({(u//l**sh)%l for u,_ in rows})==1: sh+=1
        print('     b=%2d: v = %s;  shared unit digits: %s  (resolved: %d)'%(b,sorted(set(vs)),digs(rows[0][0],l,sh),max(g,0)))
