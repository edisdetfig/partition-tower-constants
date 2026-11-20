import cypari2, sys, json
pari=cypari2.Pari(); pari.allocatemem(3*10**9, silent=True); pari('read("sysB.gp")')
l=int(sys.argv[1]); PREC=90
if l==7:  ks=[(x,6*x+4,1,(-1)**(x+1)) for x in range(6,13)];  want={5:5,11:6,13:0}; s=4
else:     ks=[(x,10*x+8,(-1)**(x+1),1) for x in range(4,12)]; want={5:0,7:4,13:1}; s=3
def toint(a,prec):   # p-adic -> integer mod l^prec (a must be integral)
    return int(pari.lift(pari.Mod(pari.truncate(a),pari(l)**prec)))
res=[]
for x,k,e2,e3 in ks:
    d,r,Q,sy,deg=pari("Bsys")(l,k,e2,e3,PREC)
    Q=[int(q) for q in Q]; hits=[]
    for ev,al in sy:
        resid={q:(toint(ev[i],1) if pari.valuation(ev[i],l)>=0 else None) for i,q in enumerate(Q) if q in want}
        if resid==want: hits.append(al)
    print('x=%2d k=%3d: newforms %3d, with these U2,U3 signs %2d, of which Q_%d-rational systems %d; matching the residues: %d'%(x,k,int(d),int(r),l,len(sy),len(hits)),flush=True)
    if len(hits)!=1: res.append([x,k,None,None]); continue
    al=hits[0]; v=int(pari.valuation(al,l)); u=al/pari(l)**s
    pr=int(pari.padicprec(u,l)); prec=min(pr,PREC-10)
    res.append([x,k,v,toint(u,prec),prec])
    print('      v_%d(a_%d) = %d;  a_l/l^%d known to %d digits'%(l,l,v,s,prec),flush=True)
json.dump(res,open('B%d.json'%l,'w'))
