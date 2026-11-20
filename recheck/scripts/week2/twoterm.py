# Cross-check (not the requested procedure): fit H_(b+4) = sigma H_(b+2) - pi H_b from pairs (n1,n2); beta = the non-unit
# root of x^2 - sigma x + pi. In a two-term tower beta = E_(b+2)/E_b exactly.
import json, sys
from itertools import combinations
P={}
for f in sys.argv[1:]:
    for l,b,n,N,v,t in json.load(open(f)): P[(l,b,n)]=v
D=60
def vl(x,l,cap):
    if x==0: return cap
    k=0
    while x%l==0 and k<cap: x//=l; k+=1
    return k
def digs(v,l,n): return [(v//l**i)%l for i in range(n)]
def div(a,b,l,prec):   # a/b in Z_l mod l^(prec - v(b)), or None
    vb=vl(b,l,prec)
    if vb>=prec or vl(a,l,prec)<vb: return None,0
    pr=prec-vb; m=l**pr
    return (a//l**vb)*pow((b//l**vb)%m,-1,m)%m, pr
for l in (7,11):
    e=(lambda b:b//2+1) if l==7 else (lambda b:b)
    H={(b,n):(v//l**e(b)) for (L,b,n),v in P.items() if L==l}
    bs=sorted({b for b,n in H}); prec=D-e(bs[-1])
    print('l = %d'%l)
    for b in bs[:-2]:
        ns=[n for n in range(8) if (b+4,n) in H]
        betas=[]
        for n1,n2 in combinations(ns,2):
            A=[[H[(b+2,n1)],-H[(b,n1)]],[H[(b+2,n2)],-H[(b,n2)]]]; y=[H[(b+4,n1)],H[(b+4,n2)]]
            Dt=A[0][0]*A[1][1]-A[0][1]*A[1][0]
            s,ps=div(y[0]*A[1][1]-A[0][1]*y[1],Dt,l,prec); p_,pp=div(A[0][0]*y[1]-A[1][0]*y[0],Dt,l,prec)
            if s is None or p_ is None: continue
            pr=min(ps,pp); m=l**pr
            if s%l==0: continue
            a=s
            for _ in range(200): a=(s-p_*pow(a,-1,m))%m
            be=(s-a)%m; vb=vl(be,l,pr)
            if vb>=pr: continue
            betas.append((vb,(be//l**vb)%l**(pr-vb),pr-vb))
        if not betas: print('   b=%2d: no usable pair'%b); continue
        vs=sorted({x[0] for x in betas}); g=min(x[2] for x in betas)
        sh=0
        while sh<g and len({(x[1]//l**sh)%l for x in betas})==1: sh+=1
        print('   b=%2d: %d pairs; v(beta) = %s; unit digits shared by all pairs: %s  (precision %d)'%(b,len(betas),vs,digs(betas[0][1],l,sh),g))
