import json, sys
from fractions import Fraction
from itertools import combinations
def leg(a,l): return 1 if pow(a,(l-1)//2,l)==1 else -1
def digs(v,l,n): return [(v//l**i)%l for i in range(n)]
def lag(pts,l,m):
    tot=0
    for xa,ya in pts:
        nu=de=1
        for xb,_ in pts:
            if xb!=xa: nu*=(-1-xb); de*=(xa-xb)
        tot+=ya*nu*pow(de%m,-1,m)
    return tot%m
def agree(a,b,l,n):
    k=0
    while k<n and (a//l**k)%l==(b//l**k)%l: k+=1
    return k
for l in (7,11):
    R=json.load(open('B%d.json'%l)); s=4 if l==7 else 3
    pr=min(r[4] for r in R); m=l**pr
    print('l = %d, s = %d'%(l,s))
    for x,k,v,u,p in R: print('   x=%2d k=%3d  v(a_l)=%d   r(x)=a_l/l^%d digits: %s'%(x,k,v,s,digs(u,l,8)))
    pts=[(r[0],r[3]%m) for r in R]
    full=leg(3,l)*lag(pts,l,m)%m
    print('   (3|%d) = %d;  (3|l)*P(-1) digits: %s'%(l,leg(3,l),digs(full,l,12)))
    ag=[agree(leg(3,l)*lag([q for q in pts if q[0]!=x],l,m)%m,full,l,40) for x,_ in pts]
    print('   leaving out one point, digits unchanged:',dict(zip([q[0] for q in pts],ag)),' -> all agree on',min(ag))
    # also: fit through the first j points, to see convergence
    for j in range(2,len(pts)+1):
        print('     through x=%d..%d: %s'%(pts[0][0],pts[j-1][0],digs(leg(3,l)*lag(pts[:j],l,m)%m,l,10)))
