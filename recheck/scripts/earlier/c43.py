from flint import nmod_poly, nmod_mat
import sys
exec(open('col.py').read().split('# 24 colours at 17')[0])
def sig(k,n): return sum(d**k for d in range(1,n+1) if n%d==0)
def forms(l,k,N):
    M=l
    E4=nmod_poly([1]+[240*sig(3,n)%M for n in range(1,N)],M); E6=nmod_poly([1]+[(-504*sig(5,n))%M for n in range(1,N)],M)
    t=lambda p: trunc(p,N,M)
    D=t((t(t(E4*E4)*E4)-t(E6*E6))*nmod_poly([pow(1728,-1,M)],M))
    # monomials E4^a E6^b of weight k-12
    basis=[]
    for a in range(0,(k-12)//4+1):
        rem=k-12-4*a
        if rem%6==0:
            b=rem//6; f=D
            for _ in range(a): f=t(f*E4)
            for _ in range(b): f=t(f*E6)
            basis.append(f)
    return basis
def Tl(l,k):
    d=None; N=l*8+5; B=forms(l,k,N); d=len(B)
    co=lambda f,i: int(f.coeffs()[i]) if i<len(f.coeffs()) else 0
    A=nmod_mat(d,d,[co(B[i],j+1) for i in range(d) for j in range(d)],l)      # rows: basis coeffs q^1..q^d
    U=nmod_mat(d,d,[co(B[i],l*(j+1)) for i in range(d) for j in range(d)],l)  # U_l images
    T=U*A.inv()   # matrix of U in basis (row convention)
    return d,T
for l,k in ((43,30),(41,28),(13,12),(29,16)):
    d,T=Tl(l,k); cp=T.charpoly() if hasattr(T,'charpoly') else None
    print(l,'weight',k,'dim S',d,'T_l mod l charpoly',[int(c) for c in cp.coeffs()] if cp is not None else T.entries(), 'factors',[( [int(c) for c in f.coeffs()],e) for f,e in cp.factor()[1]])
d,T=Tl(43,30); print('T-30I =',[(int(T[i,j])-(30 if i==j else 0))%43 for i in range(2) for j in range(2)])
l=43; N=l**3*4+5; p=pr(24,N,l)
H=lambda b,n: p[l**b*n+1]
for b0 in (0,1):
    ok=sum((H(b0+2,n)-31*H(b0+1,n)+30*H(b0,n))%l==0 for n in range(1,5)); print('43 two-term law from level',b0,ok,'of 4')
l=41; N=l**3*4+5; p=pr(24,N,l)
# any 2-term law at 41? test the three candidate pairs of roots among {1,34,33}
for r1,r2 in ((1,34),(1,33),(34,33)):
    print('41 law roots',r1,r2,sum((H(3,n)-(r1+r2)*H(2,n)+r1*r2*H(1,n))%l==0 for n in range(1,5)),'of 4')
