from flint import nmod_poly
def euler(N,M):
    c=[0]*N; k=0
    while True:
        a=k*(3*k-1)//2; b=k*(3*k+1)//2
        if a>=N: break
        s=1 if k%2==0 else M-1
        c[a]=(c[a]+s)%M
        if b<N and k>0: c[b]=(c[b]+s)%M
        k+=1
    return nmod_poly(c,M)
def trunc(p,N,M): return nmod_poly(p.coeffs()[:N],M)
def pr(r,N,M):
    inv=euler(N,M).inverse_series_trunc(N); res=nmod_poly([1],M); base=inv; e=r
    while e:
        if e&1: res=trunc(res*base,N,M)
        base=trunc(base*base,N,M); e>>=1
    c=res.coeffs(); return [int(x) for x in c]+[0]*(N-len(c))
def val(x,l,cap):
    if x==0: return cap
    k=0
    while x%l==0: x//=l;k+=1
    return k
# 24 colours at 17: towers stop changing
l=17; M=l**3; N=l**3*6+5; p=pr(24,N,M)
for b in (1,2):
    ok=sum((p[l**(b+1)*n+1]-p[l**b*n+1])%l**b==0 for n in range(0,6)); print('24col 17 stop b',b,ok,'of 6',' (and mod l^(b+1):',sum((p[l**(b+1)*n+1]-p[l**b*n+1])%l**(b+1)==0 for n in range(6)),')')
# control at 13 (should NOT stop)
l=13; M=l**2; N=l**3*4+5; q=pr(24,N,M)
print('24col 13 control',sum((q[l**2*n+1]-q[l*n+1])%l==0 for n in range(1,12)),'of 11')
# 3 colours at 17: 17^ceil(b/2); 5 colours at 11: 11^(floor(b/2)+1); 5 colours at 23: factor 23 each level
for r,l,B in ((3,17,3),(5,11,3),(5,23,2),(3,11,3)):
    M=l**4; N=l**B*5+l**B; p=pr(r,N,M)
    for b in range(1,B+1):
        m=l**b; d=(r*pow(24,-1,m))%m
        vs=[val(p[m*n+d],l,4) for n in range(5)]
        print(r,'col',l,'level',b,'delta',d,'valuations',vs)
