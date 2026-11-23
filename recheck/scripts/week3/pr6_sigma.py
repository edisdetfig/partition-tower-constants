# NOTE.md (R2): sigma_p(T) = p^((l-5)/2) (1+T)^(c_p/2), where p = omega(p)<p> and <p> = (1+l)^(c_p).
# Check l-adically (to l^12) that at weight K = l-3+(l-1)x:  sigma_p(T_K) = (p|l)^x p^((K-2)/2),  square = p^(K-2).
from sympy import jacobi_symbol
def teich(p,l,N):
    w=p%l**N
    for _ in range(N+2): w=pow(w,l,l**N)
    return w
def logbase(u,l,N):                      # c with (1+l)^c = u mod l^N, for u = 1 mod l
    c=0; cur=1; g=1+l
    for i in range(N-1):
        gi=pow(g,l**i,l**N)
        for dgt in range(l):
            if (u*pow(cur,-1,l**N)*pow(pow(gi,dgt,l**N),-1,l**N)-1)%(l**(i+2))==0:
                c+=dgt*l**i; cur=cur*pow(gi,dgt,l**N)%l**N; break
    return c
bad=0; n=0
for l in (13,17,29,37,43,499):
    N=12; m=l**N
    for p in (2,3):
        w=teich(p,l,N); u=p*pow(w,-1,m)%m; c=logbase(u,l,N)
        assert pow(1+l,c,m)==u
        half=c*pow(2,-1,l**(N-1))%l**(N-1)
        for x in range(0,7):
            K=l-3+(l-1)*x
            sig=pow(p,(l-5)//2,m)*pow(pow(1+l,K-(l-3),m),half,m)%m
            want=(int(jacobi_symbol(p,l))**x*pow(p,(K-2)//2,m))%m
            n+=1
            if sig!=want or sig*sig%m!=pow(p,K-2,m): bad+=1
print('sigma_p checks:',n,'failures:',bad)
