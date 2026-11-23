from flint import fmpz
from sympy import jacobi_symbol, factorint
def p(n): return fmpz(n).partitions_p()
# 1. lambda_b parity and K_b formula
bad=0
for l in [13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,97,101,499,607]:
    for b in range(1,9):
        m=l**b; beta=pow(24,-1,m); g=(24*beta-1)//m
        assert 0<g<24 and (g+pow(l,-b,24))%24==0
        lam=((l**b-l**(b-1))//2 if b%2 else (l**b-l**(b-1))) - (g+1)//2
        K=g+2*lam-1
        Kf=(l-1)*l**(b-1)-2 if b%2 else 2*(l-1)*l**(b-1)-2
        s2=-int(jacobi_symbol(8,g)); s3=-int(jacobi_symbol(12,g))
        e2=-int(jacobi_symbol(2,l))**b; e3=-int(jacobi_symbol(3,l))**b
        if lam%2 or lam<0 or K!=Kf or K<10 or (K//2-1)<b or (s2,s3)!=(e2,e3) or (K+2)%((l-1)*l**(b-1)): bad+=1; print('BAD',l,b,g,lam,K)
print('1. lambda_b even, K_b formula, K_b>=10, K_b/2-1>=b, signs: failures =',bad)
# 2. BW relations at 29
d={b:pow(24,-1,29**b) for b in range(1,8)}
r1=[(n,(p(29*n+23)-10*p(29**3*n+X))%29) for X in (806,d[3]) for n in range(3)]
print('2a. p(29n+23)-10p(29^3n+X) mod 29, X=806 then 19308:',[int(x[1]) for x in r1])
r2=[int((p(29**3*n+X)-329*p(29**5*n+Y))%841) for (X,Y) in ((806,19308),(d[3],d[5])) for n in range(3)]
print('2b. p(29^3n+X)-329p(29^5n+Y) mod 29^2, (806,19308) then (19308,16237993):',r2)
m3=29**3
r3a=[int((p(29**5*n+Y)-14706*p(29**7*n+d[7]))%m3) for Y in (19308,d[5]) for n in range(2)]
r3b=[int((p(29**7*n+d[7])-14706*p(29**5*n+d[5]))%m3) for n in range(2)]
print('2c. as printed p(29^5n+Y)-14706p(29^7n+d7) mod 29^3, Y=19308 then 16237993:',r3a)
print('2d. reversed p(29^7n+d7)-14706p(29^5n+d5) mod 29^3:',r3b)
# 3. arithmetic of the sample congruences
print('3. 73*11367817^3 ==',73*11367817**3==107239327321929851108449,'; 67*665279^3 ==',67*665279**3==19728144816847346813)
for l,Q,N0 in ((67,665279,6177904563614),(73,11367817,9040523964878168)):
    M=(24*N0-1); print('   ',l,Q,'24N0-1 divisible by Q^2:',M%(Q*Q)==0,'by l:',M%l==0,'(M/Q^2|Q)=',int(jacobi_symbol((M//(Q*Q))%Q,Q)),'(l|Q)=',int(jacobi_symbol(l,Q)))
print('   37*2999^3*13+1 == 24*540584055355:',37*2999**3*13+1==24*540584055355)
for Q in (4049,8941): print('   (43|%d) ='%Q,int(jacobi_symbol(43,Q)),' Q^-2 mod 43 =',pow(Q,-2,43),' (-1|Q)=',int(jacobi_symbol(-1,Q)),' (-3|Q)=',int(jacobi_symbol(-3,Q)),' (3|Q)=',int(jacobi_symbol(3,Q)))
print('4. p(540584055355) mod 37 =',int(p(540584055355)%37))
