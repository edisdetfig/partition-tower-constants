import cypari2
pari=cypari2.Pari(); pari.allocatemem(10**9)
mf=pari.mfinit([6,6],0); eb=pari.mfeigenbasis(mf); print('newforms at level 6 weight 6:',len(eb))
a=[int(x) for x in pari.mfcoefs(eb[0],870)]
print('h6 = q %+dq^2 %+dq^3 %+dq^4 %+dq^5 %+dq^6'%tuple(a[2:7]),'; Atkin-Lehner (2,3):',[int(pari.round(pari.real(pari.mfatkineigenvalues(mf,p)[0][0]))) for p in (2,3)])
print('a_859 mod 859 =',a[859]%859,'; a_p mod p at 853, 857, 863 =',[a[p]%p for p in (853,857,863)])
print('primes 5 <= p <= 863 with p | a_p(h6):',[int(p) for p in pari.primes([5,863]) if a[int(p)]%int(p)==0])
h3=[int(x) for x in pari.mfcoefs(pari.mfeigenbasis(pari.mfinit([3,6],0))[0],4)]
print('h3 = q %+dq^2 %+dq^3'%(h3[2],h3[3]),'; dim M_6(Gamma0(6)) =',int(pari.mfdim([6,6],4)),'; dim S_6 at levels 1, 2 =',int(pari.mfdim([1,6],1)),int(pari.mfdim([2,6],1)),'; Eisenstein dim =',int(pari.mfdim([6,6],3)))
