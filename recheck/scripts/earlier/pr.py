import json,re,sys
from flint import fmpz
P=lambda n:int(fmpz(n).partitions_p())
D=json.loads(open('prism.json').read())
L=37; ok=dict(sum=0,A=0,B=0,lawA=0,lawB=0,c=0,c25=0,n=0,nl=0,np=0)
for r in D['rows']:
    lv=r['levels']
    for j,x in enumerate(lv):
        ok['n']+=1
        ok['sum']+=(x['cA']+x['cB']-x['c'])%L==0
        ok['A']+=(x['c25']+(x['kappa']-32)*x['c']-6*x['cA'])%L==0
        ok['B']+=(x['c25']+(x['kappa']-1)*x['c']+6*x['cB'])%L==0
        if j+1<len(lv):
            ok['nl']+=1; ok['lawA']+=(15*x['cA']-lv[j+1]['cA'])%L==0; ok['lawB']+=(34*x['cB']-lv[j+1]['cB'])%L==0
        b=x['b']
        if b<=5:
            m=L**b; N=m*r['n']+pow(24,-1,m); M=24*N-1
            ok['np']+=1; ok['c']+=P(N)%L==x['c']; ok['c25']+=P((25*M+1)//24)%L==x['c25']
print(ok)
