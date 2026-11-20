import re,sys
from fractions import Fraction as Fr
F=sys.argv[1]
L=None; cells=[]
for line in open(F,encoding='utf-8'):
    m=re.match(r'ℓ = (\d+) ',line)
    if m: L=int(m.group(1))
    m=re.match(r'\s+(Ono|class ε=.1)\s+where χ = ([+-]\d):\s+(\d+) of\s+(\d+) primes',line)
    if m: cells.append((L,m.group(1),int(m.group(2)),int(m.group(3)),int(m.group(4))))
print(len(cells),'cells; primes',sorted(set(c[0] for c in cells)))
def stat(model,var=False):
    s=0
    for L,ty,ch,O,n in cells:
        if model=='E2': p=(L+ch)/(L*L-1)
        else: p=1/L
        E=n*p; s+=(O-E)**2/(E*(1-p) if var else E)
    return s
print('E2  sum (O-E)^2/E      = %.2f'%stat('E2'))
print('flat sum (O-E)^2/E      = %.2f'%stat('flat'))
print('E2  with binomial var   = %.2f'%stat('E2',True))
print('flat with binomial var  = %.2f'%stat('flat',True))
# 36-count statistic
tot=[];L=None
for line in open(F,encoding='utf-8'):
    m=re.match(r'\s+(class ε=.1|Ono)\s+expected\s+([\d.]+) \(derived\)\s+counted\s+(\d+)',line)
    if m: tot.append((float(m.group(2)),int(m.group(3))))
big=[(E,O) for E,O in tot if E>1]
print(len(tot),'totals,',len(big),'with E>1; chi2 (rounded E) = %.2f'%sum((O-E)**2/E for E,O in big))
# density at 37
def leg(a,p):
    a%=p
    return 0 if a==0 else (1 if pow(a,(p-1)//2,p)==1 else -1)
for l,degs in ((37,[1,1]),(41,[1,1]),(43,[1,1]),(53,[1,1])):
    tot=Fr(0); per=[]
    for kind in ('eps','eps','ono'):
        s=Fr(0)
        for Q in range(1,l):
            Qi=pow(Q,-1,l); D=pow(Qi,3,l)
            c=0 if kind=='ono' else pow(Qi,2,l)
            x=(c*c-4*D)%l
            p=Fr(1)
            for f in degs:
                q=l**f; ch=0 if x==0 else (1 if f%2==0 else leg(x,l)); p*=Fr(q+ch,q*q-1)
            s+=p
        per.append(s/(l-1)); tot+=s/(l-1)
    print(l,'density',tot,'= 1 in %.1f'%(1/tot),'; x78495 = %.1f'%(78495*tot),'x78498 = %.1f'%(78498*tot),'per kind',[ '%.2f'%(78495*x) for x in per])
