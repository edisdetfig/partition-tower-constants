from flint import nmod_poly
p=607
f=nmod_poly([458,176,331,412,360,410,588,195,1],p)
print("factor degrees:",sorted(g.degree() for g,e in f.factor()[1] for _ in range(e)))
def mulmod(a,b): return (a*b)%f
def powmod(a,n):
    r=nmod_poly([1],p)
    while n:
        if n&1: r=mulmod(r,a)
        a=mulmod(a,a); n>>=1
    return r
given=[[0,1,0,0,0,0,0,0],[445,159,605,585,381,90,31,467],[186,254,264,494,242,520,82,595],[604,548,315,605,490,284,248,77],[563,562,533,364,255,124,357,37],[438,567,88,444,563,196,196,131],[469,281,484,211,220,591,380,283],[135,56,139,332,277,16,527,231]]
s=nmod_poly([0,1],p); ok=True; tot=nmod_poly([0],p); prod=nmod_poly([1],p)
for i in range(8):
    c=[int(x) for x in s.coeffs()]+[0]*8; c=c[:8]
    print(i+1,c==given[i]); ok&=c==given[i]
    tot=tot+s; prod=mulmod(prod,s); s=powmod(s,p)
print("back to t:",[int(x) for x in s.coeffs()]==[0,1], "sum",tot.coeffs(),"prod",prod.coeffs())
