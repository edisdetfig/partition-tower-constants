# F2 recheck: orders of the constants (companion matrices of P_l) modulo l and l^2
REF={5:[4,2],7:[5,6],11:[9,3],13:[6,3],17:[6,12],19:[10,2],23:[5,20],29:[3,14],31:[30,6]}
POL={37:[13715640031306785055086320854257,12718641364044840723509473763747,1],
41:[102634213897242228709985159985173,33113711977509355733826871033800,1],
43:[250089566708442925599949223613440,456265882141526569232945444006852,1],
47:[2745333515965138293747935832325324,451629498878395392752190439239349,1],
53:[17579347254392346885913618476050148,17989064678683986962258330439435519,1],
59:[90287726400164071080001336355303676,82896254441341078936689312668002993,1],
61:[216620822251427698494782631324931280,214772260455710589155236120068070578,367333559700715629968018781117360526,1]}
def mmul(A,B,M):
    n=len(A); return [[sum(A[i][k]*B[k][j] for k in range(n))%M for j in range(n)] for i in range(n)]
def mpow(A,e,M):
    n=len(A); R=[[int(i==j) for j in range(n)] for i in range(n)]
    while e:
        if e&1: R=mmul(R,A,M)
        A=mmul(A,A,M); e>>=1
    return R
def comp(pol,M):
    n=len(pol)-1; C=[[0]*n for _ in range(n)]
    for i in range(1,n): C[i][i-1]=1
    for i in range(n): C[i][n-1]=(-pol[i])%M
    return C
def order(C,l,M):
    n=len(C); I=[[int(i==j) for j in range(n)] for i in range(n)]
    # order mod l by trial over divisors of l^n! style bound: brute force
    X=[row[:] for row in C]; X=[[x%l for x in r] for r in X]; A=[r[:] for r in X]; k=1
    while A!=I:
        A=mmul(A,X,l); k+=1
        if k>l**n+5: return None
    return k
out=[]
for l,d in REF.items():
    a=d[0]+d[1]*l; c1=next(k for k in range(1,l) if pow(a,k,l)==1); c2=next(k for k in range(1,l*l) if pow(a,k,l*l)==1)
    out.append((l,c1,c2,c2==c1*l))
for l,pol in POL.items():
    C1=comp(pol,l); c1=order(C1,l,l); C2=comp(pol,l*l); n=len(C2); I=[[int(i==j) for j in range(n)] for i in range(n)]
    back=(mpow(C2,c1,l*l)==I); full=(mpow(C2,c1*l,l*l)==I)
    out.append((l,c1,c1*l if (not back and full) else ('early' if back else '?'),(not back) and full))
for o in out: print('l=%d: order mod l = %s, order mod l^2 = %s, equals l times the first: %s'%o)
