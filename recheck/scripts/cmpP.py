# compare P_l mod l (project) with (3|l)-twisted charpoly of my T_l mod l on the class
import subprocess, re, os
PROJ={37:[29,25,1],41:[16,28,1],43:[14,6,1],47:[19,18,1],53:[36,3,1],59:[57,53,1],61:[34,15,31,1],
      67:[15,18,1,1],71:[2,2,18,1],73:[29,71,65,1],79:[30,55,52,1],83:[64,18,76,1],97:[46,3,63,22,1]}
def leg(a,p): return 1 if pow(a%p,(p-1)//2,p)==1 else -1
def pmul(a,b,p):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]=(r[i+j]+x*y)%p
    return r
ok=0
for l,P in PROJ.items():
    out=subprocess.run(['python3','-I',os.path.join(os.path.dirname(os.path.abspath(__file__)),'vr2.py'),str(l)],capture_output=True,text=True,env=dict(os.environ,OMP_NUM_THREADS='1')).stdout
    cp=[1]
    for m in re.finditer(r'^FAC \[([0-9, ]+)\] (\d+)',out,re.M):
        f=[int(x) for x in m.group(1).split(',')]
        for _ in range(int(m.group(2))): cp=pmul(cp,f,l)
    s=leg(3,l); d=len(cp)-1
    tw=[(c*(s**(d-i)))%l for i,c in enumerate(cp)]      # roots multiplied by s
    same=(tw==[c%l for c in P]); ok+=same
    print(l,'mine (twisted):',tw,' project:',P,' equal:',same)
print(ok,'of',len(PROJ))
