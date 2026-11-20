import time
from flint import fmpz
t=time.time()
N=9040523964878168
v=fmpz(N).partitions_p()
print('p(%d) mod 73 ='%N, int(v%73), ' seconds', round(time.time()-t), flush=True)
