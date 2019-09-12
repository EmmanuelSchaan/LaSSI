from headers import *

def function(i):
   return 0., 1.


with sharedmem.MapReduce(np=nProc) as pool:
   result = np.array(pool.map(function, range(4)))
