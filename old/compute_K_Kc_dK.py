from numpy import *

##################################################################################
#

# number of bins
n= 101

# min and max for bin edges
lmin = 1.e-3
lmax = 1.e2

# geometric factor
alpha = (lmax/lmin)**(1./n)

# bin edges and centers
L = zeros(n+1)
Lc = zeros(n)
dL = zeros(n)

# compute L, Lc, dL
L[0] = lmin
for i in range(n):
   L[i+1] = L[i] * alpha
   # bin center, "log average of l", ie exp(<ln(l)>)
   a = L[i]
   b = L[i+1]
   Lc[i] = exp( (b*log(b) - a*log(a))/(b-a) - 1. )
   # bin width
   dL[i] = L[i+1] - L[i]


savetxt("./input/Kc.txt", Lc)
savetxt("./input/dK.txt", dL)
savetxt("./input/K.txt", L)   # this array has one more element than the two others
