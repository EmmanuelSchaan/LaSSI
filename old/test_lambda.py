import numpy as np
import copy

nBins = 3


def createBins():
   print "dict"
   w = np.empty(nBins, dtype=object)


#   def f0(x): return 0
#   def f1(x): return 1
#   def f2(x): return 2


   for iBin in range(nBins):
      w[iBin] = lambda x: iBin*1.
      
#      if iBin==0:
#         w[iBin] = f0
#      if iBin==1:
#         w[iBin] = f1
#      if iBin==2:
#         w[iBin] = f2
#
#      print w[0], w[1], w[2]


   for iBin in range(nBins):
      print w[iBin](0.)

   

   return w

w = createBins()

#w = createBins()
for iBin in range(nBins):
   print w[iBin](0.)
for iBin in range(nBins):
   print w[iBin]
