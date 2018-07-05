import universe
reload(universe)
from universe import *

import projection_kernel
reload(projection_kernel)
from projection_kernel import *

import pn_2d
reload(pn_2d)
from pn_2d import *

import covp_2d
reload(covp_2d)
from covp_2d import *

##################################################################################

nBins = 2
nG = nBins
nS = nBins



##################################################################################

#u = Universe(params=None)
u = UnivPlanck15()

# LSST souce sample
w_glsst = WeightTracerLSSTSources(u, name='glsst')

# split it into bins
zBounds = w_glsst.splitBins(nBins)

# generate the corresponding tracer and shear bins
print "Generate tracer and shear bins"
w_g = {}
w_s = {}
for iBin in range(nBins):
   zMin = zBounds[iBin]
   zMax = zBounds[iBin+1]
   w_g[iBin] = WeightTracerLSSTSources(u, zMin=zMin, zMax=zMax, ngal=w_glsst.ngal_per_arcmin2/nBins, name='g'+str(iBin))
   w_s[iBin] = WeightLensCustom(u, w_glsst.dndz, zMinG=zMin, zMaxG=zMax, name='s'+str(iBin))
   print "- done "+str(iBin+1)+" of "+str(nBins)
print "total ngal="+str(np.sum([w_g[i].ngal_per_arcmin2 for i in range(nBins)]))+"/arcmin2, should be "+str(w_glsst.ngal_per_arcmin2)

# gg
print "Compute p_gg"
tStart = time()
p2d_gg = {}
for iBin in range(nBins):
   p2d_gg[iBin] = P2d(u, u, w_g[iBin], fPnoise=lambda l:1./w_g[iBin].ngal, doT=False, name='', nProc=1, save=True)
tStop = time()
print "took "+str(tStop-tStart)+" sec"

# ss
print "Compute p_ss"
tStart = time()
p2d_ss = {}
for iBin1 in range(nBins):
   # cross correlation: different shear bins
   for iBin2 in range(iBin1):
      p2d_ss[iBin1, iBin2] = P2d(u, u, w_s[iBin1], w_s[iBin2], doT=False, name='', nProc=1, save=True)
   # auto-correlation: same shear bin
   p2d_ss[iBin1, iBin1] = P2d(u, u, w_s[iBin1], fPnoise=lambda l:0.26**2/w_s[iBin1].ngal, doT=False, nProc=1, save=True)
tStop = time()
print "took "+str(tStop-tStart)+" sec"

# gs
print "Compute p_gs"
tStart = time()
p2d_gs = {}
for iBin1 in range(nBins):
   # shear bin should be higher z than galaxy bin
   for iBin2 in range(iBin1,nBins):
      p2d_gs[iBin1, iBin2] = P2d(u, u, w_g[iBin1], w_s[iBin2], doT=False, name='', nProc=1, save=True)
tStop = time()
print "took "+str(tStop-tStart)+" sec"



'''
p2d_gg = P2d(u, u, w_g, fPnoise=lambda l:1./w_g.ngal, doT=False, nProc=1, save=False)
p2d_ss = P2d(u, u, w_s, fPnoise=lambda l:0.26**2/w_s.ngal, doT=False, nProc=1, save=False)
p2d_gs = P2d(u, u, w_g, w_s, doT=False, nProc=1, save=False)


cov_gg_gg = CovP2d(p2d_gg, p2d_gg, p2d_gg, p2d_gg, T2d=None, HSVP2d=None, fsky=1., npairs_id='log')
#cov_gg_gg.plot(P2d=p2d_gg)
print "total SNR ss =", cov_gg_gg.snrG(P2d=p2d_gg)
#
cov_ss_ss = CovP2d(p2d_ss, p2d_ss, p2d_ss, p2d_ss, T2d=None, HSVP2d=None, fsky=1., npairs_id='log')
#cov_ss_ss.plot(P2d=p2d_ss)
print "total SNR gg =", cov_ss_ss.snrG(P2d=p2d_ss)
#
cov_gs_gs = CovP2d(p2d_gg, p2d_ss, p2d_gs, p2d_gs, T2d=None, HSVP2d=None, fsky=1., npairs_id='log')
#cov_gs_gs.plot(P2d=p2d_gs)
print "total SNR gs =", cov_gs_gs.snrG(P2d=p2d_gs)
'''
