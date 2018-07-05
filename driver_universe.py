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

#u = Universe(params=None)
u = UnivPlanck15()

# LSST souce sample
w_g = WeightTracerLSSTSources(u, name='glsst')
# lensing kernel for the LSSt sources
w_s = WeightLensCustom(u, w_g.dndz, zMin=w_g.zMin, zMax=w_g.zMax, name='slsst')

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
