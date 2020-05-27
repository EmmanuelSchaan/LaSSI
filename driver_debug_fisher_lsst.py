import cmb
reload(cmb)
from cmb import *

import cmb_lensing_rec
reload(cmb_lensing_rec)
from cmb_lensing_rec import *

import fisher_lsst
reload(fisher_lsst)
from fisher_lsst import *

##################################################################################
# Forecast parameters

## for debugging
#nBins = 3  # 10
#nL = 10  # 50, 20, 100
#fsky = 0.4

# actual params
nBins =  3   #10
nL = 10  #50
fsky = 0.35


# cosmological parameters to include
massiveNu = False  #True
wCDM = False #True
curvature = False #True

# priors to include
PlanckPrior = True

# include a known magnification bias
magBias = True

# forecast name
#name = "lcdm"
#name = "lcdm_mnu_curv_w0wa"
#name = "lcdm_mnu_curv_w0wa_newellsandunits"
#name = "lcdm_mnu_curv_w0wa_newellsandunits_perfectm"
#name = "gphotoz_lmaxmask"
#name = "gphotoz"
name = None

# Parallel evaluations
nProc = 4   # not actually used

##################################################################################
# Parameter classes

cosmoPar = CosmoParams(massiveNu=massiveNu, wCDM=wCDM, curvature=curvature, PlanckPrior=PlanckPrior)
#cosmoPar.plotParams()
galaxyBiasPar = GalaxyBiasParams(nBins=nBins)
#galaxyBiasPar.plotParams()
shearMultBiasPar = ShearMultBiasParams(nBins=nBins)
#shearMultBiasPar = ShearMultBiasParams(nBins=nBins, mStd=1.e-5)   # perfect photo-z priors
#shearMultBiasPar.plotParams()

# Gaussian photo-z only:
#photoZPar = PhotoZParams(nBins=nBins)
# Photo-z with Gaussian core and outliers:
photoZPar = PhotoZParams(nBins=nBins, outliers=0.1)
#photoZPar.plotParams()

#cosmoPar.plotContours()

#pat = PatPlanckParams()
#pat.printParams()

#u = Universe(cosmoPar.paramsClassy)
#u.plotDistances()


##################################################################################
# CMB lensing noise

# CMB S4
cmb = CMB(beam=1., noise=1., nu1=143.e9, nu2=143.e9, lMin=1., lMaxT=3.e3, lMaxP=5.e3, atm=False, name="cmbs4")
cmbLensRec = CMBLensRec(cmb, save=False, nProc=nProc)
fNk = cmbLensRec.fN_k_mv


##################################################################################
# Fisher calculation

import fisher_lsst
reload(fisher_lsst)
from fisher_lsst import *

fsky = 0.35 # 0.4

# same tomo bins for g and s
fish = FisherLsst(cosmoPar, galaxyBiasPar, shearMultBiasPar, photoZPar, nBins=nBins, nL=nL, fsky=fsky, fNk=fNk, magBias=magBias, name=name, nProc=nProc, save=False)

fish.plotDiagCov()

#fish.plotGPhotozRequirements(cosmoPar.ILCDM, name="lcdm")
#fish.plotOutlierPhotozRequirements(cosmoPar.ILCDM, name="lcdm")
#
#
#par, _ = fish.computePosterior(fish.lMaxMask+fish.gsOnlyMask, cosmoPar.ILCDM)
#par.plotContours(IPar=cosmoPar.ILCDM, marg=True, lim=4., color='#E10014', path=None)
#
