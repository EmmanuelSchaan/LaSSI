import fisher_lsst_dndz
reload(fisher_lsst_dndz)
from fisher_lsst_dndz import *

##################################################################################
# Forecast parameters

nBins = 10  #5   #10
nZ = 10
nL = 20
fsky = 0.4

# cosmological parameters
massiveNu = True  #False
wCDM = True #False
curvature = True #False

# priors to include
PlanckPrior = True

# include null crosses
fullCross = True #False # True

# include a known magnification bias
magBias = False

# forecast name
name = "lsst_dndz"

# Parallel evaluations
nProc = 3   # not actually used, because CLASS won't be pickled...

##################################################################################
# Parameter classes

cosmoPar = CosmoParams(massiveNu=massiveNu, wCDM=wCDM, curvature=curvature, PlanckPrior=PlanckPrior)
#cosmoPar.plotParams()
galaxyBiasPar = GalaxyBiasParams(nBins=nBins)
#galaxyBiasPar.plotParams()
shearMultBiasPar = ShearMultBiasParams(nBins=nBins)
#shearMultBiasPar = ShearMultBiasParams(nBins=nBins, mStd=1.e-5)   # perfect photo-z priors
#shearMultBiasPar.plotParams()
#photoZPar = PhotoZParams(nBins=nBins)
#photoZPar.plotParams()


dndzParams = DndzParams(nBins=nBins, nZ=nZ, sNgal=1.e-4)
dndzParams.plotParams()
