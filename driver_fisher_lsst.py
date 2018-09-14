import fisher_lsst
reload(fisher_lsst)
from fisher_lsst import *

##################################################################################
# Forecast parameters

nBins = 10
nL = 20
fsky = 0.4

# cosmological parameters
massiveNu = True  #False
wCDM = True #False
curvature = True #False

# priors to include
PlanckPrior = True

# forecast name
#name = "lcdm"
name = "lcdm_mnu_curv_w0wa"

##################################################################################
# Parameter classes

cosmoPar = CosmoParams(massiveNu=massiveNu, wCDM=wCDM, curvature=curvature, PlanckPrior=PlanckPrior)
#cosmoPar.plotParams()
galaxyBiasPar = GalaxyBiasParams(nBins=nBins)
#galaxyBiasPar.plotParams()
shearMultBiasPar = ShearMultBiasParams(nBins=nBins)
#shearMultBiasPar.plotParams()
photoZPar = PhotoZParams(nBins=nBins)
#photoZPar.plotParams()


#pat = PatPlanckParams()
#pat.printParams()

#u = Universe(cosmoPar.paramsClassy)
#u.plotDistances()

##################################################################################
# Fisher calculation

fisherLsst = FisherLsst(cosmoPar, galaxyBiasPar, shearMultBiasPar, photoZPar, nBins=nBins, nL=nL, fsky=0.4, name=name, magBias=False, save=True)



#fisherLsst.hackCurvature()


#fisherLsst.plotDndz()
#fisherLsst.plotPowerSpectra()
#fisherLsst.plotUncertaintyPowerSpectra()
#fisherLsst.plotCovMat()
#fisherLsst.plotDerivativeDataVector()

#fisherLsst.fullPar.printParams()
#fisherLsst.posteriorPar.printParams()

#fisherLsst.posteriorPar.plotParams(IPar=range(cosmoPar.nPar))
#fisherLsst.posteriorPar.plotParams(IPar=range(cosmoPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar))
#fisherLsst.posteriorPar.plotParams(IPar=range(cosmoPar.nPar+galaxyBiasPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar))
#fisherLsst.posteriorPar.plotParams(IPar=range(cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar+photoZPar.nPar))


#fisherLsst.posteriorPar.plotParamStd(IPar=range(cosmoPar.nPar))


##################################################################################
# Photo-z requirements


#fisherLsst.photoZRequirements()
#fisherLsst.shearBiasRequirements()


















