import fisher_lsst_dndz
reload(fisher_lsst_dndz)
from fisher_lsst_dndz import *

##################################################################################
# Forecast parameters

nBins = 5  #5   #10
nZ = 5
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


dndzPar = DndzParams(nBins=nBins, nZ=nZ, sNgal=1.e-4)
#dndzPar.plotParams()



##################################################################################
# Fisher calculation

fisherLsst = FisherLsst(cosmoPar, galaxyBiasPar, shearMultBiasPar, dndzPar, nBins=nBins, nL=nL, fsky=0.4, magBias=magBias, fullCross=fullCross, name=name, nProc=nProc, save=True)



fisherLsst.plotDndz()
fisherLsst.plotPowerSpectra()
fisherLsst.plotUncertaintyPowerSpectra()
fisherLsst.plotCovMat()
#fisherLsst.plotInvCovMat()

fisherLsst.printSnrPowerSpectra(path=fisherLsst.figurePath+"/snr.txt")





fisherLsst.plotDerivativeDataVectorCosmo()


#fisherLsst.plotSingleDerivative("gg", 0, 0)
#fisherLsst.plotSingleDerivative("ss", 0, 15)
#fisherLsst.plotSingleDerivative("gg", 0, 20)

#cosmoPar.printParams()
#fisherLsst.posteriorPar.printParams()
fisherLsst.fullPar.printParams(path=fisherLsst.figurePath+"/prior_uncertainties.txt")
fisherLsst.posteriorPar.printParams(path=fisherLsst.figurePath+"/posterior_uncertainties.txt")


#cosmoPar.plotContours(path=fisherLsst.figurePath+"/contours_cosmo_prior.pdf")
#fisherLsst.posteriorPar.plotContours(IPar=range(cosmoPar.nPar), path=fisherLsst.figurePath+"/contours_cosmo_posterior.pdf")


fisherLsst.checkConditionNumbers()

#fisherLsst.posteriorPar.plotParams(IPar=range(cosmoPar.nPar))
#fisherLsst.posteriorPar.plotParams(IPar=range(cosmoPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar))
#fisherLsst.posteriorPar.plotParams(IPar=range(cosmoPar.nPar+galaxyBiasPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar))
#fisherLsst.posteriorPar.plotParams(IPar=range(cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar+photoZPar.nPar))

















