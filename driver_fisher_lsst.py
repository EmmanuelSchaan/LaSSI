import fisher_lsst
reload(fisher_lsst)
from fisher_lsst import *

##################################################################################
# Forecast parameters

nBins = 5  #5   #10
nL = 20
fsky = 0.4

# cosmological parameters
massiveNu = True  #False
wCDM = True #False
curvature = True #False

# priors to include
PlanckPrior = True

# include null crosses
fullCross = False #False # True

# include a known magnification bias
magBias = False

# forecast name
#name = "lcdm"
#name = "lcdm_mnu_curv_w0wa"
name = "lcdm_mnu_curv_w0wa_smallerphotozstep"

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

fisherLsst = FisherLsst(cosmoPar, galaxyBiasPar, shearMultBiasPar, photoZPar, nBins=nBins, nL=nL, fsky=0.4, magBias=magBias, fullCross=fullCross, name=name, save=False)


#fisherLsst.plotDndz()
#fisherLsst.plotPowerSpectra()
#fisherLsst.plotUncertaintyPowerSpectra()
#fisherLsst.plotCovMat()
#fisherLsst.plotDerivativeDataVectorCosmo()


#fisherLsst.plotSingleDerivative("gg", 0, 0)
#fisherLsst.plotSingleDerivative("ss", 0, 15)
#fisherLsst.plotSingleDerivative("gg", 0, 20)

cosmoPar.printParams()
fisherLsst.fullPar.printParams()
fisherLsst.posteriorPar.printParams()

#fisherLsst.posteriorPar.plotParams(IPar=range(cosmoPar.nPar))
#fisherLsst.posteriorPar.plotParams(IPar=range(cosmoPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar))
#fisherLsst.posteriorPar.plotParams(IPar=range(cosmoPar.nPar+galaxyBiasPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar))
#fisherLsst.posteriorPar.plotParams(IPar=range(cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar+photoZPar.nPar))




# Check eigenvalues of Fisher matrix
eigenValPar, eigenVecPar = np.linalg.eigh(fisherLsst.covMat)
plt.semilogy(eigenValPar, '.')
plt.show()

eigenValPar, eigenVecPar = np.linalg.eigh(fisherLsst.fisherData)
plt.semilogy(1./np.sqrt(eigenValPar), '.')
#plt.semilogy(eigenValPar, '.')
plt.show()




##################################################################################
# Photo-z requirements


#fisherLsst.photoZRequirements()
#fisherLsst.shearBiasRequirements()


















