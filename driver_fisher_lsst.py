import fisher_lsst
reload(fisher_lsst)
from fisher_lsst import *

##################################################################################
# Forecast parameters

nBins = 10  #5   #10
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

fisherLsst = FisherLsst(cosmoPar, galaxyBiasPar, shearMultBiasPar, photoZPar, nBins=nBins, nL=nL, fsky=0.4, magBias=magBias, fullCross=fullCross, name=name, save=True)


fisherLsst.plotDndz()
fisherLsst.plotPowerSpectra()
fisherLsst.plotUncertaintyPowerSpectra()
fisherLsst.plotCovMat()
fisherLsst.plotDerivativeDataVectorCosmo()


#fisherLsst.plotSingleDerivative("gg", 0, 0)
#fisherLsst.plotSingleDerivative("ss", 0, 15)
#fisherLsst.plotSingleDerivative("gg", 0, 20)

#cosmoPar.printParams()
#fisherLsst.fullPar.printParams()
#fisherLsst.posteriorPar.printParams()

#fisherLsst.posteriorPar.plotParams(IPar=range(cosmoPar.nPar))
#fisherLsst.posteriorPar.plotParams(IPar=range(cosmoPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar))
#fisherLsst.posteriorPar.plotParams(IPar=range(cosmoPar.nPar+galaxyBiasPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar))
#fisherLsst.posteriorPar.plotParams(IPar=range(cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar+photoZPar.nPar))


##################################################################################
# Photo-z requirements


fisherLsst.photoZRequirements()
#fisherLsst.shearBiasRequirements()


##################################################################################
# Test inversions of cov and Fisher matrices

'''
import scipy
from scipy.sparse import csc_matrix

# Check eigenvalues of cov matrix
#
# try with numpy inversion
eigenValPar, eigenVecPar = np.linalg.eigh(np.linalg.inv(fisherLsst.covMat))
plt.semilogy(eigenValPar, 'b')
#
# Try with scipy's sparse class
sa = csc_matrix(fisherLsst.covMat)
eigenValPar2, eigenVecPar2 = np.linalg.eigh(scipy.sparse.linalg.inv(sa).todense())
plt.semilogy(eigenValPar2, 'r.')
#
# show the absolute difference, to compare
plt.semilogy(np.abs(eigenValPar2 - eigenValPar), 'g')

plt.show()



# Check eigenvalues of Fisher matrix
#
# try with numpy inversion
eigenValPar, eigenVecPar = np.linalg.eigh(np.linalg.inv(fisherLsst.fisherPosterior))
plt.semilogy(eigenValPar, 'b')
#
# Try with scipy's sparse class
sa = csc_matrix(fisherLsst.fisherPosterior)
eigenValPar2, eigenVecPar2 = np.linalg.eigh(scipy.sparse.linalg.inv(sa).todense())
plt.semilogy(eigenValPar2, 'r.')
#
# show the absolute difference, to compare
plt.semilogy(np.abs(eigenValPar2 - eigenValPar), 'g')

plt.show()
'''











