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
name = "lcdm_mnu_curv_w0wa_newellsandunits"

# Parallel evaluations
nProc = 3   # not actually used, because CLASS won't be pickled...

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

#cosmoPar.plotContours()

#pat = PatPlanckParams()
#pat.printParams()

#u = Universe(cosmoPar.paramsClassy)
#u.plotDistances()

##################################################################################
# Fisher calculation

fisherLsst = FisherLsst(cosmoPar, galaxyBiasPar, shearMultBiasPar, photoZPar, nBins=nBins, nL=nL, fsky=0.4, magBias=magBias, fullCross=fullCross, name=name, nProc=nProc, save=True)


fisherLsst.plotDndz()
fisherLsst.plotPowerSpectra()
fisherLsst.plotUncertaintyPowerSpectra()
fisherLsst.plotCovMat()
fisherLsst.plotInvCovMat()

fisherLsst.printSnrPowerSpectra(path=fisherLsst.figurePath+"/snr.txt")





fisherLsst.plotDerivativeDataVectorCosmo()


#fisherLsst.plotSingleDerivative("gg", 0, 0)
#fisherLsst.plotSingleDerivative("ss", 0, 15)
#fisherLsst.plotSingleDerivative("gg", 0, 20)

#cosmoPar.printParams()
fisherLsst.fullPar.printParams(path=fisherLsst.figurePath+"/prior_uncertainties.txt")
fisherLsst.posteriorPar.printParams(path=fisherLsst.figurePath+"/posterior_uncertainties.txt")
#fisherLsst.posteriorPar.printParams()

cosmoPar.plotContours(path=fisherLsst.figurePath+"/contours_cosmo_prior.pdf")
fisherLsst.posteriorPar.plotContours(IPar=range(cosmoPar.nPar), path=fisherLsst.figurePath+"/contours_cosmo_posterior.pdf")


fisherLsst.checkConditionNumbers()

fisherLsst.posteriorPar.plotParams(IPar=range(cosmoPar.nPar))
fisherLsst.posteriorPar.plotParams(IPar=range(cosmoPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar))
fisherLsst.posteriorPar.plotParams(IPar=range(cosmoPar.nPar+galaxyBiasPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar))
fisherLsst.posteriorPar.plotParams(IPar=range(cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar+photoZPar.nPar))


##################################################################################
# Photo-z requirements


fisherLsst.photoZRequirements()
#fisherLsst.shearBiasRequirements()


##################################################################################
# show data vector and uncertainties


fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# data vector
ax.semilogy(fisherLsst.dataVector, 'b', label=r'data')
ax.semilogy(-fisherLsst.dataVector, 'b--')
#
# cov matrix
ax.semilogy(np.sqrt(np.diag(fisherLsst.covMat)), 'r', label=r'uncertainty')
#
# relative uncertainty
ax.semilogy(np.sqrt(np.diag(fisherLsst.covMat)) / fisherLsst.dataVector, 'g', label=r'relat. uncertainty')
#
# observables
ax.axvline(fisherLsst.nL*fisherLsst.nGG, c='k', lw=1.5)
ax.axvline(fisherLsst.nL*(fisherLsst.nGG+fisherLsst.nGS), c='k', lw=1.5)
#
ax.legend()

plt.show()







import scipy
from scipy.sparse import csc_matrix



'''
##################################################################################
# Test inversions of cov and Fisher matrices


# Check eigenvalues of cov matrix
#
eigenValPar, eigenVecPar = np.linalg.eigh(fisherLsst.covMat)
plt.semilogy(eigenValPar, 'b')
#
## relative cov matrix
#invDataVector = 1./fisherLsst.dataVector
#invDataVector[np.where(np.isfinite(invDataVector)==False)] = 0.
##
#matInvDataVector = np.diag(invDataVector)
##
#relatCovMat = np.dot(matInvDataVector, np.dot(fisherLsst.covMat, matInvDataVector))
#invRelatCovMat = np.linalg.inv(relatCovMat)
##
#eigenValPar, eigenVecPar = np.linalg.eigh(relatCovMat)
#plt.semilogy(eigenValPar, 'g')


#
# Try with SVD
U, s, Vh = scipy.linalg.svd(fisherLsst.covMat)
plt.semilogy(s[::-1], 'r.')
V = np.conj(Vh.transpose())
Uh = np.conj(U.transpose())
#sInv = np.linalg.inv(np.diag(s))



plt.show()




# try with numpy inversion
inv1 = np.linalg.inv(fisherLsst.covMat)
eigenValPar, eigenVecPar = np.linalg.eigh(inv1)
plt.semilogy(eigenValPar, 'b')
#
# try with the SVD, truncated
sInv = 1./s
sInvMax = np.max(sInv)
# remove the super poorly constrained modes, that lead to numerical instabilities
sInv[sInv<=sInvMax*1.e-5] = 0.
sInv = np.diag(sInv)

inv3 = np.dot(V, np.dot(sInv, Uh))
eigenValPar3, eigenVecPar3 = np.linalg.eigh(inv3)
plt.semilogy(eigenValPar3, 'm.')


plt.show()



# matrix condition number
print "cov matrix"
print "inverse condition number:", 1./np.linalg.cond(fisherLsst.covMat)
print "number numerical precision:", np.finfo(fisherLsst.covMat.dtype).eps

#
#print "relative cov matrix"
#print "inverse condition number:", 1./np.linalg.cond(relatCovMat)
#print "number numerical precision:", np.finfo(relatCovMat.dtype).eps









# Check eigenvalues of inverse cov matrix
#
# try with numpy inversion
inv1 = np.linalg.inv(fisherLsst.covMat)
eigenValPar, eigenVecPar = np.linalg.eigh(inv1)
plt.semilogy(eigenValPar, 'b')
#
# Try with scipy's sparse class
sa = csc_matrix(fisherLsst.covMat)
inv2 = scipy.sparse.linalg.inv(sa).todense()
eigenValPar2, eigenVecPar2 = np.linalg.eigh(inv2)
plt.semilogy(eigenValPar2, 'r.')
#
# Try with SVD
U, s, Vh = scipy.linalg.svd(fisherLsst.covMat)
V = np.conj(Vh.transpose())
Uh = np.conj(U.transpose())
sInv = np.linalg.inv(np.diag(s))
inv3 = np.dot(V, np.dot(sInv, Uh))
eigenValPar3, eigenVecPar3 = np.linalg.eigh(inv3)
plt.semilogy(eigenValPar3, 'm.')
#
# show the absolute difference, to compare
plt.semilogy(np.abs(eigenValPar2 - eigenValPar), 'g')

plt.show()


#print np.std(fisherLsst.covMat), np.max(fisherLsst.covMat)
#print np.std(inv1), np.max(inv1)
#
#inv4 = np.linalg.inv(fisherLsst.covMat + 1.e-27*np.diag(np.ones(fisherLsst.nData)))
#print np.std(inv4-inv1)/np.std(inv1), np.std(inv4-inv1)/np.std(inv4)

'''



##################################################################################


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


# matrix condition number
print "Fisher matrix"
print "inverse condition number:", 1./np.linalg.cond(fisherLsst.fisherPosterior)
print "number numerical precision:", np.finfo(fisherLsst.fisherPosterior.dtype).eps




# look at specific eigenvectors

fig=plt.figure(0)
ax=fig.add_subplot(111)
#

#ax.plot(eigenVecPar[:,-1], 'b.')  # best constrained mode

ax.plot(eigenVecPar[:,0], 'r.')  # worst constrained mode

#
ax.set_xticks(range(fisherLsst.fullPar.nPar))
ax.set_xticklabels(fisherLsst.fullPar.namesLatex, fontsize=24)
[l.set_rotation(45) for l in ax.get_xticklabels()]


plt.show()





