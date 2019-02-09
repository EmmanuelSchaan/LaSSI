import fisher_lsst
reload(fisher_lsst)
from fisher_lsst import *

##################################################################################
# Forecast parameters

nBins = 10  #5   #10
nL = 20
fsky = 0.4

# cosmological parameters to include
massiveNu = True  #False
wCDM = True #False
curvature = True #False

# priors to include
PlanckPrior = True

# include a known magnification bias
magBias = False

# forecast name
#name = "lcdm"
#name = "lcdm_mnu_curv_w0wa"
#name = "lcdm_mnu_curv_w0wa_newellsandunits"
#name = "lcdm_mnu_curv_w0wa_newellsandunits_perfectm"
#name = "gphotoz_lmaxmask"
#name = "gphotoz"
name = None

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
# Fisher calculation

fish = FisherLsst(cosmoPar, galaxyBiasPar, shearMultBiasPar, photoZPar, nBins=nBins, nL=nL, fsky=0.4, magBias=magBias, name=name, nProc=nProc, save=True)



#fish.SampleDndz(photoZPar, nSamples=2)#, path=fish.figurePath+"/sampling_dndz_prior.pdf")




# Show observables and uncertainties
fish.plotDndz()
fish.plotPowerSpectra()
fish.plotUncertaintyPowerSpectra()
fish.plotCovMat()
#fish.plotInvCovMat()
fish.printSnrPowerSpectra(path=fish.figurePath+"/snr.txt")
fish.plotDerivativeDataVectorCosmo()
#fish.plotSingleDerivative("gg", 0, 0)
#fish.plotSingleDerivative("ss", 0, 15)
#fish.plotSingleDerivative("gg", 0, 20)


# Condition numbers,
# for various data combinations
fish.checkConditionNumbers(mask=fish.lMaxMask) # default
fish.checkConditionNumbers(mask=fish.lMaxMask+fish.noNullMask) # no null 2-pt functions
fish.checkConditionNumbers(mask=fish.lMaxMask+fish.gOnlyMask)   # g-only
fish.checkConditionNumbers(mask=fish.lMaxMask+fish.sOnlyMask)   # s-only



# Parameter uncertainties, and dn/dz visualizations
# for various data combinations
#
# Prior
fish.fullPar.printParams(path=fish.figurePath+"/prior_uncertainties.txt")
#cosmoPar.plotContours(path=fish.figurePath+"/contours_cosmo_prior.pdf")
# visualize photo-z uncertainties
fish.SampleDndz(photoZPar, nSamples=10, path=fish.figurePath+"/sampling_dndz_prior.pdf")
#
# Fiducial
par = fish.fullPar.copy()
fisherData, fisherPosterior = fish.generateFisher(mask=fish.lMaxMask)  # default
par.fisher = fisherPosterior
par.printParams(path=fish.figurePath+"/posterior_uncertainties.txt")
#par.plotContours(IPar=range(cosmoPar.nPar), path=fish.figurePath+"/contours_cosmo_posterior.pdf")
# visualize photo-z uncertainties
I = range(-photoZPar.nPar, 0)
newPhotoZPar = par.extractParams(I, marg=True)
fish.SampleDndz(newPhotoZPar, nSamples=10, path=fish.figurePath+"/sampling_dndz_posterior.pdf")
#
# No null 2-pt functions
par = fish.fullPar.copy()
fisherData, fisherPosterior = fish.generateFisher(mask=fish.lMaxMask+fish.noNullMask)  # no null 2-pt functions
par.fisher = fisherPosterior
par.printParams(path=fish.figurePath+"/posterior_uncertainties_nonull.txt")
# visualize photo-z uncertainties
I = range(-photoZPar.nPar, 0)
newPhotoZPar = par.extractParams(I, marg=True)
fish.SampleDndz(newPhotoZPar, nSamples=10, path=fish.figurePath+"/sampling_dndz_posterior_nonull.pdf")
#
# g-only
par = fish.fullPar.copy()
fisherData, fisherPosterior = fish.generateFisher(mask=fish.lMaxMask+fish.gOnlyMask)  # g-only
par.fisher = fisherPosterior
par.printParams(path=fish.figurePath+"/posterior_uncertainties_gonly.txt")
# visualize photo-z uncertainties
I = range(-photoZPar.nPar, 0)
newPhotoZPar = par.extractParams(I, marg=True)
fish.SampleDndz(newPhotoZPar, nSamples=10, path=fish.figurePath+"/sampling_dndz_posterior_gonly.pdf")
#
# s-only
par = fish.fullPar.copy()
fisherData, fisherPosterior = fish.generateFisher(mask=fish.lMaxMask+fish.sOnlyMask)  # s-only
par.fisher = fisherPosterior
## !!! no need to fix the galaxy bias, since I added an uninformative prior on it,
## just to be able to invert the Fisher matrix
## With s-only, bg cannot be constrained: we fix it
#I = range(cosmoPar.nPar) + range(cosmoPar.nPar+galaxyBiasPar.nPar, fish.fullPar.nPar)
#par = par.extractParams(I, marg=False)
par.printParams(path=fish.figurePath+"/posterior_uncertainties_sonly.txt")
# visualize photo-z uncertainties
I = range(-photoZPar.nPar, 0)
newPhotoZPar = par.extractParams(I, marg=True)
fish.SampleDndz(newPhotoZPar, nSamples=10, path=fish.figurePath+"/sampling_dndz_posterior_sonly.pdf")


#fish.posteriorPar.plotParams(IPar=range(cosmoPar.nPar))
#fish.posteriorPar.plotParams(IPar=range(cosmoPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar))
#fish.posteriorPar.plotParams(IPar=range(cosmoPar.nPar+galaxyBiasPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar))
#fish.posteriorPar.plotParams(IPar=range(cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar+photoZPar.nPar))


if fish.photoZPar.outliers==0.:
   # Gaussian Photo-z requirements
   fish.photoZRequirements(mask=fish.lMaxMask, name="")  # default
   fish.photoZRequirements(mask=fish.lMaxMask+fish.noNullMask, name="nonull")  # no null 2-pt functions
   fish.photoZRequirements(mask=fish.lMaxMask+fish.gOnlyMask, name="gonly")  # g-only
   #fish.photoZRequirements(mask=fish.lMaxMask+fish.sOnlyMask, name="sonly")  # s-only

else:
   # Outlier requirements
   #
   # G photo-z prior: fixed at level of LSST req
   fish.photoZOutliersRequirements(mask=fish.lMaxMask, name="", Gphotoz='req')  # default
   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.noNullMask, name="nonull", Gphotoz='req')  # no null 2-pt functions
   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.gOnlyMask, name="gonly", Gphotoz='req')  # g-only
   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.sOnlyMask, name="sonly", Gphotoz='req')  # s-only
   #
   # G photo-z prior: fixed with perfect prior
   fish.photoZOutliersRequirements(mask=fish.lMaxMask, name="", Gphotoz='perfect')  # default
   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.noNullMask, name="nonull", Gphotoz='perfect')  # no null 2-pt functions
   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.gOnlyMask, name="gonly", Gphotoz='perfect')  # g-only
   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.sOnlyMask, name="sonly", Gphotoz='perfect')  # s-only
   #
   # G photo-z prior: fixed with no prior
   fish.photoZOutliersRequirements(mask=fish.lMaxMask, name="", Gphotoz='noprior')  # default
   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.noNullMask, name="nonull", Gphotoz='noprior')  # no null 2-pt functions
   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.gOnlyMask, name="gonly", Gphotoz='noprior')  # g-only
   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.sOnlyMask, name="sonly", Gphotoz='noprior')  # s-only



#fish.shearBiasRequirements()




##################################################################################
# Debug: show data vector and uncertainties

'''
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# data vector
ax.semilogy(fish.dataVector, 'b', label=r'data')
ax.semilogy(-fish.dataVector, 'b--')
#
# cov matrix
ax.semilogy(np.sqrt(np.diag(fish.covMat)), 'r', label=r'uncertainty')
#
# relative uncertainty
ax.semilogy(np.sqrt(np.diag(fish.covMat)) / fish.dataVector, 'g', label=r'relat. uncertainty')
#
# observables
ax.axvline(fish.nL*fish.nGG, c='k', lw=1.5)
ax.axvline(fish.nL*(fish.nGG+fish.nGS), c='k', lw=1.5)
#
ax.legend()

plt.show()

'''





import scipy
from scipy.sparse import csc_matrix



'''
##################################################################################
# Test inversions of cov and Fisher matrices


# Check eigenvalues of cov matrix
#
eigenValPar, eigenVecPar = np.linalg.eigh(fish.covMat)
plt.semilogy(eigenValPar, 'b')
#
## relative cov matrix
#invDataVector = 1./fish.dataVector
#invDataVector[np.where(np.isfinite(invDataVector)==False)] = 0.
##
#matInvDataVector = np.diag(invDataVector)
##
#relatCovMat = np.dot(matInvDataVector, np.dot(fish.covMat, matInvDataVector))
#invRelatCovMat = np.linalg.inv(relatCovMat)
##
#eigenValPar, eigenVecPar = np.linalg.eigh(relatCovMat)
#plt.semilogy(eigenValPar, 'g')


#
# Try with SVD
U, s, Vh = scipy.linalg.svd(fish.covMat)
plt.semilogy(s[::-1], 'r.')
V = np.conj(Vh.transpose())
Uh = np.conj(U.transpose())
#sInv = np.linalg.inv(np.diag(s))



plt.show()




# try with numpy inversion
inv1 = np.linalg.inv(fish.covMat)
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
print "inverse condition number:", 1./np.linalg.cond(fish.covMat)
print "number numerical precision:", np.finfo(fish.covMat.dtype).eps

#
#print "relative cov matrix"
#print "inverse condition number:", 1./np.linalg.cond(relatCovMat)
#print "number numerical precision:", np.finfo(relatCovMat.dtype).eps









# Check eigenvalues of inverse cov matrix
#
# try with numpy inversion
inv1 = np.linalg.inv(fish.covMat)
eigenValPar, eigenVecPar = np.linalg.eigh(inv1)
plt.semilogy(eigenValPar, 'b')
#
# Try with scipy's sparse class
sa = csc_matrix(fish.covMat)
inv2 = scipy.sparse.linalg.inv(sa).todense()
eigenValPar2, eigenVecPar2 = np.linalg.eigh(inv2)
plt.semilogy(eigenValPar2, 'r.')
#
# Try with SVD
U, s, Vh = scipy.linalg.svd(fish.covMat)
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


#print np.std(fish.covMat), np.max(fish.covMat)
#print np.std(inv1), np.max(inv1)
#
#inv4 = np.linalg.inv(fish.covMat + 1.e-27*np.diag(np.ones(fish.nData)))
#print np.std(inv4-inv1)/np.std(inv1), np.std(inv4-inv1)/np.std(inv4)

'''



##################################################################################

'''
# Check eigenvalues of Fisher matrix
#
# try with numpy inversion
eigenValPar, eigenVecPar = np.linalg.eigh(np.linalg.inv(fish.fisherPosterior))
plt.semilogy(eigenValPar, 'b')
#
# Try with scipy's sparse class
sa = csc_matrix(fish.fisherPosterior)
eigenValPar2, eigenVecPar2 = np.linalg.eigh(scipy.sparse.linalg.inv(sa).todense())
plt.semilogy(eigenValPar2, 'r.')
#
# show the absolute difference, to compare
plt.semilogy(np.abs(eigenValPar2 - eigenValPar), 'g')

plt.show()


# matrix condition number
print "Fisher matrix"
print "inverse condition number:", 1./np.linalg.cond(fish.fisherPosterior)
print "number numerical precision:", np.finfo(fish.fisherPosterior.dtype).eps




# look at specific eigenvectors

fig=plt.figure(0)
ax=fig.add_subplot(111)
#

#ax.plot(eigenVecPar[:,-1], 'b.')  # best constrained mode

ax.plot(eigenVecPar[:,0], 'r.')  # worst constrained mode

#
ax.set_xticks(range(fish.fullPar.nPar))
ax.set_xticklabels(fish.fullPar.namesLatex, fontsize=24)
[l.set_rotation(45) for l in ax.get_xticklabels()]


plt.show()

'''



