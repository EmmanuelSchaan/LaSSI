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
nBins =  10  
nL = 50  
fsky = 0.35


# cosmological parameters to include
massiveNu = True  
wCDM = True 
curvature = True 

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
fish = FisherLsst(cosmoPar, galaxyBiasPar, shearMultBiasPar, photoZPar, nBins=nBins, nL=nL, fsky=fsky, fNk=fNk, magBias=magBias, name=name, nProc=nProc, save=True)


# different tomo bins for g and s
#fish = FisherLsst(cosmoPar, galaxyBiasPar, shearMultBiasPar, photoZPar, photoZSPar=photoZPar, nBins=nBins, nL=nL, fsky=fsky, fNk=fNk, magBias=magBias, name=name, nProc=nProc, save=True)


##################################################################################
# basic plots and checks
'''
#fish.SampleDndz(photoZPar, nSamples=2)#, path=fish.figurePath+"/sampling_dndz_prior.pdf")

# Show observables and uncertainties
fish.plotEllBins(show=False)
fish.plotDndz(show=False)
fish.plotPowerSpectra(show=False)
fish.plotUncertaintyPowerSpectra(show=False)
#fish.plotCovMat(show=False)
#fish.plotInvCovMat(show=False)
fish.printSnrPowerSpectra(path=fish.figurePath+"/snr.txt")
fish.plotDerivativeDataVectorCosmo(show=False)
#fish.plotSingleDerivative("gg", 0, 0)
#fish.plotSingleDerivative("ss", 0, 15)
#fish.plotSingleDerivative("gg", 0, 20)


# Condition numbers,
# for various data combinations
fish.checkConditionNumbers(mask=fish.lMaxMask) # default
fish.checkConditionNumbers(mask=fish.lMaxMask+fish.noNullMask) # no null 2-pt functions
fish.checkConditionNumbers(mask=fish.lMaxMask+fish.gsOnlyMask) # g,s-only
fish.checkConditionNumbers(mask=fish.lMaxMask+fish.gsOnlyMask+fish.noNullMask) # g,s-only, no null
fish.checkConditionNumbers(mask=fish.lMaxMask+fish.gOnlyMask)   # g-only
fish.checkConditionNumbers(mask=fish.lMaxMask+fish.sOnlyMask)   # s-only
'''

##################################################################################
# Parameter uncertainties, and dn/dz visualizations
# for various data combinations
'''
# Prior
fish.fullPar.printParams(path=fish.figurePath+"/prior_uncertainties.txt")
#cosmoPar.plotContours(path=fish.figurePath+"/contours_cosmo_prior.pdf")
# visualize photo-z uncertainties
#fish.SampleDndz(photoZPar, nSamples=10, path=fish.figurePath+"/sampling_dndz_prior.pdf")


# k,g,s
par = fish.fullPar.copy()
fisherData, fisherPosterior = fish.generateFisher(mask=fish.lMaxMask)  # default
par.fisher = fisherPosterior
par.printParams(path=fish.figurePath+"/posterior_uncertainties_kgs.txt")
#par.plotContours(IPar=range(cosmoPar.nPar), path=fish.figurePath+"/contours_cosmo_posterior.pdf")
# visualize photo-z uncertainties
#I = range(-photoZPar.nPar, 0)
#newPhotoZPar = par.extractParams(I, marg=True)
#fish.SampleDndz(newPhotoZPar, nSamples=10, path=fish.figurePath+"/sampling_dndz_posterior_kgs.pdf")


# k,g,s no null 2-pt functions
par = fish.fullPar.copy()
fisherData, fisherPosterior = fish.generateFisher(mask=fish.lMaxMask+fish.noNullMask)  # no null 2-pt functions
par.fisher = fisherPosterior
par.printParams(path=fish.figurePath+"/posterior_uncertainties_kgs_nonull.txt")
# visualize photo-z uncertainties
#I = range(-photoZPar.nPar, 0)
#newPhotoZPar = par.extractParams(I, marg=True)
#fish.SampleDndz(newPhotoZPar, nSamples=10, path=fish.figurePath+"/sampling_dndz_posterior_kgs_nonull.pdf")


# Fiducial: g,s
par = fish.fullPar.copy()
fisherData, fisherPosterior = fish.generateFisher(mask=fish.lMaxMask+fish.gsOnlyMask)  # default
par.fisher = fisherPosterior
par.printParams(path=fish.figurePath+"/posterior_uncertainties_gs.txt")
#par.plotContours(IPar=range(cosmoPar.nPar), path=fish.figurePath+"/contours_cosmo_posterior.pdf")
# visualize photo-z uncertainties
#I = range(-photoZPar.nPar, 0)
#newPhotoZPar = par.extractParams(I, marg=True)
#fish.SampleDndz(newPhotoZPar, nSamples=10, path=fish.figurePath+"/sampling_dndz_posterior_gs.pdf")

# g,s, no null 2pt
par = fish.fullPar.copy()
fisherData, fisherPosterior = fish.generateFisher(mask=fish.lMaxMask+fish.gsOnlyMask+fish.noNullMask)  # default
par.fisher = fisherPosterior
par.printParams(path=fish.figurePath+"/posterior_uncertainties_gs_nonull.txt")
#par.plotContours(IPar=range(cosmoPar.nPar), path=fish.figurePath+"/contours_cosmo_posterior.pdf")
# visualize photo-z uncertainties
#I = range(-photoZPar.nPar, 0)
#newPhotoZPar = par.extractParams(I, marg=True)
#fish.SampleDndz(newPhotoZPar, nSamples=10, path=fish.figurePath+"/sampling_dndz_posterior_gs_nonull.pdf")

# g-only
par = fish.fullPar.copy()
fisherData, fisherPosterior = fish.generateFisher(mask=fish.lMaxMask+fish.gOnlyMask)  # g-only
par.fisher = fisherPosterior
par.printParams(path=fish.figurePath+"/posterior_uncertainties_g.txt")
# visualize photo-z uncertainties
#I = range(-photoZPar.nPar, 0)
#newPhotoZPar = par.extractParams(I, marg=True)
#fish.SampleDndz(newPhotoZPar, nSamples=10, path=fish.figurePath+"/sampling_dndz_posterior_g.pdf")

# s-only
par = fish.fullPar.copy()
fisherData, fisherPosterior = fish.generateFisher(mask=fish.lMaxMask+fish.sOnlyMask)  # s-only
par.fisher = fisherPosterior
## !!! no need to fix the galaxy bias, since I added an uninformative prior on it,
## just to be able to invert the Fisher matrix
## With s-only, bg cannot be constrained: we fix it
#I = range(cosmoPar.nPar) + range(cosmoPar.nPar+galaxyBiasPar.nPar, fish.fullPar.nPar)
#par = par.extractParams(I, marg=False)
par.printParams(path=fish.figurePath+"/posterior_uncertainties_s.txt")
# visualize photo-z uncertainties
#I = range(-photoZPar.nPar, 0)
#newPhotoZPar = par.extractParams(I, marg=True)
#fish.SampleDndz(newPhotoZPar, nSamples=10, path=fish.figurePath+"/sampling_dndz_posterior_s.pdf")


#fish.posteriorPar.plotParams(IPar=range(cosmoPar.nPar))
#fish.posteriorPar.plotParams(IPar=range(cosmoPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar))
#fish.posteriorPar.plotParams(IPar=range(cosmoPar.nPar+galaxyBiasPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar))
#fish.posteriorPar.plotParams(IPar=range(cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar, cosmoPar.nPar+galaxyBiasPar.nPar+shearMultBiasPar.nPar+photoZPar.nPar))
'''

##################################################################################
##################################################################################
# Dependence on photo-z priors
# old version, to be removed

'''
# Gaussian Photo-z requirements
fish.photoZRequirements(mask=fish.lMaxMask, name="")  # k,g,s
fish.photoZRequirements(mask=fish.lMaxMask+fish.noNullMask, name="nonull")  # k,g,s, no null 2pt
fish.photoZRequirements(mask=fish.lMaxMask+fish.gsOnlyMask, name="gs")  # g,s
fish.photoZRequirements(mask=fish.lMaxMask+fish.gsOnlyMask+fish.noNullMask, name="gs_nonull")  # g,s, no null 2pt
#fish.photoZRequirements(mask=fish.lMaxMask+fish.gOnlyMask, name="g")  # g
#fish.photoZRequirements(mask=fish.lMaxMask+fish.sOnlyMask, name="s")  # s

# Outlier requirements
if fish.photoZPar.outliers==0.1:
   #
   # G photo-z prior: fixed at level of LSST req
   fish.photoZOutliersRequirements(mask=fish.lMaxMask, name="", Gphotoz='req')  # k,g,s
   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.noNullMask, name="nonull", Gphotoz='req')  # k,g,s, no null 2pt
   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.gsOnlyMask, name="gs", Gphotoz='req')  # g,s
   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.gsOnlyMask+fish.noNullMask, name="gs_nonull", Gphotoz='req')  # g,s, no null 2pt
#   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.gOnlyMask, name="g", Gphotoz='req')  # g
#   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.sOnlyMask, name="s", Gphotoz='req')  # s
   #
   # G photo-z prior: fixed with perfect prior
#!!! problem with the line below: matrix is singular
   fish.photoZOutliersRequirements(mask=fish.lMaxMask, name="", Gphotoz='perfect')  # k,g,s
#!!! line below: matrix is singular
   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.noNullMask, name="nonull", Gphotoz='perfect')  # k,g,s, no null 2pt
   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.gsOnlyMask, name="gs", Gphotoz='perfect')  # g,s
   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.gsOnlyMask+fish.noNullMask, name="gs_nonull", Gphotoz='perfect')  # g,s, no null 2pt
#   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.gOnlyMask, name="g", Gphotoz='perfect')  # g
#   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.sOnlyMask, name="s", Gphotoz='perfect')  # s
   #
   # G photo-z prior: fixed with no prior
   fish.photoZOutliersRequirements(mask=fish.lMaxMask, name="", Gphotoz='noprior')  # k,g,s
   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.noNullMask, name="nonull", Gphotoz='noprior')  # k,g,s, no null 2pt
   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.gsOnlyMask, name="gs", Gphotoz='noprior')  # g,s
   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.gsOnlyMask+fish.noNullMask, name="gs_nonull", Gphotoz='noprior')  # g,s, no null 2pt
#   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.gOnlyMask, name="g", Gphotoz='noprior')  # g
#   fish.photoZOutliersRequirements(mask=fish.lMaxMask+fish.sOnlyMask, name="s", Gphotoz='noprior')  # s



#fish.shearBiasRequirements()
'''


##################################################################################
##################################################################################
# Dependence on photo-z priors
# new version
'''
fish.plotGPhotozRequirements(cosmoPar.ILCDMW0Wa, name="lcdmw0wa")
fish.plotOutlierPhotozRequirements(cosmoPar.ILCDMW0Wa, name="lcdmw0wa")

fish.plotGPhotozRequirements(cosmoPar.ILCDMW0Wa, name="lcdmmnu")
fish.plotOutlierPhotozRequirements(cosmoPar.ILCDMW0Wa, name="lcdmmnu")
'''

##################################################################################
##################################################################################
# Contour plots

'''
# LCDMW0Wa
# GS
par, _ = fish.computePosterior(fish.lMaxMask+fish.gsOnlyMask, cosmoPar.ILCDMW0Wa)
par.plotContours(IPar=cosmoPar.ILCDMW0Wa, marg=True, lim=4., color='#E10014', path=fish.figurePath+"/contours_gs_lcdmw0wa.pdf")
# GS no null
par, _ = fish.computePosterior(fish.lMaxMask+fish.gsOnlyMask+fish.noNullMask, cosmoPar.ILCDMW0Wa)
par.plotContours(IPar=cosmoPar.ILCDMW0Wa, marg=True, lim=4., color='#E10014', path=fish.figurePath+"/contours_gsnonull_lcdmw0wa.pdf")
# GKS
par, _ = fish.computePosterior(fish.lMaxMask, cosmoPar.ILCDMW0Wa)
par.plotContours(IPar=cosmoPar.ILCDMW0Wa, marg=True, lim=4., color='#E10014', path=fish.figurePath+"/contours_gks_lcdmw0wa.pdf")

# LCDMMnu
# GS
par, _ = fish.computePosterior(fish.lMaxMask+fish.gsOnlyMask, cosmoPar.ILCDMMnu)
par.plotContours(IPar=cosmoPar.ILCDMMnu, marg=True, lim=4., color='#E10014', path=fish.figurePath+"/contours_gs_lcdmmnu.pdf")
# GS no null
par, _ = fish.computePosterior(fish.lMaxMask+fish.gsOnlyMask+fish.noNullMask, cosmoPar.ILCDMMnu)
par.plotContours(IPar=cosmoPar.ILCDMMnu, marg=True, lim=4., color='#E10014', path=fish.figurePath+"/contours_gsnonull_lcdmmnu.pdf")
# GKS
par, _ = fish.computePosterior(fish.lMaxMask, cosmoPar.ILCDMMnu)
par.plotContours(IPar=cosmoPar.ILCDMMnu, marg=True, lim=4., color='#E10014', path=fish.figurePath+"/contours_gks_lcdmmnu.pdf")
'''

