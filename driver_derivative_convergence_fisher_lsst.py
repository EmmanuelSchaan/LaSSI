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


# Parallel evaluations
nProc = 4   # not actually used

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

##################################################################################

save = False

#derivStepSizes = np.array([0.01, 0.05, 0.1, 0.5, 1., 1.5, 2.])
derivStepSizes = np.array([0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 1., 1.5, 2.])#, 1.75]
#derivStepSizes = np.array([0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 1.])#, 1.5, 2.])#, 1.75]
#derivStepSizes = np.array([2.])
nStep = len(derivStepSizes)

iStepFid = 6#-1 #4


nData = 11550  # read from the output of the fisher class
nPar = 140  # read from the output of the fisher class
# data vector derivatives
deriv  = np.zeros((nStep, nPar, nData))
# Posterior uncertainties
par = np.empty(nStep, dtype=object)
s = np.zeros((nStep, nPar))


for iStep in range(nStep):
   derivStepSize = derivStepSizes[iStep]
   # forecast name
   name = 'derivstepsize'+floatExpForm(derivStepSize, round=2)

   # Parameter classes
   cosmoPar = CosmoParams(massiveNu=massiveNu, wCDM=wCDM, curvature=curvature, PlanckPrior=PlanckPrior, derivStepSize=derivStepSize)
   #print 'cosmo par'
   #print cosmoPar.high - cosmoPar.low
   galaxyBiasPar = GalaxyBiasParams(nBins=nBins, derivStepSize=derivStepSize)
   #print 'gal bias'
   #print galaxyBiasPar.high - galaxyBiasPar.low
   shearMultBiasPar = ShearMultBiasParams(nBins=nBins, derivStepSize=derivStepSize)
   #print 'shear mult bias'
   #print shearMultBiasPar.high - shearMultBiasPar.low
   photoZPar = PhotoZParams(nBins=nBins, outliers=0.1, derivStepSize=derivStepSize)
   #print 'photo z'
   #print photoZPar.high - photoZPar.low

   # Perform the Fisher forecast
   fish = FisherLsst(cosmoPar, galaxyBiasPar, shearMultBiasPar, photoZPar, nBins=nBins, nL=nL, fsky=fsky, fNk=fNk, magBias=magBias, name=name, nProc=nProc, save=save)

   # read the derivatives
   deriv[iStep,:,:] = fish.derivativeDataVector
   
   # extract the posterior uncertainties
   par[iStep], s[iStep,:] = fish.computePosterior()


##################################################################################
# Fiducial analysis, just to get the figure path...

cosmoPar = CosmoParams(massiveNu=massiveNu, wCDM=wCDM, curvature=curvature, PlanckPrior=PlanckPrior)
galaxyBiasPar = GalaxyBiasParams(nBins=nBins)
shearMultBiasPar = ShearMultBiasParams(nBins=nBins)
photoZPar = PhotoZParams(nBins=nBins, outliers=0.1)

fish = FisherLsst(cosmoPar, galaxyBiasPar, shearMultBiasPar, photoZPar, nBins=nBins, nL=nL, fsky=fsky, fNk=fNk, magBias=magBias, name=None, nProc=nProc, save=False)



##################################################################################
##################################################################################

'''
# Check parameters
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.fill_between(derivStepSizes, -0.1*np.ones(nStep), 0.1*np.ones(nStep), facecolor='gray', edgecolor='', alpha=0.5)
#
for iPar in range(nPar):
   ax.plot(derivStepSizes, s[:,iPar] / s[iStepFid,iPar]-1., label=par[0].namesLatex[iPar], alpha=0.1)
#
#ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'Derivative step size')
ax.set_ylabel(r'$\sigma(\text{step size}) / \sigma_\text{fiducial} - 1$')

plt.show()
'''


##################################################################################
# Check derivatives

'''
# relative difference with fiducial step size
derivRelDiff = deriv / deriv[iStepFid, :, :][np.newaxis,:,:] - 1.
derivRelDiff *= fish.lMaxMask[np.newaxis,np.newaxis,:]

##################################################################################

# find what C_ell has a problem, for a given cosmo parameter
#fish.plotErrorDerivativeDataVectorCosmo(derivRelDiff[iStepFid-1,:,:], show=False)

##################################################################################


# for each parameter, combine all the data,
# to find if the parameter is ok or not
minRelDiff = np.min(derivRelDiff, axis=-1)
meanRelDiff = np.mean(derivRelDiff, axis=-1) * len(fish.lMaxMask) / np.sum(fish.lMaxMask) 
maxRelDiff = np.max(derivRelDiff, axis=-1)
#sRelDiff = np.std(derivRelDiff, axis=-1) * len(fish.lMaxMask) / np.sum(fish.lMaxMask) 


fig=plt.figure(0)
ax=fig.add_subplot(111)
#
for iPar in range(cosmoPar.nPar):
#for iPar in [3]:
#for iPar in [1]:
   ax.fill_between(derivStepSizes, 0.5*iPar - 0.1*np.ones(nStep), 0.5*iPar + 0.1*np.ones(nStep), facecolor='gray', edgecolor='', alpha=0.5)
   #
   ax.plot(derivStepSizes, 0.5*iPar +  meanRelDiff[:,iPar], label=par[0].namesLatex[iPar], alpha=1.)
   #
   ax.fill_between(derivStepSizes, 0.5*iPar + minRelDiff[:,iPar], 0.5*iPar + maxRelDiff[:,iPar], alpha=0.5)
   #ax.fill_between(derivStepSizes, 0.5*iPar + meanRelDiff[:,iPar] -  sRelDiff[:,iPar], 0.5*iPar + meanRelDiff[:,iPar] + sRelDiff[:,iPar], alpha=0.5)
#
ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'Derivative step size')
ax.set_ylabel(r'Relat change in derivative')

plt.show()
'''



##################################################################################
##################################################################################
# Check diagonal Fisher matrix


diagFisher = np.zeros((nStep, nPar))
for iStep in range(nStep):
   for iPar in range(nPar):
      diagFisher[iStep, iPar] = par[iStep].fisher[iPar, iPar]

diagFisherRelDiff = np.sqrt(diagFisher[:,:]) / np.sqrt(diagFisher[iStepFid,:][np.newaxis,:]) - 1.


fig=plt.figure(0, figsize=(8,6))
ax=fig.add_subplot(111)
#
for iPar in range(cosmoPar.nPar):
   ax.fill_between(derivStepSizes, 0.25*iPar - 0.1*np.ones(nStep), 0.25*iPar + 0.1*np.ones(nStep), facecolor='gray', edgecolor='', alpha=0.5)
   ax.axhline(0.25*iPar, color='gray', linestyle='--')
   #
   ax.plot(derivStepSizes[iStepFid], 0.25*iPar +  diagFisherRelDiff[iStepFid,iPar], 'ko')
   ax.plot(derivStepSizes, 0.25*iPar +  diagFisherRelDiff[:,iPar], label=par[0].namesLatex[iPar], alpha=1.)
#
ax.set_yticks(0.25*np.arange(cosmoPar.nPar))
ax.set_yticklabels(par[0].namesLatex[range(cosmoPar.nPar)])
ax.set_xlim((derivStepSizes[0], derivStepSizes[-1]))
ax.set_ylim((-0.25, (cosmoPar.nPar)*0.25))
#ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'Derivative step size')
ax.set_title(r'Relative change in $\sigma_\alpha^\text{unmarg} = \left(F_{\alpha, \alpha}\right)^{-1/2}$')
#
fig.savefig(fish.figurePath+'/convergence_uncertainty_unmarg.pdf', bbox_inches='tight')

plt.show()


##################################################################################
# Check marginalized cosmo param constraints


sRelDiff = s[:,:] / s[iStepFid,:][np.newaxis,:] - 1.


fig=plt.figure(0, figsize=(8,6))
ax=fig.add_subplot(111)
#
for iPar in range(cosmoPar.nPar):
   ax.fill_between(derivStepSizes, 0.25*iPar - 0.1*np.ones(nStep), 0.25*iPar + 0.1*np.ones(nStep), facecolor='gray', edgecolor='', alpha=0.5)
   ax.axhline(0.25*iPar, color='gray', linestyle='--')
   #
   ax.plot(derivStepSizes[iStepFid], 0.25*iPar +  sRelDiff[iStepFid,iPar], 'ko')
   ax.plot(derivStepSizes, 0.25*iPar +  sRelDiff[:,iPar], label=par[0].namesLatex[iPar], alpha=1.)
#
ax.set_yticks(0.25*np.arange(cosmoPar.nPar))
ax.set_yticklabels(par[0].namesLatex[range(cosmoPar.nPar)])
ax.set_xlim((derivStepSizes[0], derivStepSizes[-1]))
ax.set_ylim((-0.25, (cosmoPar.nPar)*0.25))
#ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'Derivative step size')
ax.set_title(r'Relative change in $\sigma_\alpha^\text{marg} = F^{-1}_{\alpha, \alpha}$')
#
fig.savefig(fish.figurePath+'/convergence_uncertainty_marg.pdf', bbox_inches='tight')

plt.show()


##################################################################################
##################################################################################
# For the marginalized posteriors, vary only the corresponding parameter step,
# not all of them at once

# first param is the one we will forecast
derivMixed  = np.zeros((nPar, nStep, nPar, nData))
fisherData  = np.zeros((nPar, nStep, nPar, nPar))
sMixed = np.zeros((nPar, nStep))

for iPar in range(cosmoPar.nPar):
   print("Computing parameter "+str(iPar))
   for iStep in range(nStep):

      # generate the derivative where only the iPar par is rescaled
      derivMixed[iPar,iStep,:,:] = fish.derivativeDataVector.copy()
      derivMixed[iPar,iStep,iPar,:] = deriv[iStep,iPar,:]

      # get the corresponding Fisher
      fisherData[iPar,iStep,:,:] = fish.generateFisher(mask=None, deriv=derivMixed[iPar,iStep,:,:])

      # get the posterior constraints
      _, sMixed[iPar,iStep] = fish.computePosterior(fisherData=fisherData[iPar,iStep,:,:])






sMixedRelDiff = sMixed[:,:] / s[:,iStepFid][:,np.newaxis] - 1.


fig=plt.figure(0, figsize=(8,6))
ax=fig.add_subplot(111)
#
for iPar in range(cosmoPar.nPar):
   ax.fill_between(derivStepSizes, 0.25*iPar - 0.1*np.ones(nStep), 0.25*iPar + 0.1*np.ones(nStep), facecolor='gray', edgecolor='', alpha=0.5)
   ax.axhline(0.25*iPar, color='gray', linestyle='--')
   #
   ax.plot(derivStepSizes[iStepFid], 0.25*iPar +  sMixedRelDiff[iStepFid,iPar], 'ko')
   ax.plot(derivStepSizes, 0.25*iPar +  sMixedRelDiff[:,iPar], label=par[0].namesLatex[iPar], alpha=1.)
#
ax.set_yticks(0.25*np.arange(cosmoPar.nPar))
ax.set_yticklabels(par[0].namesLatex[range(cosmoPar.nPar)])
ax.set_xlim((derivStepSizes[0], derivStepSizes[-1]))
ax.set_ylim((-0.25, (cosmoPar.nPar)*0.25))
#ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'Derivative step size')
ax.set_title(r'Relative change in $\sigma_\alpha^\text{marg} = F^{-1}_{\alpha, \alpha}$')
#
fig.savefig(fish.figurePath+'/convergence_uncertainty_marg_onebyone.pdf', bbox_inches='tight')

plt.show()


