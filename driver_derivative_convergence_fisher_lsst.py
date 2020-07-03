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

save = True

#derivStepSizes = np.array([0.5, 1., 2.])
derivStepSizes = np.array([0.01, 0.05, 0.1, 0.5, 1., 1.5,  5., 10.])
nStep = len(derivStepSizes)

# store all the derivatives, to compare them
nData = 11550  # read from the output of the fisher class
nPar = 140  # read from the output of the fisher class
deriv  = np.zeros((nStep, nPar, nData))

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

   # same tomo bins for g and s
   fish = FisherLsst(cosmoPar, galaxyBiasPar, shearMultBiasPar, photoZPar, nBins=nBins, nL=nL, fsky=fsky, fNk=fNk, magBias=magBias, name=name, nProc=nProc, save=save)

   deriv[iStep,:,:] = fish.derivativeDataVector


##################################################################################
# basic plots and checks
'''
fish.plotDerivativeDataVectorCosmo(show=False)
#fish.plotSingleDerivative("gg", 0, 0)
#fish.plotSingleDerivative("ss", 0, 15)
#fish.plotSingleDerivative("gg", 0, 20)
'''
