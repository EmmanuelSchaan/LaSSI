from headers import *

##################################################################################

class Parameters(object):
   
   def __init__(self):
      '''Requires:
      self.nPar: number of parameters
      self.names: array of strings for plain text parameter names
      self.namesLatex: array of strings for latex names
      self.fiducial: array of fiducial parameter values
      self.high: high values for the numerical derivative
      self.low: low values for the numerical derivative
      self.priorFisher: nPar*nPar Fisher matrix of priors on parameters
      '''
      pass


   def addParams(self, newParams):
      '''Adds new parameters to the parameter set.
      '''
      # concatenate parameter names and values
      self.names = np.concatenate((self.names, newParams.names))
      self.namesLatex = np.concatenate((self.namesLatex, newParams.namesLatex))
      self.fiducial = np.concatenate((self.fiducial, newParams.fiducial))
      self.high = np.concatenate((self.high, newParams.high))
      self.low = np.concatenate((self.low, newParams.low))
      # combine the Fisher priors
      combined = np.zeros((self.nPar+newParams.nPar, self.nPar+newParams.nPar))
#      print self.nPar, newParams.nPar
#      print shape(combined[:self.nPar, :self.nPar])
      combined[:self.nPar, :self.nPar] = self.priorFisher[:,:]
      combined[self.nPar:, self.nPar:] = newParams.priorFisher[:,:]
      self.priorFisher = combined
      # increase the number of parameters
      self.nPar += newParams.nPar

   def extractParams(self):
      '''Create new parameter class
      with only a subset of the current parameters.
      '''
      pass


   def plotParams(self):
      '''Show the parameter names, fiducial/high/low values, priors.
      '''
      invFisher = np.linalg.inv(self.priorFisher)
      std = np.sqrt(np.diag(invFisher))
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.errorbar(range(self.nPar), self.fiducial, yerr=std, fmt='o')
      #
      ax.set_xticks(range(self.nPar))
      ax.set_xticklabels(self.namesLatex, fontsize=24)
      [l.set_rotation(45) for l in ax.get_xticklabels()]

      plt.show()


##################################################################################

class ShearMultBiasParams(Parameters):

   def __init__(self, nBins=2, mStd=0.005):
      self.nPar = nBins

      self.names = np.array(['m'+str(iBin) for iBin in range(self.nPar)])
      self.namesLatex = np.array([r'$m_{'+str(iBin)+'}$' for iBin in range(self.nPar)])
      self.fiducial = np.array([0. for iBin in range(self.nPar)])
      self.high = np.array([0.05 for iBin in range(self.nPar)])
      self.low = np.array([-0.05 for iBin in range(self.nPar)])
      self.priorStd = np.array([mStd for iBin in range(self.nPar)])
      self.priorFisher = np.diagflat(1./self.priorStd**2)

##################################################################################

class PhotoZParams(Parameters):

   def __init__(self, nBins=2, dzFid=0., szFid=0.05, dzStd=0.002, szStd=0.003):
      self.nPar = 2 * nBins

      # bias and std dev of photo-z
      dz = np.array(['dz'+str(iBin) for iBin in range(nBins)])
      sz = np.array(['sz'+str(iBin) for iBin in range(nBins)])
      self.names = np.concatenate((dz, sz))
      
      # latex strings for the names
      dz = np.array([r'$\delta z_{'+str(iBin)+'}$' for iBin in range(nBins)])
      sz = np.array([r'$\sigma z_{'+str(iBin)+'}$' for iBin in range(nBins)])
      self.namesLatex = np.concatenate((dz, sz))

      # fiducial values
      dz = np.array([dzFid for iBin in range(nBins)])
      sz = np.array([szFid for iBin in range(nBins)])
      self.fiducial = np.concatenate((dz, sz))

      # high/low values for derivative
      self.high = self.fiducial * 1.05
      self.low = self.fiducial * 0.95

      # std dev of parameter prior
      dz = np.array([dzStd for iBin in range(nBins)])
      sz = np.array([szStd for iBin in range(nBins)])
      self.priorStd = np.concatenate((dz, sz))
      # corresponding Fisher matrix of priors
      self.priorFisher = np.diagflat(1./self.priorStd**2)



##################################################################################
#!!!!! incomplete

class GalaxyBiasParams(Parameters):

   def __init__(self, nBins, w_g):
      '''w_g: dictionary of projection kernel classes
      for each z-bin
      '''
      self.nPar = nBins
      self.names = np.array(['bg'+str(iBin) for iBin in range(self.nPar)])
      self.namesLatex = np.array([r'$b_{g\;'+str(iBin)+'}$' for iBin in range(self.nPar)])
      
      # tracer bias for the LSST gold sample,
      # from the LSST Science book, chapter 3 and 13.
      self.fiducial = np.array([w_g[iBin].bMean() for iBin in range(self.nPar)])

      # high/low values for derivative
      self.high = self.fiducial * 1.05
      self.low = self.fiducial * 0.95

      self.priorFisher = np.zeros((self.nPar, self.nPar))


##################################################################################
#!!!!! incomplete

class CosmoParams(Parameters):

   def __init__(self):
      self.nPar = 8

      self.names = np.array(['Omega_m',
                            'Omega_b',
#                            'Omega_k',
                            'A_S',
                            'n_S',
                            'w_0',
#                            'w_a',
                            'h',
                            'M_nu',
                            'tau'])

      self.namesLatex = np.array([r'$\Omega^0_m$',
                                  r'$\Omega^0_b$',
#                                  r'$\Omega^0_k$',
                                  r'$A_S$',
                                  r'$n_S$',
                                  r'$w_0$',
#                                  r'$w_a$',
                                  r'$h_0$',
                                  r'$M_\nu$',
                                  r'$\tau$'])





      










