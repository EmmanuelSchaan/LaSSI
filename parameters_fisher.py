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


   def copy(self):
      newPar = Parameters()
      newPar.nPar = self.nPar
      newPar.names = self.names
      newPar.namesLatex = self.namesLatex
      newPar.fiducial = self.fiducial
      newPar.high = self.high
      newPar.low = self.low
      newPar.priorFisher = self.priorFisher
      return newPar

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


   def plotParams(self, IPar=None):
      '''Show the parameter names, fiducial values and priors.
      IPar (optional): indices of parameters to show
      '''
      if IPar is None:
         IPar = range(self.nPar)

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # if a Fisher prior is available, plot it
      try:
         invFisher = np.linalg.inv(self.priorFisher)
         std = np.sqrt(np.diag(invFisher))
         ax.errorbar(range(len(IPar)), self.fiducial[IPar], yerr=std[IPar], fmt='o')
      # otherwise, just show the fiducial values
      except:
         ax.errorbar(range(len(IPar)), self.fiducial[IPar], fmt='o')
      #
      ax.set_xticks(range(len(IPar)))
      ax.set_xticklabels(self.namesLatex[IPar], fontsize=24)
      [l.set_rotation(45) for l in ax.get_xticklabels()]

      plt.show()

   def printParams(self, IPar=None):
      '''Show the parameter names, fiducial values and priors.
      IPar (optional): indices of parameters to show
      '''
      if IPar is None:
         IPar = range(self.nPar)

      # if a Fisher prior is available, print the uncertainties
      try:
         invFisher = np.linalg.inv(self.priorFisher)
         std = np.sqrt(np.diag(invFisher))
         for iPar in IPar:
            print self.names[iPar]+" = "+str(self.fiducial[iPar])+" +/- "+str(std[iPar])
      # otherwise, just print the fiducial values
      except:
         for iPar in IPar:
            print self.names[iPar]+" = "+str(self.fiducial[iPar])



#   def plotParamStd(self, IPar=None):
#      '''Show the parameter names and priors.
#      IPar (optional): indices of parameters to show
#      '''
#      if IPar is None:
#         IPar = range(self.nPar)
#
#      fig=plt.figure(0)
#      ax=fig.add_subplot(111)
#      #
#      # if a Fisher prior is available, plot it
#      try:
#         invFisher = np.linalg.inv(self.priorFisher)
#         std = np.sqrt(np.diag(invFisher))
#         ax.errorbar(range(len(IPar)), np.zeros_like(IPar), yerr=std[IPar], fmt='o')
#      # otherwise, just show the fiducial values
#      except:
#         print "Problem with the Fisher priors"
#         fig.clf()
#         return
#      #
#      ax.set_xticks(range(len(IPar)))
#      ax.set_xticklabels(self.namesLatex[IPar], fontsize=24)
#      [l.set_rotation(45) for l in ax.get_xticklabels()]
#
#      plt.show()



   def plotContours(self, IPar=None):
      '''Show confidence ellipses.
      IPar (optional): indices of parameters to show
      '''
      if IPar is None:
         IPar = range(self.nPar)

      pass


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
      sz = np.array([r'$\sigma z_{'+str(iBin)+'}/(1+z)$' for iBin in range(nBins)])
      self.namesLatex = np.concatenate((dz, sz))

      # fiducial values
      dz = np.array([dzFid for iBin in range(nBins)])
      sz = np.array([szFid for iBin in range(nBins)])
      self.fiducial = np.concatenate((dz, sz))

      # high values
      dz = np.array([dzFid+0.05 for iBin in range(nBins)])
      sz = np.array([szFid+0.01 for iBin in range(nBins)])
      self.high = np.concatenate((dz, sz))
      # low values
      dz = np.array([dzFid-0.05 for iBin in range(nBins)])
      sz = np.array([szFid-0.01 for iBin in range(nBins)])
      self.low = np.concatenate((dz, sz))

      # std dev of parameter prior
      dz = np.array([dzStd for iBin in range(nBins)])
      sz = np.array([szStd for iBin in range(nBins)])
      self.priorStd = np.concatenate((dz, sz))
      # corresponding Fisher matrix of priors
      self.priorFisher = np.diagflat(1./self.priorStd**2)



##################################################################################

class GalaxyBiasParams(Parameters):

   def __init__(self, nBins):
      '''Galaxy biases, relative to their fiducial values
      '''
      self.nPar = nBins
      self.names = np.array(['bg'+str(iBin) for iBin in range(self.nPar)])
      self.namesLatex = np.array([r'$b_{g\;'+str(iBin)+'} / b_{g\;'+str(iBin)+r'}^\text{ fid}$' for iBin in range(self.nPar)])
      
      self.fiducial = np.array([1. for iBin in range(self.nPar)])
      # high/low values for derivative
      self.high = self.fiducial * 1.05
      self.low = self.fiducial * 0.95

      self.priorFisher = np.zeros((self.nPar, self.nPar))


##################################################################################

class CosmoParams(Parameters):

   def __init__(self, massiveNu=False, wCDM=False, curvature=False):
      '''Step sizes inspired from Allison+15
      '''
      # base cosmology
      self.nPar = 5#6
      self.names = np.array(['Omega_cdm', 'Omega_b', 'A_s', 'n_s', 'h'])#, 'tau_reio'])
      self.namesLatex = np.array([r'$\Omega^0_\text{CDM}$', r'$\Omega^0_\text{b}$', r'$A_\text{s}$', r'$n_\text{s}$', r'$h_0$'])#, r'$\tau$'])
      self.fiducial = np.array([0.267, 0.0493, 2.3e-9, 0.9624, 0.6712])#, 0.06])
      self.high = np.array([0.267 + 0.0066, 0.0493 + 0.0018, 2.3e-9 + 1.e-10, 0.9624 + 0.01, 0.6712 + 0.067])#, 0.06 + 0.02])
      self.low = np.array([0.267 - 0.0066, 0.0493 - 0.0018, 2.3e-9 - 1.e-10, 0.9624 - 0.01, 0.6712 - 0.067])#, 0.06 - 0.02])
      self.paramsClassy = {
                           # Cosmological parameters
                           'Omega_cdm': 0.267,
                           'Omega_b': 0.0493,
                           'A_s': 2.3e-9,
                           'n_s': 0.9624,
                           'tau_reio': 0.06,
                           'h': 0.6712,
                           # parameters
                           'reio_parametrization': 'reio_camb',
                           'output': 'mPk dTk vTk',#'lCl tCl pCl mPk',
                           'P_k_max_1/Mpc': 10.,
                           'non linear': 'halofit',
                           'z_max_pk': 100.
                           }
      self.paramsClassyHigh = {
                           # Cosmological parameters
                           'Omega_cdm': 0.267 + 0.0066,
                           'Omega_b': 0.0493 + 0.0018,
                           'A_s': 2.3e-9 + 1.e-10,
                           'n_s': 0.9624 + 0.01,
                           'tau_reio': 0.06 + 0.02,
                           'h': 0.6712 + 0.067,
                           # parameters
                           'reio_parametrization': 'reio_camb',
                           'output': 'mPk dTk vTk',#'lCl tCl pCl mPk',
                           'P_k_max_1/Mpc': 10.,
                           'non linear': 'halofit',
                           'z_max_pk': 100.
                           }
      self.paramsClassyLow = {
                           # Cosmological parameters
                           'Omega_cdm': 0.267 - 0.0066,
                           'Omega_b': 0.0493 - 0.0018,
                           'A_s': 2.3e-9 - 1.e-10,
                           'n_s': 0.9624 - 0.01,
                           'tau_reio': 0.06 - 0.02,
                           'h': 0.6712 - 0.067,
                           # parameters
                           'reio_parametrization': 'reio_camb',
                           'output': 'mPk dTk vTk',#'lCl tCl pCl mPk',
                           'P_k_max_1/Mpc': 10.,
                           'non linear': 'halofit',
                           'z_max_pk': 100.
                           }

      # neutrino masses
      if massiveNu:
         self.nPar += 1
         self.names = np.concatenate((self.names, np.array(['m_ncdm'])))
         self.namesLatex = np.concatenate((self.namesLatex, np.array([r'$M_\nu$'])))
         #
         Mnu = 0.06 # eV, minimum possible masses
         normalHierarchy = True
         # compute neutrino masses
         self.fiducial = np.concatenate((self.fiducial, np.array([Mnu])))
         nuMasses = self.computeNuMasses(Mnu, normal=normalHierarchy)
         self.paramsClassy.update({
                                 # Massive neutrinos
                                 'N_ncdm': 3,
                                 'm_ncdm': str(nuMasses[0])+','+str(nuMasses[1])+','+str(nuMasses[2]),
                                 'deg_ncdm': '1, 1, 1',
                                 })
         self.low = np.concatenate((self.low, np.array([Mnu])))
         self.paramsClassyLow.update({
                                 # Massive neutrinos
                                 'N_ncdm': 3,
                                 'm_ncdm': str(nuMasses[0])+','+str(nuMasses[1])+','+str(nuMasses[2]),
                                 'deg_ncdm': '1, 1, 1',
                                 })
         self.high = np.concatenate((self.high, np.array([Mnu+0.02])))
         nuMasses = self.computeNuMasses(Mnu+0.02, normal=normalHierarchy)
         self.paramsClassyHigh.update({
                                 # Massive neutrinos
                                 'N_ncdm': 3,
                                 'm_ncdm': str(nuMasses[0])+','+str(nuMasses[1])+','+str(nuMasses[2]),
                                 'deg_ncdm': '1, 1, 1',
                                 })


      # wCDM
      if wCDM:
         self.nPar += 2
         self.names = np.concatenate((self.names, np.array(['w0_fld', 'wa_fld'])))
         self.namesLatex = np.concatenate((self.namesLatex, np.array([r'$w_0$', r'$w_a$'])))
         self.fiducial = np.concatenate((self.fiducial, np.array([-1., 0.])))
         
         self.high = np.concatenate((self.high, np.array([-1.+0.3, 0.+0.6])))
         self.low = np.concatenate((self.low, np.array([-1., 0.])))
         self.paramsClassy.update({
                                 # w0 and wa
                                 'Omega_Lambda': 0.,
                                 'w0_fld': -1.,
                                 'wa_fld': 0.,
                                 })
         self.paramsClassyHigh.update({
                                 # w0 and wa
                                 'Omega_Lambda': 0.,
                                 'w0_fld': -1.+0.3,
                                 'wa_fld': 0.+0.6,
                                 })
         self.paramsClassyLow.update({
                                 # w0 and wa
                                 'Omega_Lambda': 0.,
                                 'w0_fld': -1.,
                                 'wa_fld': 0.,
                                 })

      # curvature
      if curvature:
         self.nPar += 1
         self.names = np.concatenate((self.names, np.array(['Omega_k'])))
         self.namesLatex = np.concatenate((self.namesLatex, np.array([r'$\Omega_\text{k}$'])))
         self.fiducial = np.concatenate((self.fiducial, np.array([0.])))

         self.high = np.concatenate((self.high, np.array([0.+0.01])))
         self.low = np.concatenate((self.low, np.array([0.-0.01])))
         self.paramsClassy.update({
                                 # Curvature
                                 'Omega_k': 0.,
                                 })
         self.paramsClassyHigh.update({
                                 # Curvature
                                 'Omega_k': 0.+0.01,
                                 })
         self.paramsClassyLow.update({
                                 # Curvature
                                 'Omega_k': 0.-0.01,
                                 })

      self.priorFisher = np.zeros((self.nPar, self.nPar))
      return


   def computeNuMasses(self, mSum, normal=True):
      '''mSum: sum of neutrino masses in eV
      normal=True for normal hierarchy
      output: masses in eV
      '''
      dmsq_atm = 2.5e-3 # eV^2
      dmsq_solar = 7.6e-5 # eV^2
      if normal:
         f = lambda m0: m0 + np.sqrt(m0**2+dmsq_solar) + np.sqrt(m0**2+dmsq_solar+dmsq_atm) - mSum
         m0 = optimize.brentq(f , 0., mSum)
         result = np.array([m0, np.sqrt(m0**2+dmsq_solar), np.sqrt(m0**2+dmsq_solar+dmsq_atm)])
      else:
         f = lambda m0: m0 + np.sqrt(m0**2+dmsq_atm) + np.sqrt(m0**2+dmsq_atm+dmsq_solar) - mSum
         m0 = optimize.brentq(f , 0., mSum)
         result = np.array([m0, np.sqrt(m0**2+dmsq_atm), np.sqrt(m0**2+dmsq_atm+dmsq_solar)])
      return result


##################################################################################




