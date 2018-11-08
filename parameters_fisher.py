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
      self.fisher: nPar*nPar Fisher matrix of parameters
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
      newPar.fisher = self.fisher
      return newPar



   def addParams(self, newParams, pos=None):
      '''Adds new parameters to the parameter set.
      The new parameters are inserted at index position pos:
      pos=0 puts the new params at the start,
      pos=self.nPar puts them at the end
      '''
      if pos is None:
         pos = self.nPar
      # concatenate parameter names and values
      self.names = np.concatenate((self.names[:pos], newParams.names, self.names[pos:]))
      self.namesLatex = np.concatenate((self.namesLatex[:pos], newParams.namesLatex, self.namesLatex[pos:]))
      self.fiducial = np.concatenate((self.fiducial[:pos], newParams.fiducial, self.fiducial[pos:]))
      self.high = np.concatenate((self.high[:pos], newParams.high, self.high[pos:]))
      self.low = np.concatenate((self.low[:pos], newParams.low, self.low[pos:]))
      # combine the Fisher priors
      combined = np.zeros((self.nPar+newParams.nPar, self.nPar+newParams.nPar))
      combined[:pos, :pos] = self.fisher[:pos, :pos]
      combined[pos:pos+newParams.nPar, pos:pos+newParams.nPar] = newParams.fisher[:,:]
      combined[pos+newParams.nPar:, pos+newParams.nPar:] = self.fisher[pos:, pos:]
      self.fisher = combined
      # increase the number of parameters
      self.nPar += newParams.nPar



   def extractParams(self, I, marg=True):
      '''Create new parameter class
      with only a subset of the current parameters.
      I: indices of the parameters to keep
      If marg is True, marginalize over the parameters that are not kept.
      If marg is False, fix the parameters that are not kept
      '''
      newPar = Parameters()
      newPar.nPar = len(I)
      newPar.names = self.names[I]
      newPar.namesLatex = self.namesLatex[I]
      newPar.fiducial = self.fiducial[I]
      newPar.high = self.high[I]
      newPar.low = self.low[I]
      # extract Fisher matrix, by marginalizing/fixing the other parameters
      J = np.ix_(I,I)   # to extract the corresponding rows and columns
      if marg:
         inv = np.linalg.inv(self.fisher)
         newPar.fisher = np.linalg.inv(inv[J])
      else:
         newPar.fisher = self.fisher[J]
      return newPar


   def reorderParams(self, I):
      '''Reorders the parameters.
      I should be a permutation of range(self.nPar).
      '''
      newPar = Parameters()
      newPar.nPar = self.nPar
      newPar.names = self.names[I]
      newPar.namesLatex = self.namesLatex[I]
      newPar.fiducial = self.fiducial[I]
      newPar.high = self.high[I]
      newPar.low = self.low[I]
      J = np.ix_(I,I)
      newPar.fisher = self.fisher[J]
      return newPar


   def convertParamsFisher(self, D):
      '''D should be the nPar*nPar Jacobian matrix:
      D_{i,j} = \partial old[i] / \partial new[j]
      '''
      newFisher = np.dot(self.fisher, D)
      newFisher = np.dot(D.transpose(), newFisher)
      return newFisher
   
   
   def paramUncertainties(self, marg=True):
      """Returns the marginalized 1-sigma uncertainties of the parameters.
      """
      if marg:
         invFisher = np.linalg.inv(self.fisher)
         std = np.sqrt(np.diag(invFisher))
      else:
         std = 1. / np.sqrt(np.diag(invFisher))
      return std
   
   

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
         invFisher = np.linalg.inv(self.fisher)
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

   def printParams(self, IPar=None, path=None):
      '''Show the parameter names, fiducial values and priors.
      IPar (optional): indices of parameters to show
      '''
      if IPar is None:
         IPar = range(self.nPar)
      
      # if no path, print to standard output
      if path is None:
         # if a Fisher prior is available, print the uncertainties
         try:
            invFisher = np.linalg.inv(self.fisher)
            std = np.sqrt(np.diag(invFisher))
            for iPar in IPar:
               print self.names[iPar]+" = "+str(self.fiducial[iPar])+" +/- "+str(std[iPar])
         # otherwise, just print the fiducial values
         except:
            for iPar in IPar:
               print self.names[iPar]+" = "+str(self.fiducial[iPar])
         
      # if a path is provided, save to file
      else:
         with open(path, 'w') as f:
            # if a Fisher prior is available, print the uncertainties
            try:
               invFisher = np.linalg.inv(self.fisher)
               std = np.sqrt(np.diag(invFisher))
               for iPar in IPar:
                  f.write(self.names[iPar]+" = "+str(self.fiducial[iPar])+" +/- "+str(std[iPar])+"\n")
            # otherwise, just print the fiducial values
            except:
               for iPar in IPar:
                  f.write(self.names[iPar]+" = "+str(self.fiducial[iPar])+"\n")



   def plotContours(self, IPar=None, marg=True, lim=4., color='#E10014', path=None):
      '''Show confidence ellipses.
      IPar (optional): indices of parameters to show
      '''
      if IPar is None:
         par = self.copy()
      else:
         par = self.extractParams(IPar, marg=marg)
      invFisher = np.linalg.inv(par.fisher)
      

      # generate the contour plot
      fig=plt.figure(0, figsize=(18, 16))
      gs = gridspec.GridSpec(par.nPar, par.nPar)#, height_ratios=[1, 1, 1])
      gs.update(hspace=0.)
      gs.update(wspace=0.)

      # loop over each column
      for j in range(par.nPar):
         # j=i case: just plot the Gaussian pdf
         ax=plt.subplot(gs[j,j])
         
         # define Gaussian with correct mean and variance
         mean = par.fiducial[j]
         sigma = np.sqrt(invFisher[j,j])
         X = np.linspace(mean - lim*sigma, mean + lim*sigma, 201)
         Y = np.exp(-0.5*(X-mean)**2 / sigma**2)
         
         # plot it
         ax.plot(X, Y, color)
         
         # 68-95% confidence regions
         ax.plot([mean,mean], [0,1], color='b')
         ax.fill_between(X, 0., Y, where=np.abs(X-mean)<sigma, color=color, alpha=0.7)
         ax.fill_between(X, 0., Y, where=np.abs(X-mean)<2.*sigma, color=color, alpha=0.3)
         
         # print the marginalized constraint
         titleStr = par.namesLatex[j]+r'$ = '+str(np.round(par.fiducial[j], 4))+r' \pm'+str(np.round(sigma, 4))+r'$'
         ax.set_title(titleStr, fontsize=12)
         
         # always remove y ticks
         ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
         
         # remove x ticks, except if last parameter
         if j<par.nPar-1:
            ax.tick_params(axis='x', which='both', top=False, bottom=False, labelbottom=False)

         # limits
         ax.set_xlim((mean - lim*sigma, mean + lim*sigma))
         ax.set_ylim((0., 1.1))


         # loop over each row i>j
         for i in range(j+1, par.nPar):
            ax=plt.subplot(gs[i,j], sharex=ax)
            
            # means and cov
            meanx = par.fiducial[j]
            meany = par.fiducial[i]
            #
            sx = np.sqrt(invFisher[j,j])
            sy = np.sqrt(invFisher[i,i])
            sxy = invFisher[i,j]

            # plot the fiducial value
            ax.plot([meanx], [meany], 'b.')

            # covariance contour
            x = np.linspace(meanx - lim*sx, meanx + lim*sx, 501)
            y = np.linspace(meany - lim*sy, meany + lim*sy, 501)
            X, Y = np.meshgrid(x, y)
            Z = bivariate_normal(X, Y, sigmax=sx, sigmay=sy, mux=meanx, muy=meany, sigmaxy=sxy)
            norm = bivariate_normal(0., 0., sigmax=sx, sigmay=sy, mux=0., muy=0., sigmaxy=sxy)
            Z = np.log(norm/Z)   # returns the chi squared

            # confidence level contours
            conf_level = np.array([0.68, 0.95])
            chi2 = -2. * np.log(1. - conf_level)
            ax.contourf(X,Y,Z, [0., chi2[0]], colors=color, alpha=0.7)
            ax.contourf(X,Y,Z, [0., chi2[1]], colors=color, alpha=0.3)
            
            # limits
            mean = par.fiducial[i]
            sigma = np.sqrt(invFisher[i,i])
            ax.set_ylim((meany - lim*sy, meany + lim*sy))
            
            # tick labels
            if j==0:
               ax.set_ylabel(par.namesLatex[i], fontsize=14)
               plt.setp(ax.get_yticklabels(), fontsize=12)
            else:
               ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
            #
            plt.setp(ax.get_xticklabels(), visible=False)

         # tick labels
         plt.setp(ax.get_xticklabels(), visible=True, rotation=45)
         ax.set_xlabel(par.namesLatex[j], fontsize=14)
         plt.setp(ax.get_xticklabels(), fontsize=12)
      
      if path is None:
         plt.show()
      else:
         fig.savefig(path)
         fig.clf()



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
      self.fisher = np.diagflat(1./self.priorStd**2)

##################################################################################

class PhotoZParams(Parameters):

   def __init__(self, nBins=2, dzFid=0., szFid=0.05, dzStd=0.002, szStd=0.003, outliers=0., outliersStd=0.05):
      self.nPar = 2 * nBins
      if outliers<>0.:
         self.nPar += nBins*(nBins-1)
      self.outliers = outliers

      # bias and std dev of photo-z
      dz = np.array(['dz'+str(iBin) for iBin in range(nBins)])
      sz = np.array(['sz'+str(iBin) for iBin in range(nBins)])
      if outliers==0.:
         self.names = np.concatenate((dz, sz))
      else:
         cij = np.array(['c_'+str(iBin)+','+str(jBin) for iBin in range(nBins) for jBin in list(set(range(nBins)) - set([iBin]))])
         self.names = np.concatenate((dz, sz, cij))
      
      # latex strings for the names
      dz = np.array([r'$\delta z_{'+str(iBin)+'}$' for iBin in range(nBins)])
      sz = np.array([r'$\sigma z_{'+str(iBin)+'}/(1+z)$' for iBin in range(nBins)])
      if outliers==0.:
         self.namesLatex = np.concatenate((dz, sz))
      else:
         cij = np.array([r'$c_{'+str(iBin)+','+str(jBin)+'}$' for iBin in range(nBins) for jBin in list(set(range(nBins)) - set([iBin]))])
         self.namesLatex = np.concatenate((dz, sz, cij))

      # fiducial values
      dz = np.array([dzFid for iBin in range(nBins)])
      sz = np.array([szFid for iBin in range(nBins)])
      if outliers==0.:
         self.fiducial = np.concatenate((dz, sz))
      else:
         cij = np.array([outliers/(nBins-1.) for iBin in range(nBins) for jBin in list(set(range(nBins)) - set([iBin]))])
         self.fiducial = np.concatenate((dz, sz, cij))

#      # high values
#      dz = np.array([dzFid+0.05 for iBin in range(nBins)])
#      sz = np.array([szFid+0.01 for iBin in range(nBins)])
#      self.high = np.concatenate((dz, sz))
#      # low values
#      dz = np.array([dzFid-0.05 for iBin in range(nBins)])
#      sz = np.array([szFid-0.01 for iBin in range(nBins)])
#      self.low = np.concatenate((dz, sz))

      # high values
      dz = np.array([dzFid+0.002 for iBin in range(nBins)])
      sz = np.array([szFid+0.003 for iBin in range(nBins)])
      if outliers==0.:
         self.high = np.concatenate((dz, sz))
      else:
         cij = np.array([(1.+0.5)*outliers/(nBins-1.) for iBin in range(nBins) for jBin in list(set(range(nBins)) - set([iBin]))])
         self.high = np.concatenate((dz, sz, cij))

      # low values
      dz = np.array([dzFid-0.002 for iBin in range(nBins)])
      sz = np.array([szFid-0.003 for iBin in range(nBins)])
      if outliers==0.:
         self.low = np.concatenate((dz, sz))
      else:
         cij = np.array([(1.-0.5)*outliers/(nBins-1.) for iBin in range(nBins) for jBin in list(set(range(nBins)) - set([iBin]))])
         self.low = np.concatenate((dz, sz, cij))

      # std dev of parameter prior
      dz = np.array([dzStd for iBin in range(nBins)])
      sz = np.array([szStd for iBin in range(nBins)])
      if outliers==0.:
         self.priorStd = np.concatenate((dz, sz))
      else:
         cij = np.array([outliersStd/(nBins-1.) for iBin in range(nBins) for jBin in list(set(range(nBins)) - set([iBin]))])
         self.priorStd = np.concatenate((dz, sz, cij))
      # corresponding Fisher matrix of priors
      self.fisher = np.diagflat(1./self.priorStd**2)


##################################################################################

class GalaxyBiasParams(Parameters):

   def __init__(self, nBins):
      '''Galaxy biases, relative to their fiducial values
      '''
      self.nPar = nBins
      self.names = np.array(['bg'+str(iBin) for iBin in range(self.nPar)])
      self.namesLatex = np.array([r'$b_{g\;'+str(iBin)+'} / b_{g\;'+str(iBin)+r'}^\text{fid}$' for iBin in range(self.nPar)])
      
      self.fiducial = np.array([1. for iBin in range(self.nPar)])
      # high/low values for derivative
      self.high = self.fiducial * 1.05
      self.low = self.fiducial * 0.95

      self.fisher = np.zeros((self.nPar, self.nPar))


##################################################################################

class CosmoParams(Parameters):

   def __init__(self, massiveNu=False, wCDM=False, curvature=False, PlanckPrior=False):
      '''Step sizes inspired from Allison+15
      '''
      # indices to keep for various relevant combinations:
      self.IFull = range(10)
      self.ILCDMMnuW0Wa = range(9)
      self.ILCDMMnu = range(7)
      self.ILCDM = range(6)
      self.ILCDMW0Wa = [0,1,2,3,4,5,7,8]
      self.ILCDMMnuCurv = [0,1,2,3,4,5,6,9]
      self.ILCDMW0WaCurv = [0,1,2,3,4,5,7,8,9]
      # base cosmology
      self.nPar = 6
      self.names = np.array(['Omega_cdm', 'Omega_b', 'A_s', 'n_s', 'h', 'tau_reio'])
      self.namesLatex = np.array([r'$\Omega^0_\text{CDM}$', r'$\Omega^0_\text{b}$', r'$A_\text{s} / A_\text{s}^\text{fid}$', r'$n_\text{s}$', r'$h_0$', r'$\tau$'])
      self.fiducial = np.array([0.267, 0.0493, 1., 0.9624, 0.6712, 0.06])
      self.high = np.array([0.267 + 0.0066, 0.0493 + 0.0018, 1. + 1.e-10/2.3e-9, 0.9624 + 0.01, 0.6712 + 0.067, 0.06 + 0.02])
      self.low = np.array([0.267 - 0.0066, 0.0493 - 0.0018, 1. - 1.e-10/2.3e-9, 0.9624 - 0.01, 0.6712 - 0.067, 0.06 - 0.02])
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
                           'output': 'mPk',#'dTk vTk lCl tCl pCl mPk',
                           'P_k_max_1/Mpc': 1.,  #10.,
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
                           'output': 'mPk',#'dTk vTk lCl tCl pCl mPk',
                           'P_k_max_1/Mpc': 1.,  #10.,
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
                           'output': 'mPk',#'dTk vTk lCl tCl pCl mPk',
                           'P_k_max_1/Mpc': 1.,  #10.,
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

         self.high = np.concatenate((self.high, np.array([0.+0.02])))
         self.low = np.concatenate((self.low, np.array([0.-0.02])))
         self.paramsClassy.update({
                                 # Curvature
                                 'Omega_k': 0.,
                                 })
         self.paramsClassyHigh.update({
                                 # Curvature
                                 'Omega_k': 0.+0.02,
                                 })
         self.paramsClassyLow.update({
                                 # Curvature
                                 'Omega_k': 0.-0.02,
                                 })

      # load Planck priors if requested
      if PlanckPrior:
         self.loadPlanckPrior()
      else:
         self.fisher = np.zeros((self.nPar, self.nPar))
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


#   def loadPlanckPrior(self, test=True):
#      # read Planck priors from Pat McDonald
#      patPlanck = PatPlanckParams()
#      if test:
#         print "initial parameters"
#         print patPlanck.names
#
#      # reorder the params to more closely match my order
#      #I = [0, 1, 8, 9, where_h_should_be, 13, 6, 3, 4, 5, 2, 7, 10, 11, 12]
#      I = [0, 1, 8, 9, 13, 6, 3, 4, 5, 2, 7, 10, 11, 12]
#      newPar = patPlanck.reorderParams(I)
#      if test:
#         print "reordered parameters"
#         print newPar.names
#
#      # add hubble parameter
#      H = Parameters()
#      H.nPar = 1
#      H.names = np.array(['h'])
#      H.namesLatex = np.array([r'$h$'])
#      H.fiducial = np.array([0.067])
#      H.low = np.array([0.])
#      H.high = np.array([0.])
##!!! what to put for the Planck h uncertainty?
#      H.fisher = np.array([[1.]])
#
#      # insert it at the right spot
#      newPar.addParams(H, pos=4)
#      if test:
#         print "adding h"
#         print newPar.names
#
#      # remove the extra parameters
#      I = range(10)
#      newPar = newPar.extractParams(I, marg=True)
#      if test:
#         print "remove extra parameters"
#         print newPar.names
#
#      # convert units and logs in the Fisher matrix
#      D = np.diagflat(np.ones(10))
#      D[0,0] = 0.067**2 # dOmh2/dOm = h^2
#      D[0,1] = 0.067**2 # dOmh2/dOb = h^2
#      D[0,2] = 1./93.# dOmh2 / dMnu = h^2 dOnu/dMnu = 1/93., in 1/eV
#      D[1,1] = 0.067**2 # dObh2/dOb = h^2
#      D[2,2] = 1./(np.log(10.) * 2.3e-9)  # dlog10As/dAs = 1/(ln10 * A_S)
#
#      newPar.fisher = newPar.convertParamsFisher(D)
#      self.fisher = newPar.fisher.copy()
#      return


   def loadPlanckPrior(self):
      path = "./input/parameters/Fisher_Planck_TT_TE_EE_lowP.txt"
      self.fisher = np.genfromtxt(path)





##################################################################################

class PatPlanckParams(Parameters):
   '''Fisher matrix for Planck parameters from Pat McDonald, private communication.
   '''
   
   def __init__(self):
      # read Planck fisher matrix from Pat
      path = "./input/parameters/fullP.tau0.066.r1.lT2_1400.lP8_2500_3_1_0.2_fish.txt"
      self.fisher = np.genfromtxt(path)
      
      self.nPar = 14
      self.names = np.array(['omega_m_h2', 'omega_b_h2', 'theta_s', 'w_0', 'w_a', 'Omega_k', 'M_nu', 'N_nueff', 'log10A_S', 'n_S', 'alpha_S', 'beta_S', 'r', 'tau'])
      self.namesLatex = np.array([r'$\omega_m h^2$', r'$\omega_b h^2$', r'$\theta_s$', r'$w_0$', r'$w_a$', r'$\Omega_k$', r'$M_\nu$', r'$N_{\nu, \text{eff}}$', r'$\log_{10}A_S$', r'$n_S$', r'$\alpha_S$', r'$\beta_S$', r'$r$', r'$\tau$'])
      self.fiducial = np.array([0.141745, 0.0223, 0.595417, -1., 0., 0., 0.06, 3.046, -8.66918, 0.9667, 0., 0., 0., 0.066])
      self.high = np.zeros_like(self.fiducial)
      self.low = np.zeros_like(self.fiducial)


##################################################################################

class DndzParams(Parameters):

   def __init__(self, nBins=2, nZ=2, sNgal=1.e-4):
      self.nBins = nBins
      self.nZ = nZ
      self.nPar = nBins * nZ
      
      # the parameters are defined as an additive correction to dn/dz
      # these are the flattened version of (nBins, nZ)
      self.fiducial = np.zeros(self.nPar)
      self.high = 1.e-2 * np.ones(self.nPar)
      self.low = np.zeros(self.nPar)
      self.names = np.empty(self.nPar, dtype=object)
      self.namesLatex = np.empty(self.nPar, dtype=object)
      i = 0
      for iBin in range(nBins):
         for iZ in range(nZ):
            self.names[i] = 'dn_'+str(iBin)+'dz_'+str(iZ)
            self.namesLatex[i] = r'$dn_{'+str(iBin)+'}/dz_{'+str(iZ)+'}$'
            i += 1
   
      # Prior: the dn/dz have to integrate to n_gal, the known nb of galaxies per bin
      self.fisher = np.zeros((self.nPar, self.nPar))
      # the Fisher matrix is block-diagonal, and the blocks on the diagonal are all ones
      for iBin in range(nBins):
         I = range(iBin*nZ, (iBin+1)*nZ)
         J = np.ix_(I,I)
         self.fisher[J] = np.ones((nZ, nZ))
      # set the uncertainty on the measured ngal in each bin
      self.fisher /= sNgal**2





















