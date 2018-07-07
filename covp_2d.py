from headers import *

##################################################################################
##################################################################################

class CovP2d(object):
   
   
   def __init__(self, Pac, Pbd, Pad, Pbc, T2d=None, HSVP2d=None, fsky=1., npairs_id='log', L=None, dL=None):
      '''Cov[Pab, Pcd] = (PacPbd + PadPbc)^2/Nmodes + T/V + HSV
      '''
      # copy classes
      self.Pac = Pac
      self.Pbd = Pbd
      self.Pad = Pad
      self.Pbc = Pbc
      #
      self.HSVP2d = HSVP2d
      self.T2d = T2d
      #
      self.fsky =  fsky
      self.OmS = 4.*np.pi*fsky
      
      # values of ell to evaluate
      if L is None or dL is None:
         self.L = np.genfromtxt("./input/Lc.txt") # center of the bins for l
         self.Nl = len(self.L)
         self.dL = np.genfromtxt("./input/dL.txt")   # width of the bins for l
      else:
         self.L = L
         self.Nl = len(self.L)
         self.dL = dL

      # load Npairs
      self.LoadNpairs(npairs_id)
      self.LoadAll()
   
   
   ##################################################################################
   
   def LoadNpairs(self, npairs_id):
      if npairs_id == "log":
         # nb of pairs in bin for l, with width dl
         Al = 2.*np.pi* self.L * self.dL + np.pi * self.dL**2
         self.Npairs = self.OmS * Al / (2.*np.pi)**2
      elif npairs_id == "lin":
         self.Npairs = self.OmS/(4.*np.pi) * (2.*self.L+1.)
      return
   
#   def LoadAll(self):
#      # Gaussian covariance contributions
#      # cosmic variance only
#      f = lambda l: self.Pac.fPinterp(l) * self.Pbd.fPinterp(l) + self.Pad.fPinterp(l) * self.Pbc.fPinterp(l)
#      self.P = np.array(map(f, self.L))
#      self.P /= self.Npairs
#      # noise only
#      f = lambda l: self.Pac.fPnoise(l) * self.Pbd.fPnoise(l) + self.Pad.fPnoise(l) * self.Pbc.fPnoise(l)
#      self.N = np.array(map(f, self.L))
#      self.N /= self.Npairs
#       noise + cosmic variance
#      f = lambda l: self.Pac.fPtotinterp(l) * self.Pbd.fPtotinterp(l) + self.Pad.fPtotinterp(l) * self.Pbc.fPtotinterp(l)
#      self.Gauss = np.array(map(f, self.L))
#      self.Gauss /= self.Npairs
#      # trispectrum
#      if self.T2d is not None:
#         self.T = self.T2d.Ttot / self.OmS
#      else:
#         self.T = np.zeros_like(self.L)
#      # supersample variance
#      if self.HSVP2d is None:
#         self.HSV = np.zeros_like(self.L)
#      else:
#         self.HSV = self.HSVP2d.P
#      # total
#      self.Total = self.Gauss + self.T + self.HSV
#
#      # covariance matrix
#      # Gaussian covariance
#      self.covMat = np.diagflat(self.Gauss)
#      return

   def LoadAll(self):
      # Gaussian covariance
      self.Gauss = self.Pac.Ptot * self.Pbd.Ptot + self.Pad.Ptot * self.Pbc.Ptot
      self.Gauss /= self.Npairs
 
      # covariance matrix
      # Gaussian covariance
      self.covMat = np.diagflat(self.Gauss)
      return
   
   
   ##################################################################################

   def snrG(self, P2d, lMax=1.e5, cv=True):
      '''Total signal to noise ratio,
      including ell bins up to lMax,
      with or without cosmic variance
      '''
      I = np.where(self.L<=lMax)[0]
      if cv:
         result = np.sqrt(np.sum((P2d.P[I]/np.sqrt(self.P[I]))**2))
      else:
         result = np.sqrt(np.sum((P2d.P[I]/np.sqrt(self.Total[I]))**2))
      return result


   def plot(self, P2d=None):
      
      if P2d is not None:
         fig=plt.figure(0)
         #
         ax=fig.add_subplot(211)
         ax.loglog(self.L, np.sqrt(self.Total)/P2d.P, 'k', lw=3, label=r'total')
         ax.loglog(self.L, np.sqrt(self.Gauss)/P2d.P, 'g-', lw=2, label=r'Gaussian')
         ax.loglog(self.L, np.sqrt(self.P)/P2d.P, 'g--', lw=2, label=r'PP')
         ax.loglog(self.L, np.sqrt(self.N)/P2d.P, 'g:', lw=2, label=r'noise')
         if hasattr(P2d, 'Ttot'):
            ax.loglog(self.L, np.sqrt(self.T)/P2d.P, 'r-', lw=2, label=r'T')
         if self.HSVP2d is not None:
            ax.loglog(self.L, np.sqrt(self.HSV)/P2d.P, 'b-', lw=2, label=r'HSV')
            ax.loglog(self.L, np.sqrt(self.HSVP2d.P1h)/P2d.P, 'b:', lw=1, label=r'HSV 1h')
            ax.loglog(self.L, np.sqrt(self.HSVP2d.P2h)/P2d.P, 'b-.', lw=1, label=r'HSV 2h')
         #
         ax.grid()
         ax.set_ylabel(r'$\sigma_P/P$')
         plt.setp(ax.get_xticklabels(), visible=False)
         #
         ax=fig.add_subplot(212)
         ax.semilogx(self.L, np.sqrt(self.Gauss/self.Total), 'g-', lw=2, label=r'Gaussian')
         ax.loglog(self.L, np.sqrt(self.P/self.Total), 'g--', lw=2, label=r'PP')
         ax.loglog(self.L, np.sqrt(self.N/self.Total), 'g:', lw=2, label=r'noise')
         if hasattr(P2d, 'Ttot'):
            ax.semilogx(self.L, np.sqrt(self.T/self.Total), 'r-', lw=2, label=r'T')
         if self.HSVP2d is not None:
            ax.semilogx(self.L, np.sqrt(self.HSV/self.Total), 'b-', lw=2, label=r'HSV')
            ax.semilogx(self.L, np.sqrt(self.HSVP2d.P1h/self.Total), 'b:', lw=1, label=r'HSV 1h')
            ax.semilogx(self.L, np.sqrt(self.HSVP2d.P2h/self.Total), 'b-.', lw=1, label=r'HSV 2h')
         #
         ax.grid()
         ax.legend(loc='center left', labelspacing=0.1, handletextpad=0.1)
         ax.set_xlabel(r'$\ell$')
         ax.set_ylabel(r'relat. contr. to error budget')
         #
         #
         #path = "./figures/covpn2d/covp_"+str(P2d.Weight)+"_fsky"+str(round(self.fsky,2))+".pdf"
         #fig.savefig(path, bbox_inches='tight')
      

      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      factor = 1. #self.L*(self.L+1.) / (2.*np.pi)
      #
      if P2d is not None:
         ax.loglog(self.L, factor * P2d.P, 'k--', lw=3, label=r'signal')
      #
      ax.loglog(self.L, factor * np.sqrt(self.Total), 'k', lw=3, label=r'total')
      ax.loglog(self.L, factor * np.sqrt(self.Gauss), 'g-', lw=2, label=r'Gauss')
      ax.loglog(self.L, factor * np.sqrt(self.P), 'g--', lw=2, label=r'PP')
      ax.loglog(self.L, factor * np.sqrt(self.N), 'g:', lw=2, label=r'noise')
#      if hasattr(self.P2d, 'Ttot'):
      ax.loglog(self.L, factor * np.sqrt(self.T), 'r-', lw=2, label=r'T')
#      if self.HSVP2d is not None:
      ax.loglog(self.L, factor * np.sqrt(self.HSV), 'b-', lw=2, label=r'HSV')
#         ax.loglog(self.L, factor * np.sqrt(self.HSVP2d.P1h), 'b:', lw=1, label=r'HSV 1h')
#         ax.loglog(self.L, factor * np.sqrt(self.HSVP2d.P2h), 'b-.', lw=1, label=r'HSV 2h')
      #
      ax.grid()
      ax.legend(loc=4)
      ax.set_xlabel(r'$\ell$')
      #ax.set_ylabel(r'$\ell(\ell+1)\, \sigma_P/(2\pi)$')
      ax.set_ylabel(r'$\sigma_P$')
      #
      #path = "./figures/covpn2d/covp_"+str(P2d.Weight)+"_absolute_fsky"+str(round(self.fsky,2))+".pdf"
      #fig.savefig(path, bbox_inches='tight')

      plt.show()










