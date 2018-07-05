import universe
reload(universe)
from universe import *

import projection_kernel
reload(projection_kernel)
from projection_kernel import *

import pn_2d
reload(pn_2d)
from pn_2d import *

import covp_2d
reload(covp_2d)
from covp_2d import *

##################################################################################
##################################################################################

class FisherLsst(object):
   
   def __init__(self, nBins=2, nL=20, fsky=1.):
      
      # number of bins
      self.nBins = nBins
      self.nG = self.nBins
      self.nS = self.nBins
      self.nGG = self.nG
      self.nGS = self.nBins*(self.nBins+1)/2
      self.nSS = self.nBins**2

      # ell bins
      self.nL = nL
      self.lMin = 20.
      self.lMax = 1.e3
      self.Le = np.logspace(np.log10(self.lMin), np.log10(self.lMax), self.nL+1, 10.)
      self.dL = self.Le[1:] - self.Le[:-1]
      # use the average ell in the bin as bin center
      f = lambda i: 2./3. * (self.Le[i+1]**3 - self.Le[i]**3) / (self.Le[i+1]**2 - self.Le[i]**2)
      self.L = np.array(map(f, range(self.nL)))
      
      # size of data vector
      self.nData = (self.nGG + self.nGS + self.nSS) * self.nL
      print "data vector has "+str(self.nData)+" elements"

      # define base cosmology
      self.u = UnivPlanck15()

      # define the tracer/shear bins
      self.w_g, self.w_s = self.generateBins(self.u)
      
      # compute the 2-point functions
      self.p2d_gg, self.p2d_gs, self.p2d_ss = self.generatePowerSpectra(self.u, self.w_g, self.w_s)
      
      # generate data vector
      self.dataVector = self.generateDataVector(self.p2d_gg, self.p2d_gs, self.p2d_ss)
      
      # generate derivatives of the data vector
      
      # generate covariance matrix
      cov = self.generateCov(self.p2d_gg, self.p2d_gs, self.p2d_ss)



   def generateBins(self, u):
      # LSST souce sample
      w_glsst = WeightTracerLSSTSources(u, name='glsst')
      # split it into bins
      zBounds = w_glsst.splitBins(self.nBins)
      # generate the corresponding tracer and shear bins
      print "Generate tracer and shear bins"
      w_g = {}
      w_s = {}
      for iBin in range(self.nBins):
         zMin = zBounds[iBin]
         zMax = zBounds[iBin+1]
         w_g[iBin] = WeightTracerLSSTSources(u, zMin=zMin, zMax=zMax, ngal=w_glsst.ngal_per_arcmin2/self.nBins, name='g'+str(iBin))
         w_s[iBin] = WeightLensCustom(u, w_glsst.dndz, zMinG=zMin, zMaxG=zMax, name='s'+str(iBin))
         print "- done "+str(iBin+1)+" of "+str(self.nBins)
      print "total ngal="+str(np.sum([w_g[i].ngal_per_arcmin2 for i in range(self.nBins)]))+"/arcmin2, should be "+str(w_glsst.ngal_per_arcmin2)
      return w_g, w_s


   def generatePowerSpectra(self, u, w_g, w_s):
      # gg
      print "Compute p_gg"
      tStart = time()
      p2d_gg = {}
      for iBin in range(self.nBins):
         p2d_gg[iBin] = P2d(u, u, w_g[iBin], fPnoise=lambda l:1./w_g[iBin].ngal, doT=False, name='', nProc=1, save=True)
      tStop = time()
      print "took "+str(tStop-tStart)+" sec"

      # gs
      print "Compute p_gs"
      tStart = time()
      p2d_gs = {}
      for iBin1 in range(self.nBins):
         # shear bin should be higher z than galaxy bin
         for iBin2 in range(iBin1,self.nBins):
            p2d_gs[iBin1, iBin2] = P2d(u, u, w_g[iBin1], w_s[iBin2], doT=False, name='', nProc=1, save=True)
      tStop = time()
      print "took "+str(tStop-tStart)+" sec"
      
      # ss
      print "Compute p_ss"
      tStart = time()
      p2d_ss = {}
      for iBin1 in range(self.nBins):
         # cross correlation: different shear bins
         for iBin2 in range(iBin1):
            p2d_ss[iBin1, iBin2] = P2d(u, u, w_s[iBin1], w_s[iBin2], doT=False, name='', nProc=1, save=True)
         # auto-correlation: same shear bin
         p2d_ss[iBin1, iBin1] = P2d(u, u, w_s[iBin1], fPnoise=lambda l:0.26**2/w_s[iBin1].ngal, doT=False, nProc=1, save=True)
      tStop = time()
      print "took "+str(tStop-tStart)+" sec"

      return p2d_gg, p2d_gs, p2d_ss






   def generateDataVector(self, p2d_gg, p2d_ss, p2d_gs):
      dataVector = np.zeros(self.nData)
      iData = 0
      # gg
      for iBin in range(self.nBins):
         dataVector[iData*self.nL:(iData+1)*self.nL] = np.array(map(p2d_gg[iBin].fPinterp, self.L))
         iData += 1
      # gs
      for iBin1 in range(self.nBins):
         # shear bin should be higher z than galaxy bin
         for iBin2 in range(iBin1,self.nBins):
            dataVector[iData*self.nL:(iData+1)*self.nL] = np.array(map(p2d_gs[iBin1, iBin2].fPinterp, self.L))
            iData += 1
      # ss
      for iBin1 in range(self.nBins):
         for iBin2 in range(iBin1+1):
            dataVector[iData*self.nL:(iData+1)*self.nL] = np.array(map(p2d_ss[iBin1, iBin2].fPinterp, self.L))
            iData += 1
      return dataVector


   def generateCov(self, p2d_gg, p2d_ss, p2d_gs):

#cov = CovP2d(p2d_gg, p2d_gg, p2d_gg, p2d_gg, T2d=None, HSVP2d=None, fsky=1., npairs_id='log')
#covMat = cov.covMat()

      # gg-gg


      # gg-gs
      # gg-ss
      # gs-gs
      # gs-ss
      # ss-ss























