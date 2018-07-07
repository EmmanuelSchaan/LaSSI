import universe
reload(universe)
from universe import *

import parameters_fisher
reload(parameters_fisher)
from parameters_fisher import *

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

class FisherLsst(object):
   
   def __init__(self, cosmoPar, photoZPar, shearMultBiasPar, nBins=2, nL=20, fsky=1., name="", save=True):
      
      self.name = ""
      self.save = save
      
      # sky fraction
      self.fsky = fsky
      
      # number of bins
      self.nBins = nBins
      self.nG = self.nBins
      self.nS = self.nBins
      self.nGG = self.nG * (self.nG+1) / 2   # not just gg from same z-bin
      self.nGS = self.nG * self.nS # not just higher z s than g
      self.nSS = self.nS * (self.nS+1) / 2

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
      print "Data vector has "+str(self.nData)+" elements"

      # define cosmology parameters
      self.cosmoPar = cosmoPar#CosmoParams(massiveNu=False, wCDM=False, curvature=False)
      # photo-z and shear bias parameters
      self.photoZPar = photoZPar#PhotoZParams(nBins=self.nBins)
      self.shearMultBiasPar = shearMultBiasPar#ShearMultBiasParams(nBins=self.nBins)
      
      ##################################################################################

      # define base cosmology
      self.u = UnivFisher(self.cosmoPar)

      # define tracer/shear bins
      self.w_g, self.w_s, self.zBounds = self.generateBins(self.u, self.photoZPar, self.shearMultBiasPar)
      # update the galaxy bias parameters
      self.galaxyBiasPar = GalaxyBiasParams(self.nBins, self.w_g)
      
      # compute the 2-point functions
      self.p2d_gg, self.p2d_gs, self.p2d_ss = self.generatePowerSpectra(self.u, self.w_g, self.w_s, save=self.save)
      
      # generate data vector
      self.dataVector = self.generateDataVector(self.p2d_gg, self.p2d_gs, self.p2d_ss)
      
      ##################################################################################
      
      # generate covariance matrix
      self.covMat = self.generateCov(self.p2d_gg, self.p2d_gs, self.p2d_ss)
      self.invCov = np.linalg.inv(self.covMat)
      
      ##################################################################################
      
      # define complete parameter vector
      self.params = self.cosmoPar.copy()
      self.params.addParams(self.galaxyBiasPar)
      self.params.addParams(self.shearMultBiasPar)
      self.params.addParams(self.photoZPar)
      
      ##################################################################################
      
      # generate derivatives of the data vector:
      # matrix of size self.params.nPar x self.nData
      self.derivativeDataVector = self.generateDerivativeDataVector()
      



   ##################################################################################


   def generateBins(self, u, photoZPar, shearMultBiasPar):
      # LSST souce sample
      w_glsst = WeightTracerLSSTSources(u, name='glsst')
      # split it into bins
      zBounds = w_glsst.splitBins(self.nBins)
      
      # generate the corresponding tracer and shear bins
      print "Generate tracer and shear bins"
      tStart = time()
      w_g = {}
      w_s = {}
      for iBin in range(self.nBins):
         # sharp photo-z cuts
         zMinP = zBounds[iBin]
         zMaxP = zBounds[iBin+1]
         # true photo-z bounds
         zMin = 1./w_glsst.aMax-1.
         zMax = 1./w_glsst.aMin-1.
         # true dn/dz_true from dn/dz_phot
         p_z_given_zp = lambda zp,z: np.exp(-0.5*(z-zp-photoZPar.fiducial[iBin])**2/photoZPar.fiducial[self.nBins+iBin]**2) / np.sqrt(2.*np.pi*photoZPar.fiducial[self.nBins+iBin]**2)
         f = lambda zp,z: w_glsst.dndz(zp) * p_z_given_zp(zp,z)
         dndz_tForInterp = lambda z: integrate.quad(f, zMinP, zMaxP, args=(z), epsabs=0., epsrel=1.e-2)[0]
         # interpolate it for speed (for lensing kernel calculation)
         Z = np.linspace(zMin, zMax, 501)
         F = np.array(map(dndz_tForInterp, Z))
         dndz_t = interp1d(Z, F, kind='linear', bounds_error=False, fill_value=0.)
         
         # tracer bin
         w_g[iBin] = WeightTracerCustom(u,
                                        lambda z: w_glsst.b(z), # galaxy bias
                                        dndz_t, # dn/dz_true
                                        zMin=zMin,
                                        zMax=zMax,
                                        name='g'+str(iBin))
         
         # shear bin
         w_s[iBin] = WeightLensCustom(u,
                                      dndz_t, # dn/dz_true
                                      m=lambda z: shearMultBiasPar.fiducial[iBin], # multiplicative shear bias
                                      zMinG=zMin,
                                      zMaxG=zMax,
                                      name='s'+str(iBin))
         #print "- done "+str(iBin+1)+" of "+str(self.nBins)
      print "total ngal="+str(np.sum([w_g[i].ngal_per_arcmin2 for i in range(self.nBins)]))+"/arcmin2, should be "+str(w_glsst.ngal_per_arcmin2)
      tStop = time()
      print "took "+str(tStop-tStart)+" sec"
      return w_g, w_s, zBounds


   def generatePowerSpectra(self, u, w_g, w_s, save=True):
      # gg: do not impose same bin
      print "Compute p_gg"
      tStart = time()
      p2d_gg = np.empty((self.nBins, self.nBins), dtype=object)
      for iBin1 in range(self.nBins):
         # auto-correlation: same bin
         p2d_gg[iBin1, iBin1] = P2d(u, u, w_g[iBin1], fPnoise=lambda l:1./w_g[iBin1].ngal, doT=False, name='', L=self.L, nProc=1, save=save)
         # cross-correlation: different bins
         for iBin2 in range(iBin1+1, self.nBins):
            p2d_gg[iBin1, iBin2] = P2d(u, u, w_g[iBin1], w_g[iBin2], doT=False, name='', L=self.L, nProc=1, save=save)
            # so that the order doesn't matter
            p2d_gg[iBin2, iBin1] = p2d_gg[iBin1, iBin2]
      tStop = time()
      print "took "+str(tStop-tStart)+" sec"

      # gs: do not impose higher z s than g
      print "Compute p_gs"
      tStart = time()
      p2d_gs = np.empty((self.nBins, self.nBins), dtype=object)
      for iBin1 in range(self.nBins):
         for iBin2 in range(self.nBins):
            p2d_gs[iBin1, iBin2] = P2d(u, u, w_g[iBin1], w_s[iBin2], doT=False, name='', L=self.L, nProc=1, save=save)
      tStop = time()
      print "took "+str(tStop-tStart)+" sec"
      
      # ss
      print "Compute p_ss"
      tStart = time()
      p2d_ss = np.empty((self.nBins, self.nBins), dtype=object)
      for iBin1 in range(self.nBins):
         # auto-correlation: same bin
         p2d_ss[iBin1, iBin1] = P2d(u, u, w_s[iBin1], fPnoise=lambda l:0.26**2/w_s[iBin1].ngal, doT=False, L=self.L, nProc=1, save=save)
         # cross correlation: different bins
         for iBin2 in range(iBin1+1, self.nBins):
            p2d_ss[iBin1, iBin2] = P2d(u, u, w_s[iBin1], w_s[iBin2], doT=False, name='', L=self.L, nProc=1, save=save)
            # so that the order doesn't matter
            p2d_ss[iBin2, iBin1] = p2d_ss[iBin1, iBin2]
      tStop = time()
      print "took "+str(tStop-tStart)+" sec"

      return p2d_gg, p2d_gs, p2d_ss


   def generateDataVector(self, p2d_gg, p2d_gs, p2d_ss):
      dataVector = np.zeros(self.nData)
      print "Generating data vector"
      tStart = time()
      iData = 0
      # gg
      for iBin1 in range(self.nBins):
         for iBin2 in range(iBin1, self.nBins):
            dataVector[iData*self.nL:(iData+1)*self.nL] = np.array(map(p2d_gg[iBin1, iBin2].fPinterp, self.L))
            iData += 1
      # gs
      for iBin1 in range(self.nBins):
         for iBin2 in range(self.nBins):
            dataVector[iData*self.nL:(iData+1)*self.nL] = np.array(map(p2d_gs[iBin1, iBin2].fPinterp, self.L))
            iData += 1
      # ss
      for iBin1 in range(self.nBins):
         for iBin2 in range(iBin1, self.nBins):
            dataVector[iData*self.nL:(iData+1)*self.nL] = np.array(map(p2d_ss[iBin1, iBin2].fPinterp, self.L))
            iData += 1

      tStop = time()
      print "took "+str(tStop-tStart)+" sec"
      return dataVector


   def generateCov(self, p2d_gg, p2d_gs, p2d_ss):
      covMat = np.zeros((self.nData, self.nData))
      tStart = time()
      print "Computing covariance matrix"
      # below, i1 and i2 define the row and column of the nL*nL blocks for each pair of 2-point function
      # i1, i2 \in [0, nGG+nGS+nSS]
      
      print "gg-gg"
      # considering gg[i1,i2]
      i1 = 0
      for iBin1 in range(self.nBins):
         for iBin2 in range(iBin1, self.nBins):
            # considering gg[j1,j2]
            i2 = 0
            for jBin1 in range(self.nBins):
               for jBin2 in range(jBin1, self.nBins):
                  # compute only upper diagonal
                  if i2>=i1:
                     covBlock = CovP2d(p2d_gg[iBin1,jBin1], p2d_gg[iBin2,jBin2], p2d_gg[iBin1,jBin2], p2d_gg[iBin2,jBin1], fsky=self.fsky, L=self.L, dL=self.dL)
                     covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = covBlock.covMat
                  # move to next column
                  i2 += 1
            # move to next row
            i1 += 1

      print "gg-gs"
      # considering gg[i1,i2]
      i1 = 0
      for iBin1 in range(self.nBins):
         for iBin2 in range(iBin1, self.nBins):
            # considering gs[j1,j2]
            i2 = self.nGG
            for jBin1 in range(self.nBins):
               for jBin2 in range(self.nBins):
                  # compute only upper diagonal
                  if i2>=i1:
                     covBlock = CovP2d(p2d_gg[iBin1,jBin1], p2d_gs[iBin2,jBin2], p2d_gs[iBin1,jBin2], p2d_gg[iBin2,jBin1], fsky=self.fsky, L=self.L, dL=self.dL)
                     covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = covBlock.covMat
                  # move to next column
                  i2 += 1
            # move to next row
            i1 += 1

      print "gg-ss"
      # considering gg[i1,i2]
      i1 = 0
      for iBin1 in range(self.nBins):
         for iBin2 in range(iBin1, self.nBins):
            # considering ss[j1,j2]
            i2 = self.nGG + self.nGS
            for jBin1 in range(self.nBins):
               for jBin2 in range(jBin1, self.nBins):
                  # compute only upper diagonal
                  if i2>=i1:
                     covBlock = CovP2d(p2d_gs[iBin1,jBin1], p2d_gs[iBin2,jBin2], p2d_gs[iBin1,jBin2], p2d_gs[iBin2,jBin1], fsky=self.fsky, L=self.L, dL=self.dL)
                     covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = covBlock.covMat
                  # move to next column
                  i2 += 1
            # move to next row
            i1 += 1

      print "gs-gs"
      # considering gs[i1,i2]
      i1 = self.nGG
      for iBin1 in range(self.nBins):
         for iBin2 in range(self.nBins):
            # considering gs[j1,j2]
            i2 = self.nGG
            for jBin1 in range(self.nBins):
               for jBin2 in range(self.nBins):
                  # compute only upper diagonal
                  if i2>=i1:
                     # watch the order for gs
                     covBlock = CovP2d(p2d_gg[iBin1,jBin1], p2d_ss[iBin2,jBin2], p2d_gs[iBin1,jBin2], p2d_gs[jBin1,iBin2], fsky=self.fsky, L=self.L, dL=self.dL)
                     covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = covBlock.covMat
                  # move to next column
                  i2 += 1
            # move to next row
            i1 += 1

      print "gs-ss"
      # considering gs[i1,i2]
      i1 = self.nGG
      for iBin1 in range(self.nBins):
         for iBin2 in range(self.nBins):
            # considering ss[j1,j2]
            i2 = self.nGG + self.nGS
            for jBin1 in range(self.nBins):
               for jBin2 in range(jBin1, self.nBins):
                  # compute only upper diagonal
                  if i2>=i1:
                     covBlock = CovP2d(p2d_gs[iBin1,jBin1], p2d_ss[iBin2,jBin2], p2d_gs[iBin1,jBin2], p2d_ss[iBin2,jBin1], fsky=self.fsky, L=self.L, dL=self.dL)
                     covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = covBlock.covMat
                  # move to next column
                  i2 += 1
            # move to next row
            i1 += 1

      print "ss-ss"
      # considering ss[i1,i2]
      i1 = self.nGG + self.nGS
      for iBin1 in range(self.nBins):
         for iBin2 in range(iBin1, self.nBins):
            # considering ss[j1,j2]
            i2 = self.nGG + self.nGS
            for jBin1 in range(self.nBins):
               for jBin2 in range(jBin1, self.nBins):
                  # compute only upper diagonal
                  if i2>=i1:
                     covBlock = CovP2d(p2d_ss[iBin1,jBin1], p2d_ss[iBin2,jBin2], p2d_ss[iBin1,jBin2], p2d_ss[iBin2,jBin1], fsky=self.fsky, L=self.L, dL=self.dL)
                     covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = covBlock.covMat
                  # move to next column
                  i2 += 1
            # move to next row
            i1 += 1

      # fill lower diagonal by symmetry
      # here i1 and i2 don't index the matrix blocks, but the matrix elements
      for i1 in range(self.nData):
         for i2 in range(i1):
            covMat[i1, i2] = covMat[i2, i1]

   
      tStop = time()
      print "took "+str(tStop-tStart)+" sec"
      return covMat



   def generateDerivativeDataVector(self):
      # Derivatives of the data vector:
      # matrix of size self.params.nPar x self.nData
      result = np.zeros((self.params.nPar, self.nData))
   
      # derivatives wrt cosmology
      
      # derivatives wrt galaxy bias:
      # divide by bMean[iBin] every 2-pt function involving the tracer bin iBin
      for iBin in range(self.nBins):
         # derivative for gg
         p2d_gg = self.p2d_gg.copy()
         
         
         # derivative for gs
         p2d_gs = self.p2d_gs.copy()
#         p2d_gs[iBin,:] /= self.galaxyBiasPar.fiducial[iBin]

   
#      # compute the 2-point functions
#      self.p2d_gg, self.p2d_gs, self.p2d_ss = self.generatePowerSpectra(self.u, self.w_g, self.w_s, save=self.save)
#
#      # generate data vector
#      self.dataVector = self.generateDataVector(self.p2d_gg, self.p2d_gs, self.p2d_ss)

   
   
   
   
   
   
   
   
   
   
      # derivatives wrt shear multiplicative bias
      
      # derivatives wrt photo-z parameters
      
      

   ##################################################################################
   
   def plotDndz(self):
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for iBin in range(self.nBins):
         # define z range of bin
         zMin = 1./self.w_g[iBin].aMax-1.
         zMax = 1./self.w_g[iBin].aMin-1.
         Z = np.linspace(zMin, zMax, 501)
         # evaluate dn/dz
         dndz = np.array(map(self.w_g[iBin].dndz, Z))
         dndz /= (180.*60./np.pi)**2 # convert from 1/sr to 1/arcmin^2
         # plot it
         ax.fill_between(Z, 0., dndz, facecolor=plt.cm.autumn(1.*iBin/self.nBins), edgecolor='', alpha=0.5)
      #
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$dN / d\Omega\; dz$ [arcmin$^{-2}$]')
      
      plt.show()
   
   
   def plotCovMat(self):
      fig=plt.figure(0, figsize=(12,8))
      ax=fig.add_subplot(111)
      #
      # compute correlation matrix
      corMat = np.zeros_like(self.covMat)
      for i in range(self.nData):
         for j in range(self.nData):
            corMat[i,j] = self.covMat[i,j] / np.sqrt(self.covMat[i,i] * self.covMat[j,j])
      upperDiag = np.triu(np.ones(self.nData))
      plt.imshow(corMat * upperDiag, interpolation='nearest', norm=LogNorm(vmin=1.e-4, vmax=1), cmap=cmaps.viridis_r)
      #
      ax.plot(np.arange(self.nData+1)-0.5, np.arange(self.nData+1)-0.5, 'k', lw=1)
      #
      # block delimiters
      ax.axhline(self.nL*self.nGG-0.5, xmin=(self.nL*self.nGG-0.5)/self.nData, c='k', lw=1.5)
      ax.axhline(self.nL*(self.nGG+self.nGS)-0.5, xmin=(self.nL*(self.nGG+self.nGS)-0.5)/self.nData, c='k', lw=1.5)
      #
      ax.axvline(self.nL*self.nGG-0.5, ymin=1.-(self.nL*self.nGG-0.5)/self.nData, c='k', lw=1.5)
      ax.axvline(self.nL*(self.nGG+self.nGS)-0.5, ymin=1.-(self.nL*(self.nGG+self.nGS)-0.5)/self.nData, c='k', lw=1.5)
      #
      # 2-pt function delimiters
      for i in range(1, self.nGG+self.nGS+self.nSS):
         ax.axhline(self.nL*i-0.5, xmin=(self.nL*i-0.5)/self.nData, c='gray', lw=0.5, ls='--')
         ax.axvline(self.nL*i-0.5, ymin=1.-(self.nL*i-0.5)/self.nData, c='gray', lw=0.5, ls='--')
      #
      plt.colorbar()
      ax.set_xlim((-0.5, (self.nData-1)+0.5))
      ax.set_ylim((-0.5, (self.nData-1)+0.5))
      ax.invert_yaxis()
      #ax.xaxis.tick_top()
      ax.xaxis.set_ticks([])
      ax.yaxis.set_ticks([])
      #ax.grid(True)
      #ax.set_title(r'Full cor: '+infile)
      #
      #fig.savefig(pathFig+"cor_small.pdf", bbox_inches='tight', format='pdf', dpi=1000)
      #fig.clf()

      plt.show()


   def plotPowerSpectra(self):

      # gg
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      i1 = 0
      for iBin1 in range(self.nBins):
         for iBin2 in range(iBin1, self.nBins):
            d = self.dataVector[i1*self.nL:(i1+1)*self.nL]
            std = np.sqrt(np.diag(self.covMat[i1*self.nL:(i1+1)*self.nL, i1*self.nL:(i1+1)*self.nL]))
            ax.errorbar(self.L, d, yerr=std, ls='-', lw=1, elinewidth=1.5, marker='.', markersize=2, label=r'$g_{'+str(iBin1)+'} g_{'+str(iBin2)+'}$')
            # move to next row
            i1 += 1
      #
      ax.legend(loc=1)
      ax.set_xscale('log')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$\ell$')
      ax.set_ylabel(r'$C_\ell^{gg}$')
      
      plt.show()


      # gs
      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      i1 = self.nGG
      for iBin1 in range(self.nBins):
         for iBin2 in range(self.nBins):
            d = self.dataVector[i1*self.nL:(i1+1)*self.nL]
            std = np.sqrt(np.diag(self.covMat[i1*self.nL:(i1+1)*self.nL, i1*self.nL:(i1+1)*self.nL]))
            ax.errorbar(self.L, d, yerr=std, ls='-', lw=1, elinewidth=1.5, marker='.', markersize=2, label=r'$g_{'+str(iBin1)+'} \gamma_{'+str(iBin2)+'}$')
            # move to next row
            i1 += 1
      #
      ax.legend(loc=1)
      ax.set_xscale('log')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$\ell$')
      ax.set_ylabel(r'$C_\ell^{g\gamma}$')

      plt.show()


      # ss
      fig=plt.figure(2)
      ax=fig.add_subplot(111)
      #
      i1 = self.nGG + self.nGS
      for iBin1 in range(self.nBins):
         for iBin2 in range(iBin1, self.nBins):
            d = self.dataVector[i1*self.nL:(i1+1)*self.nL]
            std = np.sqrt(np.diag(self.covMat[i1*self.nL:(i1+1)*self.nL, i1*self.nL:(i1+1)*self.nL]))
            ax.errorbar(self.L, d, yerr=std, ls='-', lw=1, elinewidth=1.5, marker='.', markersize=2, label=r'$\gamma_{'+str(iBin1)+'} \gamma_{'+str(iBin2)+'}$')
            # move to next row
            i1 += 1
      #
      ax.legend(loc=1)
      ax.set_xscale('log')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$\ell$')
      ax.set_ylabel(r'$C_\ell^{\gamma\gamma}$')

      plt.show()









