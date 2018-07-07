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
   
   def __init__(self, cosmoPar, galaxyBiasPar, shearMultBiasPar, photoZPar, nBins=2, nL=20, fsky=1., name=None, save=True):
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
      print "Tomographic bins: "+str(self.nBins)

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
      
      if name is None:
         self.name = "lsst_gg_gs_ss_nBins"+str(self.nBins)+"_nL"+str(self.nL)
      else:
         self.name = name

      ##################################################################################

      # cosmology parameters
      self.cosmoPar = cosmoPar
      
      # nuisance parameters
      self.galaxyBiasPar = galaxyBiasPar
      self.shearMultBiasPar = shearMultBiasPar
      self.photoZPar = photoZPar
      # combined nuisance parameters
      self.nuisancePar = self.galaxyBiasPar.copy()
      self.nuisancePar.addParams(self.shearMultBiasPar)
      self.nuisancePar.addParams(self.photoZPar)
      
      # all parameters
      self.fullPar = self.cosmoPar.copy()
      self.fullPar.addParams(self.nuisancePar)
      print "Params: "+str(self.fullPar.nPar)+" total = "+str(self.cosmoPar.nPar)+" cosmo + "+str(self.nuisancePar.nPar)+" nuisance"

      ##################################################################################
      # Fiducial data vector and covariance
      
      tStartFisher = time()
      
      print "Run CLASS",
      tStart = time()
      self.u = Universe(self.cosmoPar.paramsClassy)
      tStop = time()
      print "("+str(np.round(tStop-tStart,1))+" sec)"

      print "Tracer and shear bins",
      tStart = time()
      self.w_g, self.w_s, self.zBounds = self.generateBins(self.u, self.nuisancePar.fiducial)
      tStop = time()
      print "("+str(np.round(tStop-tStart,1))+" sec)"
      
      print "Power spectra",
      tStart = time()
      self.p2d_gg, self.p2d_gs, self.p2d_ss = self.generatePowerSpectra(self.u, self.w_g, self.w_s, save=self.save)
      tStop = time()
      print "("+str(np.round(tStop-tStart,1))+" sec)"

      print "Data vector",
      tStart = time()
      self.dataVector = self.generateDataVector(self.p2d_gg, self.p2d_gs, self.p2d_ss)
      tStop = time()
      print "("+str(np.round(tStop-tStart,1))+" sec)"

      print "Covariance matrix",
      tStart = time()
      self.covMat = self.generateCov(self.p2d_gg, self.p2d_gs, self.p2d_ss)
      self.invCov = np.linalg.inv(self.covMat)
      tStop = time()
      print "("+str(np.round(tStop-tStart,1))+" sec)"
      
      ##################################################################################
      
      # Derivatives of the data vector
      if self.save:
         self.saveDerivativeDataVector()
      self.loadDerivativeDataVector()
      
      tStopFisher = time()
      print "Full calculation took "+str(np.round((tStopFisher-tStartFisher)/60.,1))+" min"

   ##################################################################################


   def generateBins(self, u, nuisancePar):
      # split the nuisance parameters
      galaxyBiasPar = nuisancePar[:self.galaxyBiasPar.nPar]
      shearMultBiasPar = nuisancePar[self.galaxyBiasPar.nPar:self.galaxyBiasPar.nPar+self.shearMultBiasPar.nPar]
      photoZPar = nuisancePar[self.galaxyBiasPar.nPar+self.shearMultBiasPar.nPar:]
      # LSST source sample
      w_glsst = WeightTracerLSSTSources(u, name='glsst')
      # split it into bins
      zBounds = w_glsst.splitBins(self.nBins)
      
      # generate the corresponding tracer and shear bins

      w_g = np.empty(self.nBins, dtype=object)
      w_s = np.empty(self.nBins, dtype=object)
      for iBin in range(self.nBins):
         # sharp photo-z cuts
         zMinP = zBounds[iBin]
         zMaxP = zBounds[iBin+1]
         # photo-z bias and uncertainty for this bin
#!!! here the param is just sz, not sz/(1+z). I may want to change this.
         dz = photoZPar[iBin]
         sz = photoZPar[self.nBins+iBin]
         # true z bounds: truncate at 3 sigma
         # careful for the first and last bin
         zMin = max(zMinP - 3.*sz, 1./w_glsst.aMax-1.)   # 1./w_glsst.aMax-1.
         zMax = min(zMaxP + 3.*sz, 1./w_glsst.aMin-1.)   # 1./w_glsst.aMin-1.
         # true dn/dz_true from dn/dz_phot

         p_z_given_zp = lambda zp,z: np.exp(-0.5*(z-zp-dz)**2/sz**2) / np.sqrt(2.*np.pi*sz**2)
         f = lambda zp,z: w_glsst.dndz(zp) * p_z_given_zp(zp,z)
         dndz_tForInterp = lambda z: integrate.quad(f, zMinP, zMaxP, args=(z), epsabs=0., epsrel=1.e-2)[0]
         # interpolate it for speed (for lensing kernel calculation)
         Z = np.linspace(zMin, zMax, 101)
         F = np.array(map(dndz_tForInterp, Z))
         dndz_t = interp1d(Z, F, kind='linear', bounds_error=False, fill_value=0.)
         
         # tracer bin
         w_g[iBin] = WeightTracerCustom(u,
                                        lambda z: galaxyBiasPar[iBin] * w_glsst.b(z), # galaxy bias
                                        dndz_t, # dn/dz_true
                                        zMin=zMin,
                                        zMax=zMax,
                                        name='g'+str(iBin))
         
         # shear bin
         w_s[iBin] = WeightLensCustom(u,
                                      dndz_t, # dn/dz_true
                                      m=lambda z: shearMultBiasPar[iBin], # multiplicative shear bias
                                      zMinG=zMin,
                                      zMaxG=zMax,
                                      name='s'+str(iBin))

         #print "- done "+str(iBin+1)+" of "+str(self.nBins)
      #print "total ngal="+str(np.sum([w_g[i].ngal_per_arcmin2 for i in range(self.nBins)]))+"/arcmin2, should be "+str(w_glsst.ngal_per_arcmin2)
      return w_g, w_s, zBounds


   ##################################################################################

   def generatePowerSpectra(self, u, w_g, w_s, name=None, save=True):
      if name is None:
         name = "_"+self.name
      # gg: do not impose same bin
      p2d_gg = np.empty((self.nBins, self.nBins), dtype=object)
      for iBin1 in range(self.nBins):
         # auto-correlation: same bin
         p2d_gg[iBin1, iBin1] = P2d(u, u, w_g[iBin1], fPnoise=lambda l:1./w_g[iBin1].ngal, doT=False, name=name, L=self.L, nProc=1, save=save)
         # cross-correlation: different bins
         for iBin2 in range(iBin1+1, self.nBins):
            p2d_gg[iBin1, iBin2] = P2d(u, u, w_g[iBin1], w_g[iBin2], doT=False, name=name, L=self.L, nProc=1, save=save)
            # so that the order doesn't matter
            p2d_gg[iBin2, iBin1] = p2d_gg[iBin1, iBin2]

      # gs: do not impose higher z s than g
      p2d_gs = np.empty((self.nBins, self.nBins), dtype=object)
      for iBin1 in range(self.nBins):
         for iBin2 in range(self.nBins):
            p2d_gs[iBin1, iBin2] = P2d(u, u, w_g[iBin1], w_s[iBin2], doT=False, name=name, L=self.L, nProc=1, save=save)
      
      # ss
      p2d_ss = np.empty((self.nBins, self.nBins), dtype=object)
      for iBin1 in range(self.nBins):
         # auto-correlation: same bin
         p2d_ss[iBin1, iBin1] = P2d(u, u, w_s[iBin1], fPnoise=lambda l:0.26**2/w_s[iBin1].ngal, doT=False, name=name, L=self.L, nProc=1, save=save)
         # cross correlation: different bins
         for iBin2 in range(iBin1+1, self.nBins):
            p2d_ss[iBin1, iBin2] = P2d(u, u, w_s[iBin1], w_s[iBin2], doT=False, name=name, L=self.L, nProc=1, save=save)
            # so that the order doesn't matter
            p2d_ss[iBin2, iBin1] = p2d_ss[iBin1, iBin2]

      return p2d_gg, p2d_gs, p2d_ss


   ##################################################################################

   def generateDataVector(self, p2d_gg, p2d_gs, p2d_ss):
      dataVector = np.zeros(self.nData)
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
      return dataVector


   ##################################################################################

   def generateCov(self, p2d_gg, p2d_gs, p2d_ss):
      covMat = np.zeros((self.nData, self.nData))
      # below, i1 and i2 define the row and column of the nL*nL blocks for each pair of 2-point function
      # i1, i2 \in [0, nGG+nGS+nSS]
      
      #print "gg-gg"
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

      #print "gg-gs"
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

      #print "gg-ss"
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

      #print "gs-gs"
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

      #print "gs-ss"
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

      #print "ss-ss"
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

      return covMat


   ##################################################################################

   def saveDerivativeDataVector(self):
      # Derivatives of the data vector:
      # matrix of size self.params.nPar x self.nData
      derivative = np.zeros((self.fullPar.nPar, self.nData))
      
      for iPar in range(self.cosmoPar.nPar):
         print "Derivative wrt "+self.cosmoPar.names[iPar],
         tStart = time()
         # high
         name = self.name+self.cosmoPar.names[iPar]+"high"
         cosmoParClassy = self.cosmoPar.paramsClassy.copy()
         print cosmoParClassy
         print "#"
         cosmoParClassy[self.cosmoPar.names[iPar]] = self.cosmoPar.paramsClassyHigh[self.cosmoPar.names[iPar]]
         print cosmoParClassy
         u = Universe(cosmoParClassy)
         w_g, w_s, zBounds = self.generateBins(u, self.nuisancePar.fiducial)
         p2d_gg, p2d_gs, p2d_ss = self.generatePowerSpectra(u, w_g, w_s, name=name, save=True)
         dataVectorHigh = self.generateDataVector(p2d_gg, p2d_gs, p2d_ss)
#         # low
#         name = self.name+self.cosmoPar.names[iPar]+"low"
#         cosmoParClassy = self.cosmoPar.paramsClassy.copy()
#         cosmoParClassy[self.cosmoPar.names[iPar]] = self.cosmoPar.paramsClassyLow[self.cosmoPar.names[iPar]]
#         u = Universe(cosmoParClassy)
#         w_g, w_s, zBounds = self.generateBins(u, self.nuisancePar.fiducial)
#         p2d_gg, p2d_gs, p2d_ss = self.generatePowerSpectra(u, w_g, w_s, name=name, save=True)
#         dataVectorHigh = self.generateDataVector(p2d_gg, p2d_gs, p2d_ss)
#         # derivative
#         derivative[iPar,:] = (dataVectorHigh-dataVectorLow) / (self.cosmoPar.high[iPar]-self.cosmoPar.low[iPar])
         derivative[iPar,:] = (dataVectorHigh-self.dataVector) / (self.cosmoPar.high[iPar]-self.cosmoPar.fiducial[iPar])
         tStop = time()
         print "("+str(np.round(tStop-tStart,1))+" sec)"
      
      # Nuisance parameters
      for iPar in range(self.nuisancePar.nPar):
         print "Derivative wrt "+self.nuisancePar.names[iPar],
         tStart = time()
         params = self.nuisancePar.fiducial.copy()
         # high
         name = "_"+self.name+self.nuisancePar.names[iPar]+"high"
         params[iPar] = self.nuisancePar.high[iPar]
         w_g, w_s, zBounds = self.generateBins(self.u, params)
         p2d_gg, p2d_gs, p2d_ss = self.generatePowerSpectra(self.u, w_g, w_s, name=name, save=True)
         dataVectorHigh = self.generateDataVector(p2d_gg, p2d_gs, p2d_ss)
#         # low
#         name = self.name+self.nuisancePar.names[iPar]+"low"
#         params[iPar] = self.nuisancePar.low[iPar]
#         w_g, w_s, zBounds = self.generateBins(self.u, params)
#         p2d_gg, p2d_gs, p2d_ss = self.generatePowerSpectra(self.u, w_g, w_s, name=name, save=True)
#         dataVectorLow = self.generateDataVector(p2d_gg, p2d_gs, p2d_ss)
#         # derivative
#         derivative[self.cosmoPar.nPar+iPar,:] = (dataVectorHigh-dataVectorLow) / (self.nuisancePar.high[iPar]-self.nuisancePar.low[iPar])
         derivative[self.cosmoPar.nPar+iPar,:] = (dataVectorHigh-self.dataVector) / (self.nuisancePar.high[iPar]-self.nuisancePar.fiducial[iPar])
         tStop = time()
         print "("+str(np.round(tStop-tStart,1))+" sec)"
      
      path = "dDatadPar_"+self.name
      np.savetxt(path, derivative)

   def loadDerivativeDataVector(self):
      path = "dDatadPar_"+self.name
      self.derivativeDataVector = np.genfromtxt(path)


   ##################################################################################
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
            ax.errorbar(self.L*(1.+0.01*i1/self.nGG), d, yerr=std, ls='-', lw=1, elinewidth=1.5, marker='.', markersize=2, label=r'$g_{'+str(iBin1)+'} g_{'+str(iBin2)+'}$')
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
            ax.errorbar(self.L*(1.+0.01*i1/self.nGS), d, yerr=std, ls='-', lw=1, elinewidth=1.5, marker='.', markersize=2, label=r'$g_{'+str(iBin1)+'} \gamma_{'+str(iBin2)+'}$')
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
            ax.errorbar(self.L*(1.+0.01*i1/self.nSS), d, yerr=std, ls='-', lw=1, elinewidth=1.5, marker='.', markersize=2, label=r'$\gamma_{'+str(iBin1)+'} \gamma_{'+str(iBin2)+'}$')
            # move to next row
            i1 += 1
      #
      ax.legend(loc=1)
      ax.set_xscale('log')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$\ell$')
      ax.set_ylabel(r'$C_\ell^{\gamma\gamma}$')

      plt.show()


   def plotDerivativeDataVector(self):
      # one color per cosmo param
      Colors = plt.cm.jet(1.*np.arange(self.cosmoPar.nPar)/self.cosmoPar.nPar)
      
      # gg
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # for each cosmo parameter
      for iPar in range(self.cosmoPar.nPar):
         # plot all the 2pt functions
         for i2pt in range(self.nGG):
            dlnDdlnP = self.derivativeDataVector[iPar, i2pt*self.nL:(i2pt+1)*self.nL] * self.cosmoPar.fiducial[iPar] / self.dataVector[i2pt*self.nL:(i2pt+1)*self.nL]
            color = Colors[iPar]
            color = darkerLighter(color, amount=-0.8*i2pt/self.nGG)
            ax.plot(self.L, dlnDdlnP, c=color, lw=3)
         ax.plot([],[], c=Colors[iPar], label=self.cosmoPar.namesLatex[iPar])
      #
      ax.grid()
      ax.legend(loc=4, ncol=4, labelspacing=0.05, frameon=False, handlelength=0.4, borderaxespad=0.01)
      ax.set_xscale('log', nonposx='clip')
      ax.set_ylim((-3., 4.))
      ax.set_xlabel(r'$\ell$')
      ax.set_ylabel(r'$d\ln C_\ell^{gg} / d\ln \text{Param.}$')


      # gs
      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      # for each cosmo parameter
      for iPar in range(self.cosmoPar.nPar):
         # plot all the 2pt functions
         for i2pt in range(self.nGG, self.nGG+self.nGS):
            dlnDdlnP = self.derivativeDataVector[iPar, i2pt*self.nL:(i2pt+1)*self.nL] * self.cosmoPar.fiducial[iPar] / self.dataVector[i2pt*self.nL:(i2pt+1)*self.nL]
            color = Colors[iPar]
            color = darkerLighter(color, amount=-0.8*(i2pt-self.nGG)/self.nGS)
            ax.plot(self.L, dlnDdlnP, c=color, lw=3)
         ax.plot([],[], c=Colors[iPar], label=self.cosmoPar.namesLatex[iPar])
      #
      ax.grid()
      ax.legend(loc=4, ncol=4, labelspacing=0.05, frameon=False, handlelength=0.4, borderaxespad=0.01)
      ax.set_xscale('log', nonposx='clip')
      ax.set_ylim((-3., 4.))
      ax.set_xlabel(r'$\ell$')
      ax.set_ylabel(r'$d\ln C_\ell^{gs} / d\ln \text{Param.}$')

      # ss
      fig=plt.figure(2)
      ax=fig.add_subplot(111)
      #
      # for each cosmo parameter
      for iPar in range(self.cosmoPar.nPar):
         # plot all the 2pt functions
         for i2pt in range(self.nGG+self.nGS, self.nGG+self.nGS+self.nSS):
            dlnDdlnP = self.derivativeDataVector[iPar, i2pt*self.nL:(i2pt+1)*self.nL] * self.cosmoPar.fiducial[iPar] / self.dataVector[i2pt*self.nL:(i2pt+1)*self.nL]
            color = Colors[iPar]
            color = darkerLighter(color, amount=-0.8*(i2pt-(self.nGG+self.nGS))/self.nSS)
            ax.plot(self.L, dlnDdlnP, c=color, lw=3)
         ax.plot([],[], c=Colors[iPar], label=self.cosmoPar.namesLatex[iPar])
      #
      ax.grid()
      ax.legend(loc=4, ncol=4, labelspacing=0.05, frameon=False, handlelength=0.4, borderaxespad=0.01)
      ax.set_xscale('log', nonposx='clip')
      ax.set_ylim((-3., 4.))
      ax.set_xlabel(r'$\ell$')
      ax.set_ylabel(r'$d\ln C_\ell^{ss} / d\ln \text{Param.}$')

      plt.show()








