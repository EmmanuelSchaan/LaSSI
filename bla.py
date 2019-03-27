   def generateTomoBins(self, u, nuisancePar, save=True, doS=True):
      '''The option doS=False is only used to save time when sampling dn/dz,
      since the s bins take much longer to generate (by a factor ~100).
      '''
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
      
      # extract the outlier contamination matrix, if needed
      if len(photoZPar)==self.nBins*(self.nBins+1):
         cij = photoZPar[2*self.nBins:]











      def createTomoBin(iBin:)

#      for iBin in range(self.nBins):
         # sharp photo-z cuts
         zMinP = zBounds[iBin]
         zMaxP = zBounds[iBin+1]
         # photo-z bias and uncertainty for this bin:
#!!!! I am using a loose mean redshift
         dz = photoZPar[iBin]
         sz = photoZPar[self.nBins+iBin] * (1.+0.5*(zMinP+zMaxP))
         # true z bounds: truncate at 5 sigma
         # careful for the first and last bin
         zMin = 1./w_glsst.aMax-1.  #max(zMinP - 5.*sz, 1./w_glsst.aMax-1.)   # 1./w_glsst.aMax-1.
         zMax = 1./w_glsst.aMin-1.  #min(zMaxP + 5.*sz, 1./w_glsst.aMin-1.)   # 1./w_glsst.aMin-1.
         
         
         # true dn/dz_true from dn/dz_phot
         p_z_given_zp = lambda zp,z: np.exp(-0.5*(z-zp-dz)**2/sz**2) / np.sqrt(2.*np.pi*sz**2)
         
         
         # If outliers are not included
         if len(photoZPar)==2*self.nBins:
            f = lambda zp,z: w_glsst.dndz(zp) * p_z_given_zp(zp,z)
            dndz_tForInterp = lambda z: integrate.quad(f, zMinP, zMaxP, args=(z), epsabs=0., epsrel=1.e-3)[0]
      
         # If outliers are included
         elif len(photoZPar)==self.nBins*(self.nBins+1):
            
            def dndzp_outliers(zp):
               ''' This is the dn/dz_p, such that for bin i:
               n_i^new = n_i * (1-sum_{j \neq i} c_ij) + sum_{j \neq i} c_ji n_j.
               '''
               result = w_glsst.dndz(zp)
               # if zp is in the bin
               if zp>=zBounds[iBin] and zp<zBounds[iBin+1]:
                  result *= 1. - np.sum([cij[iBin*(self.nBins-1)+j] for j in range(self.nBins-1)])
               # if zp is in another bin
               else:
                  # find which bin this is
                  jBin = np.where(np.array([(zp>=zBounds[j])*(zp<zBounds[j+1]) for j in range(self.nBins)])==1)[0][0]
                  # since the diagonal c_ii is not encoded, make sure to skip it if iBin > jBin
                  i = iBin - (iBin>jBin)
                  result *= cij[jBin*(self.nBins-1)+i]
               return result

            f = lambda zp,z: dndzp_outliers(zp) * p_z_given_zp(zp,z)
            dndz_tForInterp = lambda z: integrate.quad(f, zMin, zMax, args=(z), epsabs=0., epsrel=1.e-3)[0]
      
         else:
            print "Error: PhotoZPar does not have the right size"

         
         # interpolate it for speed (for lensing kernel calculation)
#         Z = np.linspace(zMin, zMax, 101)
         Z = np.linspace(zMin, zMax, 501)
         F = np.array(map(dndz_tForInterp, Z))
         dndz_t = interp1d(Z, F, kind='linear', bounds_error=False, fill_value=0.)



#!!!!!! Bottleneck is the shear bin, by a factor 100 compared to lens bin and getting dndz_t
#         tStart = time()
         # shear bin
         if doS:
            w_s_current = WeightLensCustom(u,
                                         dndz_t, # dn/dz_true
                                         m=lambda z: shearMultBiasPar[iBin], # multiplicative shear bias
                                         zMinG=zMin,
                                         zMaxG=zMax,
                                         name='s'+str(iBin))
#         tStop = time()
#         print "-- shear bin took", tStop-tStart, "sec"


#         tStart = time()
         # tracer bin
         w_g_current = WeightTracerCustom(u,
                                        lambda z: galaxyBiasPar[iBin] * w_glsst.b(z), # galaxy bias
                                        dndz_t, # dn/dz_true
                                        zMin=zMin,
                                        zMax=zMax,
                                        name='g'+str(iBin))
#         tStop = time()
#         print "-- clustering bin took", tStop-tStart, "sec"


         # add magnification bias, if requested
         if self.magBias:
   #!!!! I am using a loose mean redshift
            alpha = w_glsst.magnificationBias(0.5*(zMinP+zMaxP))
            print "bin "+str(iBin)+": mag bias alpha="+str(alpha)
            w_g_nomagbias_current = WeightTracerCustom(u,
                                           lambda z: galaxyBiasPar[iBin] * w_glsst.b(z), # galaxy bias
                                           dndz_t, # dn/dz_true
                                           zMin=zMin,
                                           zMax=zMax,
                                           name='g'+str(iBin))
            w_g_current.f = lambda a: w_g_nomagbias[iBin].f(a) + 2.*(alpha-1.)*w_s[iBin].f(a)


         #print "- done "+str(iBin+1)+" of "+str(self.nBins)
      #print "total ngal="+str(np.sum([w_g[i].ngal_per_arcmin2 for i in range(self.nBins)]))+"/arcmin2, should be "+str(w_glsst.ngal_per_arcmin2)
      
      
         return w_g_current, w_s_current
      
      
      
      
      
      with sharedmem.MapReduce(np=nProc) as pool:
         result = np.array(pool.map(createTomoBin, range(self.nBins)))
      w_g = result[:,0]
      w_s = result[:,1]
      
      
      
      if doS:
         return w_g, w_s, zBounds
      else:
         return w_g, zBounds
