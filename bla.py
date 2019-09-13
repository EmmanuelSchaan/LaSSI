      # Gaussian photo-z
      dndzG = {}
      
      tStart = time()
      # Gaussian photo-z bias and uncertainty for all bins
      dz = photoZPar[:self.nBins]
      ddz = dz[:,None,None]
      sz = photoZPar[self.nBins:2*self.nBins] * (1.+0.5*(zBounds[:-1]+zBounds[1:]))
      ssz = sz[:,None,None]
      
      # axes: [iBin, z, zp]
      z = np.linspace(zMin, zMax, 201)
      zz = z[None,:,None]
      zp = z.copy()
      zzp = zp[None,None,:]
      
      # integrand
      # sharp bins in zp
      integrand = (zBounds[:-1,None,None]<zzp) * (zzp<zBounds[1:,None,None])
      integrand *= np.exp(-0.5*(zz-zzp-ddz)**2/ssz**2) / np.sqrt(2.*np.pi*ssz**2)
      integrand *= w_glsst.dndz(zzp)
      integrand *= zzp[:,:,1:] - zzp[:,:,:-1]
      
      # integrals
      # axes: [iBin, z]
      result = np.sum(integrand, axis=-1)
      
      # interpolate
      for iBin in range(self.nBins):
         dndzG[iBin] = interp1d(z, result[iBin,:], kind='linear', bounds_error=False, fill_value=0.)
      tStop = time()
      if test:
         print "-- getting dn/dz took", tStop-tStart, "sec"
         
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
         # Gaussian photo-z error
         p_z_given_zp = lambda zp,z: np.exp(-0.5*(z-zp-dz)**2/sz**2) / np.sqrt(2.*np.pi*sz**2)
         f = lambda zp,z: w_glsst.dndz(zp) * p_z_given_zp(zp,z)
         dndz_GForInterp = lambda z: integrate.quad(f, zMinP, zMaxP, args=(z), epsabs=0., epsrel=1.e-2)[0]

         tStop = time()
         if test:
            print "-- before dn/dz took", tStop-tStart, "sec"

         # interpolate it for speed (for lensing kernel calculation)
         tStart = time()
         Z = np.linspace(zMin, zMax, 201)
         with sharedmem.MapReduce(np=self.nProc) as pool:
            F = np.array(pool.map(dndz_GForInterp, Z))
         dndzG[iBin] = interp1d(Z, F, kind='linear', bounds_error=False, fill_value=0.)
         tStop = time()
         if test:
            print "-- getting dn/dz took", tStop-tStart, "sec"





########################################

      # Gaussian photo-z
      dndzG = {}
      for iBin in range(self.nBins):
         tStart = time()
         # sharp photo-z cuts
         zMinP = zBounds[iBin]
         zMaxP = zBounds[iBin+1]
         # photo-z bias and uncertainty for this bin:
         dz = photoZPar[iBin]
         sz = photoZPar[self.nBins+iBin] * (1.+0.5*(zMinP+zMaxP))

         # Gaussian photo-z error
         p_z_given_zp = lambda zp,z: np.exp(-0.5*(z-zp-dz)**2/sz**2) / np.sqrt(2.*np.pi*sz**2)
         f = lambda zp,z: w_glsst.dndz(zp) * p_z_given_zp(zp,z)
         dndz_GForInterp = lambda z: integrate.quad(f, zMinP, zMaxP, args=(z), epsabs=0., epsrel=1.e-2)[0]

         tStop = time()
         if test:
            print "-- before dn/dz took", tStop-tStart, "sec"

         # interpolate it for speed (for lensing kernel calculation)
         tStart = time()
         Z = np.linspace(zMin, zMax, 201)
         with sharedmem.MapReduce(np=self.nProc) as pool:
            F = np.array(pool.map(dndz_GForInterp, Z))
         dndzG[iBin] = interp1d(Z, F, kind='linear', bounds_error=False, fill_value=0.)
         tStop = time()
         if test:
            print "-- getting dn/dz took", tStop-tStart, "sec"





