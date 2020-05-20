
   def generatelMaxMask(self, kMaxG=0.3):
      '''Creates a mask to discard the clustering modes
      with ell >= kMax * chi - 0.5,
      where chi is the mean comoving distance to the bin,
      and kMax=0.3 h/Mpc,
      as in the DESC SRD 2018.
      mask is 1 for modes to mask out, 0 otherwise.
      Should be called after the tomo bins have been generated.
      '''
      lMaxMask = np.zeros(self.nData)
      iData = 0
      
      # kk
      iData += self.nKK # ie +=1
      
      # kg
      for iBin1 in range(self.nBins):
         z1 = self.w_g[iBin1].zMean()
         chi1 = self.u.bg.comoving_distance(z1)
         lMax = kMaxG * chi1 - 0.5
         lMaxMask[iData*self.nL:(iData+1)*self.nL] = self.L > lMax
         iData += 1
         
      # ks
      iData += self.nKS
      
      # gg
      for iBin1 in range(self.nBins):
         z1 = self.w_g[iBin1].zMean()
         chi1 = self.u.bg.comoving_distance(z1)
         for iBin2 in range(iBin1, self.nBins):
            z2 = self.w_g[iBin2].zMean()
            chi2 = self.u.bg.comoving_distance(z2)
            chi = min(chi1, chi2)
            lMax = kMaxG * chi - 0.5
            lMaxMask[iData*self.nL:(iData+1)*self.nL] = self.L > lMax
            iData += 1
      # gs
      for iBin1 in range(self.nBins):
         z1 = self.w_g[iBin1].zMean()
         chi1 = self.u.bg.comoving_distance(z1)
         for iBin2 in range(self.nBins):
            lMax = kMaxG * chi1 - 0.5
            lMaxMask[iData*self.nL:(iData+1)*self.nL] = self.L > lMax
            iData += 1
      return lMaxMask
