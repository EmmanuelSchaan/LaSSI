   def generateDataVector(self, p2d_gg, p2d_gs, p2d_ss):
      '''The data vector is made of the various power spectra,
      with a choice of units that makes the covariance matrix for gg and ss more similar,
      and with an ell^alpha factor that makes the covariance matrix more ell-independent.
      '''
      dataVector = np.zeros(self.nData)
      iData = 0
      # gg
      for iBin1 in range(self.nBins):
         for iBin2 in range(iBin1, self.nBins):
            if (iBin2==iBin1) or self.fullCross:
               dataVector[iData*self.nL:(iData+1)*self.nL] = np.array(map(p2d_gg[iBin1, iBin2].fPinterp, self.L))
               dataVector[iData*self.nL:(iData+1)*self.nL] *= self.L**self.alpha
               iData += 1
      # gs
      for iBin1 in range(self.nBins):
         for iBin2 in range(self.nBins):
            if (iBin2>=iBin1) or self.fullCross:
               dataVector[iData*self.nL:(iData+1)*self.nL] = np.array(map(p2d_gs[iBin1, iBin2].fPinterp, self.L))
               dataVector[iData*self.nL:(iData+1)*self.nL] *= self.sUnit * self.L**self.alpha
               iData += 1
      # ss
      for iBin1 in range(self.nBins):
         for iBin2 in range(iBin1, self.nBins):
            dataVector[iData*self.nL:(iData+1)*self.nL] = np.array(map(p2d_ss[iBin1, iBin2].fPinterp, self.L))
            dataVector[iData*self.nL:(iData+1)*self.nL] *= self.sUnit**2 * self.L**self.alpha
            iData += 1
      return dataVector
