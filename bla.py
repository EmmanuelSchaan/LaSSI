
   
   
   def generateCov(self, p_kk, p_kg, p_ks, p_gg, p_gs, p_ss, p_kk_shot, p_gg_shot, p_ss_shot):
      covMat = np.zeros((self.nData, self.nData))
      # below, i1 and i2 define the row and column of the nL*nL blocks for each pair of 2-point function
      # i1, i2 \in [0, n2pt]
   
      # include the shot noises
      p_kk += p_kk_shot
      p_gg += p_gg_shot
      p_ss += p_ss_shot
      # generic Gaussian cov
      cov = lambda Pac, Pbd, Pad, Pbc, Npairs: np.diagflat((Pac * Pbd + Pad * Pbc) / self.Npairs)
      
      # "kk-kk"
      i1 = 0
      i2 = 0
      covBlock = cov(p_kk, p_kk, p_kk, p_kk, self.Nmodes)
      covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.L**(2*self.alpha) * covBlock
      
      # "kk-kg"
      i1 = 0
      i2 = self.nKK
      # considering kg[i1]
      for iBin1 in range(self.nBins):
         covBlock = cov(p_kk, p_kg[iBin1], p_kg[iBin1], p_kk, self.Nmodes)
         covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.L**(2*self.alpha) * covBlock
         # move to next column
         i2 += 1
      
      # "kk-ks"
      i1 = 0
      i2 = self.nKK + self.nKG
      # considering ks[i1]
      for iBin1 in range(self.nBins):
         covBlock = cov(p_kk, p_ks[iBin1], p_ks[iBin1], p_kk, self.Nmodes)
         covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.sUnit * self.L**(2*self.alpha) * covBlock
         # move to next column
         i2 += 1

      # "kk-gg"
      i1 = 0
      i2 = self.nKK + self.nKG + self.nKS
      # considering gg[i1, j1]
      for iBin1 in range(self.nBins):
         for iBin2 in range(iBin1, self.nBins):
            covBlock = cov(p_kg[iBin1], p_kg[iBin2], p_kg[iBin2], p_kg[iBin1], self.Nmodes)
            covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.L**(2*self.alpha) * covBlock
            # move to next column
            i2 += 1

      # "kk-gs"
      i1 = 0
      i2 = self.nKK + self.nKG + self.nKS + self.nGG
      # considering gs[i1, j1]
      for iBin1 in range(self.nBins):
         for iBin2 in range(self.nBins):
            covBlock = cov(p_kg[iBin1], p_ks[iBin2], p_ks[iBin2], p_kg[iBin1], self.Nmodes)
            covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.sUnit * self.L**(2*self.alpha) * covBlock
            # move to next column
            i2 += 1

      # "kk-ss"
      i1 = 0
      i2 = self.nKK + self.nKG + self.nKS + self.nGG + self.nGS
      # considering ss[i1, j1]
      for iBin1 in range(self.nBins):
         for iBin2 in range(iBin1, self.nBins):
            covBlock = cov(p_ks[iBin1], p_ks[iBin2], p_ks[iBin2], p_ks[iBin1], self.Nmodes)
            covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.sUnit**2 * self.L**(2*self.alpha) * covBlock
            # move to next column
            i2 += 1

      # "kg-kg"
      i1 = self.nKK
      # considering kg[i1]
      for iBin1 in range(self.nBins):
         i2 = self.nKK
         # considering kg[j1]
         for jBin1 in range(self.nBins):
            # compute only upper diagonal
            if i2>=i1:
               covBlock = cov(p_kk, p_gg[iBin1, jBin1], p_kg[jBin1], p_kg[iBin1], self.Nmodes)
               covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.L**(2*self.alpha) * covBlock
            # move to next column
            i2 += 1
         # move to next row
         i1 += 1
      
      # "kg-ks"
      i1 = self.nKK
      # considering kg[i1]
      for iBin1 in range(self.nBins):
         i2 = self.nKK + self.nKG
         # considering ks[j1]
         for jBin1 in range(self.nBins):
            # compute only upper diagonal
            if i2>=i1:
               covBlock = cov(p_kk, p_gs[iBin1, jBin1], p_ks[jBin1], p_kg[iBin1], self.Nmodes)
               covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.sUnit * self.L**(2*self.alpha) * covBlock
            # move to next column
            i2 += 1
         # move to next row
         i1 += 1

      # "kg-gg"
      i1 = self.nKK
      # considering kg[i1]
      for iBin1 in range(self.nBins):
         i2 = self.nKK + self.nKG + self.nKS
         # considering gg[j1, j2]
         for jBin1 in range(self.nBins):
            for jBin2 in range(jBin1, self.nBins):
               # compute only upper diagonal
               if i2>=i1:
                  covBlock = cov(p_kg[jBin1], p_gg[iBin1, jBin2], p_kg[jBin2], p_gg[iBin1, jBin1], self.Nmodes)
                  covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.L**(2*self.alpha) * covBlock
               # move to next column
               i2 += 1
         # move to next row
         i1 += 1

      # "kg-gs"
      i1 = self.nKK
      # considering kg[i1]
      for iBin1 in range(self.nBins):
         i2 = self.nKK + self.nKG + self.nKS + self.nGG
         # considering gs[j1, j2]
         for jBin1 in range(self.nBins):
            for jBin2 in range(self.nBins):
               # compute only upper diagonal
               if i2>=i1:
                  covBlock = cov(p_kg[jBin1], p_gs[iBin1, jBin2], p_ks[jBin2], p_gg[iBin1, jBin1], self.Nmodes)
                  covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.sUnit * self.L**(2*self.alpha) * covBlock
               # move to next column
               i2 += 1
         # move to next row
         i1 += 1

      # "kg-ss"
      # considering kg[i1]
      i1 = self.nKK
      for iBin1 in range(self.nBins):
         i2 = self.nKK + self.nKG + self.nKS + self.nGG + self.nGS
         # considering ss[j1, j2]
         for jBin1 in range(self.nBins):
            for jBin2 in range(jBin1, self.nBins):
               # compute only upper diagonal
               if i2>=i1:
                  covBlock = cov(p_ks[jBin1], p_gs[iBin1, jBin2], p_ks[jBin2], p_gs[iBin1, jBin1], self.Nmodes)
                  covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.sUnit**2 * self.L**(2*self.alpha) * covBlock
               # move to next column
               i2 += 1
         # move to next row
         i1 += 1

      # "ks-ks"
      i1 = self.nKK + self.nKG
      # considering ks[i1]
      for iBin1 in range(self.nBins):
         i2 = self.nKK + self.nKG
         # considering ks[j1]
         for jBin1 in range(self.nBins):
            # compute only upper diagonal
            if i2>=i1:
               covBlock = cov(p_kk, p_ss[iBin1, jBin1], p_ks[jBin1], p_ks[iBin1], self.Nmodes)
               covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.sUnit**2 * self.L**(2*self.alpha) * covBlock
            # move to next column
            i2 += 1
         # move to next row
         i1 += 1

      # "ks-gg"
      i1 = self.nKK + self.nKG
      # considering ks[i1]
      for iBin1 in range(self.nBins):
         i2 = self.nKK + self.nKG + self.nKS
         # considering gg[j1, j2]
         for jBin1 in range(self.nBins):
            for jBin2 in range(jBin1, self.nBins):
               # compute only upper diagonal
               if i2>=i1:
                  # watch the order for gs
                  covBlock = cov(p_kg[jBin1], p_gs[jBin2, iBin1], p_kg[jBin2], p_gs[jBin1, iBin1], self.Nmodes)
                  covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.sUnit * self.L**(2*self.alpha) * covBlock
               # move to next column
               i2 += 1
         # move to next row
         i1 += 1

      # "ks-gs"
      i1 = self.nKK + self.nKG
      # considering ks[i1]
      for iBin1 in range(self.nBins):
         i2 = self.nKK + self.nKG + self.nKS + self.nGG
         # considering gs[j1, j2]
         for jBin1 in range(self.nBins):
            for jBin2 in range(self.nBins):
               # compute only upper diagonal
               if i2>=i1:
                  # watch the order for gs
                  covBlock = cov(p_kg[jBin1], p_ss[iBin1, jBin2], p_ks[jBin2], p_gs[jBin1, iBin1], self.Nmodes)
                  covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.sUnit**2 * self.L**(2*self.alpha) * covBlock
               # move to next column
               i2 += 1
         # move to next row
         i1 += 1

      # "ks-ss"
      i1 = self.nKK + self.nKG
      # considering ks[i1]
      for iBin1 in range(self.nBins):
         i2 = self.nKK + self.nKG + self.nKS + self.nGG + self.nGS
         # considering ss[j1, j2]
         for jBin1 in range(self.nBins):
            for jBin2 in range(jBin1, self.nBins):
               # compute only upper diagonal
               if i2>=i1:
                  covBlock = cov(p_ks[jBin1], p_ss[iBin1, jBin2], p_ks[jBin2], p_ss[iBin1, jBin1], self.Nmodes)
                  covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.sUnit**3 * self.L**(2*self.alpha) * covBlock
               # move to next column
               i2 += 1
         # move to next row
         i1 += 1
      
      #print "gg-gg"
      # considering gg[i1,i2]
      i1 = self.nKK + self.nKG + self.nKS
      for iBin1 in range(self.nBins):
         for iBin2 in range(iBin1, self.nBins):
            # considering gg[j1,j2]
            i2 = self.nKK + self.nKG + self.nKS
            for jBin1 in range(self.nBins):
               for jBin2 in range(jBin1, self.nBins):
                  # compute only upper diagonal
                  if i2>=i1:
                     covBlock = cov(p_gg[iBin1,jBin1], p_gg[iBin2,jBin2], p_gg[iBin1,jBin2], p_gg[iBin2,jBin1], self.Nmodes)
                     covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.L**(2*self.alpha) * covBlock
                  # move to next column
                  i2 += 1
            # move to next row
            i1 += 1

      #print "gg-gs"
      # considering gg[i1,i2]
      i1 = self.nKK + self.nKG + self.nKS
      for iBin1 in range(self.nBins):
         for iBin2 in range(iBin1, self.nBins):
            # considering gs[j1,j2]
            i2 = self.nKK + self.nKG + self.nKS + self.nGG
            for jBin1 in range(self.nBins):
               for jBin2 in range(self.nBins):
                  # compute only upper diagonal
                  if i2>=i1:
                     covBlock = cov(p_gg[iBin1,jBin1], p_gs[iBin2,jBin2], p_gs[iBin1,jBin2], p_gg[iBin2,jBin1], self.Nmodes)
                     covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.sUnit * self.L**(2*self.alpha) *  covBlock
                  # move to next column
                  i2 += 1
            # move to next row
            i1 += 1

      #print "gg-ss"
      # considering gg[i1,i2]
      i1 = self.nKK + self.nKG + self.nKS
      for iBin1 in range(self.nBins):
         for iBin2 in range(iBin1, self.nBins):
            # considering ss[j1,j2]
            i2 = self.nKK + self.nKG + self.nKS + self.nGG + self.nGS
            for jBin1 in range(self.nBins):
               for jBin2 in range(jBin1, self.nBins):
                  # compute only upper diagonal
                  if i2>=i1:
                     covBlock = cov(p_gs[iBin1,jBin1], p_gs[iBin2,jBin2], p_gs[iBin1,jBin2], p_gs[iBin2,jBin1], self.Nmodes)
                     covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.sUnit**2 * self.L**(2*self.alpha) * covBlock
                  # move to next column
                  i2 += 1
            # move to next row
            i1 += 1

      #print "gs-gs"
      # considering gs[i1,i2]
      i1 = self.nKK + self.nKG + self.nKS + self.nGG
      for iBin1 in range(self.nBins):
         for iBin2 in range(self.nBins):
            # considering gs[j1,j2]
            i2 = self.nKK + self.nKG + self.nKS + self.nGG
            for jBin1 in range(self.nBins):
               for jBin2 in range(self.nBins):
                  # compute only upper diagonal
                  if i2>=i1:
                     # watch the order for gs
                     covBlock = cov(p_gg[iBin1,jBin1], p_ss[iBin2,jBin2], p_gs[iBin1,jBin2], p_gs[jBin1,iBin2], self.Nmodes)
                     covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.sUnit**2 * self.L**(2*self.alpha) * covBlock
                  # move to next column
                  i2 += 1
            # move to next row
            i1 += 1

      #print "gs-ss"
      # considering gs[i1,i2]
      i1 = self.nKK + self.nKG + self.nKS + self.nGG
      for iBin1 in range(self.nBins):
         for iBin2 in range(self.nBins):
            # considering ss[j1,j2]
            i2 = self.nKK + self.nKG + self.nKS + self.nGG + self.nGS
            for jBin1 in range(self.nBins):
               for jBin2 in range(jBin1, self.nBins):
                  # compute only upper diagonal
                  if i2>=i1:
                     covBlock = cov(p_gs[iBin1,jBin1], p_ss[iBin2,jBin2], p_gs[iBin1,jBin2], p_ss[iBin2,jBin1], self.Nmodes)
                     covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.sUnit**3 * self.L**(2*self.alpha) * covBlock
                  # move to next column
                  i2 += 1
            # move to next row
            i1 += 1

      #print "ss-ss"
      # considering ss[i1,i2]
      i1 = self.nKK + self.nKG + self.nKS + self.nGG + self.nGS
      for iBin1 in range(self.nBins):
         for iBin2 in range(iBin1, self.nBins):
            # considering ss[j1,j2]
            i2 = self.nKK + self.nKG + self.nKS + self.nGG + self.nGS
            for jBin1 in range(self.nBins):
               for jBin2 in range(jBin1, self.nBins):
                  # compute only upper diagonal
                  if i2>=i1:
                     covBlock = cov(p_ss[iBin1,jBin1], p_ss[iBin2,jBin2], p_ss[iBin1,jBin2], p_ss[iBin2,jBin1], self.Nmodes)
                     covMat[i1*self.nL:(i1+1)*self.nL, i2*self.nL:(i2+1)*self.nL] = self.sUnit**4 * self.L**(2*self.alpha) * covBlock
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
