   def printSnrPowerSpectra(self, path):
      with open(path, 'w') as f:
         f.write("SNR\n\n")
         
         ###########################################################
         # gg
         
         # gg: auto
         f.write("GG\n")
         f.write("auto\n")
         i1 = 0
         Itotal = []
         for iBin1 in range(self.nBins):
            I = range(i1*self.nL, (i1+1)*self.nL)
            d = self.dataVector[I]
            J = np.ix_(I,I)
            cov = self.covMat[J]
            invCov = np.linalg.inv(cov)
            snr = np.dot(d.transpose(), np.dot(invCov, d))
            snr = np.sqrt(snr)
            f.write("   "+str(iBin1)+","+str(iBin1)+": "+str(snr)+"\n")
            i1 += self.nBins - iBin1
            Itotal += I
         # gg: total auto
         d = self.dataVector[Itotal]
         J = np.ix_(Itotal,Itotal)
         cov = self.covMat[J]
         invCov = np.linalg.inv(cov)
         snr = np.dot(d.transpose(), np.dot(invCov, d))
         snr = np.sqrt(snr)
         f.write("total auto: "+str(snr))
         
         
         # gg: cross i,i+1
         f.write("cross i,i+1\n")
         i1 = 1
         Itotal = []
         for iBin1 in range(self.nBins-1):
            I = range(i1*self.nL, (i1+1)*self.nL)
            d = self.dataVector[I]
            J = np.ix_(I,I)
            cov = self.covMat[J]
            invCov = np.linalg.inv(cov)
            snr = np.dot(d.transpose(), np.dot(invCov, d))
            snr = np.sqrt(snr)
            f.write("   "+str(iBin1)+","+str(iBin1+i1)+": "+str(snr)+"\n")
            i1 += self.nBins - iBin1
            Itotal += I
         # gg: total i,i+1
         d = self.dataVector[Itotal]
         J = np.ix_(Itotal,Itotal)
         cov = self.covMat[J]
         invCov = np.linalg.inv(cov)
         snr = np.dot(d.transpose(), np.dot(invCov, d))
         snr = np.sqrt(snr)
         f.write("total i,i+1: "+str(snr))


         # gg: cross i,i+2
         f.write("cross i,i+2\n")
         i1 = 2
         for iBin1 in range(self.nBins-2):
            I = range(i1*self.nL, (i1+1)*self.nL)
            d = self.dataVector[I]
            J = np.ix_(I,I)
            cov = self.covMat[J]
            invCov = np.linalg.inv(cov)
            snr = np.dot(d.transpose(), np.dot(invCov, d))
            snr = np.sqrt(snr)
            f.write("   "+str(iBin1)+","+str(iBin1+i1)+": "+str(snr)+"\n")
            i1 += self.nBins - iBin1
         # gg: total i,i+2
         d = self.dataVector[Itotal]
         J = np.ix_(Itotal,Itotal)
         cov = self.covMat[J]
         invCov = np.linalg.inv(cov)
         snr = np.dot(d.transpose(), np.dot(invCov, d))
         snr = np.sqrt(snr)
         f.write("total i,i+2: "+str(snr))
         

         # gg: total
         I = range(0., self.nGG*self.nL)
         d = self.dataVector[I]
         J = np.ix_(I,I)
         cov = self.covMat[J]
         invCov = np.linalg.inv(cov)
         snr = np.dot(d.transpose(), np.dot(invCov, d))
         snr = np.sqrt(snr)
         f.write("total gg: "+str(snr)+"\n\n")


         ###########################################################
         # gs

         # gs: all
         f.write("GS\n")
         i1 = self.nGG
         for iBin1 in range(self.nBins):
            for iBin2 in range(self.nBins):
               I = range(i1*self.nL, (i1+1)*self.nL)
               d = self.dataVector[I]
               J = np.ix_(I,I)
               cov = self.covMat[J]
               invCov = np.linalg.inv(cov)
               snr = np.dot(d.transpose(), np.dot(invCov, d))
               snr = np.sqrt(snr)
               f.write("   "+str(iBin1)+","+str(iBin2)+": "+str(snr)+"\n")
               i1 += 1
         # gs: total
         I = range(self.nGG*self.nL, (self.nGG+self.nGS)*self.nL)
         d = self.dataVector[I]
         J = np.ix_(I,I)
         cov = self.covMat[J]
         invCov = np.linalg.inv(cov)
         snr = np.dot(d.transpose(), np.dot(invCov, d))
         snr = np.sqrt(snr)
         f.write("total gs: "+str(snr)+"\n\n")


         ###########################################################
         # ss

         # ss: all
         f.write("SS\n")
         i1 = self.nGG + self.nGS
         for iBin1 in range(self.nBins):
            for iBin2 in range(iBin1, self.nBins):
               I = range(i1*self.nL, (i1+1)*self.nL)
               d = self.dataVector[I]
               J = np.ix_(I,I)
               cov = self.covMat[J]
               invCov = np.linalg.inv(cov)
               snr = np.dot(d.transpose(), np.dot(invCov, d))
               snr = np.sqrt(snr)
               f.write("   "+str(iBin1)+","+str(iBin2)+": "+str(snr)+"\n")
               i1 += 1

         # ss: total
         I = range((self.nGG+self.nGS)*self.nL, (self.nGG+self.nGS+self.nSS)*self.nL)
         d = self.dataVector[I]
         J = np.ix_(I,I)
         cov = self.covMat[J]
         invCov = np.linalg.inv(cov)
         snr = np.dot(d.transpose(), np.dot(invCov, d))
         snr = np.sqrt(snr)
         f.write("total ss: "+str(snr)+"\n\n")











