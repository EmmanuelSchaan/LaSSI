
   def plotUncertaintyPowerSpectra(self, show=False):

      # kk
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      i1 = 0
      I = range(i1*self.nL, (i1+1)*self.nL)
      L = extractMaskedVec(self.L, mask=self.lMaxMask[I])
      d = extractMaskedVec(self.dataVector, mask=self.lMaxMask, I=I)
      shot = extractMaskedVec(self.shotNoiseVector, mask=self.lMaxMask, I=I)
      cov = extractMaskedMat(self.covMat, mask=self.lMaxMask, I=I)
      std = np.sqrt(np.diag(cov))
      #
      ax.plot(L, std / d, '.-', lw=2, color='b')
      #
      ax.set_xlim((90., 2.5e3))
      ax.set_xscale('log')
      ax.set_yscale('log', nonposy='clip')
      #
      ax.set_title(r'CMB lensing')
      ax.set_xlabel(r'$\ell$')
      ax.set_ylabel(r'$\sigma \left( C_\ell^{\kappa_\text{CMB} \kappa_\text{CMB}} \right) / C_\ell^{\kappa_\text{CMB} \kappa_\text{CMB}}$', fontsize=18)
      #
      fig.savefig(self.figurePath+"/sp2d_kk.pdf", bbox_inches='tight')
      if show:
         plt.show()
      fig.clf()


      # kg
      Colors = plt.cm.autumn(1.*np.arange(self.nBins)/(self.nBins-1.))
      #
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      i1 = self.nKK
      for iBin1 in range(self.nBins):
         I = range(i1*self.nL, (i1+1)*self.nL)
         L = extractMaskedVec(self.L, mask=self.lMaxMask[I])
         d = extractMaskedVec(self.dataVector, mask=self.lMaxMask, I=I)
         shot = extractMaskedVec(self.shotNoiseVector, mask=self.lMaxMask, I=I)
         cov = extractMaskedMat(self.covMat, mask=self.lMaxMask, I=I)
         std = np.sqrt(np.diag(cov))
         #
         color = Colors[iBin1]
         ax.plot(L, std/d, '.-', lw=2, color=color, label=r'$\langle \kappa_\text{CMB}\, g_{'+str(iBin1)+r'}  \rangle$')
         # advance counter in data vector
         i1 += 1
      #
      ax.set_xlim((90., 2.5e3))
      ax.legend(loc=1, fontsize='x-small', labelspacing=0., handlelength=0.4, borderaxespad=0.05)
      ax.set_xscale('log')
      ax.set_yscale('log', nonposy='clip')
      #
      ax.set_title(r'Galaxy $\times$ CMB lensing')
      ax.set_xlabel(r'$\ell$')
      ax.set_ylabel(r'$\sigma \left( C_\ell^{g \kappa_\text{CMB}} \right) / C_\ell^{g \kappa_\text{CMB}}$', fontsize=18)
      #
      fig.savefig(self.figurePath+"/sp2d_kg.pdf")
      if show:
         plt.show()
      fig.clf()
      

      # ks
      Colors = plt.cm.autumn(1.*np.arange(self.nBins)/(self.nBins-1.))
      #
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      i1 = self.nKK + self.nKG
      for iBin1 in range(self.nBins):
         I = range(i1*self.nL, (i1+1)*self.nL)
         L = extractMaskedVec(self.L, mask=self.lMaxMask[I])
         d = extractMaskedVec(self.dataVector, mask=self.lMaxMask, I=I) / self.sUnit
         shot = extractMaskedVec(self.shotNoiseVector, mask=self.lMaxMask, I=I) / self.sUnit
         cov = extractMaskedMat(self.covMat, mask=self.lMaxMask, I=I)
         std = np.sqrt(np.diag(cov)) / self.sUnit
         #
         color = Colors[iBin1]
         ax.plot(L, std/d, '.-', lw=2, color=color, label=r'$\langle \kappa_\text{CMB}\, \kappa_{g_{'+str(iBin1)+r'}}  \rangle$')
         ax.plot(L, shot, ls='--', lw=1, color='grey')  # same color for all tomo bins, since they have the same n_gal
         # advance counter in data vector
         i1 += 1
      #
      ax.set_xlim((90., 2.5e3))
      ax.legend(loc=1, fontsize='x-small', labelspacing=0., handlelength=0.4, borderaxespad=0.05)
      ax.set_xscale('log')
      ax.set_yscale('log', nonposy='clip')
      #
      ax.set_title(r'CMB lensing $\times$ galaxy lensing')
      ax.set_xlabel(r'$\ell$')
      ax.set_ylabel(r'$\sigma \left( C_\ell^{\kappa_\text{CMB} \kappa_g} \right) / C_\ell^{\kappa_\text{CMB} \kappa_g}$', fontsize=18)
      #
      fig.savefig(self.figurePath+"/sp2d_ks.pdf")
      if show:
         plt.show()
      fig.clf()


      # gg: panels
      Colors = plt.cm.autumn(1.*np.arange(self.nBins)/(self.nBins-1.))
      #
      fig=plt.figure(0)
      gs = gridspec.GridSpec(3, 1)#, height_ratios=[1, 1, 1])
      gs.update(hspace=0.)
      
      # auto
      ax0=plt.subplot(gs[0])
      i1 = self.nKK + self.nKG + self.nKS
      for iBin1 in [0,1,2,3]:
         color = Colors[iBin1]
         ax0.plot([], [], c=color, label=r'$i='+str(iBin1)+'$')
      for iBin1 in range(self.nBins):
         I = range(i1*self.nL, (i1+1)*self.nL)
         L = extractMaskedVec(self.L, mask=self.lMaxMask[I])
         d = extractMaskedVec(self.dataVector, mask=self.lMaxMask, I=I)
         shot = extractMaskedVec(self.shotNoiseVector, mask=self.lMaxMask, I=I)
         cov = extractMaskedMat(self.covMat, mask=self.lMaxMask, I=I)
         std = np.sqrt(np.diag(cov))
         #
         color = Colors[iBin1]
         ax0.plot(L, std / d, '.-', lw=2, color=color)
         # advance counter in data vector
         i1 += self.nBins - iBin1
      #
      ax0.legend(loc=1, fontsize='x-small', labelspacing=0., handlelength=0.4, borderaxespad=0.05)
      ax0.set_xscale('log')
      ax0.set_yscale('log', nonposy='clip')
      plt.setp(ax0.get_xticklabels(), visible=False)
      #
      ax0.set_title(r'Clustering: $\sigma\left( C_\ell \right) / C_\ell$')
      ax0.set_ylabel(r'$g_i g_i$', fontsize=18)

      # cross i,i+1
      ax1=plt.subplot(gs[1])
      i1 = self.nKK + self.nKG + self.nKS + 1
      for iBin1 in [4,5,6,7]:
         color = Colors[iBin1]
         ax1.plot([], [], c=color, label=r'$i='+str(iBin1)+'$')
      for iBin1 in range(self.nBins-1):
         I = range(i1*self.nL, (i1+1)*self.nL)
         L = extractMaskedVec(self.L, mask=self.lMaxMask[I])
         d = extractMaskedVec(self.dataVector, mask=self.lMaxMask, I=I)
         cov = extractMaskedMat(self.covMat, mask=self.lMaxMask, I=I)
         std = np.sqrt(np.diag(cov))
         #
         color = Colors[iBin1]
         ax1.plot(L, std / d, '-', lw=2, color=color)
         # advance counter in data vector
         i1 += self.nBins - iBin1
      #
      ax1.legend(loc=1, fontsize='x-small', labelspacing=0., handlelength=0.4, borderaxespad=0.05)
      ax1.set_xscale('log')
      ax1.set_yscale('log', nonposy='clip')
      plt.setp(ax1.get_xticklabels(), visible=False)
      #
      ax1.set_ylabel(r'$g_i g_{i+1}$', fontsize=18)

      # cross i,i+2
      ax2=plt.subplot(gs[2])
      i1 = self.nKK + self.nKG + self.nKS + 2
      for iBin1 in [8,9]:
         color = Colors[iBin1]
         ax2.plot([], [], c=color, label=r'$i='+str(iBin1)+'$')
      for iBin1 in range(self.nBins-2):
         I = range(i1*self.nL, (i1+1)*self.nL)
         L = extractMaskedVec(self.L, mask=self.lMaxMask[I])
         d = extractMaskedVec(self.dataVector, mask=self.lMaxMask, I=I)
         cov = extractMaskedMat(self.covMat, mask=self.lMaxMask, I=I)
         std = np.sqrt(np.diag(cov))
         #
         color = Colors[iBin1]
         ax2.plot(L, std / d, '.-', lw=2, color=color)
         # advance counter in data vector
         i1 += self.nBins - iBin1
      #
      ax2.legend(loc=1, fontsize='x-small', labelspacing=0., handlelength=0.4, borderaxespad=0.05)
      #
      ax2.set_xlim((90., 2.5e3))
      ax2.set_xscale('log')
      ax2.set_yscale('log', nonposy='clip')
      #
      ax2.set_ylabel(r'$g_i g_{i+2}$', fontsize=18)
      ax2.set_xlabel(r'$\ell$')
      #
      fig.savefig(self.figurePath+"/sp2d_gg.pdf")
      if show:
         plt.show()
      fig.clf()



      # gs
      Colors = plt.cm.winter(1.*np.arange(self.nBins)/(self.nBins-1.))
      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      i1 = self.nKK + self.nKG + self.nKS + self.nGG
      for iBin1 in range(self.nBins):
         color = Colors[iBin1]
         ax.plot([], [], c=color, label=r'$\langle g_{'+str(iBin1)+r'}\, \kappa_{g_j} \rangle$')

#         # show only gs for the highest s bin 
#         i1 += self.nBins - 1
#         #
#         I = range(i1*self.nL, (i1+1)*self.nL)
#         L = extractMaskedVec(self.L, mask=self.lMaxMask[I])
#         d = extractMaskedVec(self.dataVector, mask=self.lMaxMask, I=I) / self.sUnit
#         cov = extractMaskedMat(self.covMat, mask=self.lMaxMask, I=I)
#         std = np.sqrt(np.diag(cov)) / self.sUnit
#         #
#         ax.errorbar(L*(1.+0.01*i1/self.nGS), d, yerr=std, ls='-', lw=1, elinewidth=1.5, marker='.', markersize=2, color=color)# label=r'$\ell \langle g_{'+str(iBin1)+'} \gamma_{'+str(iBin2)+r'}\rangle$')
#         # move to next row
#         i1 += 1

         for iBin2 in range(self.nBins):  # show all the cross-correlations
            I = range(i1*self.nL, (i1+1)*self.nL)
            L = extractMaskedVec(self.L, mask=self.lMaxMask[I])
            d = extractMaskedVec(self.dataVector, mask=self.lMaxMask, I=I) / self.sUnit
            cov = extractMaskedMat(self.covMat, mask=self.lMaxMask, I=I)
            std = np.sqrt(np.diag(cov)) / self.sUnit
            #
            ax.plot(L*(1.+0.01*i1/self.nGS), std / d, '.-', lw=1, color=color)
            # move to next row
            i1 += 1

      #
      ax.set_xlim((90., 2.5e3))
      ax.legend(loc=1, fontsize='x-small', labelspacing=0., handlelength=0.4, borderaxespad=0.05)
      ax.set_xscale('log')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$\ell$')
      ax.set_ylabel(r'$\sigma\left( C_\ell^{g\kappa_g} \right) / C_\ell^{g\kappa_g}$')
      ax.set_title(r'Galaxy $\times$ galaxy lensing')
      #
      fig.savefig(self.figurePath+"/sp2d_gs.pdf")
      if show:
         plt.show()
      fig.clf()


      # ss: all on same plot
      Colors = plt.cm.jet(1.*np.arange(self.nBins)/(self.nBins-1.))
      fig=plt.figure(2)
      ax=fig.add_subplot(111)
      #
      i1 = self.nKK + self.nKG + self.nKS + self.nGG + self.nGS
      for iBin1 in range(self.nBins):
         # add entry to caption
         color = Colors[iBin1]
         ax.plot([], [], c=color, label=r'$\langle\kappa_{g_i} \kappa_{g_{i+'+str(iBin1)+r'}} \rangle $')


#         # show only the auto
#         I = range(i1*self.nL, (i1+1)*self.nL)
#         L = extractMaskedVec(self.L, mask=self.lMaxMask[I])
#         d = extractMaskedVec(self.dataVector, mask=self.lMaxMask, I=I) / self.sUnit**2
#         shot = extractMaskedVec(self.shotNoiseVector, mask=self.lMaxMask, I=I) / self.sUnit**2
#         cov = extractMaskedMat(self.covMat, mask=self.lMaxMask, I=I)
#         std = np.sqrt(np.diag(cov)) / self.sUnit**2
#         #
#         ax.errorbar(L*(1.+0.01*i1/self.nSS), d, yerr=std, ls='-', lw=1, elinewidth=1.5, marker='.', markersize=2, color=color)#, label=r'$\gamma_{'+str(iBin1)+'} \gamma_{'+str(iBin2)+'}$')
##            ax.plot(L*(1.+0.01*i1/self.nSS), shot, ls='--', lw=1, color=color)   # different color for each bin
#         ax.plot(L*(1.+0.01*i1/self.nSS), shot, ls='--', lw=1, color='grey')   # same color for all bins
#         # move to next row
#         i1 += self.nBins - iBin1

         # show all the cross-correlations
         for iBin2 in range(iBin1, self.nBins):
            color = Colors[iBin2-iBin1]
            #color = Colors[iBin1]
            #
            I = range(i1*self.nL, (i1+1)*self.nL)
            L = extractMaskedVec(self.L, mask=self.lMaxMask[I])
            d = extractMaskedVec(self.dataVector, mask=self.lMaxMask, I=I) / self.sUnit**2
            shot = extractMaskedVec(self.shotNoiseVector, mask=self.lMaxMask, I=I) / self.sUnit**2
            cov = extractMaskedMat(self.covMat, mask=self.lMaxMask, I=I)
            std = np.sqrt(np.diag(cov)) / self.sUnit**2
            #
            ax.plot(L*(1.+0.01*i1/self.nSS), std / d, '.-', lw=1, color=color)
            # move to next row
            i1 += 1

      #
      ax.set_xlim((90., 2.5e3))
      ax.legend(loc=1, fontsize='x-small', labelspacing=0., handlelength=0.4, borderaxespad=0.05)
      ax.set_xscale('log')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$\ell$')
      ax.set_ylabel(r'$\sigma\left( C_\ell^{\gamma\gamma} \right) / C_\ell^{\gamma\gamma}$')
      ax.set_title(r'Shear tomography')
      #
      fig.savefig(self.figurePath+"/sp2d_ss.pdf")
      if show:
         plt.show()
      fig.clf()
