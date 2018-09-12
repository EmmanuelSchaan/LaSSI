   def plotPowerSpectra(self):

      '''
      # gg: all on same plot
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      Colors = plt.cm.jet(1.*np.arange(self.nBins)/(5.))
      i1 = 0
      for iBin1 in range(self.nBins):
         # add entry to caption
         color = Colors[iBin1]
         ax.plot([], [], c=color, label=r'$\langle g_{i} g_{i+'+str(iBin1)+r'}\rangle$')
         for iBin2 in range(iBin1, self.nBins):
            d = self.dataVector[i1*self.nL:(i1+1)*self.nL]
            std = np.sqrt(np.diag(self.covMat[i1*self.nL:(i1+1)*self.nL, i1*self.nL:(i1+1)*self.nL]))
            #
            color = Colors[iBin2-iBin1]
            ax.errorbar(self.L*(1.+0.01*i1/self.nGG), d, yerr=std, ls='-', lw=2, elinewidth=1.5, marker='.', markersize=2, color=color)
            # move to next row
            i1 += 1
      #
      ax.legend(loc=1, labelspacing=0.05, handlelength=0.4, borderaxespad=0.01)
      ax.set_xscale('log')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$\ell$')
      ax.set_ylabel(r'$C_\ell^{gg}$')
      '''
      
      
      # gg: panels
      Colors = plt.cm.autumn(1.*np.arange(self.nBins)/(self.nBins-1.))
      #
      fig=plt.figure(0)
      gs = gridspec.GridSpec(3, 1)#, height_ratios=[1, 1, 1])
      gs.update(hspace=0.)
      
      # auto
      ax0=plt.subplot(gs[0])
      i1 = 0
      for iBin1 in range(self.nBins):
         d = self.dataVector[i1*self.nL:(i1+1)*self.nL]
         std = np.sqrt(np.diag(self.covMat[i1*self.nL:(i1+1)*self.nL, i1*self.nL:(i1+1)*self.nL]))
         #
         color = Colors[iBin1]
         ax0.errorbar(self.L, d, yerr=std, ls='-', lw=2, elinewidth=1.5, marker='.', markersize=2, color=color)
         # advance counter in data vector
         i1 += self.nBins - iBin1
      #
      ax0.set_xscale('log')
      ax0.set_yscale('log', nonposy='clip')
      plt.setp(ax0.get_xticklabels(), visible=False)
      #
      ax0.set_title(r'Clustering: $C_\ell^{gg}$')
      ax0.set_ylabel(r'$\langle g_i g_i\rangle$', fontsize=18)

      # cross i,i+1
      ax1=plt.subplot(gs[1])
      i1 = 1
      for iBin1 in range(self.nBins-1):
         d = self.dataVector[i1*self.nL:(i1+1)*self.nL]
         std = np.sqrt(np.diag(self.covMat[i1*self.nL:(i1+1)*self.nL, i1*self.nL:(i1+1)*self.nL]))
         #
         color = Colors[iBin1]
         ax1.errorbar(self.L, d, yerr=std, ls='-', lw=2, elinewidth=1.5, marker='.', markersize=2, color=color)
         # advance counter in data vector
         i1 += self.nBins - iBin1
      #
      ax1.set_xscale('log')
      ax1.set_yscale('log', nonposy='clip')
      plt.setp(ax1.get_xticklabels(), visible=False)
      #
      ax1.set_ylabel(r'$\langle g_i g_{i+1}\rangle$', fontsize=18)

      # cross i,i+2
      ax2=plt.subplot(gs[2])
      i1 = 2
      for iBin1 in range(self.nBins-2):
         d = self.dataVector[i1*self.nL:(i1+1)*self.nL]
         std = np.sqrt(np.diag(self.covMat[i1*self.nL:(i1+1)*self.nL, i1*self.nL:(i1+1)*self.nL]))
         #
         color = Colors[iBin1]
         ax2.errorbar(self.L, d, yerr=std, ls='-', lw=2, elinewidth=1.5, marker='.', markersize=2, color=color)
         # advance counter in data vector
         i1 += self.nBins - iBin1
      #
      ax2.set_xscale('log')
      ax2.set_yscale('log', nonposy='clip')
      #
      ax2.set_ylabel(r'$\langle g_i g_{i+2}\rangle$', fontsize=18)
      ax2.set_xlabel(r'$\ell$')
      


      # gs
      Colors = plt.cm.winter(1.*np.arange(self.nBins)/(self.nBins-1.))
      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      i1 = self.nGG
      for iBin1 in range(self.nBins):
         color = Colors[iBin1]
         for iBin2 in range(self.nBins):
            d = self.dataVector[i1*self.nL:(i1+1)*self.nL]
            std = np.sqrt(np.diag(self.covMat[i1*self.nL:(i1+1)*self.nL, i1*self.nL:(i1+1)*self.nL]))
            ax.errorbar(self.L*(1.+0.01*i1/self.nGS), d, yerr=std, ls='-', lw=1, elinewidth=1.5, marker='.', markersize=2, color=color)# label=r'$\langle g_{'+str(iBin1)+'} \gamma_{'+str(iBin2)+r'}\rangle$')
            # move to next row
            i1 += 1
      #
      ax.legend(loc=1)
      ax.set_xscale('log')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$\ell$')
      ax.set_ylabel(r'$C_\ell^{g\gamma}$')
      ax.set_title(r'Galaxy - galaxy lensing')

      
      # ss: all on same plot
      Colors = plt.cm.jet(1.*np.arange(self.nBins)/(self.nBins-1.))
      fig=plt.figure(2)
      ax=fig.add_subplot(111)
      #
      i1 = self.nGG + self.nGS
      for iBin1 in range(self.nBins):
         # add entry to caption
         color = Colors[iBin1]
         ax.plot([], [], c=color, label=r'$\langle\gamma_{i} \gamma_{i+'+str(iBin1)+r'} \rangle $')
         for iBin2 in range(iBin1, self.nBins):
            d = self.dataVector[i1*self.nL:(i1+1)*self.nL]
            #
            color = Colors[iBin2-iBin1]
            #
            std = np.sqrt(np.diag(self.covMat[i1*self.nL:(i1+1)*self.nL, i1*self.nL:(i1+1)*self.nL]))
            ax.errorbar(self.L*(1.+0.01*i1/self.nSS), d, yerr=std, ls='-', lw=1, elinewidth=1.5, marker='.', markersize=2, color=color)#, label=r'$\gamma_{'+str(iBin1)+'} \gamma_{'+str(iBin2)+'}$')
            # move to next row
            i1 += 1
      #
      ax.legend(loc=1, labelspacing=0.05, handlelength=0.4, borderaxespad=0.01)
      ax.set_xscale('log')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$\ell$')
      ax.set_ylabel(r'$C_\ell^{\gamma\gamma}$')
      ax.set_title(r'Shear tomography')
      
      
      plt.show()
