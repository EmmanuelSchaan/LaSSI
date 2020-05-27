
   def plotOutlierPhotozRequirements(self, ICosmoPar, name):
      # Self-calibration of Gaussian photo-z: compare data sets
      Photoz, sCosmoGks, sGPhotozGks, sOutlierPhotozGks = self.varyOutlierPhotozPrior(self.lMaxMask, ICosmoPar, dzStd=0.002, szStd=0.003)
      Photoz, sCosmoGs, sGPhotozGs, sOutlierPhotozGs = self.varyOutlierPhotozPrior(self.lMaxMask+self.gsOnlyMask, ICosmoPar, dzStd=0.002, szStd=0.003)
      Photoz, sCosmoGsnonull, sGPhotozGsnonull, sOutlierPhotozGsnonull = self.varyOutlierPhotozPrior(self.lMaxMask+self.gsOnlyMask+self.noNullMask, ICosmoPar, dzStd=0.002. szStd=0.003)


      ###############################################################
      ###############################################################
      # Degradation in cosmological parameters

      # get the names of the cosmo params
      parCosmo = self.cosmoPar.extractParams(ICosmoPar, marg=False)

      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # fiducial prior
      ax.axvline(0.002, ls='--', lw=1, c='gray', label=r'Requirement')
      prop_cycle = plt.rcParams['axes.prop_cycle']
      colors = prop_cycle.by_key()['color']
      for iPar in range(parCosmo.nPar):
         ax.plot(Photoz/(self.nBins-1), sCosmoGs[iPar, :] / sCosmoGs[iPar, 0], c=colors[iPar], label=parCosmo.namesLatex[iPar])
      #
      ax.set_xscale('log', nonposx='clip')
      ax.legend(loc=2, labelspacing=0.1, frameon=False, handlelength=1.)
      #ax.set_ylabel(r'$\sigma_\text{Param} / \sigma_\text{Perfect photo-z}$')
      ax.set_ylabel(r'Degradation in posterior uncertainty')
      ax.set_xlabel(r'Prior on outlier fraction $c_{ij}$')
      #
      path = "/outlierphotozreq_deg_cosmo_"+name+".pdf"
      fig.savefig(self.figurePath+path, bbox_inches='tight')
      fig.clf()


      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # fiducial prior
      ax.axvline(0.002, ls='--', lw=1, c='gray', label=r'Requirement')
      prop_cycle = plt.rcParams['axes.prop_cycle']
      colors = prop_cycle.by_key()['color']
      for iPar in range(parCosmo.nPar):
         ax.plot(Photoz/(self.nBins-1), sCosmoGks[iPar,:] / sCosmoGs[iPar,:], c=colors[iPar], label=parCosmo.namesLatex[iPar])
      #
      ax.set_xscale('log', nonposx='clip')
      ax.legend(loc=2, labelspacing=0.1, frameon=False, handlelength=1.)
      ax.set_ylabel(r'$\sigma_\text{Param, gks} / \sigma_\text{Param}$')
      ax.set_xlabel(r'Prior on outlier fraction $c_{ij}$')
      #
      path = "/outlierphotozreq_cosmo_vs_gks_"+name+".pdf"
      fig.savefig(self.figurePath+path, bbox_inches='tight')
      fig.clf()


      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # fiducial prior
      ax.axvline(0.002, ls='--', lw=1, c='gray', label=r'Requirement')
      prop_cycle = plt.rcParams['axes.prop_cycle']
      colors = prop_cycle.by_key()['color']
      for iPar in range(parCosmo.nPar):
         ax.plot(Photoz/(self.nBins-1), sCosmoGsnonull[iPar,:] / sCosmoGs[iPar,:], c=colors[iPar], label=parCosmo.namesLatex[iPar])
      #
      ax.set_xscale('log', nonposx='clip')
      ax.legend(loc=2, labelspacing=0.1, frameon=False, handlelength=1.)
      ax.set_ylabel(r'$\sigma_\text{Param, no null} / \sigma_\text{Param}$')
      ax.set_xlabel(r'Prior on outlier fraction $c_{ij}$')
      #
      path = "/outlierphotozreq_cosmo_vs_gsnonull_"+name+".pdf"
      fig.savefig(self.figurePath+path, bbox_inches='tight')
      fig.clf()




      ###############################################################
      ###############################################################
      # Degradation in Gaussian photo-z
      
      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      # fiducial prior
      ax.axvline(0.002, ls='-', lw=1, color='gray', label=r'Requirement')
      ax.axhline(0.002, ls='-', lw=1, color='gray')
      #
      # photo-z shifts
      # add legend entry
      color = 'r'
      ax.plot([], [], color=color, label=r'$\delta z$')
      darkLight = 0.
      for iPar in range(self.nBins):
         color = darkerLighter(color, amount=darkLight)
         darkLight += -0.75 * 1./self.nBins
         ax.plot(Photoz/(self.nBins-1), sGPhotozGs[iPar,:]/sGPhotozGs[iPar,0], color=color)
         #ax.plot(Photoz, sGPhotozGks[iPar, :], color=color, ls='--')
         #ax.plot(Photoz, sGPhotozGsnonull[iPar, :], color=color, ls=':')
      #
      # photo-z scatter
      # add legend entry
      color = 'b'
      ax.plot([], [], color=color, label=r'$\sigma_z / (1+z)$')
      darkLight = 0.
      for iPar in range(self.nBins, 2*self.nBins):
         color = darkerLighter(color, amount=darkLight)
         darkLight += -0.75 * 1./self.nBins
         ax.plot(Photoz/(self.nBins-1), sGPhotozGs[iPar,:]/sGPhotozGs[iPar,0], color=color)
         #ax.plot(Photoz, sGPhotozGks[iPar, :], color=color, ls='--')
         #ax.plot(Photoz, sGPhotozGsnonull[iPar, :], color=color, ls=':')
      #
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.legend(loc=2, labelspacing=0.1, frameon=False, handlelength=1.)
      ax.set_ylabel(r'Degradation in Gaussian photo-z posterior')
      ax.set_xlabel(r'Prior on outlier fraction $c_{ij}$')
      #
      path = "/outlierphotozreq_deg_gphotoz_"+name+".pdf"
      fig.savefig(self.figurePath+path, bbox_inches='tight')
      fig.clf()


      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      # fiducial prior
      ax.axvline(0.002, ls='-', lw=1, color='gray', label=r'Requirement')
      #ax.axhline(0.002, ls='-', lw=1, color='gray')
      #
      # photo-z shifts
      # add legend entry
      color = 'r'
      ax.plot([], [], color=color, label=r'$\delta z$')
      darkLight = 0.
      for iPar in range(self.nBins):
         color = darkerLighter(color, amount=darkLight)
         darkLight += -0.75 * 1./self.nBins
         ax.plot(Photoz/(self.nBins-1), sGPhotozGks[iPar,:]/sGPhotozGs[iPar,:], color=color)
      #
      # photo-z scatter
      # add legend entry
      color = 'b'
      ax.plot([], [], color=color, label=r'$\sigma_z / (1+z)$')
      darkLight = 0.
      for iPar in range(self.nBins, 2*self.nBins):
         color = darkerLighter(color, amount=darkLight)
         darkLight += -0.75 * 1./self.nBins
         ax.plot(Photoz/(self.nBins-1), sGPhotozGks[iPar,:]/sGPhotozGs[iPar,:], color=color)
      #
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.legend(loc=2, labelspacing=0.1, frameon=False, handlelength=1.)
      ax.set_ylabel(r'$\sigma_\text{G photo-z, gks} / \sigma_\text{G photo-z}$')
      ax.set_xlabel(r'Prior on outlier fraction $c_{ij}$')
      #
      path = "/outlierphotozreq_gphotoz_vs_gks_"+name+".pdf"
      fig.savefig(self.figurePath+path, bbox_inches='tight')
      fig.clf()


      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      # fiducial prior
      ax.axvline(0.002, ls='-', lw=1, color='gray', label=r'Requirement')
      #ax.axhline(0.002, ls='-', lw=1, color='gray')
      #
      # photo-z shifts
      # add legend entry
      color = 'r'
      ax.plot([], [], color=color, label=r'$\delta z$')
      darkLight = 0.
      for iPar in range(self.nBins):
         color = darkerLighter(color, amount=darkLight)
         darkLight += -0.75 * 1./self.nBins
         ax.plot(Photoz/(self.nBins-1), sGPhotozGsnonull[iPar,:]/sGPhotozGs[iPar,:], color=color)
      #
      # photo-z scatter
      # add legend entry
      color = 'b'
      ax.plot([], [], color=color, label=r'$\sigma_z / (1+z)$')
      darkLight = 0.
      for iPar in range(self.nBins, 2*self.nBins):
         color = darkerLighter(color, amount=darkLight)
         darkLight += -0.75 * 1./self.nBins
         ax.plot(Photoz/(self.nBins-1), sGPhotozGsnonull[iPar,:]/sGPhotozGs[iPar, :], color=color)
      #
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.legend(loc=2, labelspacing=0.1, frameon=False, handlelength=1.)
      ax.set_ylabel(r'$\sigma_\text{G photo-z, no null} / \sigma_\text{G photo-z}$')
      ax.set_xlabel(r'Prior on outlier fraction $c_{ij}$')
      #
      path = "/outlierphotozreq_gphotoz_vs_gsnonull_"+name+".pdf"
      fig.savefig(self.figurePath+path, bbox_inches='tight')
      fig.clf()



      ###############################################################
      ###############################################################
      # Outlier photo-z posteriors
      
      nCij = self.nBins * (self.nBins - 1)

      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      ax.axvline(self.photoZPar.priorStd[-1], ls='-', lw=1, color='gray', label=r'Fiducial')
      ax.axhline(self.photoZPar.priorStd[-1], ls='-', lw=1, color='gray')
      #ax.plot(Photoz/(self.nBins-1), Photoz/(self.nBins-1), c='gray', ls=':', lw=1, label=r'Posterior=Prior')
      #
      color = 'r'
      darkLight = 0.
      for iPar in range(nCij):
         color = 'r'
         color = darkerLighter(color, amount=darkLight)
         darkLight += -0.75 * 1./nCij
         ax.plot(Photoz/(self.nBins-1), sOutlierPhotozGs[iPar, :], color=color, alpha=0.3)
         #ax.plot(Photoz/(self.nBins-1), sOutlierPhotozGks[iPar, :], color=color, alpha=0.3, ls='--')
         #ax.plot(Photoz/(self.nBins-1), sOutlierPhotozGsnonull[iPar, :], color=color, alpha=0.3, ls=':')
      #
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.legend(loc=2)
      ax.set_ylabel(r'Posterior on outlier fraction $c_{ij}$')
      ax.set_xlabel(r'Prior on outlier fraction $c_{ij}$')
      #
      path = "/outlierphotozreq_outlierphotoz_"+name+".pdf"
      fig.savefig(self.figurePath+path, bbox_inches='tight')
      fig.clf()



      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      ax.axvline(self.photoZPar.priorStd[-1], ls='-', lw=1, color='gray', label=r'Fiducial')
      ax.axhline(self.photoZPar.priorStd[-1], ls='-', lw=1, color='gray')
      #ax.plot(Photoz/(self.nBins-1), Photoz/(self.nBins-1), c='gray', ls=':', lw=1, label=r'Posterior=Prior')
      #
      color = 'r'
      darkLight = 0.
      for iPar in range(nCij):
         color = 'r'
         color = darkerLighter(color, amount=darkLight)
         darkLight += -0.75 * 1./nCij
         ax.plot(Photoz/(self.nBins-1), sOutlierPhotozGks[iPar,:]/ sOutlierPhotozGs[iPar,:], color=color, alpha=0.3)
      #
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.legend(loc=2)
      ax.set_ylabel(r'Posterior on outlier fraction $c_{ij}$')
      ax.set_ylabel(r'$\sigma_{c_{ij}, \text{gks}} / \sigma_{c_{ij}}$')
      ax.set_xlabel(r'Prior on outlier fraction $c_{ij}$')
      #
      path = "/gphotozreq_outlierphotoz_vs_gks_"+name+".pdf"
      fig.savefig(self.figurePath+path, bbox_inches='tight')
      fig.clf()


      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      ax.axvline(self.photoZPar.priorStd[-1], ls='-', lw=1, color='gray', label=r'Fiducial')
      ax.axhline(self.photoZPar.priorStd[-1], ls='-', lw=1, color='gray')
      #ax.plot(Photoz/(self.nBins-1), Photoz/(self.nBins-1), c='gray', ls=':', lw=1, label=r'Posterior=Prior')
      #
      color = 'r'
      darkLight = 0.
      for iPar in range(nCij):
         color = 'r'
         color = darkerLighter(color, amount=darkLight)
         darkLight += -0.75 * 1./nCij
         ax.plot(Photoz/(self.nBins-1), sOutlierPhotozGsnonull[iPar,:]/ sOutlierPhotozGs[iPar,:], color=color, alpha=0.3)
      #
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.legend(loc=2)
      ax.set_ylabel(r'Posterior on outlier fraction $c_{ij}$')
      ax.set_ylabel(r'$\sigma_{c_{ij}, \text{no null}} / \sigma_{c_{ij}}$')
      ax.set_xlabel(r'Prior on outlier fraction $c_{ij}$')
      #
      path = "/gphotozreq_outlierphotoz_vs_gsnonull_"+name+".pdf"
      fig.savefig(self.figurePath+path, bbox_inches='tight')
      fig.clf()






