

def plotEllBins(self):

   Z = np.linspace(1.e-4, 4., 201)

   def fLMax(z, kMaxG=0.3):
      '''Enforce that k < kMaxG = 0.3h/Mpc
      for clustering, following DESC SRD v1,
      where nonlin gal bias is a 10% effect on clustering.
      '''
      return kMaxG * self.u.bg.comoving_distance(z) - 0.5

   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   #
   ax.plot(Z, fLMax(Z, kMaxG=0.3), label=r'$k_\text{max}=0.3$ h/Mpc')
   #ax.plot(Z, fLMax(Z, kMaxG=0.5))
   #
   # show the actual ell and z for the various bins
   for iBin in range(self.nBins):
      zBin = self.w_g[iBin].zMean()
      lMaxG = fLMax(zBin, kMaxG=0.3)
      for iL in range(self.nL):
         if self.L[iL]<lMaxG:
            ax.plot(zBin, self.L[iL], 'ko')
         else:
            ax.plot(zBin, self.L[iL], 'ko', alpha=0.3)
   #
   ax.set_xlim((0., 4.))
   ax.set_ylim((0., 1100.))
   ax.set_xlabel(r'$z_\text{mean}$')
   ax.set_ylabel(r'$\ell$')
   #
   fig.savefig(self.figurePath+"/ell_bins.pdf", bbox_inches='tight')
   fig.clf()
