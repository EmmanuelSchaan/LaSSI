   def fSingleSource(self, a, aSource):
      """lensing projection kernel
      for single source
      a is dimless, Wlensing(a) in (h Mpc^-1)
      """
      # Expression valid for curved cosmology
      z = 1./a - 1.
      zSource = 1./aSource - 1.
      d_L = self.U.bg.comoving_transverse_distance(z)
      d_S = self.U.bg.comoving_transverse_distance(zSource)
      chi_L = self.U.bg.comoving_distance(z)
      chi_S = self.U.bg.comoving_distance(zSource)
      d_LS = self.U.fK(chi_S - chi_L)
      
      wlensing = 1.5 * (100./3.e5)**2 * self.U.bg.Omega0_m
      wlensing *= d_L/a * d_LS / d_S
      return wlensing


   def fForInterp(self, a):
      integrand = lambda a_s: self.fdpdz(1./a_s-1.) /a_s**2 * self.fSingleSource(a, a_s)
      result *= integrate.quad(integrand, self.aMin, a, epsabs=0, epsrel=1.e-2)[0]
      return result
