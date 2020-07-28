
import universe
reload(universe)
from universe import *


u = UnivPlanck15()

def dndz(z, zMean=1., zWidth=0.15):
   result = np.exp(-0.5 * (z-zMean)**2 / zWidth**2)
   result /= np.sqrt(2. * np.pi * zWidth**2)
   return result

Z = np.linspace(0., 4., 101)
Dndz = np.array(map(dndz, Z))

plt.plot(Z, Dndz)
plt.show()

def b(z, a=1.):
   return 1. / u.bg.scale_independent_growth_factor(z)**a


def zEff(a=1.):

   def integrand(z): return b(z, a=a) * dndz(z) * z
   result = integrate.quad(integrand, 0., 4., epsabs=0., epsrel=1.e-3)[0]

   def integrand(z): return b(z, a=a) * dndz(z)
   result /= integrate.quad(integrand, 0., 4., epsabs=0., epsrel=1.e-3)[0]

   return result

print "a, sEff:", 1., zEff(a=1.)
print "a, zEff:",2., zEff(a=2.)
print "Relative difference in zEff:", zEff(a=2) / zEff(a=1) - 1.
