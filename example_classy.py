# neutrino masses: https://github.com/ThomasTram/iCLASS/blob/master/neutrinohierarchy.ipynb
# cosmo doc: https://lesgourg.github.io/class-tour/Narbonne.pdf

from classy import Class
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt





def computeNuMasses(mSum, normal=True):
   '''mSum: sum of neutrino masses in eV
   normal=True for normal hierarchy
   output: masses in eV
   '''
   dmsq_atm = 2.5e-3 # eV^2
   dmsq_solar = 7.6e-5 # eV^2
   if normal:
      f = lambda m0: m0 + np.sqrt(m0**2+dmsq_solar) + np.sqrt(m0**2+dmsq_solar+dmsq_atm) - mSum
      m0 = optimize.brentq(f , 0., mSum)
      result = np.array([m0, np.sqrt(m0**2+dmsq_solar), np.sqrt(m0**2+dmsq_solar+dmsq_atm)])
   else:
      f = lambda m0: m0 + np.sqrt(m0**2+dmsq_atm) + np.sqrt(m0**2+dmsq_atm+dmsq_solar) - mSum
      m0 = optimize.brentq(f , 0., mSum)
      result = np.array([m0, np.sqrt(m0**2+dmsq_atm), np.sqrt(m0**2+dmsq_atm+dmsq_solar)])
   return result



nuMasses = computeNuMasses(0.2, normal=True)
print "nu masses, sum= ", nuMasses, np.sum(nuMasses)


params = {
    'output': 'lCl tCl pCl mPk',
    'l_max_scalars': 2000,
    'lensing': 'yes',
    'A_s': 2.3e-9,
    'n_s': 0.9624,
    'h': 0.6711,
    'N_ur': 3.046,
    'omega_b': 0.022068,
    'omega_cdm': 0.12029,
    'P_k_max_1/Mpc': 3.0,
    'N_ncdm': 3,
    'm_ncdm': str(nuMasses[0])+','+str(nuMasses[1])+','+str(nuMasses[2]),
    'deg_ncdm': '1, 1, 1',
    'non linear': 'none,halofit',
    'z_pk': '0.,5.'
}


# Create an instance of the CLASS wrapper
cosmo = Class()

# Set the parameters to the cosmological code
cosmo.set(params)

# Run the whole code. Depending on your output, it will call the
# CLASS modules more or less fast. For instance, without any
# output asked, CLASS will only compute background quantities,
# thus running almost instantaneously.
# This is equivalent to the beginning of the `main` routine of CLASS,
# with all the struct_init() methods called.
cosmo.compute()

# Access the lensed cl
cls = cosmo.lensed_cl()

# Print on screen to see the output
print "cls contains:", cls.viewkeys()
#print cls
# It is a dictionnary that contains the fields: tt, te, ee, bb, pp, tp


#cosmo.get_current_derived_parameters(['Omega_Lambda'])

#cosmo.pk(k, z)   # function that returns P(k,z). Watch that there is no h factor anywhere
K = np.logspace(np.log10(1.e-4), np.log10(3.), 501, 10.)
f = lambda k: cosmo.pk(k, 0.)
Plin = np.array(map(f, K))

plt.loglog(K, Plin)
plt.show()









# Clean CLASS (the equivalent of the struct_free() in the `main`
# of CLASS. This step is primordial when running in a loop over different
# cosmologies, as you will saturate your memory very fast if you ommit
# it.
#cosmo.struct_cleanup()

# If you want to change completely the cosmology, you should also
# clean the arguments, otherwise, if you are simply running on a loop
# of different values for the same parameters, this step is not needed
#cosmo.empty()
