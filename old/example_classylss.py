# taken from: http://classylss.readthedocs.io/en/stable/examples.html
import classylss
import classylss.binding as CLASS
import numpy as np
from matplotlib import pyplot as plt



print "classylss version = ", classylss.__version__
print "CLASS version = ", classylss.class_version

# various ways of inputting parameters
#engine = CLASS.ClassEngine({'H0':70, 'Omega_m':0.31})
#default_engine = CLASS.ClassEngine() # default parameters
#params = classylss.load_ini('concise.ini')
#print(params)




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
#    'Omega_cdm': 0.32,
#    'Omega_k': 0.,
    'P_k_max_1/Mpc': 3.0,
#    'N_ncdm': 3,
#    'm_ncdm': str(nuMasses[0])+','+str(nuMasses[1])+','+str(nuMasses[2]),
#    'deg_ncdm': '1, 1, 1',
    'non linear': 'none,halofit',
    'z_max_pk': 100.,
#    'z_pk': '0.,5.'
}


engine = CLASS.ClassEngine(params)

#pre_params = classylss.load_precision('pk_ref.pre')
#print(pre_params)

#high_pre_engine = CLASS.ClassEngine(pre_params)

# initialize the background module
bg = CLASS.Background(engine)


# print out some cosmological parameters
print("h = ", bg.h)
print("Omega0_m = ", bg.Omega0_m)
print("Omega0_lambda = ", bg.Omega0_lambda)
print("Omega0_r = ", bg.Omega0_r)
print("Omega0_k = ", bg.Omega0_k)




plt.plot(z, bg.scale_independent_growth_factor(z), label=r"$D(z)$")
plt.plot(z, bg.scale_independent_growth_rate(z), label=r"$f(z)$")
plt.legend()
plt.xlabel(r"$z$")
plt.ylabel("growth factor/rate")
plt.show()



plt.loglog(1+z, bg.Omega_m(z), label=r"$\Omega_m$")
plt.loglog(1+z, bg.Omega_r(z), label=r"$\Omega_r$")
plt.legend()
plt.xlabel(r"$1+z$")
plt.ylabel("density parameter")
plt.show()


# initialize with proper output
engine = CLASS.ClassEngine({'output': 'dTk vTk mPk', 'non linear': 'halofit', 'P_k_max_h/Mpc' : 20., "z_max_pk" : 100.0})
sp = CLASS.Spectra(engine)
k = np.logspace(-2, 0, 100)
pk_nl = sp.get_pk(k=k, z=0.5)
pk_lin = sp.get_pklin(k=k, z=0.5)

plt.loglog(k, pk_nl, label='nonlinear')
plt.loglog(k, pk_lin, label='linear')
plt.legend()
plt.xlabel(r"$k$ $[h\mathrm{Mpc}^{-1}]$")
plt.ylabel(r"$P$ $[h^{-3} \mathrm{Mpc}^3]$")
plt.show()

z = np.linspace(0., 2., 1024)
plt.plot(1. + z, sp.sigma8_z(z))
plt.xlabel(r"$1+z$")
plt.ylabel(r"$\sigma_8(z)$")
plt.show()


transfer = sp.get_transfer(z=0)
print(transfer.dtype.names)
plt.subplot(211)

plt.plot(transfer['k'], transfer['d_tot'])
plt.ylabel("total density transfer")
plt.subplot(212)
plt.plot(transfer['k'], transfer['t_tot'])
plt.xlabel(r"$k$ $[h\mathrm{Mpc}^{-1}]$")
plt.ylabel("total velocity transfer")
plt.show()














