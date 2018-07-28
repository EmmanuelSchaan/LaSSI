# Parameter combinations to look at:
base LCDM
base LCDM + massive neutrinos
base LCDM + massive neutrinos + curvature
base LCDM + wCDM
base LCDM + massive neutrinos + wCDM (agreed, this is a bit stupid, but LSST people want it)
base LCDM + massive neutrinos + curvature + wCDM (also a bit stupid)




# base LCDM
names = np.array(['Omega_cdm', 'Omega_b', 'A_s', 'n_s', 'h', 'tau_reio'])
namesLatex = np.array([r'$\Omega^0_\text{CDM}$', r'$\Omega^0_\text{b}$', r'$A_\text{s}$', r'$n_\text{s}$', r'$h_0$', r'$\tau$'])
fiducial = np.array([0.267, 0.0493, 2.3e-9, 0.9624, 0.6712], 0.06])

# massive neutrinos
names = np.concatenate((names, np.array(['m_ncdm'])))
namesLatex = np.concatenate((namesLatex, np.array([r'$M_\nu$'])))
fiducial = np.concatenate((fiducial, np.array([0.06])))

# wCDM
names = np.concatenate((names, np.array(['w0_fld', 'wa_fld'])))
namesLatex = np.concatenate((namesLatex, np.array([r'$w_0$', r'$w_a$'])))
fiducial = np.concatenate((fiducial, np.array([-1., 0.])))

# curvature
names = np.concatenate((names, np.array(['Omega_k'])))
namesLatex = np.concatenate((namesLatex, np.array([r'$\Omega_\text{k}$'])))
fiducial = np.concatenate((fiducial, np.array([0.])))
