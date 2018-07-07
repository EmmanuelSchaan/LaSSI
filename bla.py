# clustering
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
#for iCosmo in reversed(range(7)):
for iCosmo in [0,1,2,6,3,5,4, 37]:
   #
   # gg
   I = datav_gg[0,:]>0.
   ax.fill_between(Ell[I], np.min(dlndatav_dlnparam_gg[:,I, iCosmo], axis=0), np.max(dlndatav_dlnparam_gg[:,I, iCosmo], axis=0), edgecolor='', facecolor=Colors[iCosmo], alpha=0.8)
   ax.plot([], [], color=Colors[iCosmo], alpha=0.8, linewidth=20, label=Cosmo[iCosmo])
   #
   ax.plot([], [], color=Colors[iCosmo], alpha=0.8, linewidth=20, label=Cosmo[iCosmo])
#
# error bars
for i in range(Ng):
   I = datav_gg[i,:]>0.
   ax.errorbar(Ell[I]*1.1**(i/(Ng-1.)), 0.*Ell[I], yerr=np.sqrt(np.diagonal(cov_gg[i*Nell:(i+1)*Nell,i*Nell:(i+1)*Nell])[I])/datav_gg[i,I], fmt='.-', c='k')
#
ax.axhline(0., c='k', ls='-')
#
#ax.legend(loc=4, ncol=2, labelspacing=0.1, frameon=False)
ax.set_ylim((-3., 3.))
ax.set_xlim((2.e1, 5.e2))
ax.set_xscale('log')
ax.set_xlabel(r'$\ell$', fontsize=24)
ax.set_ylabel(r'$d \ln C_\ell^{\; gg}/ d \ln \# $', fontsize=24)
