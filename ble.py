   def plotContours(self, IPar=None):
      '''Show confidence ellipses.
      IPar (optional): indices of parameters to show
      '''
      if IPar is None:
         IPar = range(self.nPar)
      nPar = len(IPar)


#ax1=plt.subplot(gs[1], sharex=ax0)




from matplotlib import patches











invFisher = np.array([[1., 0.7, 0.0],
               [0.7, 2., -0.99*np.sqrt(2.*3.)],
               [0.0, -0.99*np.sqrt(2.*3.), 3.]])

fiducial = np.array([-10., -20., -30.])

nPar = 3

#invFisher = np.linalg.inv(Fisher)

namesLatex = np.array(['a', 'b', 'c'])


# plot range, in sigmas
lim = 4.


#colors = ['r', 'g', 'b']

color='#E10014'













def plotContours

   fig=plt.figure(0)
   gs = gridspec.GridSpec(nPar, nPar)#, height_ratios=[1, 1, 1])
   gs.update(hspace=0.)
   gs.update(wspace=0.)

   # loop over each column
   for j in range(nPar):
      # j=i case: just plot the Gaussian pdf
      ax=plt.subplot(gs[j,j])
      
      # define Gaussian with correct mean and variance
      mean = fiducial[j]
      sigma = np.sqrt(invFisher[j,j])
      X = np.linspace(mean - lim*sigma, mean + lim*sigma, 201)
      Y = np.exp(-0.5*(X-mean)**2 / sigma**2)
      
      # plot it
      ax.plot(X, Y, color)
      
      # 68-95% confidence regions
      ax.plot([mean,mean], [0,1], color='b')
      ax.fill_between(X, 0., Y, where=np.abs(X-mean)<sigma, color=color, alpha=0.7)
      ax.fill_between(X, 0., Y, where=np.abs(X-mean)<2.*sigma, color=color, alpha=0.3)
      
      # print the marginalized constraint
      titleStr = namesLatex[j]+r'$ = '+str(np.round(fiducial[j], 2))+r' \pm'+str(np.round(sigma, 2))+r'$'
      ax.set_title(titleStr, fontsize=12)
      
      # always remove y ticks
      ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
      
      # remove x ticks, except if last parameter
      if j<nPar-1:
         ax.tick_params(axis='x', which='both', top=False, bottom=False, labelbottom=False)

      # limits
      ax.set_xlim((mean - lim*sigma, mean + lim*sigma))
      ax.set_ylim((0., 1.1))


      # loop over each row i>j
      for i in range(j+1, nPar):
         ax=plt.subplot(gs[i,j], sharex=ax)
         
         # plot the fiducial value
         ax.plot([fiducial[j]], [fiducial[i]], 'b.')

         # means and cov
         meanx = fiducial[j]
         meany = fiducial[i]
         #
         sx = np.sqrt(invFisher[j,j])
         sy = np.sqrt(invFisher[i,i])
         sxy = invFisher[i,j]
         

         # covariance contour
         x = np.linspace(meanx - lim*sx, meanx + lim*sx, 501)
         y = np.linspace(meany - lim*sy, meany + lim*sy, 501)
         X, Y = np.meshgrid(x, y)
         Z = bivariate_normal(X, Y, sigmax=sx, sigmay=sy, mux=meanx, muy=meany, sigmaxy=sxy)
         norm = bivariate_normal(0., 0., sigmax=sx, sigmay=sy, mux=0., muy=0., sigmaxy=sxy)
         Z = np.log(norm/Z)   # returns the chi squared

         # confidence level contours
         conf_level = np.array([0.68, 0.95])
         chi2 = -2. * np.log(1. - conf_level)
         ax.contourf(X,Y,Z, [0., chi2[0]], colors=color, alpha=0.7)
         ax.contourf(X,Y,Z, [0., chi2[1]], colors=color, alpha=0.3)
         
         # limits
         mean = fiducial[i]
         sigma = np.sqrt(invFisher[i,i])
         ax.set_ylim((mean - lim*sigma, mean + lim*sigma))
         
         # tick labels
         if j==0:
            ax.set_ylabel(namesLatex[i])
         else:
            ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
         #
         plt.setp(ax.get_xticklabels(), visible=False)

      # tick labels
      plt.setp(ax.get_xticklabels(), visible=True, rotation=45)
      ax.set_xlabel(namesLatex[j])

   plt.show()





