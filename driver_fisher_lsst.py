import fisher_lsst
reload(fisher_lsst)
from fisher_lsst import *

##################################################################################
# Forecast parameters

nBins = 2
nL = 20
fsky = 0.4

# cosmological parameters
massiveNu = False
wCDM = False
curvature = False

##################################################################################

# Parameter classes
cosmoPar = CosmoParams(massiveNu=massiveNu, wCDM=wCDM, curvature=curvature)
#cosmoPar.plotParams()
galaxyBiasPar = GalaxyBiasParams(nBins=nBins)
#galaxyBiasPar.plotParams()
shearMultBiasPar = ShearMultBiasParams(nBins=nBins)
#shearMultBiasPar.plotParams()
photoZPar = PhotoZParams(nBins=nBins)
#photoZPar.plotParams()


#params.addParams(ShearMultBiasParams(nBins=1))
#params.addParams(PhotoZParams(nBins=1))
#params.plotParams()

##################################################################################
# Fisher calculation

fisherLsst = FisherLsst(cosmoPar, galaxyBiasPar, shearMultBiasPar, photoZPar, nBins=nBins, nL=nL, fsky=0.4, name=None, save=False)


#fisherLsst.plotDndz()
#fisherLsst.plotPowerSpectra()
#fisherLsst.plotCovMat()
fisherLsst.plotDerivativeDataVector()


##################################################################################

#u = UnivPlanck15()
#
#tStart = time()
#N = int(1.e3)
#for i in range(N):
#   result = u.fPinterp(1.e-2, 5.)
#tStop = time()
#print (tStop-tStart)/N, "sec"

#X = np.linspace(0., 2.*np.pi, 101)
#Y = np.sin(X)
#f = interp1d(X, Y, kind='linear', bounds_error=False, fill_value=0.)
#
#tStart = time()
#N = int(1.e3)
#for i in range(N):
#   result = f(1.987)
#tStop = time()
#print (tStop-tStart)/N, "sec"
#
#print "wowowowow"
#tStart = time()
#N = int(1.e3)
#for i in range(1):
#   result = self.f(0.5*np.ones(N))
#tStop = time()
#print (tStop-tStart)/N, "sec"



#def darkerLighter(color, amount=0.):
#   """
#   Adapted from https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
#   Input can be matplotlib color string, hex string, or RGB tuple.
#   amount=0: color unchanged
#   amount=-1: returns white
#   amount=1: returns black
#
#   Examples:
#   >> lighten_color('g', 0.3)
#   >> lighten_color('#F034A3', 0.6)
#   >> lighten_color((.3,.55,.1), 0.5)
#   """
#   # force amount between 0. and 1.
#   amount = min(amount, 1.)
#   amount = max(amount, -1.)
#   # read the color
#   try:
#      c = mc.cnames[color]
#   except:
#      c = color
#   c = colorsys.rgb_to_hls(*mc.to_rgb(c))
#   # my Lagrange interpolation polynomial
#   newC1 = 0.5*amount*(amount-1.) - c[1]*(amount+1.)*(amount-1.)
#   return colorsys.hls_to_rgb(c[0], newC1, c[2])




#color = 'b'
#
#X = np.linspace(0., 1., 101)
#Y = np.sin(2.*np.pi*X)
#plt.plot(X, Y, lw=3, c=darkerLighter('b', amount = 1.))
#plt.show()


#X = np.linspace(-1., 1., 101)
#c = 0.3
#Y = 0.5*X*(X-1.) - c * (X+1.)*(X-1.)
#plt.plot(X, Y)
#plt.show()


