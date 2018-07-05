import fisher_lsst
reload(fisher_lsst)
from fisher_lsst import *

##################################################################################
# define parameters

#params = GalaxyBiasParams(nBins=1)
#params.addParams(ShearMultBiasParams(nBins=1))
#params.addParams(PhotoZParams(nBins=1))
#
#params.plotParams()


##################################################################################

#u = UnivPlanck15()

fisherLsst = FisherLsst(nBins=2, nL=20, fsky=1., save=False)


#fisherLsst.plotDndz()
#fisherLsst.plotPowerSpectra()
#fisherLsst.plotCovMat()

