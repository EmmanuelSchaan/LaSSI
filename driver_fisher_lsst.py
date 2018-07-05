import fisher_lsst
reload(fisher_lsst)
from fisher_lsst import *

import parameters_fisher
reload(parameters_fisher)
from parameters_fisher import *

##################################################################################


#u = UnivPlanck15()

fisherLsst = FisherLsst(nBins=2, nL=20, fsky=1., save=False)

# show covariance matrix
#fisherLsst.plotCovMat()
#fisherLsst.plotPowerSpectra()
