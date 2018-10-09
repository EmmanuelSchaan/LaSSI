      for iPhotoz in range(nPhotoz):
         photoz = Photoz[iPhotoz]
         # update the photo-z priors
         newPhotoZPar = PhotoZParams(nBins=self.nBins, dzFid=0., szFid=0.05, dzStd=photoz, szStd=photoz*1.5)
         # update the full parameter object
         newPar = self.cosmoPar.copy()
         newPar.addParams(self.galaxyBiasPar)
         newPar.addParams(self.shearMultBiasPar)
         newPar.addParams(newPhotoZPar)
         # get the new posterior Fisher matrix, including the prior
         newPar.fisher += self.fisherData
         
         # Extract full uncertainties
         sFull[:,iPhotoz] = newPar.paramUncertainties(marg=True)
         
         # Extract parameter combinations:
         #
         # Full: LCDM + Mnu + curv + w0,wa
         # reject unwanted cosmo params
         I = self.cosmoPar.IFull + range(self.cosmoPar.nPar, self.fullPar.nPar)
         par = newPar.extractParams(I, marg=False)
         if iPhotoz==0 or iPhotoz==nPhotoz-1:
            par.printParams(path=self.figurePath+"/posterior_full_photozprior_"+floatExpForm(photoz)+".txt")
         # cosmology
         parCosmoFull = par.extractParams(range(len(self.cosmoPar.IFull)), marg=True)
         sCosmoFull[:, iPhotoz] = parCosmoFull.paramUncertainties(marg=True)
         # photo-z
         parPhotozFull = par.extractParams(range(-self.photoZPar.nPar, 0), marg=True)
         sPhotozFull[:, iPhotoz] = parPhotozFull.paramUncertainties(marg=True)
         #
         # LCDM + Mnu
         # reject unwanted cosmo params
         I = self.cosmoPar.ILCDMMnu + range(self.cosmoPar.nPar, self.fullPar.nPar)
         par = newPar.extractParams(I, marg=False)
         if iPhotoz==0 or iPhotoz==nPhotoz-1:
            par.printParams(path=self.figurePath+"/posterior_ldcmmnu_photozprior_"+floatExpForm(photoz)+".txt")
         # cosmology
         parCosmoLCDMMnu = par.extractParams(range(len(self.cosmoPar.ILCDMMnu)), marg=True)
         sCosmoLCDMMnu[:, iPhotoz] = parCosmoLCDMMnu.paramUncertainties(marg=True)
         # photo-z
         parPhotozLCDMMnu = par.extractParams(range(-self.photoZPar.nPar, 0), marg=True)
         sPhotozLCDMMnu[:, iPhotoz] = parPhotozLCDMMnu.paramUncertainties(marg=True)
         #
         # LCDM + Mnu + w0,wa
         # reject unwanted cosmo params
         I = self.cosmoPar.ILCDMMnuW0Wa + range(self.cosmoPar.nPar, self.fullPar.nPar)
         par = newPar.extractParams(I, marg=False)
         if iPhotoz==0 or iPhotoz==nPhotoz-1:
            par.printParams(path=self.figurePath+"/posterior_lcdmmnuw0wa_photozprior_"+floatExpForm(photoz)+".txt")
         # cosmology
         parCosmoLCDMMnuW0Wa = par.extractParams(range(len(self.cosmoPar.ILCDMMnuW0Wa)), marg=True)
         sCosmoLCDMMnuW0Wa[:, iPhotoz] = parCosmoLCDMMnuW0Wa.paramUncertainties(marg=True)
         # photo-z
         parPhotozLCDMMnuW0Wa = par.extractParams(range(-self.photoZPar.nPar, 0), marg=True)
         sPhotozLCDMMnuW0Wa[:, iPhotoz] = parPhotozLCDMMnuW0Wa.paramUncertainties(marg=True)
         #
         # LCDM + Mnu + curv
         # reject unwanted cosmo params
         I = self.cosmoPar.ILCDMMnuCurv + range(self.cosmoPar.nPar, self.fullPar.nPar)
         par = newPar.extractParams(I, marg=False)
         if iPhotoz==0 or iPhotoz==nPhotoz-1:
            par.printParams(path=self.figurePath+"/posterior_lcdmmnucurv_photozprior_"+floatExpForm(photoz)+".txt")
         # cosmology
         parCosmoLCDMMnuCurv = par.extractParams(range(len(self.cosmoPar.ILCDMMnuCurv)), marg=True)
         sCosmoLCDMMnuCurv[:, iPhotoz] = parCosmoLCDMMnuCurv.paramUncertainties(marg=True)
         # photo-z
         parPhotozLCDMMnuCurv = par.extractParams(range(-self.photoZPar.nPar, 0), marg=True)
         sPhotozLCDMMnuCurv[:, iPhotoz] = parPhotozLCDMMnuCurv.paramUncertainties(marg=True)
         #
         # LCDM + w0,wa
         # reject unwanted cosmo params
         I = self.cosmoPar.ILCDMW0Wa + range(self.cosmoPar.nPar, self.fullPar.nPar)
         par = newPar.extractParams(I, marg=False)
         if iPhotoz==0 or iPhotoz==nPhotoz-1:
            par.printParams(path=self.figurePath+"/posterior_lcdmw0wa_photozprior_"+floatExpForm(photoz)+".txt")
         # cosmology
         parCosmoLCDMW0Wa = par.extractParams(range(len(self.cosmoPar.ILCDMW0Wa)), marg=True)
         sCosmoLCDMW0Wa[:, iPhotoz] = parCosmoLCDMW0Wa.paramUncertainties(marg=True)
         # photo-z
         parPhotozLCDMW0Wa = par.extractParams(range(-self.photoZPar.nPar, 0), marg=True)
         sPhotozLCDMW0Wa[:, iPhotoz] = parPhotozLCDMW0Wa.paramUncertainties(marg=True)
         #
         # LCDM + w0,wa + curvature
         # reject unwanted cosmo params
         I = self.cosmoPar.ILCDMW0WaCurv + range(self.cosmoPar.nPar, self.fullPar.nPar)
         par = newPar.extractParams(I, marg=False)
         if iPhotoz==0 or iPhotoz==nPhotoz-1:
            par.printParams(path=self.figurePath+"/posterior_lcdmw0wacurv_photozprior_"+floatExpForm(photoz)+".txt")
         # cosmology
         parCosmoLCDMW0WaCurv = par.extractParams(range(len(self.cosmoPar.ILCDMW0WaCurv)), marg=True)
         sCosmoLCDMW0WaCurv[:, iPhotoz] = parCosmoLCDMW0WaCurv.paramUncertainties(marg=True)
         # photo-z
         parPhotozLCDMW0WaCurv = par.extractParams(range(-self.photoZPar.nPar, 0), marg=True)
         sPhotozLCDMW0WaCurv[:, iPhotoz] = parPhotozLCDMW0WaCurv.paramUncertainties(marg=True)
