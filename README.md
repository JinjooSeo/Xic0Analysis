# Xic0Analysis
Xic0 analysis code
macros are not ready yet (need to confirm)

//-----------------------------------------------

Updates at Sep. 9 (kimc)

	1. Enforced EventTree sync to eXiTree
		a. EventTree is filled only if eXiTree is filled
		b. Added variable " fNeXiPair " under EventTree to enable event by event distinction:
		   if ( fNeXiPair != 0 ) it means those entries are duplication of same event
		c. Added " fVtxZ " under EventTree
			c-1. for data, fVtxZ = AliVVertex->GetZ()
			c-2. for MC, fVtxZ = AliAODMCHeader->GetVtxZ()

	2. Updated fWeightFit to get the parameters via functions

//-----------------------------------------------

Updates at Aug. 28 (kimc):

	1. Added high multiplicity triggers feature
		a. Total 3 triggers (kINT7, kHighMultV0, and kHighMultSPD) can be toggled on/off in steering macro
		b. If more than 1 trigger (> 1) enabled simultaneously, triggerbit will be saved under EventTree
		c. To proceed at least one trigger (ex. kINT7) should be on
		d. "hEevntNumbers" histogram is updated accordingly

	2. Added SPD centrality / # of SPDTracklets branches under EventTree
		a. EventTree->fCentrality : centrality calculated by V0M
		b. EventTree->fCentralSPD : centrality calculated by SPD
		c. EventTree->fNSPDTracklets : # of SPD tracklets of the event
		* Checked strong correlation in 'fCentrality:fCentralSPD' and in 'fCentralSPD:fNSPDTracklets'

	3. Removed fRunOffset feature from relevant histograms
		a. Affected histograms: NumOfEvtperRun, NumOfe, and NumOfXi
		b. Now absolute run # saved in those histograms
		c. Double checked working properly in later step (MakeROOTResultsXic0.C)

	4. Bug check: confirmed working properly in local/Grid, for both data and MC

