# Xic0Analysis
Xic0 analysis code

//-----------------------------------------------

Update at Jan.25, 2021 (kimc)

	- Added 2nd production macro (Xic0AnaMakeRoot.C) for analysis
		a. A few modifications compared to original, especially adding > 0 on several denominators
		b. Tested on Lego train output available at Jan. 25

//-----------------------------------------------

Update at Dec. 18 (kimc)

	- Added TH2 histograms (fired trigger vs. mult) for nomalization purpose for data
		a. Histograms:
			a-1. hNorm_multV0: uses V0M centrality (fCentrality)
			a-2. hNorm_multSPD: uses SPD centrality (fCentrality)
		b. Conditions:
			b-1. Triggers: all (any trigger, 0), kINT7 (1), kHMV0 (2), kHMSPD (3), and 'kHMV0 || kHMSPD' (4)
			b-2. 1,000 bins assigned for multiplicity axis (bin width 0.1)
			b-3. Histograms are being filled after trigger check finished

		* 2nd update: histogram fillup moment changed (before: trigger only -> now: trigger + zVtx + PV)

//-----------------------------------------------

Updates at Oct. 8 (kimc)

	1. AliAnalysisTaskSEXi0c... is updated by Dr. Bok for compatability with pPb dataset
		a. pp: Get multiplicity from V0M
		b. PPb: Get multiplicity from V0A
	2. Macro for LEGO train running: AddTaskXic0Semileptonic.C
		a. Added High multiplicity trigger option (bool UseTrifHM)
			a-1. By default it is disabled: only kINY7 will be used
			a-2. Add kHMV0 and kHMSPD triggers f enabled
		b. Tested at Grid: no problem

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

