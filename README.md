# Xic0Analysis
Xic0 analysis code

//-----------------------------------------------

Update at Mar. 24, 2021 (kimc)

	- AliNormalizationCounter update
		a. Moved "StoreEvent" invoking point:
		   after "pileup rejection and Vtx->GetNContributors() > 1" and before |zVtx| < 10
		b. To match the events being analyzed and normalization factor obatained via ANC

//-----------------------------------------------

Update at Mar. 8, 2021 (kimc)

	- Save a triggerbit for data only: for MC, 0 will be assigned
	- Cleaned up some obsolete objects (items commented out, SPD histograms, etc)
	- Bug (?) found: all trees' float type were "/f (24 bit truncated float)" -> updated to "/F (32 bit float)"
	- Added hard coded run number cut to accept valid runs only (RUN2: 252000 - 295000)
	- Added a boolean variable "fINEL" under EventTree: true for "INEL > 0"
		a. data: AliPPVsMultUtils::IsINELgtZERO(event) = kINT7 + SPD tracklets >= 1 + |eta| < 1"
		b. MC: require "IsPhysicalPrimary + IsCharged + |eta| < 1"
		   * Make sure to run "general purpose MC" during the train run, when you deal with "INEL > 0" !
	- AliNormalizationCounter update
		a. For AliAnalysisTaseSE...,
		   a-1. Added 4 more counter object: MB_0to100, MB_0p1to30, MB_30to100, and HMV0_0to0p1
		        * Cannot use one common counter since official method accepts only integer multiplicity
		   a-2. New AliRDHFCutsXictoeleXifromAODtracks object (fEvtCuts_HMV0) added for HMV0
		b. For AddTask..., added above counter objects as container 8-11
	- macro/kimc/Xic0AnaMakeRoot.C also updated accordingly

//-----------------------------------------------

Update at Feb. 10, 2021 (kimc)

	- Added analysis macros for final cross section estimation
		a. Xi0cAnaMakeRoot.C (updated Xi topology cut after WDK)
		b. Xi0cAnaFunction.C (collection of xSec analysis subroutines)
		c. Xi0cAnaExecute.C (steering macro)

//-----------------------------------------------

Update at Jan. 25, 2021 (kimc)

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

