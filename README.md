# Xic0Analysis
Xic0 analysis code

Oct. 21 (kimc)

	- AliAnalysistask... .cxx
		a. Now updated pileup cut (MV) works properly for both MV and HMV0
		b. DEBUG: there was a bug in saving EventTree info -
		   it's likely substantial eXi pairs have been lost since last September (2020)

	- AddTaskXic0Semileptonic.C
		a. By default "UseTrigHM" is on
		b. Masked out HMSPD trigger

	- macro_kimc/Xi0cAnaMakeRoot.C
		now check pileup status by trigger by using variable EventTree->fPileup

//-----------------------------------------------

 Oct. 1 (kimc)

	- AliAnalysisTask...
		
		a.	Newly added pileup rejection cut (based on ' AliRDHFCutsXictoeleXifromAODtracks ') is problematic
			a-1. It works as adding a flag on the cut object:
			     e.g., fEvtCuts->SetOptPileup(AliRDHFCuts::kRejectMVPileupEvent);
			a-2. No return statement was invoked - I'm afraid pileup events will be accepted in this state
			a-3. A bool type member function exists to check the pileup status:
			     AliRDHFCuts->IsEventRejectedDueToPileup()

		b.	Two cut objects exist by triggers: fEvtCuts (kINT7) and fEvtCuts_HMV0 (kHighMultV0)
			b-1. Pileup status depends on the trigger used? YES, IT DOES
				 * I see many cases of " MB is ok, but HMV0 is pileup "
			b-2. Added new tree variable to judge later:
			     fEventTree->fPileup 0 (no pileup in MB/HMV0), 1 (MB ok, HMV0 pileup) and 2 (MB pileup, HMV0 ok)

		c.	Wrapped AliNormalizationCounter objects by pileup status boolean variable:
			the events being counted by ANC are the ones " NOT piled up "

	- /macro_kimc/Xi0cAnaMakeRoot.C

		a. New variable to deal with pileup status by trigger: EventTree->fPileup (L246)
		b. Invoke the condition during the event selection (L847)

//-----------------------------------------------

 Sep. 27 (kimc)

	- macro_kimc/Xic0AnaMakeRoot.C
		a. Added new TH1 histogram " hRunNumber ", for data only (!IsMC)
		b. Fill each eXi candidate's run number after following cuts:
		   runnumber validity, trigger, multiplicity, INEL>0 (if invoked), Xi mass tolerance,
		   pair mass low limit, and dummy tree rejection
		   -> L874
		c. WARNING: this fill point is BEFORE the quality cut (VLoose, Loose, ..., VTight)

	- macro_kimc/Xi0cAnaFunction.C
		a. FUnction " GetNormFac " : returns ANC->GetNeventsForNorm() for specific setup
		   * CAVEAT: require the train output (which has AliNormalizationCounter object)
		b. Function " GetRunFrac " : returns certain fraction of runs from given hRunNumber
		c. Function " GetV0xSec " : contains hard-coded V0 x-section values
		* To check how to use above a. - c., refer macro_kimc/Xi0cAnaExecute.C, L39 - L54

//-----------------------------------------------

 Sep. 17 (kimc)

	- AliAnalysisTaskSEXic0Semileptonic.cxx
		a. Updated pileup rejection configurations
			* https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGtoolsPileup
		b. Implemented # of contributors cut:
			b-1. For pp, nSPDTracklet < 20 : 3 (conventional, for low multiplicity)
			b-2. For pp, 20 <= nSPD... < 50 : 4
			b-3. For pp, nSPD... > 50 : 5
			b-4. For pA (pPb) : 5 (* NOTE exists in above link: PLEASE CHECK)

//-----------------------------------------------

 July 7 (kimc)

	- AliAnalysisTask...
		a. 4 new AliNormalizationCounter objects are added to address INEL>0 separately
		b. Moved "INEL>0" judge point before applying cuts (pileup, vertex < 10)
		c. Not pushed to AliPhysics yet

	- macro_kimc/Xic0AnaMakeRoot.C
		a. Modified function arguments to use a script (Xi0cAnaRun.sh)
		b. INEL>0 condition added for later use: currently disabled by default (conventional)
		c. Modified pT weighting function for MC: "1" will be applied if given parameters are 1, otherwise "expo"
		c. Unified all output histograms' binning
		d. Removed SPD trigger related part (some may still linger)

	- macro_kimc/Xic0AnaFunction.C : collection of functions to be used in the XS analysis
		a. NOT finished yet
		b. Prepared: get norm factor from ANC, prefilter eff corr, unfolding, Xic0 eff corr, XS, and pT weighting
		c. No hard-corded parameters in each function: especially get pT binning from histogram

	- macro_kimc/Xic0AnaExecute.C : steering macro for analysis
	- macro_kimc/Xic0AnaRun.sh : script for overall analysis

//-----------------------------------------------

 June 23 (kimc)

	- Added 4 new AliNormalizationCounter objects to account INEL>0 :
		a. In code, fCntINEL0_MB_0to100, fCntINEL0_MB_0p1to30, fCntINEL0_MB_30to100, and fCntINEL0_HMV0_0to0p1
		b. In Tree, ANCINEL0_MB_0to100, ANCINEL0_MB_0p1to30, ANCINEL0_MB_30to100, ANCINEL0_HMV0_0to0p1

//-----------------------------------------------

 June 16 (kimc)

	- Recovering analysis chain (still ongoing)
	- Xi0cAnaMakeRoot.C
		a. Modified arguments interface for better script based handling
		b. Argument " INELLgt0 " is by default on (being applied) - watch out
		c. In " eXiPairTree() ", enforced all histograms' binning to be same (1st one, 0 < pT < 20)
		d. In " eXiPairTree() ", modified fit function definition again (expo only when parameters are not 1.0)
	- Newly added Xi0cAnaRun.sh

//-----------------------------------------------

 Mar. 24 (kimc)

	- AliNormalizationCounter update
		a. Moved "StoreEvent" invoking point:
		   after "pileup rejection and Vtx->GetNContributors() > 1" and before |zVtx| < 10
		b. To match the events being analyzed and normalization factor obatained via ANC

//-----------------------------------------------

 Mar. 8 (kimc)

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

 Feb. 10 (kimc)

	- Added analysis macros for final cross section estimation
		a. Xi0cAnaMakeRoot.C (updated Xi topology cut after WDK)
		b. Xi0cAnaFunction.C (collection of xSec analysis subroutines)
		c. Xi0cAnaExecute.C (steering macro)

//-----------------------------------------------

 Jan. 25, 2021 (kimc)

	- Added 2nd production macro (Xic0AnaMakeRoot.C) for analysis
		a. A few modifications compared to original, especially adding > 0 on several denominators
		b. Tested on Lego train output available at Jan. 25

//-----------------------------------------------

 Dec. 18 (kimc)

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

 Oct. 8 (kimc)

	1. AliAnalysisTaskSEXi0c... is updated by Dr. Bok for compatability with pPb dataset
		a. pp: Get multiplicity from V0M
		b. PPb: Get multiplicity from V0A
	2. Macro for LEGO train running: AddTaskXic0Semileptonic.C
		a. Added High multiplicity trigger option (bool UseTrifHM)
			a-1. By default it is disabled: only kINY7 will be used
			a-2. Add kHMV0 and kHMSPD triggers f enabled
		b. Tested at Grid: no problem

//-----------------------------------------------

 Sep. 9 (kimc)

	1. Enforced EventTree sync to eXiTree
		a. EventTree is filled only if eXiTree is filled
		b. Added variable " fNeXiPair " under EventTree to enable event by event distinction:
		   if ( fNeXiPair != 0 ) it means those entries are duplication of same event
		c. Added " fVtxZ " under EventTree
			c-1. for data, fVtxZ = AliVVertex->GetZ()
			c-2. for MC, fVtxZ = AliAODMCHeader->GetVtxZ()

	2. Updated fWeightFit to get the parameters via functions

//-----------------------------------------------

 Aug. 28 (kimc):

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

