void ReadRunInfo(const char* listPath, const char* taskOpt, string& datPath, string& datPatt, vector<int>& vecRuns)
{
	datPath = "";
	datPatt = "";
	vecRuns.clear();

	int counter = 0;
	int tempRun = 0;

	ifstream in;
	in.open(Form("%s/%s.txt", listPath, taskOpt));
	if (!in.is_open()) { cout <<"Cannot open input runs file!\n"; return; }
	while (in.is_open())
	{
		if      (counter == 0) in >> datPath;
		else if (counter == 1) in >> datPatt;
		else
		{
			in >> tempRun;
			if (!in.good())	{ break; in.close(); }
			vecRuns.push_back(tempRun);
		}
		counter++;
	}

	return;
}//ReadRunInfo

//==============================================================================================
void Xi0c_run(const int mode = 0, const char* listPath = ".", const char* taskOpt = "LHC18pAOD")
{
	enum {local, test, full, terminate, collect}; //Modes (argument)

	const char* taskName = "Xi0cSemiL";
	const char* taskAdd  = "Xi0c_add.C";
	const char* taskSrc  = "AliAnalysisTaskSEXic0Semileptonic.cxx";
	const char* taskLib  = "AliAnalysisTaskSEXic0Semileptonic.cxx AliAnalysisTaskSEXic0Semileptonic.h";

	//-------------------------------------------

	TString taskOptStr = taskOpt;
	bool isMC = (taskOptStr.Contains("MC"))?true:false;

    //Analysis manager, Input handler
    AliAnalysisManager* mgr = new AliAnalysisManager(Form("%s_%s", taskName, taskOpt));
    AliInputEventHandler* hdr;
	if (taskOptStr.Contains("AOD") || taskOptStr.Contains("MC")) hdr = new AliAODInputHandler();
	else if (taskOptStr.Contains("ESD")) hdr = new AliESDInputHandler();
	else { cout <<"Input type?" <<endl; return; }
	hdr->SetNeedField(1); //?
    mgr->SetInputEventHandler(hdr);

    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
    gROOT->LoadMacro(Form("%s++g", taskSrc));

	//-------------------------------------------

	//Physics selection
	if (taskOptStr.Contains("AOD"))
	{
		gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
		AliPhysicsSelectionTask* taskPhys = AddTaskPhysicsSelection(isMC, 1); //0 for PbPb, 1 for pp/pPb
		if (!taskPhys) { cout <<"Cannot find AliPhysicsSelectionTask!" <<endl; return; }
	}

	//Multiplicity
	gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
	AliMultSelectionTask* taskMult = AddTaskMultSelection(true);
	if (!taskMult) { cout <<"Cannot find AliMultSelectionTask!" <<endl; return; }

	//pID
	gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTaskPIDResponse* taskPid = AddTaskPIDResponse(isMC);
	if (!taskPid) { cout <<"Cannot find AliAnalysisTaskPIDResponse!" <<endl; return; }

	//Main task
	gROOT->LoadMacro(taskAdd);
	TString taskNameStr = taskName;
	AliAnalysisTaskSEXic0Semileptonic* task = Xi0c_add(taskNameStr, taskOptStr);
	if (!task) { cout <<"Cannot find AliAnalysisTaskSEXic0Semileptonic!" <<endl; return; }

	//-------------------------------------------

    if (!mgr->InitAnalysis()) return;
    //mgr->SetDebugLevel(2);
    //mgr->PrintStatus();

	if (mode == local)
	{
        TChain* chain = new TChain("aodTree");
        chain->Add(Form("AliAOD_%s.root", isMC?"mc":"data"));
		mgr->SetUseProgressBar(1, 500);
        mgr->StartAnalysis("local", chain);
	}
	else //Grid
	{
		//Setup
        AliAnalysisAlien *alien = new AliAnalysisAlien();
		alien->AddIncludePath("-I$ALICE_ROOT/include \
				               -I$ALICE_ROOT/lib \
							   -I$ALICE_PHYSICS/include \
							   -I$ALICE_PHYSICS/lib \
							   -I$ALICE_PHYSICS/OADB/macros");
		alien->SetAPIVersion("V1.1x");
        alien->SetAdditionalLibs(taskLib);
        alien->SetAnalysisSource(taskSrc);
        alien->SetAliPhysicsVersion("vAN-20200322-1");

		//Running
        alien->SetExecutable(Form("%s_%s.sh", taskName, taskOpt));
        alien->SetFileForTestMode("TestMode.txt"); //?
        alien->SetInputFormat("xml-single"); //?
        alien->SetJDLName(Form("%s_%s.jdl", taskName, taskOpt));
        alien->SetMasterResubmitThreshold(90);
		alien->SetNrunsPerMaster(6);
        alien->SetPrice(1); //?
        alien->SetSplitMaxInputFileNumber(60);
        alien->SetSplitMode("se"); //?
        alien->SetTTL(3600 * 12);
		alien->SetUseSubmitPolicy(); //?

		//Output
		alien->SetDefaultOutputs(false); //Enables all outputs of the tasks connected to the analysis manager
		alien->SetDropToShell(0); //?
		alien->SetGridWorkingDir(Form("%s_%s", taskName, taskOpt));
		alien->SetGridOutputDir("out");
		alien->SetOutputFiles(AliAnalysisManager::GetCommonFileName()); //AnalysisResults.root
		alien->SetOutputToRunNo();
		alien->SetKeepLogs(true);

		//Merge
		alien->SetMergeViaJDL((mode==collect)?false:true);
		alien->SetMaxMergeStages(1);
		//alien->SetMaxMergeFiles(500); //?

		//***************************************

		if (isMC == false) alien->SetRunPrefix("000"); //Data

		string datPath;
		string datPatt;
		vector<int> vecRuns;
		ReadRunInfo(listPath, taskOpt, datPath, datPatt, vecRuns);
		if (vecRuns.size() == 0) { cout <<"No run info found! Stop...\n"; return; }
		else
		{
			alien->SetGridDataDir(datPath.c_str());
			alien->SetDataPattern(datPatt.c_str());
			for (unsigned int i=0; i<vecRuns.size(); i++) alien->AddRunNumber(vecRuns[i]);
		}

		//***************************************

		if (mode == test)
		{
			alien->SetRunMode("test");
			alien->SetNtestFiles(1);
		}
		else if (mode == full) 
		{
			alien->SetRunMode("full");
		}
		else if ((mode==terminate) || (mode==collect))
		{
			alien->SetRunMode("terminate");
		}

		cout <<Form("Starting analysis %s_%s...\n", taskName, taskOpt) <<endl;
		mgr->SetGridHandler(alien);
		mgr->StartAnalysis("grid");

		delete alien; //Close connection to Grid
	}//Grid

	return;
}//Main
