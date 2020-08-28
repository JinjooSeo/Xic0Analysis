AliAnalysisTaskSEXic0Semileptonic* Xi0c_add(
		TString taskName,
		TString taskOpt
		)
{
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr || !mgr->GetInputEventHandler()) return 0x0;

	AliAnalysisTaskSEXic0Semileptonic *task = new AliAnalysisTaskSEXic0Semileptonic(taskName, taskOpt);
	if (!task) return 0x0;
	mgr->AddTask(task);

	//-------------------------------------------

	//Data or MC
	task->SetMC((taskOpt.Contains("MC"))?true:false);
	cout <<Form("\nSet input type: %s\n", (taskOpt.Contains("MC"))?"MC":"DATA") <<endl;

	//Add triggers
	task->UseTrig_kINT7();
	task->UseTrig_kHMV0();
	task->UseTrig_kHMSPD();

	//Containers
	TString taskOut = AliAnalysisManager::GetCommonFileName();
	AliAnalysisDataContainer *taskI0 = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *taskO1, *taskO2, *taskO3, *taskO4, *taskO5, *taskO6, *taskO7;
	taskO1 = mgr->CreateContainer("histogram",    TDirectory::Class(), AliAnalysisManager::kOutputContainer, taskOut);
	taskO2 = mgr->CreateContainer("cut",          TList::Class(), AliAnalysisManager::kOutputContainer, taskOut);
	taskO3 = mgr->CreateContainer("MCCutTree",    TTree::Class(), AliAnalysisManager::kOutputContainer, taskOut);
	taskO4 = mgr->CreateContainer("PairTree",     TTree::Class(), AliAnalysisManager::kOutputContainer, taskOut);
	taskO5 = mgr->CreateContainer("MCPairTree",   TTree::Class(), AliAnalysisManager::kOutputContainer, taskOut);
	taskO6 = mgr->CreateContainer("EventTree",    TTree::Class(), AliAnalysisManager::kOutputContainer, taskOut);
	taskO7 = mgr->CreateContainer("eleXiCounter", AliNormalizationCounter::Class(),
			AliAnalysisManager::kOutputContainer, taskOut);

	int iCon = 0;
	mgr->ConnectInput (task, iCon, taskI0); iCon++;
	mgr->ConnectOutput(task, iCon, taskO1); iCon++;
	mgr->ConnectOutput(task, iCon, taskO2); iCon++;
	mgr->ConnectOutput(task, iCon, taskO3); iCon++;
	mgr->ConnectOutput(task, iCon, taskO4); iCon++;
	mgr->ConnectOutput(task, iCon, taskO5); iCon++;
	mgr->ConnectOutput(task, iCon, taskO6); iCon++;
	mgr->ConnectOutput(task, iCon, taskO7);

	//-------------------------------------------

	//Trigger (official, exc)
	//task->SelectCollisionCandidates(AliVEvent::kAnyINT);     //kMB | kINT7 | kINT5 | kINT8 | kSPI7
	//task->SelectCollisionCandidates(AliVEvent::kMB);         //Minimum bias in PbPb 2010-11
	//task->SelectCollisionCandidates(AliVEvent::kINT1);       //V0A | V0C | SPD minimum bias trigger
	//task->SelectCollisionCandidates(AliVEvent::kINT5);       //V0OR minimum bias trigger
	//task->SelectCollisionCandidates(AliVEvent::kINT7);       //V0AND minimum bias trigger
	//task->SelectCollisionCandidates(AliVEvent::kHighMultV0); //High-multiplicity V0 trigger
	//task->SetFilterBit(128);

	return task;
}//Xi0c_add
