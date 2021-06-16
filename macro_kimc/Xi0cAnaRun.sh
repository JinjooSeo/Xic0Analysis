#/bin/bash -f

export rlbq="root -l -b -q"

$rlbq "Xi0cAnaMakeRoot.C+( \"AnalysisResults_data.root\", \"MB\", 0, 100)"
$rlbq "Xi0cAnaMakeRoot.C+( \"AnalysisResults_data.root\", \"MB\", 0.1, 30)"
$rlbq "Xi0cAnaMakeRoot.C+( \"AnalysisResults_data.root\", \"MB\", 30, 100)"
$rlbq "Xi0cAnaMakeRoot.C+( \"AnalysisResults_data.root\", \"HMV0\", 0, 0.1)"
$rlbq "Xi0cAnaMakeRoot.C+( \"AnalysisResults_MC.root\")"
$rlbq "Xi0cAnaExecute.C+"
