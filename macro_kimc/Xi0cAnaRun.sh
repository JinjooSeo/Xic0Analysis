#/bin/bash -f

export rlbq="root -l -b -q"

# Arguments: file, trigger, multiplicity min, multiplicity max, apply weight on MC, and INEL>0

$rlbq "Xi0cAnaMakeRoot.C+(\"AnalysisResults_data.root\", \"MB\", 0, 100)" #, false, true)"
$rlbq "Xi0cAnaMakeRoot.C+(\"AnalysisResults_data.root\", \"MB\", 0.1, 30)" #, false, true)"
$rlbq "Xi0cAnaMakeRoot.C+(\"AnalysisResults_data.root\", \"MB\", 30, 100)" #, false, true)"
$rlbq "Xi0cAnaMakeRoot.C+(\"AnalysisResults_data.root\", \"HMV0\", 0, 0.1)" #, false, true)"

$rlbq "Xi0cAnaMakeRoot.C+(\"AnalysisResults_MC.root\", \"MB\", 0, 100, true)" #, true)"
$rlbq "Xi0cAnaMakeRoot.C+(\"AnalysisResults_MC.root\", \"MB\", 0.1, 30, true)" #, true)"
$rlbq "Xi0cAnaMakeRoot.C+(\"AnalysisResults_MC.root\", \"MB\", 30, 100, true)" #, true)"
$rlbq "Xi0cAnaMakeRoot.C+(\"AnalysisResults_MC.root\", \"HMV0\", 0, 0.1, true)" #, true)"
# $rlbq "Xi0cAnaMakeRoot.C+(\"AnalysisResults_MC.root\", \"MB\", 0, 100, false)" #, true)" # w/o weight on MC

#root -l "Xi0cAnaExecute.C+(false)"
