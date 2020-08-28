#!/bin/bash

#Script for bash submission to Grid
#kimc#cern.ch

MODE=$1 #2 (full), 3 (terminate), and 4 (collect merged output)
LIST=$2
echo "Submitting Grid jobs in mode" $1 "w/ list" $2...

#Cheek required files exist
if [ ! -f "./AliAnalysisTaskSEXic0Semileptonic.cxx" ] ||
   [ ! -f "./AliAnalysisTaskSEXic0Semileptonic.h" ] ||
   [ ! -f "./Xi0c_add.C" ] ||
   [ ! -f "./Xi0c_run.C" ] ; then
	echo "Cannot find required file(s): stop."
	set -e
fi

#Submit jobs
for ITEM in `cat $LIST`; do

	LSPATH="$(dirname $ITEM)"
	PREFIX="$(basename $ITEM | awk -F '.' '{print $1}')"
	echo "Submitting..." $LSPATH"/"$PREFIX".txt"

	root -l -b -q Xi0c_run.C\($MODE,\"$LSPATH\",\"$PREFIX\"\)

	if [ $MODE != 0 ] ; then
		echo "Cleaning up Grid files..."
		rm *.d *.so *.pcm *.xml std* myAnalysis.C
	fi

	#Conenction to Alien kept open even after the job submitted: it causes crach for next submit...
done
