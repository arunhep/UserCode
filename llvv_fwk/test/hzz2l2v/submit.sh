#!/usr/bin/env bash

#--------------------------------------------------
# Global Code 
#--------------------------------------------------

if [[ $# -eq 0 ]]; then 
    printf "NAME\n\tsubmit.sh - Main driver to submit jobs\n"
    printf "\nSYNOPSIS\n"
    printf "\n\t%-5s\n" "./submit.sh [OPTION]" 
    printf "\nOPTIONS\n" 
    printf "\n\t%-5s  %-40s\n"  "0"  "completely clean up the directory" 
    printf "\n\t%-5s  %-40s\n"  "1"  "run 'runHZZ2l2vAnalysis' on samples.json" 
    printf "\n\t%-5s  %-40s\n"  "1.1"  "run 'runHZZ2l2vAnalysis' on photon_samples.json" 
    printf "\n\t%-5s  %-40s\n"  "1.2"  "run 'runHZZ2l2vAnalysis' on photon_samples.json with photon re-weighting"
    printf "\n\t%-5s  %-40s\n"  "2"  "compute integrated luminosity from processed samples" 
    printf "\n\t%-5s  %-40s\n"  "2.1"  "compute integrated luminosity from processed photon samples" 
    printf "\n\t%-5s  %-40s\n"  "3"  "make plots and combine root files" 
    printf "\n\t%-5s  %-40s\n"  "3.1"  "make plots for photon_samples" 
    printf "\n\t%-5s  %-40s\n"  "3.2"  "make plots for photon_samples with photon re-weighting" 
fi

step=$1   #variable that store the analysis step to run

#Additional arguments to take into account
arguments=''; for var in "${@:2}"; do arguments=$arguments" "$var; done
if [[ $# -ge 4 ]]; then echo "Additional arguments will be considered: "$arguments ;fi 

#--------------------------------------------------
# Global Variables
#--------------------------------------------------

SUFFIX=_forE

#SUFFIX=_2016_04_19
#SUFFIX=_debug
#SUFFIX=$(date +"_%Y_%m_%d") 
MAINDIR=$CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v
JSON=$MAINDIR/samples_dataE.json
#JSON=$MAINDIR/samples_dataEFG.json
RESULTSDIR=$MAINDIR/results$SUFFIX
PLOTSDIR=$MAINDIR/plots${SUFFIX}
PLOTTER=$MAINDIR/plotter${SUFFIX}

echo "Input: " $JSON
echo "Output: " $RESULTSDIR

queue='8nh'
#IF CRAB3 is provided in argument, use crab submissiong instead of condor/lsf
if [[ $arguments == *"crab3"* ]]; then queue='crab3' ;fi 


################################################# STEPS between 0 and 1
if [[ $step == 0 ]]; then   
        #analysis cleanup
	echo "ALL DATA WILL BE LOST! [N/y]?"
	read answer
	if [[ $answer == "y" ]];
	then
	    echo "CLEANING UP..."
	    rm -rdf $RESULTSDIR $PLOTSDIR LSFJOB_* core.* *.sh.e* *.sh.o*
	fi
fi #end of step0

###  ############################################## STEPS between 1 and 2
   if [[ $step == 1  ]]; then        #submit jobs for 2l2v analysis
	echo "JOB SUBMISSION"
	runAnalysisOverSamples.py -e runEfficencyMiniAODExample_HWW -j $JSON -o $RESULTSDIR  -c $MAINDIR/../runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s $queue --report True 
   fi
