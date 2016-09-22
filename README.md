This framework contains codes for Trigger efficiency calculations.

Instructions to get the code and run it.

cmsrel CMSSW_8_0_18

cd CMSSW_8_0_18/src

cmsenv

scramv1 b

cd UserCode/llvv_fwk/test/hzz2l2v/

sh submit 1

Important Points:

1. To change the datasets and JSON file etc : modify the samples_dataE.json file.

2. submit.sh is the script where you can change the name of executables, output directories etc.

3. Main code is here : UserCode/llvv_fwk/bin , latest one is runEfficencyMiniAODExample_HWW.cc which saves a tree having branches of pt and eta of passing probes as well as total probes for each HLT leg defined inside it.

4. ReferenceMethod codes are being used to compute the efficiency using reference method where you can compute the efficiencies for the soup of triggers.
