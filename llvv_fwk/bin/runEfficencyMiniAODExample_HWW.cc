#include <iostream>
#include <math.h>
#include <boost/shared_ptr.hpp>

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

//Load here all the dataformat that we will need
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
//#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h" //for svfit

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
//#include "UserCode/llvv_fwk/interface/HiggsUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/TMVAUtils.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"
#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"
#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"

#include "UserCode/llvv_fwk/interface/PatUtils.h"
#include "UserCode/llvv_fwk/interface/LumiUtils.h"

//L1 EM particles
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"


#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TNtuple.h"
#include <Math/VectorUtil.h>
#include <TMath.h>

using namespace std;


bool passFilter(pat::TriggerObjectStandAlone& Triggerobj, std::string& FilterName)
   {   
     for (unsigned h = 0; h < Triggerobj.filterLabels().size(); ++h){
//cout << Triggerobj.filterLabels()[h] << endl;
       if(FilterName.compare(Triggerobj.filterLabels()[h])==0) return true;
   }   
       return false;

   }   

double DeltaRtrig(const pat::Electron& el, std::vector<pat::TriggerObjectStandAlone> object)
   {   

    std::vector<double> coll;
    for (unsigned int i = 0 ; i < object.size() ; i++){
        float dR = deltaR(el.superCluster()->eta(),el.phi(),object[i].eta(),object[i].phi()) ;

        coll.push_back(dR);
   }   
        std::sort(coll.begin(),coll.end());
        if(coll.size() <=  0 ) return 999.0 ;
        else return coll.at(0);
   }   




/*double DeltaR(const pat::Electron& el, std::vector<reco::GenParticle> object)
   {   
    std::vector<double> coll;
    for (unsigned int i = 0 ; i < object.size() ; i++){
        float dR = deltaR(el.eta(),el.phi(),object[i].eta(),object[i].phi()) ;

        coll.push_back(dR);
   }
        std::sort(coll.begin(),coll.end());
        if(coll.size() <=  0 ) return 999.0 ;
        else return coll.at(0);
   }


*/

bool hasZMother(const reco::GenParticle  p)
{
    bool foundZ(false);
    const reco::Candidate  *part = (p.mother());
    // loop on the mother particles to check if is has a W has mother
         while ((part->numberOfMothers()>0)) {
         const reco::Candidate  *MomPart =part->mother();
         if ((fabs(MomPart->pdgId())==23)){
             foundZ = true;
             break;
             }
            part = MomPart;
            }
            return foundZ;
          }


int main(int argc, char* argv[])
{
  //##############################################
  //########    GLOBAL INITIALIZATION     ########
  //##############################################
 
  // check arguments
  if(argc<2){ std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl; exit(0); }
  
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  
  // configure the process
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");

  bool isMC = runProcess.getParameter<bool>("isMC");
  double xsec = runProcess.getParameter<double>("xsec");
  int mctruthmode=runProcess.getParameter<int>("mctruthmode");
  TString dtag=runProcess.getParameter<std::string>("dtag");

  TString suffix=runProcess.getParameter<std::string>("suffix");
  std::vector<std::string> urls=runProcess.getUntrackedParameter<std::vector<std::string> >("input");

TString outUrl = runProcess.getParameter<std::string>("outfile");

TFile *ofile=TFile::Open(outUrl, "recreate");
TTree *ntuple     = new TTree("ntuple","Efficiency Tree");

double Nvtx;
double ZMass;
double tagPT;
double totalProbePT;
double passProbePTLeg1;
double passProbePTLeg2;
double passProbePTLeg3;
double passProbePTLeg4;
double passProbePTLeg5;
double passProbePTLeg6;
double passProbePTLeg7;
double passProbePTLeg8;
double passProbePTLeg9;
double passProbePTLeg10;
double passProbePTLeg11;
double passProbePTLeg12;
double passProbePTLeg13;
double passProbePTLeg14;
double passProbePTLeg15;
double passProbePTLeg16;
double tagEta;
double totalProbeEta;
double passProbeEtaLeg1;
double passProbeEtaLeg2;
double passProbeEtaLeg3;
double passProbeEtaLeg4;
double passProbeEtaLeg5;
double passProbeEtaLeg6;
double passProbeEtaLeg7;
double passProbeEtaLeg8;
double passProbeEtaLeg9;
double passProbeEtaLeg10;
double passProbeEtaLeg11;
double passProbeEtaLeg12;
double passProbeEtaLeg13;
double passProbeEtaLeg14;
double passProbeEtaLeg15;
double passProbeEtaLeg16;
double totalProbePhi;
double passProbePhiLeg1;
double passProbePhiLeg2;
double passProbePhiLeg3;
double passProbePhiLeg4;
double passProbePhiLeg5;
double passProbePhiLeg6;
double passProbePhiLeg7;
double passProbePhiLeg8;
double passProbePhiLeg9;
double passProbePhiLeg10;
double passProbePhiLeg11;
double passProbePhiLeg12;
double passProbePhiLeg13;
double passProbePhiLeg14;
double passProbePhiLeg15;
double passProbePhiLeg16;
double PUWeight;
double lumi_weight;
double counting;
double count_tag;
double count_probe;
double count_passprobeleg1;
double count_passprobeleg2;
double count_passprobeleg3;
double passL1PTNorm;
double passL1PTIso;
double passL1PTER;
double passL1PT;
double passL1EtaNorm;
double passL1EtaIso;
double passL1EtaER;
double passL1Eta;

TBranch *b_Nvtx = ntuple->Branch("Nvtx",&Nvtx,"Nvtx/D");
TBranch *b_ZMass = ntuple->Branch("ZMass",&ZMass,"ZMass/D");
TBranch *b_tagPT = ntuple->Branch("tagPT",&tagPT,"tagPT/D");
TBranch *b_tagEta = ntuple->Branch("tagEta",&tagEta,"tagEta/D");
TBranch *b_totalProbePT = ntuple->Branch("totalProbePT",&totalProbePT,"totalProbePT/D");
TBranch *b_totalProbeEta = ntuple->Branch("totalProbeEta",&totalProbeEta,"totalProbeEta/D");
TBranch *b_passProbePTLeg1 = ntuple->Branch("passProbePTLeg1",&passProbePTLeg1,"passProbePTLeg1/D");
TBranch *b_passProbePTLeg2 = ntuple->Branch("passProbePTLeg2",&passProbePTLeg2,"passProbePTLeg2/D");
TBranch *b_passProbePTLeg3 = ntuple->Branch("passProbePTLeg3",&passProbePTLeg3,"passProbePTLeg3/D");
TBranch *b_passProbePTLeg4 = ntuple->Branch("passProbePTLeg4",&passProbePTLeg4,"passProbePTLeg4/D");
TBranch *b_passProbePTLeg5 = ntuple->Branch("passProbePTLeg5",&passProbePTLeg5,"passProbePTLeg5/D");
TBranch *b_passProbePTLeg6 = ntuple->Branch("passProbePTLeg6",&passProbePTLeg6,"passProbePTLeg6/D");
TBranch *b_passProbePTLeg7 = ntuple->Branch("passProbePTLeg7",&passProbePTLeg7,"passProbePTLeg7/D");
TBranch *b_passProbePTLeg8 = ntuple->Branch("passProbePTLeg8",&passProbePTLeg8,"passProbePTLeg8/D");
TBranch *b_passProbePTLeg9 = ntuple->Branch("passProbePTLeg9",&passProbePTLeg9,"passProbePTLeg9/D");
TBranch *b_passProbePTLeg10 = ntuple->Branch("passProbePTLeg10",&passProbePTLeg10,"passProbePTLeg10/D");
TBranch *b_passProbePTLeg11 = ntuple->Branch("passProbePTLeg11",&passProbePTLeg11,"passProbePTLeg11/D");
TBranch *b_passProbePTLeg12 = ntuple->Branch("passProbePTLeg12",&passProbePTLeg12,"passProbePTLeg12/D");
TBranch *b_passProbePTLeg13 = ntuple->Branch("passProbePTLeg13",&passProbePTLeg13,"passProbePTLeg13/D");
TBranch *b_passProbePTLeg14 = ntuple->Branch("passProbePTLeg14",&passProbePTLeg14,"passProbePTLeg14/D");
TBranch *b_passProbePTLeg15 = ntuple->Branch("passProbePTLeg15",&passProbePTLeg15,"passProbePTLeg15/D");
TBranch *b_passProbePTLeg16 = ntuple->Branch("passProbePTLeg16",&passProbePTLeg16,"passProbePTLeg16/D");
TBranch *b_passProbeEtaLeg1 = ntuple->Branch("passProbeEtaLeg1",&passProbeEtaLeg1,"passProbeEtaLeg1/D");
TBranch *b_passProbeEtaLeg2 = ntuple->Branch("passProbeEtaLeg2",&passProbeEtaLeg2,"passProbeEtaLeg2/D");
TBranch *b_passProbeEtaLeg3 = ntuple->Branch("passProbeEtaLeg3",&passProbeEtaLeg3,"passProbeEtaLeg3/D");
TBranch *b_passProbeEtaLeg4 = ntuple->Branch("passProbeEtaLeg4",&passProbeEtaLeg4,"passProbeEtaLeg4/D");
TBranch *b_passProbeEtaLeg5 = ntuple->Branch("passProbeEtaLeg5",&passProbeEtaLeg5,"passProbeEtaLeg5/D");
TBranch *b_passProbeEtaLeg6 = ntuple->Branch("passProbeEtaLeg6",&passProbeEtaLeg6,"passProbeEtaLeg6/D");
TBranch *b_passProbeEtaLeg7 = ntuple->Branch("passProbeEtaLeg7",&passProbeEtaLeg7,"passProbeEtaLeg7/D");
TBranch *b_passProbeEtaLeg8 = ntuple->Branch("passProbeEtaLeg8",&passProbeEtaLeg8,"passProbeEtaLeg8/D");
TBranch *b_passProbeEtaLeg9 = ntuple->Branch("passProbeEtaLeg9",&passProbeEtaLeg9,"passProbeEtaLeg9/D");
TBranch *b_passProbeEtaLeg10 = ntuple->Branch("passProbeEtaLeg10",&passProbeEtaLeg10,"passProbeEtaLeg10/D");
TBranch *b_passProbeEtaLeg11 = ntuple->Branch("passProbeEtaLeg11",&passProbeEtaLeg11,"passProbeEtaLeg11/D");
TBranch *b_passProbeEtaLeg12 = ntuple->Branch("passProbeEtaLeg12",&passProbeEtaLeg12,"passProbeEtaLeg12/D");
TBranch *b_passProbeEtaLeg13 = ntuple->Branch("passProbeEtaLeg13",&passProbeEtaLeg13,"passProbeEtaLeg13/D");
TBranch *b_passProbeEtaLeg14 = ntuple->Branch("passProbeEtaLeg14",&passProbeEtaLeg14,"passProbeEtaLeg14/D");
TBranch *b_passProbeEtaLeg15 = ntuple->Branch("passProbeEtaLeg15",&passProbeEtaLeg15,"passProbeEtaLeg15/D");
TBranch *b_passProbeEtaLeg16 = ntuple->Branch("passProbeEtaLeg16",&passProbeEtaLeg16,"passProbeEtaLeg16/D");
TBranch *b_passProbePhiLeg1 = ntuple->Branch("passProbePhiLeg1",&passProbePhiLeg1,"passProbePhiLeg1/D");
TBranch *b_passProbePhiLeg2 = ntuple->Branch("passProbePhiLeg2",&passProbePhiLeg2,"passProbePhiLeg2/D");
TBranch *b_passProbePhiLeg3 = ntuple->Branch("passProbePhiLeg3",&passProbePhiLeg3,"passProbePhiLeg3/D");
TBranch *b_passProbePhiLeg4 = ntuple->Branch("passProbePhiLeg4",&passProbePhiLeg4,"passProbePhiLeg4/D");
TBranch *b_passProbePhiLeg5 = ntuple->Branch("passProbePhiLeg5",&passProbePhiLeg5,"passProbePhiLeg5/D");
TBranch *b_passProbePhiLeg6 = ntuple->Branch("passProbePhiLeg6",&passProbePhiLeg6,"passProbePhiLeg6/D");
TBranch *b_passProbePhiLeg7 = ntuple->Branch("passProbePhiLeg7",&passProbePhiLeg7,"passProbePhiLeg7/D");
TBranch *b_passProbePhiLeg8 = ntuple->Branch("passProbePhiLeg8",&passProbePhiLeg8,"passProbePhiLeg8/D");
TBranch *b_passProbePhiLeg9 = ntuple->Branch("passProbePhiLeg9",&passProbePhiLeg9,"passProbePhiLeg9/D");
TBranch *b_passProbePhiLeg10 = ntuple->Branch("passProbePhiLeg10",&passProbePhiLeg10,"passProbePhiLeg10/D");
TBranch *b_passProbePhiLeg11 = ntuple->Branch("passProbePhiLeg11",&passProbePhiLeg11,"passProbePhiLeg11/D");
TBranch *b_passProbePhiLeg12 = ntuple->Branch("passProbePhiLeg12",&passProbePhiLeg12,"passProbePhiLeg12/D");
TBranch *b_passProbePhiLeg13 = ntuple->Branch("passProbePhiLeg13",&passProbePhiLeg13,"passProbePhiLeg13/D");
TBranch *b_passProbePhiLeg14 = ntuple->Branch("passProbePhiLeg14",&passProbePhiLeg14,"passProbePhiLeg14/D");
TBranch *b_passProbePhiLeg15 = ntuple->Branch("passProbePhiLeg15",&passProbePhiLeg15,"passProbePhiLeg15/D");
TBranch *b_passProbePhiLeg16 = ntuple->Branch("passProbePhiLeg16",&passProbePhiLeg16,"passProbePhiLeg16/D");
TBranch *b_totalProbePhi = ntuple->Branch("totalProbePhi",&totalProbePhi,"totalProbePhi/D");


TBranch *b_PUWeight = ntuple->Branch("PUWeight",&PUWeight,"PUWeight/D");
TBranch *b_lumi_weight = ntuple->Branch("lumi_weight",&lumi_weight,"lumi_weight/D");
TBranch *b_counting = ntuple->Branch("counting",&counting);
TBranch *b_count_tag = ntuple->Branch("count_tag",&count_tag);
TBranch *b_count_probe = ntuple->Branch("count_probe",&count_probe);
TBranch *b_count_paasprobeleg1 = ntuple->Branch("count_passprobeleg1",&count_passprobeleg1,"count_passprobeleg1/D");
TBranch *b_count_paasprobeleg2 = ntuple->Branch("count_passprobeleg2",&count_passprobeleg2,"count_passprobeleg2/D");
TBranch *b_passL1PTNorm = ntuple->Branch("passL1PTNorm",&passL1PTNorm,"passL1PTNorm/D");
TBranch *b_passL1PTIso = ntuple->Branch("passL1PTIso",&passL1PTIso,"passL1PTIso/D");
TBranch *b_passL1PTER = ntuple->Branch("passL1PTER",&passL1PTER,"passL1PTER/D");
TBranch *b_passL1PT = ntuple->Branch("passL1PT",&passL1PT,"passL1PT/D");
TBranch *b_passL1EtaNorm = ntuple->Branch("passL1EtaNorm",&passL1EtaNorm,"passL1EtaNorm/D");
TBranch *b_passL1EtaIso = ntuple->Branch("passL1EtaIso",&passL1EtaIso,"passL1EtaIso/D");
TBranch *b_passL1EtaER = ntuple->Branch("passL1EtaER",&passL1EtaER,"passL1EtaER/D");
TBranch *b_passL1Eta = ntuple->Branch("passL1Eta",&passL1Eta,"passL1Eta/D");

// good lumi mask
  lumiUtils::GoodLumiFilter goodLumiFilter(runProcess.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> >("lumisToProcess", std::vector<edm::LuminosityBlockRange>()));

  bool filterOnlyEE(false), filterOnlyMUMU(false), filterOnlyEMU(false);
  bool isSingleMuPD(!isMC && dtag.Contains("SingleMu"));  
  bool isV0JetsMC(isMC && (dtag.Contains("DYJetsToLL_50toInf") || dtag.Contains("_WJets"))); // #FIXME should be reactivated as soon as we have exclusive jet samples
  bool isWGmc(isMC && dtag.Contains("WG"));
  bool isZGmc(isMC && dtag.Contains("ZG"));
  bool isMC_GG  = isMC && ( string(dtag.Data()).find("GG" )  != string::npos);
  bool isMC_VBF = isMC && ( string(dtag.Data()).find("VBF")  != string::npos);
  bool isMC_125OnShell = isMC && (mctruthmode==521);
  if(isMC_125OnShell) mctruthmode=125;
  bool isMC_ZZ  = isMC && ( string(dtag.Data()).find("TeV_ZZ")  != string::npos);
  bool isMC_WZ  = isMC && ( string(dtag.Data()).find("TeV_WZ")  != string::npos);

  bool isData_DoubleEle = string(dtag.Data()).find("Data13TeV_DoubleElectron");


  const edm::ParameterSet& myVidElectronIdConf 		= runProcess.getParameterSet("electronidparas");
  const edm::ParameterSet& myVidElectronTightIdWPConf 	= myVidElectronIdConf.getParameterSet("tight");
  const edm::ParameterSet& myVidElectronMediumIdWPConf 	= myVidElectronIdConf.getParameterSet("medium");
  const edm::ParameterSet& myVidElectronLooseIdWPConf 	= myVidElectronIdConf.getParameterSet("loose");

  VersionedPatElectronSelector electronVidTightId(myVidElectronTightIdWPConf);
  VersionedPatElectronSelector electronVidMediumId(myVidElectronMediumIdWPConf);
  VersionedPatElectronSelector electronVidLooseId(myVidElectronLooseIdWPConf);

  //##############################################
  //########    INITIATING HISTOGRAMS     ########
  //##############################################
  SmartSelectionMonitor mon;
  
  int evcounter = 0;
  int counter1 = 0;
  int counter2 = 0; 
 
  edm::InputTag nullTag("None");
  edm::EDGetToken m_dmxEGToken;
  bool m_doDmxEGs;

//  edm::InputTag dmxEGTag  = runProcess.getParameter<edm::InputTag>("caloStage2Digis");
//   m_dmxEGToken          = consumes<l1t::EGamma>(dmxEGTag);
//   m_doDmxEGs            = !(dmxEGTag==nullTag); 
 

  //##############################################
  //######## GET READY FOR THE EVENT LOOP ########
  //##############################################

  fwlite::ChainEvent ev(urls);
  const size_t totalEntries= ev.size();
  
  //MC normalization (to 1/pb)
  double xsecWeight = xsec/totalEntries;
  if(!isMC) xsecWeight=1.0;
/*

  //pileup weighting
  edm::LumiReWeighting* LumiWeights = NULL;
  utils::cmssw::PuShifter_t PuShifters;
  double PUNorm[] = {1,1,1};
  if(isMC){
          std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
          std::vector<float> dataPileupDistribution; for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);}
          std::vector<float> mcPileupDistribution;
          utils::getMCPileupDistributionFromMiniAOD(ev,dataPileupDistribution.size(), mcPileupDistribution);
          while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
          while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);
          gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
          LumiWeights = new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution);
          PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution,0.05);
          utils::getPileupNormalization(mcPileupDistribution, PUNorm, LumiWeights, PuShifters);
  }
 
  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE/
*/
  //##############################################
  //########           EVENT LOOP         ########
  //##############################################
  //loop on all the events
  printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
  printf("Scanning the ntuple :");
  int treeStep(totalEntries/50);
  //DuplicatesChecker duplicatesChecker;
  //int nDuplicates(0);
  //totalEntries->10
  int bothlegs = 0;


 count_probe = 0;
 count_passprobeleg1 = 0;
 count_passprobeleg2 = 0;
 count_passprobeleg3 = 0;

 count_tag = 0;
// count_Z = 0;
 counting = 0;

//if(ev%10000==0) cout << ev << endl;

  ofstream outputfile;
  outputfile.open("bothLegPassed.txt");

cout << "totalEntries = " << totalEntries << endl; 
  for( size_t iev=0; iev< totalEntries; iev++){
//  for( size_t iev=0; iev< 2; iev++){  
//cout <<"****************************" << endl;

if(iev%10000==0) {cout <<"Number of events processed =   " << iev << endl;}
    

       //##############################################   EVENT LOOP STARTS   ##############################################
       ev.to(iev); //load the event content from the EDM file
       edm::EventBase const & myEvent = ev;
       //if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) { nDuplicates++; continue; }

//       bool passOnlyMuon(false);
//       bool passOnlyEle(false); 

Nvtx = -999;
ZMass = -999.;
tagPT = -999.;
totalProbePT = -999.;
passProbePTLeg1 = -999.;
passProbePTLeg2 = -999.;
passProbePTLeg3 = -999.;
passProbePTLeg4 = -999.;
passProbePTLeg5 = -999.;
passProbePTLeg6 = -999.;
passProbePTLeg7 = -999.;
passProbePTLeg8 = -999.;
passProbePTLeg9 = -999.;
passProbePTLeg10 = -999.;
passProbePTLeg11 = -999.;
passProbePTLeg12 = -999.;
passProbePTLeg13 = -999.;
passProbePTLeg14 = -999.;
passProbePTLeg15 = -999.;
passProbePTLeg16 = -999.;
passL1PTNorm = -999.;
passL1PTIso = -999.;
passL1PTER = -999.;
passL1PT = -999.;
tagEta = -999.;
totalProbeEta = -999.;
passProbeEtaLeg1 = -999.;
passProbeEtaLeg2 = -999.;
passProbeEtaLeg3 = -999.;
passProbeEtaLeg4 = -999.;
passProbeEtaLeg5 = -999.;
passProbeEtaLeg6 = -999.;
passProbeEtaLeg7 = -999.;
passProbeEtaLeg8 = -999.;
passProbeEtaLeg9 = -999.;
passProbeEtaLeg10 = -999.;
passProbeEtaLeg11 = -999.;
passProbeEtaLeg12 = -999.;
passProbeEtaLeg13 = -999.;
passProbeEtaLeg14 = -999.;
passProbeEtaLeg15 = -999.;
passProbeEtaLeg16 = -999.;
totalProbePhi = -999.;
passProbePhiLeg1 = -999.;
passProbePhiLeg2 = -999.;
passProbePhiLeg3 = -999.;
passProbePhiLeg4 = -999.;
passProbePhiLeg5 = -999.;
passProbePhiLeg6 = -999.;
passProbePhiLeg7 = -999.;
passProbePhiLeg8 = -999.;
passProbePhiLeg9 = -999.;
passProbePhiLeg10 = -999.;
passProbePhiLeg11 = -999.;
passProbePhiLeg12 = -999.;
passProbePhiLeg13 = -999.;
passProbePhiLeg14 = -999.;
passProbePhiLeg15 = -999.;
passProbePhiLeg16 = -999.;
passL1EtaNorm = -999.;
passL1EtaIso = -999.;
passL1EtaER = -999.;
passL1Eta = -999.;
PUWeight = 1.;
lumi_weight = 1.;

double rho = 0.;

	++evcounter;
       //apply trigger and require compatibilitiy of the event with the PD
       edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
       if(!tr.isValid())return false;
      
       //std::cout << "Definition of the Booleian variable" << std::endl;
 
       bool passOnlyMuon(false);
       bool passOnlyEle(false); 

              
       //load all the objects we will need to access
       reco::VertexCollection vtx;
       fwlite::Handle< reco::VertexCollection > vtxHandle; 
       vtxHandle.getByLabel(ev, "offlineSlimmedPrimaryVertices");
       if(vtxHandle.isValid()){ vtx = *vtxHandle;}
/*
       reco::GenParticleCollection gen;
       fwlite::Handle< reco::GenParticleCollection > genHandle;
       genHandle.getByLabel(ev, "prunedGenParticles");
       if(genHandle.isValid()){ gen = *genHandle;}
*/
       pat::MuonCollection muons;
       fwlite::Handle< pat::MuonCollection > muonsHandle;
       muonsHandle.getByLabel(ev, "slimmedMuons");
       if(muonsHandle.isValid()){ muons = *muonsHandle;}

       pat::ElectronCollection electrons;
       fwlite::Handle< pat::ElectronCollection > electronsHandle;
       electronsHandle.getByLabel(ev, "slimmedElectrons");
       if(electronsHandle.isValid()){ electrons = *electronsHandle;}

       pat::JetCollection jets;
       fwlite::Handle< pat::JetCollection > jetsHandle;
       jetsHandle.getByLabel(ev, "slimmedJets");
       if(jetsHandle.isValid()){ jets = *jetsHandle;}

       fwlite::Handle<edm::TriggerResults> triggerBits;
       triggerBits.getByLabel(ev,"TriggerResults","","HLT");
       const edm::TriggerNames &names = ev.triggerNames(*triggerBits);

       fwlite::Handle< double > rhoH;
       rhoH.getByLabel(ev,"fixedGridRhoFastjetAllCalo");
       rho = *rhoH;


/*
       fwlite::Handle<reco::GenParticleCollection> genParticles;
       genParticles.getByLabel(ev, "prunedGenParticles");
       int theNbOfGenParticles = genParticles->size();
       std::vector<reco::GenParticle> genElectrons;
*/
   //check the L1 candidateds
   //    
       fwlite::Handle<l1extra::L1EmParticleCollection> theL1extraParticles;
       theL1extraParticles.getByLabel(ev,"l1extraParticles","NonIsolated","RECO");

       fwlite::Handle<l1extra::L1EmParticleCollection> theL1extraParticlesIsolated;
       theL1extraParticlesIsolated.getByLabel(ev,"l1extraParticles","Isolated","RECO");

//       fwlite::Handle<l1t::L1EGamma> theL1extraEGM;
//       theL1extraEGM.getByLabel(ev,"caloStage2Digis","EGamma","RECO");
//
       fwlite::Handle< BXVector<l1t::EGamma> > dmxegs;
       dmxegs.getByLabel(ev,"caloStage2Digis","EGamma","RECO");

   
       fwlite::Handle< pat::TriggerObjectStandAloneCollection > triggerObjects;
       triggerObjects.getByLabel(ev, "selectedPatTrigger");

       std::vector<pat::TriggerObjectStandAlone> TagtriggerObj;
       std::vector<pat::TriggerObjectStandAlone> ProbetriggerObj;
       std::vector<pat::TriggerObjectStandAlone> passLeg2triggerObj;
       std::vector<pat::TriggerObjectStandAlone> passLeg1Obj;
       std::vector<pat::TriggerObjectStandAlone> passLeg2Obj;
       std::vector<pat::TriggerObjectStandAlone> passLeg3Obj;
       std::vector<pat::TriggerObjectStandAlone> passLeg4Obj;
       std::vector<pat::TriggerObjectStandAlone> passLeg5Obj;
       std::vector<pat::TriggerObjectStandAlone> passLeg6Obj;
       std::vector<pat::TriggerObjectStandAlone> passLeg7Obj;
       std::vector<pat::TriggerObjectStandAlone> passLeg8Obj;
       std::vector<pat::TriggerObjectStandAlone> passLeg9Obj;
       std::vector<pat::TriggerObjectStandAlone> passLeg10Obj;
       std::vector<pat::TriggerObjectStandAlone> passLeg11Obj;
       std::vector<pat::TriggerObjectStandAlone> passLeg12Obj;
       std::vector<pat::TriggerObjectStandAlone> passLeg14Obj;       
       std::vector<pat::TriggerObjectStandAlone> passLeg15Obj;

       // DERIVE WEIGHTS TO APPLY TO SAMPLE
       //

  /*     //pileup weight
       float weight = 1.0;
       double TotalWeight_plus = 1.0;
       double TotalWeight_minus = 1.0;
       float puWeight(1.0);

       if(isMC){          
          int ngenITpu = 0;

          fwlite::Handle< std::vector<PileupSummaryInfo> > puInfoH;
          puInfoH.getByLabel(ev, "addPileupInfo");
       //   puInfoH.getByLabel(ev, "slimmedAddPileupInfo");   
          for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++){
             if(it->getBunchCrossing()==0)      { ngenITpu += it->getPU_NumInteractions(); }
          }

          puWeight          = LumiWeights->weight(ngenITpu) * PUNorm[0];
          weight            = xsecWeight*puWeight;
          TotalWeight_plus  = PuShifters[utils::cmssw::PUUP  ]->Eval(ngenITpu) * (PUNorm[2]/PUNorm[0]);
          TotalWeight_minus = PuShifters[utils::cmssw::PUDOWN]->Eval(ngenITpu) * (PUNorm[1]/PUNorm[0]);
       } 
*/
//	PUWeight = puWeight;
//	lumi_weight = xsecWeight;


       ///////////////////////
       ///                 ///
       /// LEPTON ANALYSIS ///
       ///                 ///
       ///////////////////////


       std::vector<patUtils::GenericLepton> leptons;
       for(size_t l=0;l<electrons.size();l++){leptons.push_back(patUtils::GenericLepton(electrons[l]));}      
       for(size_t l=0;l<muons    .size();l++){leptons.push_back(patUtils::GenericLepton(muons    [l]));}      
       std::sort(leptons.begin(),   leptons.end(), utils::sort_CandidatesByPt);

       LorentzVector muDiff(0,0,0,0);
       std::vector<patUtils::GenericLepton> selLeptons, extraLeptons;
       //Request exactly two leptons
       
       //PRE-SELECTION based to pt value
       for(unsigned int j=0; j<leptons.size(); j++){
         if(leptons[j].pt() > 8) selLeptons.push_back(leptons[j]);
       }

       std::sort(selLeptons.begin(),   selLeptons.end(), utils::sort_CandidatesByPt);
       if(selLeptons.size() != 2) continue;

	 // opposite charge condition
	 
       if (! (selLeptons[0].charge()*selLeptons[1].charge() < 0) ) continue;

       if(abs(selLeptons[0].pdgId()) == 13 && abs(selLeptons[1].pdgId()) == 13) passOnlyMuon = true;
       if(abs(selLeptons[0].pdgId()) == 11 && abs(selLeptons[1].pdgId()) == 11) passOnlyEle = true;

//cout << "Selected Lepton Pdg Id = " << selLeptons[0].pdgId() << "   " << selLeptons[1].pdgId() << endl;
     
       //std::cout << "=> Passed the selection of the events containing only muon or ele" << std::endl; 
       bool passKinMu(false),passIdMu(false),passIsoMu(false);
       bool passKinEle(false),passIdEle(false),passIsoEle(false);

       /// MUON TAG ///
       //Logical tag for Tag and Probe Muon
       bool Tag(false), Probe(false),ProbeMuKin; 
       //Logical tag for Id Muons
       bool passIdTh(false), passIdSf(false), passIdLo(false), passIsoProbeMu(false);
       //Logical tag for Pass Probe Muons  
       bool PassProbeTh(false), PassProbeSf(false), PassProbeLo(false);

       /// ELE TAG ///
       //Logical tag for Tag and Probe Ele
       bool TagEle(false), ProbeEle(false), ProbeEleKin(false);
       bool TagTightIdSTD(false),ProbeTightIdSTD(false);
       bool ProbeISO(false), TagISO(false);
       //Logical tag for Id Ele
       bool passIdEleTh(false), passIdEleMd(false), passIdEleLo(false), passIsoProbeEle(false);
       //Logical tag for Pass Probe Ele
       bool PassProbeEleTh(false), PassProbeEleMd(false), PassProbeEleLo(false);
       bool ZPick(false);

       bool passControlLeg1(false), passControlLeg2(false), passControlLeg3(false) , passProbeLeg1(false), passProbeLeg2(false), passProbeLeg3(false), passProbeLeg4(false), passProbeLeg5(false), passProbeLeg6(false), passProbeLeg7(false), passProbeLeg8(false), passProbeLeg9(false), passProbeLeg10(false), passProbeLeg11(false),passProbeLeg12(false),passProbeLeg14(false),passProbeLeg15(false);
       std::string ControlFilterLeg;
      //////////////////////////
      ///                    ///
      /// ELECTRON EFFICENCY ///
      ///                    ///
      ////////////////////////// 
         if(passOnlyEle){
         int first  = rand()%2;
         int second = (first+1)%2;
         
         Nvtx = vtx.size();
         double etaf = selLeptons[first].el.superCluster()->eta();
         double ptf  = selLeptons[first].pt();
         double etas = selLeptons[second].el.superCluster()->eta();         
         double pts  = selLeptons[second].pt();
         double phif = selLeptons[first].phi();
         double phis = selLeptons[second].phi();

	 double ecalPFIsof = selLeptons[first].el.ecalPFClusterIso();	
	 double hcalPFIsof = selLeptons[first].el.hcalPFClusterIso();	
	 double trackIsof = selLeptons[first].el.trackIso();	
	 double ecalPFIsos = selLeptons[second].el.ecalPFClusterIso();	
	 double hcalPFIsos = selLeptons[second].el.hcalPFClusterIso();	
	 double trackIsos = selLeptons[second].el.trackIso();	
         double detaseedf = selLeptons[first].el.deltaEtaSeedClusterTrackAtCalo();
         double detaseeds = selLeptons[second].el.deltaEtaSeedClusterTrackAtCalo();
         double misshit = selLeptons[second].el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
         double dphiin = selLeptons[second].el.deltaPhiSuperClusterTrackAtVtx();
         double chi2 = selLeptons[second].el.gsfTrack()->normalizedChi2();        

        //cout << "detain seed = " << selLeptons[second].el.deltaEtaSeedClusterTrackAtCalo() << endl;
        //cout <<" mising hits = " << selLeptons[second].el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) << endl;
if(fabs(etas) <= 1.479) {
//         ProbeISO = ((ecalPFIsos-rho*0.165)/pts) < 0.160 && ((hcalPFIsos-rho*0.060)/pts) < 0.120 && (trackIsos/pts) < 0.08 && abs(detaseeds) < 0.004 && misshit < 1 && abs(dphiin) < 0.020;
         ProbeISO = ((ecalPFIsos-rho*0.165)/pts) < 0.160 && ((hcalPFIsos-rho*0.060)/pts) < 0.120 && (trackIsos/pts) < 0.08 && abs(detaseeds) < 0.004 && abs(dphiin) < 0.020;         
}
else if(fabs(etas) > 1.479 && fabs(etas) < 2.5) {
 //        ProbeISO = ((ecalPFIsos-rho*0.132)/pts) < 0.120 && ((hcalPFIsos-rho*0.131)/pts) < 0.120 && (trackIsos/pts) < 0.08 && misshit < 1 && abs(chi2) < 3;
          ProbeISO = ((ecalPFIsos-rho*0.132)/pts) < 0.120 && ((hcalPFIsos-rho*0.131)/pts) < 0.120 && (trackIsos/pts) < 0.08 && abs(chi2) < 3;
}

//	 TagISO = ecalPFIsof < 0.45 && hcalPFIsof < 0.25 && trackIsof < 0.2;
//	 ProbeISO = ecalPFIsos < 0.45 && hcalPFIsos < 0.25 && trackIsos < 0.2;

         for ( pat::TriggerObjectStandAlone obj: *triggerObjects ) {
             obj.unpackPathNames(names);

//DoubleElectron Trigger


               std::string leg1Filter("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter"); // Leg 1 of HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2
               std::string leg2Filter("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter"); // Leg 2 of HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2 
               std::string leg3Filter("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter"); // Leg 1 of HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2
               std::string leg4Filter("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter"); // Leg 2 of HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2

//SingleElectron Trigger
                                     
             std::string leg5Filter("hltEle27noerWPLooseGsfTrackIsoFilter"); //HLT_Ele27_WPLoose_Gsf_v*
             std::string leg6Filter("hltEle27WPTightGsfTrackIsoFilter"); //HLT_Ele27_WPTight_Gsf_v*
             std::string leg7Filter("hltEle23WPLooseGsfTrackIsoFilter"); //HLT_Ele23_WPLoose_Gsf_v*
             std::string leg8Filter("hltEle25WPTightGsfTrackIsoFilter"); //HLT_Ele25_WPTight_Gsf_v*
             std::string leg9Filter("hltEle25erWPLooseGsfTrackIsoFilter"); //HLT_Ele25_eta2p1WPLoose_Gsf_v*
             std::string leg10Filter("hltEle25erWPTightGsfTrackIsoFilter"); //HLT_Ele25_eta2p1_WPTight_Gsf_v*
             std::string leg11Filter("hltEle35WPLooseGsfTrackIsoFilter"); //HLT_Ele35_WPLoose_Gsf_v*
             std::string leg12Filter("hltEle27erWPLooseGsfTrackIsoFilter"); //HLT_Ele27_eta2p1_WPLoose_Gsf_v*
             std::string leg14Filter("hltEle45WPLooseGsfTrackIsoFilter"); //HLT_Ele45_WPLoose_Gsf_v*
             std::string leg15Filter("hltEle27WPTightGsfTrackIsoFilter"); //HLT_Ele27_WPTight_Gsf_v

               // ======================= Triggers for Tag ================================  
         
               ControlFilterLeg = "hltEle27WPTightGsfTrackIsoFilter"; 
  
               passControlLeg1   = passFilter(obj,ControlFilterLeg);
             
             passProbeLeg1     = passFilter(obj,leg1Filter);
             passProbeLeg2     = passFilter(obj,leg2Filter);
             passProbeLeg3     = passFilter(obj,leg3Filter);
             passProbeLeg4     = passFilter(obj,leg4Filter);
             passProbeLeg5     = passFilter(obj,leg5Filter);
             passProbeLeg6     = passFilter(obj,leg6Filter);
             passProbeLeg7     = passFilter(obj,leg7Filter);
             passProbeLeg8     = passFilter(obj,leg8Filter);
             passProbeLeg9     = passFilter(obj,leg9Filter);
             passProbeLeg10     = passFilter(obj,leg10Filter);
             passProbeLeg11     = passFilter(obj,leg11Filter);
             passProbeLeg12     = passFilter(obj,leg12Filter);
             passProbeLeg14     = passFilter(obj,leg14Filter);
             passProbeLeg15     = passFilter(obj,leg15Filter);

             if  (passProbeLeg1){
             passLeg1Obj.push_back(obj);
                }

             if  (passProbeLeg2){
           passLeg2Obj.push_back(obj);
            }

            if  (passProbeLeg3){
             passLeg3Obj.push_back(obj);
                }

            if  (passProbeLeg4){
             passLeg4Obj.push_back(obj);
            }

            if  (passProbeLeg5){
             passLeg5Obj.push_back(obj);
                }

             if  (passProbeLeg6){
           passLeg6Obj.push_back(obj);
            }

            if  (passProbeLeg7){
             passLeg7Obj.push_back(obj);
                }

            if  (passProbeLeg8){
             passLeg8Obj.push_back(obj);
                }

            if  (passProbeLeg9){
             passLeg9Obj.push_back(obj);
                }

            if  (passProbeLeg10){
             passLeg10Obj.push_back(obj);
                }

            if  (passProbeLeg11){
             passLeg11Obj.push_back(obj);
                }

            if  (passProbeLeg12){
             passLeg12Obj.push_back(obj);
                }

            if  (passProbeLeg14){
             passLeg14Obj.push_back(obj);
                }

            if  (passProbeLeg15){
             passLeg15Obj.push_back(obj);
                }

             if  (passControlLeg1){             
             TagtriggerObj.push_back(obj);
                }
        }

                passKinEle = (((abs(etaf) >= 0 && abs(etaf) <= 1.4442) || (abs(etaf) >= 1.5660 && abs(etaf) <= 2.5)) && ptf > 30);


	 TagTightIdSTD		= patUtils::passId(electronVidTightId, myEvent, selLeptons[first].el);
	 ProbeTightIdSTD	= patUtils::passId(electronVidTightId, myEvent, selLeptons[second].el);
       
	 bool passTagdRCut(false);

         double dRTag =DeltaRtrig(selLeptons[first].el,TagtriggerObj);

         if (dRTag < 0.1) {
         passTagdRCut = true;}

         if(passKinEle && TagTightIdSTD && passTagdRCut) TagEle = true;

         if(!TagEle)  continue;
	tagPT = ptf;
	tagEta = etaf;

	//Selection of the Probe
         bool   ProbeEle = false;
        ProbeEleKin = ((abs(etas) >= 0 && abs(etas) <= 2.5) && pts > 0);//Fix it       
           if((ProbeTightIdSTD && ProbeEleKin && ProbeISO)) ProbeEle = true;
        //   if((ProbeTightIdSTD && ProbeEleKin)) ProbeEle = true;


	 if (ProbeEle) {
//count_probe++;
	 TLorentzVector lep1(selLeptons[first].px(),selLeptons[first].py(),selLeptons[first].pz(),selLeptons[first].energy());
	 TLorentzVector lep2(selLeptons[second].px(),selLeptons[second].py(),selLeptons[second].pz(),selLeptons[second].energy());  
         double mass = (lep1+lep2).M();
         ZMass = mass;
         ZPick = mass > 60 && mass < 120;
 	 if(!ZPick)  continue;

totalProbePT = pts;
totalProbeEta = etas;
totalProbePhi = phis;
//count_probe++;

//cout  << "Probe PT   = " << pts << "  " << "Probe Eta  = " << etas << endl;


//now loop on the L1 particles

        l1extra::L1EmParticleCollection::const_iterator itrEm;
        float maxL1matched = -1;
        for( itrEm = theL1extraParticles->begin();  itrEm != theL1extraParticles->end(); ++itrEm ) {
            float deltaR = sqrt(pow(itrEm->eta()-etas,2)+ pow(acos(cos(itrEm->phi()-phis)),2));
            if (deltaR< 0.5) {
                if (itrEm->pt()>maxL1matched) maxL1matched = itrEm->pt();
            }
        }

        l1extra::L1EmParticleCollection::const_iterator itrEmIso;
        for( itrEmIso = theL1extraParticlesIsolated->begin();  itrEmIso != theL1extraParticlesIsolated->end(); ++itrEmIso ) {
            float deltaR = sqrt(pow(itrEmIso->eta()-etas,2)+ pow(acos(cos(itrEmIso->phi()-phis)),2));
            if (deltaR< 0.5) {
                if (itrEmIso->pt()>maxL1matched) maxL1matched = itrEmIso->pt();
            }
        }


         bool    passleg1dRCut = false;
//if(maxL1matched >= 15) {
//cout << "maximum L1 Threshold = " << maxL1matched << endl;
         //Counting the Passing Prob
        
         double  dRleg1 = DeltaRtrig(selLeptons[second].el,passLeg1Obj);
         if(dRleg1 < 0.1) passleg1dRCut = true;
//} // L1 Condition

         bool    passleg2dRCut = false;
         double  dRleg2 = DeltaRtrig(selLeptons[second].el,passLeg2Obj);

         if(dRleg2 < 0.1) passleg2dRCut = true;
 
         bool    passleg3dRCut = false;
         double  dRleg3 = DeltaRtrig(selLeptons[second].el,passLeg3Obj);
         if(dRleg3 < 0.1) passleg3dRCut = true;

         bool    passleg4dRCut = false;
         double  dRleg4 = DeltaRtrig(selLeptons[second].el,passLeg4Obj);
         if(dRleg4 < 0.1) passleg4dRCut = true;

         bool    passleg5dRCut = false;
         double  dRleg5 = DeltaRtrig(selLeptons[second].el,passLeg5Obj);
         if(dRleg5 < 0.1) passleg5dRCut = true;

         bool    passleg6dRCut = false;
         double  dRleg6 = DeltaRtrig(selLeptons[second].el,passLeg6Obj);
         if(dRleg6 < 0.1) passleg6dRCut = true;

         bool    passleg7dRCut = false;
         double  dRleg7 = DeltaRtrig(selLeptons[second].el,passLeg7Obj);
         if(dRleg7 < 0.1) passleg7dRCut = true;

         bool    passleg8dRCut = false;
         double  dRleg8 = DeltaRtrig(selLeptons[second].el,passLeg8Obj);
         if(dRleg8 < 0.1) passleg8dRCut = true;

         bool    passleg9dRCut = false;
         double  dRleg9 = DeltaRtrig(selLeptons[second].el,passLeg9Obj);
         if(dRleg9 < 0.1) passleg9dRCut = true;

         bool    passleg10dRCut = false;
         double  dRleg10 = DeltaRtrig(selLeptons[second].el,passLeg10Obj);
         if(dRleg10 < 0.1) passleg10dRCut = true;

         bool    passleg11dRCut = false;
         double  dRleg11 = DeltaRtrig(selLeptons[second].el,passLeg11Obj);
         if(dRleg11 < 0.1) passleg11dRCut = true;

         bool    passleg12dRCut = false;
         double  dRleg12 = DeltaRtrig(selLeptons[second].el,passLeg12Obj);
         if(dRleg12 < 0.1) passleg12dRCut = true;

         bool    passleg14dRCut = false;
         double  dRleg14 = DeltaRtrig(selLeptons[second].el,passLeg14Obj);
         if(dRleg14 < 0.1) passleg14dRCut = true;

         bool    passleg15dRCut = false;
         double  dRleg15 = DeltaRtrig(selLeptons[second].el,passLeg15Obj);
//cout << "dRleg15 = " << dRleg15 << endl;
         if(dRleg15 < 0.1) passleg15dRCut = true;
//cout << "passleg15dRCut = " << passleg15dRCut << endl;

         	  if(passleg1dRCut){
//	count_passprobeleg1++;
	passProbePTLeg1 = pts;
	passProbeEtaLeg1 = etas;
        passProbePhiLeg1 = phis;
          }
           	if(passleg2dRCut){
	passProbePTLeg2 = pts;
	passProbeEtaLeg2 = etas;
        passProbePhiLeg2 = phis;
             }

          	 if(passleg3dRCut){
        passProbePTLeg3 = pts;
        passProbeEtaLeg3 = etas;
        passProbePhiLeg3 = phis;
                                        }

           	if(passleg4dRCut){
        passProbePTLeg4 = pts;
        passProbeEtaLeg4 = etas;
        passProbePhiLeg4 = phis;
                                        }

                if(passleg5dRCut){
  //      count_passprobeleg1++;
        passProbePTLeg5 = pts;
        passProbeEtaLeg5 = etas;
        passProbePhiLeg5 = phis;
                                        }

                if(passleg6dRCut){
        passProbePTLeg6 = pts;
        passProbeEtaLeg6 = etas;
        passProbePhiLeg6 = phis;
   //     count_passprobeleg2++;
                                        }

                if(passleg7dRCut){
        passProbePTLeg7 = pts;
        passProbeEtaLeg7 = etas;
        passProbePhiLeg7 = phis;
 //       count_passprobeleg3++;
                                        }

                if(passleg8dRCut){
        passProbePTLeg8 = pts;
        passProbeEtaLeg8 = etas;
        passProbePhiLeg8 = phis;
                                                }
 
    
                if(passleg9dRCut){
        passProbePTLeg9 = pts;
        passProbeEtaLeg9 = etas;
        passProbePhiLeg9 = phis;
                                                }

                if(passleg10dRCut){
        passProbePTLeg10 = pts;
        passProbeEtaLeg10 = etas;
        passProbePhiLeg10 = phis;
                                                }
                                           

                if(passleg11dRCut){
        passProbePTLeg11 = pts;
        passProbeEtaLeg11 = etas;
        passProbePhiLeg11 = phis;
                                                }

                if(passleg12dRCut){
        passProbePTLeg12 = pts;
        passProbeEtaLeg12 = etas;
        passProbePhiLeg12 = phis;
                                                }

                if(passleg11dRCut || passleg12dRCut) {
        passProbePTLeg13 = pts;
        passProbeEtaLeg13 = etas;
        passProbePhiLeg13 = phis;
                                                }

                if(passleg14dRCut){
        passProbePTLeg14 = pts;
        passProbeEtaLeg14 = etas;
        passProbePhiLeg14 = phis;
                                                }

                if(passleg12dRCut || passleg14dRCut) {
        passProbePTLeg15 = pts;
        passProbeEtaLeg15 = etas;
        passProbePhiLeg15 = phis;
                                                }

                if(passleg15dRCut) {
        passProbePTLeg16 = pts;
        passProbeEtaLeg16 = etas;
        passProbePhiLeg16 = phis;
                                                }

//================== L1 Trigger Efficiency
float maxL1MatchedNorm = -1;
float maxL1MatchedIso = -1;
float maxL1MatchedER = -1;
bool L1Norm(false),L1Iso(false),L1ER(false);
     for ( int ibx=dmxegs->getFirstBX(); ibx<=dmxegs->getLastBX(); ++ibx) { 
     for ( auto itr = dmxegs->begin(ibx); itr != dmxegs->end(ibx); ++itr ) {
         float deltaR = sqrt(pow(itr->eta()-etas,2)+ pow(acos(cos(itr->phi()-phis)),2));
         if (deltaR < 0.5) {
 //        cout << "ET =  : " << itr->et() << " pt=" << itr->pt() << " eta=" << itr->eta() <<  " phi=" << itr->phi() << std::endl;
  //       cout << "Iso = " << itr->hwIso() << endl; 
         if(itr->pt() > maxL1MatchedNorm) maxL1MatchedNorm = itr->pt();
         if(itr->hwIso() == 1 && itr->pt()>maxL1MatchedIso) maxL1MatchedIso = itr->pt();
         if(itr->hwIso() == 1 && fabs(itr->eta()) < 2.1 && itr->pt()>maxL1MatchedER) maxL1MatchedER = itr->pt();
            }
       }
     }

	if(maxL1MatchedNorm > 40) L1Norm = true;
	if(maxL1MatchedIso > 24) L1Iso = true;
	if(maxL1MatchedER > 22) L1ER = true;

if(L1Norm){ 
passL1PTNorm = pts;
passL1EtaNorm = etas;
}
if(L1Iso){
passL1PTIso = pts;
passL1EtaIso = etas;
}
if(L1ER){
passL1PTER = pts;
passL1EtaER = etas;
}

if(L1Norm == 1 || L1Iso == 1 || L1ER == 1){
passL1PT = pts;
passL1Eta = etas;
}
   //  cout << "maxL1MatchedNorm = " << maxL1MatchedNorm << "   " << "maxL1MatchedIso   = " << maxL1MatchedIso << "  " << "maxL1MatchedER = " << maxL1MatchedER << endl;
   //  cout << "L1Norm = " << L1Norm <<  "  "  << "L1Iso = " << L1Iso << " " << "L1ER = " << L1ER << endl;   
//cout << "*****************************" << endl;
	ntuple->Fill();
	} //If you find Probe Electron
       }  
       
  }

//cout << "Number of Total Probes = " << count_probe << endl;
//cout << "Number of passing Probes HLT_Ele27_Loose = " << count_passprobeleg1 << endl;
//cout << "Number of passing Probes HLT_Ele27_Tight = " << count_passprobeleg2 << endl;
//cout << "Number of passing Probes HLT_Ele23_Loose = " << count_passprobeleg3 << endl;


printf("Results save in %s\n", outUrl.Data());
  
  ofile->Write();
  ofile->Close();
}
