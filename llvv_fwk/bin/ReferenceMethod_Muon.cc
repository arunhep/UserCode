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




double DeltaRMC(const pat::Electron& el, std::vector<reco::GenParticle> object)
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

  double Mu1PT;
  double Mu1Eta;
  double Mu1Phi;
  double Mu2PT;
  double Mu2Eta;
  double Mu2Phi;
  int    ReferenceTriggerPass;
  int    isTriggerSoupPass;
  int    pass_Mu17_8;
  int    pass_IsoMu20;
  int    pass_IsoMu22;
  int    pass_IsoMu27;
  int    pass_IsoTkMu27;
  int    pass_Mu17_Tk8;
  int    nVtx;
  double Zmass;
  int    ReferenceTriggerPass_New;
  int    isTriggerSoupPass_New;
  int    pass_Mu17_8_New;
  int    pass_IsoMu20_New;
  int    pass_IsoMu22_New;
  int    pass_IsoMu27_New;
  int    pass_IsoTkMu27_New;
  int    pass_Mu17_Tk8_New;


  TBranch *b_Mu1PT = ntuple->Branch("Mu1PT",&Mu1PT,"Mu1PT/D");
  TBranch *b_Mu1Eta = ntuple->Branch("Mu1Eta",&Mu1Eta,"Mu1Eta/D");
  TBranch *b_Mu1Phi = ntuple->Branch("Mu1Phi",&Mu1Phi,"Mu1Phi/D");
  TBranch *b_Mu2PT = ntuple->Branch("Mu2PT",&Mu2PT,"Mu2PT/D");
  TBranch *b_Mu2Eta = ntuple->Branch("Mu2Eta",&Mu2Eta,"Mu2Eta/D");
  TBranch *b_Mu2Phi = ntuple->Branch("Mu2Phi",&Mu2Phi,"Mu2Phi/D");

  TBranch *b_ReferenceTriggerPass= ntuple->Branch("ReferenceTriggerPass",&ReferenceTriggerPass,"ReferenceTriggerPass/I");
  TBranch *b_isTriggerSoupPass= ntuple->Branch("isTriggerSoupPass",&isTriggerSoupPass,"isTriggerSoupPass/I");
  TBranch *b_pass_Mu17_8= ntuple->Branch("pass_Mu17_8",&pass_Mu17_8,"pass_Mu17_8/I");
  TBranch *b_pass_IsoMu20= ntuple->Branch("pass_IsoMu20",&pass_IsoMu20,"pass_IsoMu20/I");
  TBranch *b_pass_IsoMu22= ntuple->Branch("pass_IsoMu22",&pass_IsoMu22,"pass_IsoMu22/I");
  TBranch *b_pass_IsoMu27= ntuple->Branch("pass_IsoMu27",&pass_IsoMu27,"pass_IsoMu27/I");
  TBranch *b_pass_IsoTkMu27= ntuple->Branch("pass_IsoTkMu27",&pass_IsoTkMu27,"pass_IsoTkMu27/I");
  TBranch *b_pass_Mu17_Tk8= ntuple->Branch("pass_Mu17_Tk8",&pass_Mu17_Tk8,"pass_Mu17_Tk8/I");
  TBranch *b_Zmass = ntuple->Branch("Zmass",&Zmass,"Zmass/D");
  TBranch *b_nVtx = ntuple->Branch("nVtx",&nVtx,"nVtx/I");
  TBranch *b_pass_Mu17_8_New= ntuple->Branch("pass_Mu17_8_New",&pass_Mu17_8_New,"pass_Mu17_8_New/I");
  TBranch *b_pass_IsoMu20_New= ntuple->Branch("pass_IsoMu20_New",&pass_IsoMu20_New,"pass_IsoMu20_New/I");
  TBranch *b_pass_IsoMu22_New= ntuple->Branch("pass_IsoMu22_New",&pass_IsoMu22_New,"pass_IsoMu22_New/I");
  TBranch *b_pass_IsoMu27_New= ntuple->Branch("pass_IsoMu27_New",&pass_IsoMu27_New,"pass_IsoMu27_New/I");
  TBranch *b_pass_IsoTkMu27_New= ntuple->Branch("pass_IsoTkMu27_New",&pass_IsoTkMu27_New,"pass_IsoTkMu27_New/I");
  TBranch *b_pass_Mu17_Tk8_New= ntuple->Branch("pass_Mu17_Tk8_New",&pass_Mu17_Tk8_New,"pass_Mu17_Tk8_New/I");
  TBranch *b_ReferenceTriggerPass_New= ntuple->Branch("ReferenceTriggerPass_New",&ReferenceTriggerPass_New,"ReferenceTriggerPass_New/I");
  TBranch *b_isTriggerSoupPass_New= ntuple->Branch("isTriggerSoupPass_New",&isTriggerSoupPass_New,"isTriggerSoupPass_New/I");
  // good lumi mask
  lumiUtils::GoodLumiFilter goodLumiFilter(runProcess.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> >("lumisToProcess", std::vector<edm::LuminosityBlockRange>()));

  bool isV0JetsMC(isMC && (dtag.Contains("DYJetsToLL_50toInf") || dtag.Contains("_WJets"))); // #FIXME should be reactivated as soon as we have exclusive jet samples
  bool isData_DoubleEle = string(dtag.Data()).find("Data13TeV_DoubleElectron");

  const edm::ParameterSet& myVidElectronIdConf          = runProcess.getParameterSet("electronidparas");
  const edm::ParameterSet& myVidElectronTightIdWPConf   = myVidElectronIdConf.getParameterSet("tight");
  const edm::ParameterSet& myVidElectronMediumIdWPConf  = myVidElectronIdConf.getParameterSet("medium");
  const edm::ParameterSet& myVidElectronLooseIdWPConf   = myVidElectronIdConf.getParameterSet("loose");

  VersionedPatElectronSelector electronVidTightId(myVidElectronTightIdWPConf);
  VersionedPatElectronSelector electronVidMediumId(myVidElectronMediumIdWPConf);
  VersionedPatElectronSelector electronVidLooseId(myVidElectronLooseIdWPConf);


  // List of Triggers
  std::string TriggerList[6] = {"HLT_IsoMu20_v*","HLT_IsoMu22_v*","HLT_Mu17_TrkIsoVVL_Mu8_Trk IsoVVL_DZ_v*","HLT_IsoMu27_v*","HLT_IsoTkMu27_v3","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*"};
// std::string TriggerList[1] = {"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"};  
  
//  std::string ReferenceTrigger = "HLT_Ele27_WPLoose_Gsf_v*";
  std::string ReferenceTrigger = "HLT_Mu17_v*"; 
 
  //##############################################
  //######## GET READY FOR THE EVENT LOOP ########
  //##############################################

  fwlite::ChainEvent ev(urls);
  const size_t totalEntries= ev.size();
  
  //MC normalization (to 1/pb)
  double xsecWeight = xsec/totalEntries;
  if(!isMC) xsecWeight=1.0;

  //pileup weighting
//  edm::LumiReWeighting* LumiWeights = NULL;
//  utils::cmssw::PuShifter_t PuShifters;
  double PUNorm[] = {1,1,1};
/* DISABLE  if(isMC){
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
 */
  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE


  //##############################################
  //########           EVENT LOOP         ########
  //##############################################
  //loop on all the events
  printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
  printf("Scanning the ntuple :");
  int treeStep(totalEntries/50);

//if(ev%10000==0) cout << ev << endl;

cout << "totalEntries = " << totalEntries << endl; 
  for( size_t iev=0; iev< totalEntries; iev++){
       ev.to(iev); //load the event content from the EDM file
       edm::EventBase const & myEvent = ev;

Mu1PT	= -999.;
Mu1Eta	= -999.;
Mu1Phi = -999.; 
Mu2PT	= -999.;
Mu2Eta	= -999.;
Mu2Phi = -999.;
ReferenceTriggerPass = 0;
isTriggerSoupPass = 0; 
pass_IsoMu20 =0;
pass_IsoMu22 =0;
pass_Mu17_8 =0;
pass_IsoMu27 =0;
pass_IsoTkMu27 =0;
pass_Mu17_Tk8 =0;
nVtx    = -999;
Zmass	= -999.;
ReferenceTriggerPass_New = 0;
isTriggerSoupPass_New = 0;
pass_IsoMu20_New =0;
pass_IsoMu22_New =0;
pass_Mu17_8_New =0;
pass_IsoMu27_New =0;
pass_IsoTkMu27_New =0;
pass_Mu17_Tk8_New =0;

       //apply trigger and require compatibilitiy of the event with the PD
       edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
       edm::TriggerResultsByName tr_New = ev.triggerResultsByName("TEST");
       if(!tr.isValid())return false;
       if(!tr_New.isValid())return false;
      

       //std::cout << "Definition of the Booleian variable" << std::endl;
 
       bool passOnlyMuon(false);
       bool passOnlyEle(false); 

              
       //load all the objects we will need to access
       reco::VertexCollection vtx;
       fwlite::Handle< reco::VertexCollection > vtxHandle; 
       vtxHandle.getByLabel(ev, "offlineSlimmedPrimaryVertices");
       if(vtxHandle.isValid()){ vtx = *vtxHandle;}

       reco::GenParticleCollection gen;
       fwlite::Handle< reco::GenParticleCollection > genHandle;
       genHandle.getByLabel(ev, "prunedGenParticles");
       if(genHandle.isValid()){ gen = *genHandle;}

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

       fwlite::Handle<edm::TriggerResults> triggerBits_New;
       triggerBits_New.getByLabel(ev,"TriggerResults","","TEST");
       const edm::TriggerNames &names_New = ev.triggerNames(*triggerBits_New);

   //check the L1 candidateds
   //    
       fwlite::Handle<l1extra::L1EmParticleCollection> theL1extraParticles;
       theL1extraParticles.getByLabel(ev,"l1extraParticles","NonIsolated","RECO");

       fwlite::Handle<l1extra::L1EmParticleCollection> theL1extraParticlesIsolated;
       theL1extraParticlesIsolated.getByLabel(ev,"l1extraParticles","Isolated","RECO");

       fwlite::Handle< pat::TriggerObjectStandAloneCollection > triggerObjects;
       triggerObjects.getByLabel(ev, "selectedPatTrigger");


       //pileup weight
       float weight = 1.0;


       ///////////////////////
       ///                 ///
       /// LEPTON ANALYSIS ///
       ///                 ///
       ///////////////////////


       //start by merging electrons and muons
       //std::cout << "=>Merging Leptons" << std::endl;
       std::vector<patUtils::GenericLepton> leptons;
  //     std::vector<int> index;
       for(size_t l=0;l<electrons.size();l++){leptons.push_back(patUtils::GenericLepton(electrons[l]));}      
       for(size_t l=0;l<muons    .size();l++){leptons.push_back(patUtils::GenericLepton(muons    [l]));}      
       std::sort(leptons.begin(),   leptons.end(), utils::sort_CandidatesByPt);

       std::vector<patUtils::GenericLepton> selLeptons, extraLeptons;
       
//cout << "total number of leptons = " << leptons.size() << endl;

       //PRE-SELECTION based to pt value
       for(unsigned int j=0; j<leptons.size(); j++){
         if(leptons[j].pt() > 0) selLeptons.push_back(leptons[j]);
       }

       std::sort(selLeptons.begin(),   selLeptons.end(), utils::sort_CandidatesByPt);
       //Request atleast two leptons
       if(selLeptons.size() < 2) continue;

 //      cout << "Number of leptons  = " << selLeptons.size() << endl;

	 // opposite charge condition
	 
       if (! (selLeptons[0].charge()*selLeptons[1].charge() < 0) ) continue;

//       for(int s = 0 ; s < selLeptons.size(); s++) {
//cout << "selLeptons[s].pdgId() = " << selLeptons[s].pdgId() << endl;
//       if(abs(selLeptons[s].pdgId()) == 11) {index.push_back(s);}

//}

//cout << "Number of electrons in an event = " << index.size() << endl;

       if(abs(selLeptons[0].pdgId()) == 13 && abs(selLeptons[1].pdgId()) == 13) passOnlyMuon = true;
       if(abs(selLeptons[0].pdgId()) == 11 && abs(selLeptons[1].pdgId()) == 11) passOnlyEle = true;

//if(index.size() >= 2) passOnlyEle = true;

     
       //std::cout << "=> Passed the selection of the events containing only muon or ele" << std::endl; 
       bool passKinMu1(false),passIdMu1(false),passIsoMu1(false);
       bool passKinMu2(false),passIdMu2(false),passIsoMu2(false);
       bool ProbeISO(false), TagISO(false), passTightIdSTDMu1(false),passTightIdSTDMu2(false);
       bool ZPick(false);

      //////////////////////////
      ///                    ///
      /// MUON EFFICENCY ///
      ///                    ///
      ////////////////////////// 
      //
 //     cout << "passOnlyEle = " << passOnlyEle << endl;
         if(passOnlyMuon){

         int first  = rand()%2;
         int second = (first+1)%2;
//         int first  = 0;
//         int second = 1;
   
//for(int t = 0; t < index.size(); t++){
//cout << "Electrons PT = " << selLeptons[index[t]].el.pt() << endl;
//} 
      
         double ptf  = selLeptons[first].mu.pt();
         double etaf = selLeptons[first].mu.eta();
         double phif = selLeptons[first].mu.phi();
         double pts  = selLeptons[second].mu.pt();
         double etas = selLeptons[second].mu.eta();
         double phis = selLeptons[second].mu.phi();




         passKinMu1 = ((abs(etaf) >= 0 && abs(etaf) <= 2.4));
  //       cout << "passKinEle1  = " << passKinEle1 << endl;
     //    passTightIdSTDEle1       = patUtils::passId(electronVidTightId, myEvent, selLeptons[first].el);     
  //       passTightIdSTDEle1       = patUtils::passId(electronVidMediumId, myEvent, selLeptons[first].el);
  
//  cout << "passTightIdSTDEle1   = " << passTightIdSTDEle1 << endl;
 
         passIdMu1  = patUtils::passId(selLeptons[first].mu, vtx[0], patUtils::llvvMuonId::StdTight);
         passIsoMu1 = patUtils::passIso(selLeptons[first].mu,  patUtils::llvvMuonIso::Tight);

         passKinMu2 = ((abs(etas) >= 0 && abs(etas) <= 2.4));       
   //      cout << "passKinEle2 = " << passKinEle2 << endl;  
 //        passTightIdSTDEle2       = patUtils::passId(electronVidTightId, myEvent, selLeptons[second].el);
//         passTightIdSTDEle2       = patUtils::passId(electronVidMediumId, myEvent, selLeptons[second].el);

         passIdMu2 = patUtils::passId(selLeptons[second].mu, vtx[0], patUtils::llvvMuonId::StdTight);
         passIsoMu2 = patUtils::passIso(selLeptons[second].mu,  patUtils::llvvMuonIso::Tight);

	 bool GoodMu1 = false;
	 bool GoodMu2 = false;

	 if(passKinMu1 && passIdMu1 && passIsoMu1) GoodMu1 = true;
	 if(passKinMu2 && passIdMu2 && passIsoMu2) GoodMu2 = true;
	 
	 
	 if((GoodMu1 && GoodMu2) != true) continue;

	 TLorentzVector lep1(selLeptons[first].px(),selLeptons[first].py(),selLeptons[first].pz(),selLeptons[first].energy());
	 TLorentzVector lep2(selLeptons[second].px(),selLeptons[second].py(),selLeptons[second].pz(),selLeptons[second].energy());  
         double mass = (lep1+lep2).M();

   //      ZPick = mass > 60 && mass < 120;
//	 if(!ZPick) continue;	 
	 Mu1PT = ptf;
	 Mu1Eta = etaf;
	 Mu1Phi = phif;
	 Mu2PT = pts;
	 Mu2Eta = etas;
	 Mu2Phi = phis;
	 Zmass = mass;

	 bool refTrig(false);
	 bool soupOfTrig(false);	 
	 bool IsoMu20(false);
	 bool IsoMu22(false);
	 bool Mu17_8(false);
         bool IsoMu27(false);
         bool IsoTkMu27(false);
         bool Mu17_Tk8(false);

         bool refTrig_New(false);
         bool soupOfTrig_New(false);
         bool IsoMu20_New(false);
         bool IsoMu22_New(false);
         bool Mu17_8_New(false);
         bool IsoMu27_New(false);
         bool IsoTkMu27_New(false);
         bool Mu17_Tk8_New(false);

	 refTrig   = utils::passTriggerPatterns(tr, ReferenceTrigger);
         IsoMu20     = utils::passTriggerPatterns(tr, TriggerList[0]);
         IsoMu22  = utils::passTriggerPatterns(tr, TriggerList[1]);
         Mu17_8  = utils::passTriggerPatterns(tr, TriggerList[2]);
         IsoMu27  = utils::passTriggerPatterns(tr, TriggerList[3]);
         IsoTkMu27  = utils::passTriggerPatterns(tr, TriggerList[4]);
         Mu17_Tk8  = utils::passTriggerPatterns(tr, TriggerList[5]);

	 soupOfTrig= (IsoMu20 || Mu17_8);

	 if(refTrig)	 { ReferenceTriggerPass = 1;}
	 if(soupOfTrig)  { isTriggerSoupPass  = 1;}
	 if(IsoMu20) 	 { pass_IsoMu20  = 1;}
	 if(IsoMu22) 	 { pass_IsoMu22  = 1;}
	 if(Mu17_8)  	 { pass_Mu17_8  = 1;}
         if(IsoMu27)     { pass_IsoMu27  = 1;}
         if(IsoTkMu27)     { pass_IsoTkMu27  = 1;}
         if(Mu17_Tk8)      { pass_Mu17_Tk8  = 1;}

         refTrig_New   = utils::passTriggerPatterns(tr_New, ReferenceTrigger);
         IsoMu20_New     = utils::passTriggerPatterns(tr_New, TriggerList[0]);
         IsoMu22_New  = utils::passTriggerPatterns(tr_New, TriggerList[1]);
         Mu17_8_New  = utils::passTriggerPatterns(tr_New, TriggerList[2]);
         IsoMu27_New  = utils::passTriggerPatterns(tr_New, TriggerList[3]);
         IsoTkMu27_New  = utils::passTriggerPatterns(tr_New, TriggerList[4]);
         Mu17_Tk8_New  = utils::passTriggerPatterns(tr_New, TriggerList[5]);
         soupOfTrig_New= (IsoMu20_New || Mu17_8_New);

//cout << "Ele23_New  = " << Ele23_New << endl;

         if(refTrig_New)     { ReferenceTriggerPass_New = 1;}
         if(soupOfTrig_New)  { isTriggerSoupPass_New  = 1;}
         if(IsoMu20_New)       { pass_IsoMu20_New  = 1;}
         if(IsoMu22_New)    { pass_IsoMu22  = 1;}
         if(Mu17_8_New)    { pass_Mu17_8_New  = 1;}
         if(IsoMu27_New)     { pass_IsoMu27_New  = 1;}
         if(IsoTkMu27_New)     { pass_IsoTkMu27_New  = 1;}
         if(Mu17_Tk8_New)      { pass_Mu17_Tk8_New  = 1;}

	
	ntuple->Fill();
       } 

//cout <<"****************************************" << endl;
  }


  //##############################################
  //########     SAVING HISTO TO FILE     ########
  //##############################################
  //save control plots to file

//TFile *ofile=TFile::Open(outUrl, "recreate");
printf("Results save in %s\n", outUrl.Data());
  
  //save all to the file
//  mon.Write();
  ofile->Write();
  ofile->Close();
//std::cout<<"evcounter : "<<evcounter<<" counter1 : "<<counter1<<" counter2 : "<<counter2;
  //if(outTxtFile)fclose(outTxtFile);
}
