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

  double Ele1PT;
  double Ele1Eta;
  double Ele1EtaSC;
  double Ele1Phi;
  double Ele2PT;
  double Ele2Eta;
  double Ele2EtaSC;
  double Ele2Phi;
  int    ReferenceTriggerPass;
  int    isTriggerSoupPass;
  int    pass_Ele23_12;
  int    pass_Ele23_12_DZ;
  int    pass_Ele17_12;
  int    pass_Ele23;
  int    pass_Ele27Loose;
  int    pass_Ele27LooseEtar;
  int    pass_Ele27TightEtar;
  int    nVtx;
  double Zmass;
  int    passdRleg1;
  int    passdRleg2;
  int    passdRleg3;
  int    passdRleg1N;
  int    passdRleg2N;
  int    passdRleg3N;
  double diEle_DZ;

  TBranch *b_diEle_DZ = ntuple->Branch("diEle_DZ",&diEle_DZ,"diEle_DZ/D");
  TBranch *b_Ele1PT = ntuple->Branch("Ele1PT",&Ele1PT,"Ele1PT/D");
  TBranch *b_Ele1Eta = ntuple->Branch("Ele1Eta",&Ele1Eta,"Ele1Eta/D");
  TBranch *b_Ele1EtaSC = ntuple->Branch("Ele1EtaSC",&Ele1EtaSC,"Ele1EtaSC/D");
  TBranch *b_Ele1Phi = ntuple->Branch("Ele1Phi",&Ele1Phi,"Ele1Phi/D");
  TBranch *b_Ele2PT = ntuple->Branch("Ele2PT",&Ele2PT,"Ele2PT/D");
  TBranch *b_Ele2Eta = ntuple->Branch("Ele2Eta",&Ele2Eta,"Ele2Eta/D");
  TBranch *b_Ele2EtaSC = ntuple->Branch("Ele2EtaSC",&Ele2EtaSC,"Ele2EtaSC/D");
  TBranch *b_Ele2Phi = ntuple->Branch("Ele2Phi",&Ele2Phi,"Ele2Phi/D");
  TBranch *b_ReferenceTriggerPass= ntuple->Branch("ReferenceTriggerPass",&ReferenceTriggerPass,"ReferenceTriggerPass/I");
  TBranch *b_isTriggerSoupPass= ntuple->Branch("isTriggerSoupPass",&isTriggerSoupPass,"isTriggerSoupPass/I");
  TBranch *b_pass_Ele23_12= ntuple->Branch("pass_Ele23_12",&pass_Ele23_12,"pass_Ele23_12/I");
  TBranch *b_pass_Ele23_12_DZ= ntuple->Branch("pass_Ele23_12_DZ",&pass_Ele23_12_DZ,"pass_Ele23_12_DZ/I");
  TBranch *b_pass_Ele17_12= ntuple->Branch("pass_Ele17_12",&pass_Ele17_12,"pass_Ele17_12/I");
  TBranch *b_pass_Ele23= ntuple->Branch("pass_Ele23",&pass_Ele23,"pass_Ele23/I");
  TBranch *b_pass_Ele27Loose= ntuple->Branch("pass_Ele27Loose",&pass_Ele27Loose,"pass_Ele27Loose/I");
  TBranch *b_pass_Ele27LooseEtar= ntuple->Branch("pass_Ele27LooseEtar",&pass_Ele27LooseEtar,"pass_Ele27LooseEtar/I");
  TBranch *b_pass_Ele27TightEtar= ntuple->Branch("pass_Ele27TightEtar",&pass_Ele27TightEtar,"pass_Ele27TightEtar/I");
  TBranch *b_passdRleg1 = ntuple->Branch("passdRleg1",&passdRleg1,"passdRleg1/I");
  TBranch *b_passdRleg2 = ntuple->Branch("passdRleg2",&passdRleg2,"passdRleg2/I");
  TBranch *b_passdRleg3 = ntuple->Branch("passdRleg3",&passdRleg3,"passdRleg3/I");
  TBranch *b_passdRleg1N = ntuple->Branch("passdRleg1N",&passdRleg1N,"passdRleg1N/I");
  TBranch *b_passdRleg2N = ntuple->Branch("passdRleg2N",&passdRleg2N,"passdRleg2N/I");
  TBranch *b_passdRleg3N = ntuple->Branch("passdRleg3N",&passdRleg3N,"passdRleg3N/I");
  TBranch *b_Zmass = ntuple->Branch("Zmass",&Zmass,"Zmass/D");
  TBranch *b_nVtx = ntuple->Branch("nVtx",&nVtx,"nVtx/I");
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
  std::string TriggerList[3] = {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_Ele27_eta2p1_WPLoose_Gsf_v*"};
  
  std::string ReferenceTrigger = "HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v*"; 
 
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
  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE


  //##############################################
  //########           EVENT LOOP         ########
  //##############################################
  //loop on all the events
  printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
  printf("Scanning the ntuple :");
  int treeStep(totalEntries/50);


cout << "totalEntries = " << totalEntries << endl; 
  for( size_t iev=0; iev< totalEntries; iev++){

//if(iev%10000==0) cout <<"Number of events processed =   " << iev << endl;

       ev.to(iev); //load the event content from the EDM file
       edm::EventBase const & myEvent = ev;

//if(ev%10==0) cout << ev << endl;

double rho = 0.;

Ele1PT	= -999.;
Ele1Eta	= -999.;
Ele1EtaSC = -999.;
Ele1Phi = -999.; 
Ele2PT	= -999.;
Ele2Eta	= -999.;
Ele2EtaSC = -999.;
Ele2Phi = -999.;
ReferenceTriggerPass = 0;
isTriggerSoupPass = 0; 
pass_Ele23_12 =0;
pass_Ele17_12 =0;
pass_Ele23 =0;
nVtx    = -999;
Zmass	= -999.;
passdRleg1 = 0;

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

   //check the L1 candidateds
   //    
       fwlite::Handle<l1extra::L1EmParticleCollection> theL1extraParticles;
       theL1extraParticles.getByLabel(ev,"l1extraParticles","NonIsolated","RECO");

       fwlite::Handle<l1extra::L1EmParticleCollection> theL1extraParticlesIsolated;
       theL1extraParticlesIsolated.getByLabel(ev,"l1extraParticles","Isolated","RECO");

       fwlite::Handle< pat::TriggerObjectStandAloneCollection > triggerObjects;
       triggerObjects.getByLabel(ev, "selectedPatTrigger");

       fwlite::Handle< double > rhoH;
       rhoH.getByLabel(ev,"fixedGridRhoFastjetAllCalo");
       rho = *rhoH;

 //     fwlite::Handle< pat::TriggerObjectStandAloneCollection > triggerObjectsNew;
  //    triggerObjectsNew.getByLabel(ev, "selectedPatTriggerUpdated");


       std::vector<pat::TriggerObjectStandAlone> TagtriggerObj;
       std::vector<pat::TriggerObjectStandAlone> ProbetriggerObj;
       std::vector<pat::TriggerObjectStandAlone> passLeg2triggerObj;
       std::vector<pat::TriggerObjectStandAlone> passLeg1Obj;
       std::vector<pat::TriggerObjectStandAlone> passLeg2Obj;
       std::vector<pat::TriggerObjectStandAlone> passLeg3Obj;
       std::vector<pat::TriggerObjectStandAlone> passLeg4Obj;
       std::vector<pat::TriggerObjectStandAlone> passLeg5Obj;
       std::vector<pat::TriggerObjectStandAlone> passLeg6Obj;

       // DERIVE WEIGHTS TO APPLY TO SAMPLE

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

  //     cout << "Number of leptons  = " << selLeptons.size() << endl;

	 // opposite charge condition
	 
 //      if (! (selLeptons[0].charge()*selLeptons[1].charge() < 0) ) continue;

//       for(int s = 0 ; s < selLeptons.size(); s++) {
//cout << "selLeptons[s].pdgId() = " << selLeptons[s].pdgId() << endl;
//       if(abs(selLeptons[s].pdgId()) == 11) {index.push_back(s);}

//}

//cout << "Number of electrons in an event = " << index.size() << endl;

       if(abs(selLeptons[0].pdgId()) == 13 && abs(selLeptons[1].pdgId()) == 13) passOnlyMuon = true;
       if(abs(selLeptons[0].pdgId()) == 11 && abs(selLeptons[1].pdgId()) == 11) passOnlyEle = true;

//if(index.size() >= 2) passOnlyEle = true;

     
       //std::cout << "=> Passed the selection of the events containing only muon or ele" << std::endl; 
       bool passKinEle1(false),passIdEle1(false),passIsoEle1(false);
       bool passKinEle2(false),passIdEle2(false),passIsoEle2(false);
       bool Probe1ISO(false), Probe2ISO(false), TagISO(false), passTightIdSTDEle1(false),passTightIdSTDEle2(false);
       bool ZPick(false);

      //////////////////////////
      ///                    ///
      /// ELECTRON EFFICENCY ///
      ///                    ///
      ////////////////////////// 
      //
 //     cout << "passOnlyEle = " << passOnlyEle << endl;
         if(passOnlyEle){

   //      int first  = rand()%2;
   //      int second = (first+1)%2;
         int first  = 0;
         int second = 1;
   
//for(int t = 0; t < index.size(); t++){
//cout << "Electrons PT = " << selLeptons[index[t]].el.pt() << endl;
//} 
      
         double ptf  = selLeptons[first].el.pt();
         double etaf = selLeptons[first].el.eta();
         double etafSC = selLeptons[first].el.superCluster()->eta();
         double phif = selLeptons[first].el.phi();
         double pts  = selLeptons[second].el.pt();
         double etas = selLeptons[second].el.eta();
         double etasSC = selLeptons[second].el.superCluster()->eta();         
         double phis = selLeptons[second].el.phi();

         double ecalPFIsof = selLeptons[first].el.ecalPFClusterIso();
         double hcalPFIsof = selLeptons[first].el.hcalPFClusterIso();
         double trackIsof = selLeptons[first].el.trackIso();
         double ecalPFIsos = selLeptons[second].el.ecalPFClusterIso();
         double hcalPFIsos = selLeptons[second].el.hcalPFClusterIso();
         double trackIsos = selLeptons[second].el.trackIso();
         double detaseedf = selLeptons[first].el.deltaEtaSeedClusterTrackAtCalo();
         double detaseeds = selLeptons[second].el.deltaEtaSeedClusterTrackAtCalo();
         double misshit = selLeptons[second].el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
         double misshitf = selLeptons[first].el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
         double dphiin = selLeptons[second].el.deltaPhiSuperClusterTrackAtVtx();
         double dphiinf = selLeptons[first].el.deltaPhiSuperClusterTrackAtVtx();
         double chi2 = selLeptons[second].el.gsfTrack()->normalizedChi2();
         double chi2f = selLeptons[first].el.gsfTrack()->normalizedChi2();

                if(fabs(etasSC) <= 1.479) {
               Probe2ISO = ((ecalPFIsos-rho*0.165)/pts) < 0.160 && ((hcalPFIsos-rho*0.060)/pts) < 0.120 && (trackIsos/pts) < 0.08 && abs(detaseeds) < 0.004 && misshit < 1 && abs(dphiin) < 0.020;
                         }
               else if(fabs(etasSC) > 1.479 && fabs(etasSC) < 2.5) {
               Probe2ISO = ((ecalPFIsos-rho*0.132)/pts) < 0.120 && ((hcalPFIsos-rho*0.131)/pts) < 0.120 && (trackIsos/pts) < 0.08 && misshit < 1 && abs(chi2) < 3;
                         }

                if(fabs(etafSC) <= 1.479) {
               Probe1ISO = ((ecalPFIsof-rho*0.165)/ptf) < 0.160 && ((hcalPFIsof-rho*0.060)/ptf) < 0.120 && (trackIsof/ptf) < 0.08 && abs(detaseedf) < 0.004 && misshitf < 1 && abs(dphiinf) < 0.020;
                         }
               else if(fabs(etafSC) > 1.479 && fabs(etafSC) < 2.5) {
               Probe1ISO = ((ecalPFIsof-rho*0.132)/ptf) < 0.120 && ((hcalPFIsof-rho*0.131)/ptf) < 0.120 && (trackIsof/ptf) < 0.08 && misshitf < 1 && abs(chi2f) < 3;
                         }


         passKinEle1 = ((abs(etafSC) >= 0 && abs(etafSC) <= 2.5));
         passTightIdSTDEle1       = patUtils::passId(electronVidTightId, myEvent, selLeptons[first].el);     
         passKinEle2 = ((abs(etasSC) >= 0 && abs(etasSC) <= 2.5));       
         passTightIdSTDEle2       = patUtils::passId(electronVidTightId, myEvent, selLeptons[second].el);
	 bool GoodEle1 = false;
	 bool GoodEle2 = false;

      if(passKinEle1 && passTightIdSTDEle1 && Probe1ISO) GoodEle1 = true;
      if(passKinEle2 && passTightIdSTDEle2 && Probe2ISO) GoodEle2 = true;
	 
	 if((GoodEle1 && GoodEle2) != true) continue;

//cout << " GOOD Leptons Found " << endl;

//`cout << "Number of verticies = " << vtx.size() << endl;

	 TLorentzVector lep1(selLeptons[first].px(),selLeptons[first].py(),selLeptons[first].pz(),selLeptons[first].energy());
	 TLorentzVector lep2(selLeptons[second].px(),selLeptons[second].py(),selLeptons[second].pz(),selLeptons[second].energy());  
         double mass = (lep1+lep2).M();

         ZPick = mass > 60 && mass < 120;
	 if(!ZPick) continue;	 
	 Ele1PT = ptf;
	 Ele1Eta = etaf;
	 Ele1EtaSC = etafSC;
	 Ele1Phi = phif;
	 Ele2PT = pts;
	 Ele2Eta = etas;
	 Ele2EtaSC = etasSC;
	 Ele2Phi = phis;
	 Zmass = mass;
         nVtx = vtx.size();
	 bool refTrig(false);
	 bool soupOfTrig(false);	 
	 bool Ele23(false);
	 bool Ele23_12(false);
         bool Ele23_12_DZ(false);
	 bool Ele17_12(false);
         bool Ele27_Loose(false);
         bool Ele27_LooseEr(false);
         bool Ele27_TightEr(false);

//	 refTrig   = utils::passTriggerPatterns(tr, ReferenceTrigger);
 //        Ele23     = utils::passTriggerPatterns(tr, TriggerList[0]);
         Ele23_12  = utils::passTriggerPatterns(tr, TriggerList[0]);
         Ele23_12_DZ  = utils::passTriggerPatterns(tr, TriggerList[1]);
//         Ele17_12  = utils::passTriggerPatterns(tr, TriggerList[2]);
//	 soupOfTrig= (Ele23 || Ele17_12);
//         Ele27_Loose = utils::passTriggerPatterns(tr, TriggerList[3]);
//         Ele27_TightEr = utils::passTriggerPatterns(tr, TriggerList[4]);
//         Ele27_LooseEr = utils::passTriggerPatterns(tr, TriggerList[5]);

//	 if(refTrig)	 { ReferenceTriggerPass = 1;}
//	 if(soupOfTrig)  { isTriggerSoupPass  = 1;}
//	 if(Ele23) 	 { pass_Ele23  = 1;}
	 if(Ele23_12) 	 { pass_Ele23_12  = 1;}
         if(Ele23_12_DZ)    { pass_Ele23_12_DZ  = 1;}
//	 if(Ele17_12)  	 { pass_Ele17_12  = 1;}
 //        if(Ele27_Loose) { pass_Ele27Loose = 1;}
  //       if(Ele27_LooseEr) { pass_Ele27LooseEtar = 1;}
  //       if(Ele27_TightEr) { pass_Ele27TightEtar = 1;}




//================================================ Trigger Objects wala Natak ===============================
/*8
bool passProbeLeg1(false), passProbeLeg2(false), passProbeLeg3(false), passProbeLeg4(false),passProbeLeg5(false),passProbeLeg6(false);
bool passleg1dRCut(false),passleg2dRCut(false),passleg3dRCut(false),passleg4dRCut(false),passleg5dRCut(false),passleg6dRCut(false);

bool passProbeLeg1New(false), passProbeLeg2New(false), passProbeLeg3New(false), passProbeLeg4New(false),passProbeLeg5New(false),passProbeLeg6New(false);
bool passleg1dRCutNew(false),passleg2dRCutNew(false),passleg3dRCutNew(false),passleg4dRCutNew(false),passleg5dRCutNew(false),passleg6dRCutNew(false);


         for ( pat::TriggerObjectStandAlone obj: *triggerObjects ) {
             obj.unpackPathNames(names);

	std::string leg1Filter("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter");
	std::string leg2Filter("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter");
	std::string leg3Filter("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter");
	std::string leg4Filter("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter");
	std::string leg5Filter("hltEle23WPLooseGsfTrackIsoFilter");
	std::string leg6Filter("hltEle27noerWPLooseGsfTrackIsoFilter");

  	passProbeLeg1     = passFilter(obj,leg1Filter);
  	passProbeLeg2     = passFilter(obj,leg2Filter);
  	passProbeLeg3     = passFilter(obj,leg3Filter);
  	passProbeLeg4     = passFilter(obj,leg4Filter);
  	passProbeLeg5     = passFilter(obj,leg5Filter);
  	passProbeLeg6     = passFilter(obj,leg6Filter);

        if(passProbeLeg1) passLeg1Obj.push_back(obj);
        if(passProbeLeg2) passLeg2Obj.push_back(obj);
        if(passProbeLeg3) passLeg3Obj.push_back(obj);
        if(passProbeLeg4) passLeg4Obj.push_back(obj);
        if(passProbeLeg5) passLeg5Obj.push_back(obj);
        if(passProbeLeg6) passLeg6Obj.push_back(obj);

	} // Trigger Objects loop 

       double  dRleg1_leadele = DeltaRtrig(selLeptons[first].el,passLeg1Obj);
       double  dRleg1_subleadele = DeltaRtrig(selLeptons[second].el,passLeg1Obj);
       double  dRleg2_leadele = DeltaRtrig(selLeptons[first].el,passLeg2Obj);
       double  dRleg2_subleadele = DeltaRtrig(selLeptons[second].el,passLeg2Obj);


if(dRleg1_leadele < 0.1 || dRleg1_subleadele < 0.1) passleg1dRCut = true;
passdRleg1 = passleg1dRCut;

if(dRleg1_leadele < 0.1 || dRleg1_subleadele < 0.1) passleg2dRCut = true;
passdRleg2 = passleg2dRCut;

if((dRleg1_leadele < 0.1 && dRleg2_subleadele < 0.1) || (dRleg1_subleadele < 0.1 && dRleg2_leadele < 0.1)) passleg3dRCut = true;
passdRleg3 = passleg3dRCut;
*/
/*
//=================================== New Trigger Objects ====================================

         for ( pat::TriggerObjectStandAlone obj1: *triggerObjectsNew ) {
             obj1.unpackPathNames(names_New);

        std::string leg1FilterN("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter");
        std::string leg2FilterN("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter");
        std::string leg3FilterN("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter");
        std::string leg4FilterN("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter");
        std::string leg5FilterN("hltEle23WPLooseGsfTrackIsoFilter");
        std::string leg6FilterN("hltEle27noerWPLooseGsfTrackIsoFilter");

        passProbeLeg1New     = passFilter(obj1,leg1FilterN);
        passProbeLeg2New     = passFilter(obj1,leg2FilterN);
        passProbeLeg3New     = passFilter(obj1,leg3FilterN);
        passProbeLeg4New     = passFilter(obj1,leg4FilterN);
        passProbeLeg5New     = passFilter(obj1,leg5FilterN);
        passProbeLeg6New     = passFilter(obj1,leg6FilterN);

        if(passProbeLeg1New) passLeg1ObjNew.push_back(obj1);
        if(passProbeLeg2New) passLeg2ObjNew.push_back(obj1);
        if(passProbeLeg3New) passLeg3ObjNew.push_back(obj1);
        if(passProbeLeg4New) passLeg4ObjNew.push_back(obj1);
        if(passProbeLeg5New) passLeg5ObjNew.push_back(obj1);
        if(passProbeLeg6New) passLeg6ObjNew.push_back(obj1);
        }

       double  dRleg1_leadeleN = DeltaRtrig(selLeptons[first].el,passLeg1ObjNew);
       double  dRleg1_subleadeleN = DeltaRtrig(selLeptons[second].el,passLeg1ObjNew);
       double  dRleg2_leadeleN = DeltaRtrig(selLeptons[first].el,passLeg2ObjNew);
       double  dRleg2_subleadeleN = DeltaRtrig(selLeptons[second].el,passLeg2ObjNew);


if(dRleg1_leadeleN < 0.1 || dRleg1_subleadeleN < 0.1) passleg1dRCutNew = true;
passdRleg1N = passleg1dRCutNew;

if(dRleg1_leadeleN < 0.1 || dRleg1_subleadeleN < 0.1) passleg2dRCutNew = true;
passdRleg2N = passleg2dRCutNew;

if((dRleg1_leadeleN < 0.1 && dRleg2_subleadeleN < 0.1) || (dRleg1_subleadeleN < 0.1 && dRleg2_leadeleN < 0.1)) passleg3dRCutNew = true;
passdRleg3N = passleg3dRCutNew;
*/

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
