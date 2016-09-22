#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"

using namespace std;

//
GammaWeightsHandler::GammaWeightsHandler(const edm::ParameterSet &runProcess,TString ewkSupWgt,bool forceAllToData)
{
  //cfg
  bool isMC = runProcess.getParameter<bool>("isMC");
  if(forceAllToData) isMC=false;
  std::vector<std::string> gammaPtWeightsFiles =  runProcess.getParameter<std::vector<std::string> >("weightsFile");  
  if(gammaPtWeightsFiles.size()==0) return;
  std::vector<TString> wgtNames; 
  wgtNames.push_back("qt");
  if(ewkSupWgt!="") wgtNames.push_back(ewkSupWgt);
  TString wgtType( isMC ? "mcfitwgts" : "datafitwgtsrebin");
  TString massType( isMC ? "mczmass" : "zmass");
    
  //categories to consider, add more if needed but keep these ones 
  TString cats[]   =  {"eq0jets","eq1jets","eq2jets","geq3jets","vbf","geq1jets","novbf","mjjq100","mjjq092","mjjq083","mjjq066","mjjq049","mjjq033","mjjq016"};
  dilCats_.push_back("ee"); dilCats_.push_back("mumu"); dilCats_.push_back("ll");
  
  //retrieve from file
  TString gammaPtWeightsFile(gammaPtWeightsFiles[0].c_str());
  gSystem->ExpandPathName(gammaPtWeightsFile);
  TFile *fwgt=TFile::Open(gammaPtWeightsFile);
  if(fwgt)
    {
      cout << "[GammaWeightsHandler] retrieving weights from: " << gammaPtWeightsFile << endl;
      
      std::map<TString, TGraph*> iWgtsH;
      for(size_t ic=0; ic<sizeof(cats)/sizeof(TString); ic++)
	{
	  for(size_t id=0; id<dilCats_.size(); id++)
	    {
	      TString key = dilCats_[id] + cats[ic];
	      std::vector<TGraph *> iwgts;
	      for(size_t iw=0; iw<wgtNames.size(); iw++)
		{ 
		  //weights
		  TString hname= key + "_" + wgtNames[iw] + "_" + wgtType;
		  cout << hname << endl;
		  TGraph *h = (TGraph *) fwgt->Get(hname);
		  if(h) iwgts.push_back(h);
		  
		  //mass shape
		  if(iw>0) continue;
		  hname= key+"_"+massType;
		  TH1 *massh = (TH1 *) fwgt->Get(hname);
		  if(massh!=0) { massh->SetDirectory(0); zmassH_[key]= massh; }
		}
	      
	      if(iwgts.size()) wgtsH_[key] = iwgts;
	    }
	}
      fwgt->Close();
    }
  
  if(wgtsH_.size()==0) return; 
  std::cout << "[GammaWeightsHandler] gamma spectrum will be reweighted using distributions found in "  
	    << gammaPtWeightsFiles.size() 
	    << " files" 
	    << std::endl;
}

//
LorentzVector GammaWeightsHandler::getMassiveP4(LorentzVector &gamma,TString evCategoryLabel)
{
  //generate a mass from the line shape (0 if not available)
  float mass(0);
  if(zmassH_.find(evCategoryLabel)!=zmassH_.end())
    {
      if(zmassH_[evCategoryLabel]->Integral())
	while(fabs(mass-91)>15) 
	  mass = zmassH_[evCategoryLabel]->GetRandom();
    }
  return LorentzVector(gamma.px(),gamma.py(),gamma.pz(),sqrt(pow(mass,2)+pow(gamma.P(),2)));
}

float GammaWeightsHandler::getWeightFor(std::vector<Float_t> &vars, TString evCategoryLabel)
{
  //get the weight (1.0 if not available)
  float weight(1.0);
  if(wgtsH_.find(evCategoryLabel) != wgtsH_.end()){
      std::vector<TGraph *> &availableWeights=wgtsH_[evCategoryLabel];
      for(size_t ivar=0; ivar<availableWeights.size(); ivar++)
	{
	  TGraph *h = availableWeights[ivar];
	  if(vars.size()>ivar) weight*=h->Eval(vars[ivar]);
	  if(weight<0) weight=0;
	}
  }
    
  return weight;
}


//
GammaWeightsHandler::~GammaWeightsHandler()
{
}
