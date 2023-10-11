#define NanoSkim_cxx

#include "NanoSkim.h"
#include <TH2.h>
#include <TStyle.h>


void NanoSkim::Begin(TTree * /*tree*/)
{
   TString option = GetOption();
   nEvtTotal = 0;
   nEvtRan = 0;
   nEvtSkim = 0;
   _HstFile = new TFile(_HstFileName,"recreate");
   BookHistograms();


}

void NanoSkim::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();
}
void NanoSkim::SlaveTerminate()
{
  _HstFile->Write();
  _HstFile->Close();

  cout<<"Writing skim file .... ";
  outT->Write();
  _SkimFile->Close();
  cout<<"  Done. "<<endl;
  

  cout<<"Total events         = "<<nEvtTotal<<endl;
  cout<<"Total events ran     = "<<nEvtRan<<endl;
  cout<<"Total events skimmed = "<<nEvtSkim<<endl;

  ofstream fout(_SumFileName);
  fout<<"Total events         = "<<nEvtTotal<<endl;
  fout<<"Total events ran     = "<<nEvtRan<<endl;
  fout<<"Total events skimmed = "<<nEvtSkim<<endl;

}

void NanoSkim::Terminate()
{
  //   cout<<"Inside Terminate()"<<endl;
}

Bool_t NanoSkim::Process(Long64_t entry)
{

  // Must read in whole event
  // DO NOT EDIT NEXT LINE
  int readevent = ReadLimited(0,entry);
  if(readevent==0){ cout<<"Did not read in any branches.. quitting."<<endl; return kTRUE;}
 
  if(_verbosity==0 && nEvtTotal%10000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;      
  else if(_verbosity>0 && nEvtTotal%50000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;

  nEvtTotal++;
  h.nevt->Fill(0);

  GoodEvt2018 = (_year==2018 ? Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && (_data ? Flag_eeBadScFilter : 1) : 1);
  GoodEvt2017 = (_year==2017 ? Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && (_data ? Flag_eeBadScFilter : 1) : 1);
  GoodEvt2016 = (_year==2016 ? Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && (_data ? Flag_eeBadScFilter : 1) : 1);

  GoodEvt = GoodEvt2018 && GoodEvt2017 && GoodEvt2016;

  if(GoodEvt){
    nEvtRan++;
    h.nevt->Fill(1);
  
    // LEP SELECTION BLOCK BEGINS                                               

    //clear array from previous event
    llep.clear(); taus.clear();

    //MUONS
    for(unsigned int i=0; i< (nMuon); i++){
      Lepton temp; temp.v.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],0.105); 
      temp.id = -13*Muon_charge[i]; temp.ind = i;  temp.charge = Muon_charge[i];
      bool passCuts = temp.v.Pt()>5 && fabs(temp.v.Eta())<2.4 && Muon_mediumId[i];
      passCuts = passCuts && Muon_pfRelIso04_all[i]<0.15;
      passCuts = passCuts && fabs(Muon_dxy[i])<0.05 && fabs(Muon_dz[i])<0.1;
      if(passCuts)
	llep.push_back(temp);
    }    

    //ELECTRONS
    for(unsigned int i=0; i< (nElectron); i++){
      Lepton temp; temp.v.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],0.000511); 
      temp.id = -11*Electron_charge[i]; temp.ind = i; temp.charge = Electron_charge[i];
      bool isprompt = false;
      if(fabs(temp.v.Eta())<=1.479)
	if(fabs(Electron_dxy[i])<0.05 && fabs(Electron_dz[i])<0.1)
	  isprompt = true;
      
      if(fabs(temp.v.Eta())>1.479)
	if(fabs(Electron_dxy[i])<0.1 && fabs(Electron_dz[i])<0.2)
	  isprompt = true;

      bool passCuts = temp.v.Pt()>5 && fabs(temp.v.Eta())<2.5 && Electron_cutBased[i]>2;
      passCuts = passCuts && isprompt;
      if(passCuts)
	llep.push_back(temp);
    }
    
    Sort(0); // Sort the light leptons.

    
    //TAUS
    for(unsigned int i=0; i< (nTau); i++){
      if(Tau_idDecayMode[i]&&(Tau_decayMode[i]<3||Tau_decayMode[i]>9)){
	//Tau energy scale correction
	float tlv_corr = 1.;
	if(_year==2016){
	  if(Tau_decayMode[i]==0) tlv_corr = 0.994;
	  if(Tau_decayMode[i]==1) tlv_corr = 0.995;
	  if(Tau_decayMode[i]>9)  tlv_corr = 1;
	}
	if(_year==2017){
	  if(Tau_decayMode[i]==0) tlv_corr = 1.007;
	  if(Tau_decayMode[i]==1) tlv_corr = 0.998;
	  if(Tau_decayMode[i]==10) tlv_corr = 1.001;
	  if(Tau_decayMode[i]==11) tlv_corr = 0.999;
	}
	if(_year==2018){
	  if(Tau_decayMode[i]==0) tlv_corr = 0.987;
	  if(Tau_decayMode[i]==1) tlv_corr = 0.995;
	  if(Tau_decayMode[i]==10) tlv_corr = 0.998;
	  if(Tau_decayMode[i]==11) tlv_corr = 1;
	}
	Lepton temp; temp.v.SetPtEtaPhiM(Tau_pt[i],Tau_eta[i],Tau_phi[i],1.77);
	temp.v *= tlv_corr; //energy correction
	temp.id = -15*Tau_charge[i]; temp.ind = i; temp.charge = Tau_charge[i];
	temp.lepcleaning = TaulepCleaning(temp.v);


	bool passCuts = temp.v.Pt()>20 && fabs(temp.v.Eta()<2.3);
	bool useDeepTau = true;
	if(!useDeepTau){
	  // MVA Tau ID.
	  passCuts = passCuts && Tau_idAntiEle[i]>1 && Tau_idAntiMu[i]>0; //looseWP both
	  passCuts = passCuts && fabs(Tau_dz[i])<0.2;
	  passCuts = passCuts && Tau_idMVAoldDM2017v2[i]>1 && temp.lepcleaning;
	}
	else if(useDeepTau){
	  // Deep Tau ID.
	  passCuts = passCuts &&Tau_idDecayModeNewDMs[i]==1 && (Tau_decayMode[i]<3 || Tau_decayMode[i]>9);
	  passCuts = passCuts && fabs(Tau_dz[i])<0.2;
	  passCuts = passCuts && Tau_idDeepTau2017v2p1VSe[i] >= 15
	    && Tau_idDeepTau2017v2p1VSmu[i] >= 3; //loose WP
	  passCuts = passCuts && Tau_idDeepTau2017v2p1VSjet[i] >= 31; //medium WP
	}
	
	if(passCuts)
	  taus.push_back(temp);
      }
    }
    
    Sort(1); //Sort the Taus
    
    //MET
    if(_year == 2017){
      metpt = METFixEE2017_pt;
      metphi = METFixEE2017_phi;
    }
    else{
      metpt = MET_pt;
      metphi = MET_phi;
    }
    
    
    // LEP SELECTION DONE

    
    bool keepThisEvent=false;    

    bool is_3L_event = false;
    bool is_2L1T_event = false;
    bool is_1L2T_event = false;

    if((int)llep.size()>2){ //3L
      bool passTrigger = false;
      if(abs(llep.at(0).id)==11)
	passTrigger = llep.at(0).v.Pt()>30;
      else if(abs(llep.at(0).id)==13)
	passTrigger = llep.at(0).v.Pt()>26;

      if(passTrigger)
	is_3L_event=true;
    }
    if((int)llep.size()>1){ //2L
      bool passTrigger = false;
      if(abs(llep.at(0).id)==11)
	passTrigger = llep.at(0).v.Pt()>30;
      else if(abs(llep.at(0).id)==13)
	passTrigger = llep.at(0).v.Pt()>26;

      if(passTrigger && (int)taus.size()>0) //2L1T
	is_2L1T_event=true;
    }
    if((int)llep.size()>0){ //1L2T
      bool passTrigger = false;
      if(abs(llep.at(0).id)==11)
	passTrigger = llep.at(0).v.Pt()>30;
      else if(abs(llep.at(0).id)==13)
	passTrigger = llep.at(0).v.Pt()>26;

      if(passTrigger && (int)taus.size()>1) //2L1T
	is_1L2T_event=true;
    }

    // ---------------------------------
    //if(is_3L_event) keepThisEvent=true;
    if(is_2L1T_event) keepThisEvent=true;
    //if(is_1L2T_event) keepThisEvent=true;
    
    if(keepThisEvent){
      nEvtSkim++;
      h.nevt->Fill(2);
      outT->Fill();
    }
    
  }
  return kTRUE;
}
void NanoSkim::Sort(int opt)
{
  if(opt==0){
    for(int i=0; i<(int)llep.size()-1; i++){
      for(int j=i+1; j<(int)llep.size(); j++){
	if( llep[i].v.Pt() < llep[j].v.Pt() ) swap(llep.at(i),llep.at(j)); }}
  }
  if(opt==1){
    for(int i=0; i<(int)taus.size()-1; i++){
      for(int j=i+1; j<(int)taus.size(); j++){
	if( taus[i].v.Pt() < taus[j].v.Pt() ) swap(taus.at(i),taus.at(j)); }}
  }

}
float NanoSkim::delta_phi(float phi1, float phi2)
{
  phi1 = TVector2::Phi_0_2pi(phi1);
  phi2 = TVector2::Phi_0_2pi(phi2);
  float dphi = fabs(phi1 - phi2);
  if(dphi>TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
  return dphi;
}
bool NanoSkim::TaulepCleaning(TLorentzVector t)
{
  bool result1=false;
  if((int)llep.size()==0)
    result1=true;
  if((int)llep.size()>0){
    for(int i=0; i<(int)llep.size(); i++){
      if(t.DeltaR(llep[i].v)>0.4)
        result1=true;
      else{
        result1=false;
        break;
      }
    }
  }
  return result1;
} 

void NanoSkim::BookHistograms()
{
  h.nevt = new TH1F("nEvents","0-Total events, 1-Total events ran, 2-Total events skimmed",5,-1,4);
}

