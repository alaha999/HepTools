#include <iostream>
#include <fstream>
#include <iomanip> 
#include <vector>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include "TF1.h"
#include <Stntuple/loop/TStnAna.hh>
#include <Stntuple/obj/TStnDBManager.hh>
#include <Stntuple/obj/TStnTriggerTable.hh>
#include <Stntuple/obj/TStnGoodRunList.hh>
#include <Stntuple/obj/TStnRunSummary.hh>
#include <Stntuple/data/TStnBeamPos.hh>
#include "Stntuple/alg/TStntuple.hh"


#include "JetUser/JetEnergyCorrections.hh"
#include "TSeel.h"

#define CMX_CMUP_LUM_RATIO 0.9263
// 232.027/250.484 = 0.9263
#define FINAL_JET_CUT 15.

TSeel::TSeel(const char* name, const char* title):
        TStnModule(name,title)
{
    fVerbose = 0;
    _jetAlgo = "JetClu";
    _coneSize = 4;
    _JES_version = 5;
    _JES_level = 0;
    _JET_ET_CUT = 5.;
    _LEP_PT_CUT3[0] = 10.;
    _LEP_PT_CUT3[1] = 5.;
    _LEP_PT_CUT3[2] = 4.;
    _LEP_PT_CUT2[0] = 10.;
    _LEP_PT_CUT2[1] = 10.;
    _USE_TRACK = 1;
    L3_name[0] = "DIELECTRON_CENTRAL_4";
    L3_name[1] = "DIMUON_CMUCMU4";
    L3_name[2] = "DIMUON_CMUPCMUP4";
    L3_name[3] = "DIMUON_CMU4_CMX4";
    L3_name[4] = "DIMUON_CMUP4_CMX4";
    L3_name[5] = "CEM4_CMU4";
    L3_name[6] = "CEM4_CMUP4";
    L3_name[7] = "CEM4_CMX4";
    L3_name[8] = "CEM4_PEM8";
    L3_name[9] = "CEM8_PEM8";
    L3_name[10] = "CMU4_PEM8";
    L3_name[11] = "CMUP4_PEM8";
    L3_name[12] = "CMX4_PEM8";
    L3_name[13] = "ELECTRON_4_LOOSE";
    L3_name[14] = "MUON_CMUP4";


  //fake rate as function of Et for each category

   //central tight electrons
  _centralelectightfakerate = new TF1("centralelectight","[0]+expo(1)",4,50);
  _centralelectightfakerate->SetParameter(0,0.00013);
  _centralelectightfakerate->SetParameter(1,-7.94);
  _centralelectightfakerate->SetParameter(2,-0.194);

  //central loose electrons
  _centralelecloosefakerate = new TF1("centralelecloose","expo",4,50);
  _centralelecloosefakerate->SetParameter(0,-8.0);
  _centralelecloosefakerate->SetParameter(1,0.022);

  //CMUP muons
  _CMUPfakerate = new TF1("CMUP","pol1",4,50);
  _CMUPfakerate->SetParameter(0,0.00086);
  _CMUPfakerate->SetParameter(1,0.00017);
 
  //CMX muons
  _CMXfakerate = new TF1("CMX","pol1",4,50);
  _CMXfakerate->SetParameter(0,0.00082);
  _CMXfakerate->SetParameter(1,0.0002);

  //CMIO muons
  _CMIOfakerate = new TF1("CMIO","pol1",10,50);
  _CMIOfakerate->SetParameter(0,-0.014);
  _CMIOfakerate->SetParameter(1,0.0037);

  _smearMean = 0.89;
  _smearRMS  = 0.06;

}


TSeel::~TSeel()
{}

int TSeel::BeginJob()
{

         // register the data block
    RegisterDataBlock("HeaderBlock","TStnHeaderBlock",&fHeaderBlock);
    RegisterDataBlock("CalDataBlock","TCalDataBlock",&fCalDataBlock);
    RegisterDataBlock("DcasDataBlock","TDcasDataBlock",&fDcasDataBlock);
    RegisterDataBlock("TrackBlock","TStnTrackBlock",&fTrackBlock);
    RegisterDataBlock("TriggerBlock","TStnTriggerBlock",&fTriggerBlock);
    RegisterDataBlock("PhotonBlock","TStnPhotonBlock",&fPhotonBlock);
    RegisterDataBlock("ElectronBlock","TStnElectronBlock",&fElectronBlock);
    RegisterDataBlock("CesDataBlock","TCesDataBlock",&fCesDataBlock);
    RegisterDataBlock("XftBlock","TXftBlock",&fXftBlock);
    RegisterDataBlock("MuonBlock","TStnMuonBlock",&fMuonBlock);
    RegisterDataBlock("JetBlock","TStnJetBlock",&fJetBlock);
    RegisterDataBlock("MetBlock","TStnMetBlock",&fMetBlock); 
    RegisterDataBlock("ZVertexBlock","TStnVertexBlock",&fVertexBlock);
    RegisterDataBlock("TStnCosmicBlock","TStnCosmicBlock",&fCosmicBlock);
    RegisterDataBlock("L3SummaryBlock","TL3SummaryBlock",&fL3SummaryBlock);
    RegisterDataBlock("GenpBlock","TGenpBlock",&fGenpBlock);
 
    if(_coneSize == 7){
        if(_jetAlgo == "KtClus"){
            RegisterDataBlock("PROD@KtClusModule-vz-cone0.7","TStnJetBlock",&fJetBlock);
        }
        else if(_jetAlgo == "MidPoint"){
            RegisterDataBlock("PROD@MidPointModule-vz-cone0.7","TStnJetBlock",&fJetBlock);
        }
        else{
            RegisterDataBlock("PROD@JetCluModule-vz-cone0.7","TStnJetBlock",&fJetBlock);
        }
    }
    else if (_coneSize == 4) {
        if(_jetAlgo == "KtClus"){
            RegisterDataBlock("PROD@KtClusModule-vz-cone0.4","TStnJetBlock",&fJetBlock);
        }
        else if(_jetAlgo == "MidPoint"){
            RegisterDataBlock("PROD@MidPointModule-vz-cone0.4","TStnJetBlock",&fJetBlock);
        }
        else{
            RegisterDataBlock("JetBlock","TStnJetBlock",&fJetBlock);
        }
    }
    else if (_coneSize == 5) {
        if(_jetAlgo == "KtClus"){
            RegisterDataBlock("PROD@KtClusModule-vz-cone0.5","TStnJetBlock",&fJetBlock);
        }
        else if(_jetAlgo == "MidPoint"){
            RegisterDataBlock("PROD@MidPointModule-vz-cone0.5","TStnJetBlock",&fJetBlock);
        }
        else{
            RegisterDataBlock("PROD@JetCluModule-vz-cone0.5","TStnJetBlock",&fJetBlock);
        }
    }
    else
        cout << "check coneSize type" << endl;

       // book histograms
    BookHistograms();
//
    if ( !fHeaderBlock) {
        printf(" >>> branch *** %s *** doesn't exist \n", "HeaderBlock");
        fEnabled = 0;
    }
//
    if ( !fDcasDataBlock) {
        printf(" >>> branch *** %s *** doesn't exist \n", "DcasDataBlock");
        fEnabled = 0;
    }
//
    if ( !fCalDataBlock) {
        printf(" >>> branch *** %s *** doesn't exist \n", "CalDataBlock");
        fEnabled = 0;
    }
//
    if ( !fTrackBlock ) {
        printf(" >>> branch *** %s *** doesn't exist \n", "TrackBlock");
        fEnabled = 0;
    }
//
    if ( !fTriggerBlock ) {
        printf(" >>> branch *** %s *** doesn't exist \n", "TriggerBlock");
        fEnabled = 0;
    }
//
    if ( !fPhotonBlock ) {
        printf(" >>> branch *** %s *** doesn't exist \n", "PhotonBlock");
        fEnabled = 0;
    }
//
    if ( !fXftBlock  ) {
      printf(">>> branch *** %s *** doesn't exist \n", "XftBlock");
      fEnabled = 0;
    }
//
    if ( !fCesDataBlock  ) {
      printf(">>> branch *** %s *** doesn't exist \n", "CesDataBlock");
      fEnabled = 0;
    }
//
    if ( !fElectronBlock ) {
        printf(" >>> branch *** %s *** doesn't exist \n", "ElectronBlock");
        fEnabled = 0;
    }
//
    if ( !fMuonBlock ) {
        printf(" >>> branch *** %s *** doesn't exist \n", "MuonBlock");
        fEnabled = 0;
    }
//
    if ( !fJetBlock ) {
        printf(" >>> branch *** %s *** doesn't exist \n", "JetBlock");
        fEnabled = 0;
    }
//
    if ( !fMetBlock ) {
        printf(" >>> branch *** %s *** doesn't exist \n", "MetBlock");
        fEnabled = 0;
    }
//
    if ( !fVertexBlock ) {
        printf(" >>> branch *** %s *** doesn't exist \n", "ZVertexBlock");
        fEnabled = 0;
    }
//
    if ( !fCosmicBlock ) {
        printf(" >>> branch *** %s *** doesn't exist \n", "CosmicBlock");
        fEnabled = 0;
    }
//
    if ( !fL3SummaryBlock ) {
        printf(" >>> branch *** %s *** doesn't exist \n", "L3SummaryBlock");
        fEnabled = 0;
    }
//
    if ( !fGenpBlock ) {
        printf(" >>> branch *** %s *** doesn't exist \n", "GenpBlock");
        fEnabled = 0;
    }

  _initialized = true;
  _prevRun = -1;
  _prevGood = true;
  // Read in from the good run file
  FILE *ifp = fopen(_GoodRunList,"r");
    if ( ifp!=NULL ) {
        Int_t tmp =0;
        float xs=0.0;
        double offline_lumi=0.0;
        _good.clear();
        while (fscanf(ifp, "%i %f", &tmp, &xs) != -1) {
            offline_lumi += xs;
            _good.push_back(tmp);
        }
//        offline_lumi *= 0.951 * 1.019;  // 1-|Z<60| effect, 2-inelastic xsec
        offline_lumi *= 0.9699 * 1.019;  // 1-|Z<60| effect, 2-inelastic xsec
        std::sort(_good.begin(),_good.end());
        std::cout << "Good run list contains: " << _good.size() << " runs" << std::endl;
        std::cout << "Total offline integrated luminosity: " << offline_lumi
                  << " +/- " << offline_lumi*0.06 << " pb^-1" << std::endl;
        std::cout << "NOTE: two factors of 0.951 (|Zvtx|<60 cm) and 1.019 (inelastic xsec correction) have been applied." << std::endl;
        std::cout << "Offline luminosity uncertainty set to 6%" << std::endl;
        fclose(ifp);
    }
    else{
        std::cout << " GOODRUN FATAL ERROR opening file goodrun_em_cmup_nosi_v7.list " << std::endl;
        _initialized = false;
    }
    EvtTotalWeight2Lep=0.;
    EvtTotalWeight2LepErrorPlus=0.;
    EvtTotalWeight2LepErrorMinus=0.;

    nEvtCosmic = 0;
    Nevents_FakeIncremented = 0;
    Nevents_FakeIncrementedOpt = 0;
    Nevents_Total = 0;
    Nevents_PassAllCuts = 0;
    NeventsGoodEleReject = 0;
    EvtWeightTotal = 0.;
    EvtWeightErrorPlus = 0.;
    EvtWeightErrorMinus = 0.;
    EvtWeightOptTotal = 0.;
    EvtWeightOptErrorPlus = 0.;
    EvtWeightOptErrorMinus = 0.;
    fakeEvtWeightTotal = 0.;
    fakeEvtWeightErrorPlus = 0.;
    fakeEvtWeightErrorMinus = 0.;
    fakeEvtWeightGoodMassTotal = 0.;
    fakeEvtWeightGoodMassErrorPlus = 0.;
    fakeEvtWeightGoodMassErrorMinus = 0.;
    fakeEvtWeightOptTotal = 0.;
    fakeEvtWeightOptErrorPlus = 0.;
    fakeEvtWeightOptErrorMinus = 0.;
     
    for ( int i=0; i<10; i++ ) {
        EvtWeightControlRegion[i] = 0.;
        EvtWeightControlRegionErrorPlus[i] = 0.;
        EvtWeightControlRegionErrorMinus[i] = 0.;
        Nevents_ControlRegion[i] = 0;
    }
    Nevents_3goodlep_hepg = 0;
    Nevents_3goodele_hepg = 0;
    Nevents_3goodlep_base = 0;  // This is actually 2goodlep+track
    Nevents_3goodele_base = 0;
    Nevents_3goodlep_matched = 0;
    Nevents_BaseGoodMass = 0;
    Nevents_PassDilepOnly = 0;
    _HstFile = new TFile(_HstFileName,"recreate"); 

    nev_2ele=nev_2mu=nev_3obj=nev_3e=nev_3m=0;
    nev_sametrack=0;
    nev_eem=nev_mmm=0;

    return 0;
}

//______________________________________________________________________________

int TSeel::BeginRun()
{
//    if ( good_run(fHeaderBlock->RunNumber()) ){
//        TStnGoodRunList grl;
//        TStnRunSummary* rs;
//        grl.Init();
//        rs = grl.GetRunSummary(fHeaderBlock->RunNumber());
//        goodRunLuminosity += rs->OfflineLumiRS();
//    }

 

    if ( fVerbose>=1 ) cout << "Trying to connect to DB..." << endl;
    TStnDBManager* dbm = TStnDBManager::Instance();

    TStnBeamPos *fgCotBeamPos = (TStnBeamPos*) dbm->GetTable("CotBeamPos");
    TStnBeamPos *fgSvxBeamPos = (TStnBeamPos*) dbm->GetTable("SvxBeamPos");
    if ( fgSvxBeamPos->Status() >= 0 )
        fBeamPos = fgSvxBeamPos;
    else
        fBeamPos = fgCotBeamPos;

    h.beam_xy->Fill(fBeamPos->X0(),fBeamPos->Y0());
    fTriggerTable = (TStnTriggerTable*) dbm->GetTable("TriggerTable");

    return 0;
}

//______________________________________________________________________________

int TSeel::Event(int ientry)
{
//    if ( !fHeaderBlock->McFlag() && !good_run(fHeaderBlock->RunNumber()) ) return 1;
    if ( fHeaderBlock) fHeaderBlock->GetEntry(ientry);
    if ( fDcasDataBlock) fDcasDataBlock->GetEntry(ientry);
    if ( fCalDataBlock) fCalDataBlock->GetEntry(ientry);
    if ( fTrackBlock) fTrackBlock->GetEntry(ientry);
    if ( fTriggerBlock) fTriggerBlock->GetEntry(ientry);
    if ( fPhotonBlock) fPhotonBlock->GetEntry(ientry);
    if ( fXftBlock) fXftBlock->GetEntry(ientry);
    if ( fCesDataBlock) fCesDataBlock->GetEntry(ientry);
    if ( fElectronBlock) fElectronBlock->GetEntry(ientry);
    if ( fMuonBlock) fMuonBlock->GetEntry(ientry);
    if ( fJetBlock) fJetBlock->GetEntry(ientry);
    if ( fMetBlock) fMetBlock->GetEntry(ientry);
    if ( fVertexBlock) fVertexBlock->GetEntry(ientry);
    if ( fCosmicBlock) fCosmicBlock->GetEntry(ientry);
    if ( fL3SummaryBlock) fL3SummaryBlock->GetEntry(ientry);
    if ( _JES_IMode==0 ) if ( fGenpBlock ) fGenpBlock->GetEntry(ientry);
  
//    if ( fCosmicBlock->HasCosmicRay() ) return 1;
//    if ( fCosmicBlock->HasNotOutgoingPair() ) {
//        nEvtCosmic++;
//        return 1;
//    }
//    

    Nevents_Total++;
    if ( _Data==1 && !good_run(fHeaderBlock->RunNumber()) ) return 1;
    if ( fVertexBlock->NVertices()<=0 ) return 1;

    
     iVtx = -1;
    for ( int i=0; i<fVertexBlock->NVertices(); i++ ) {
        TStnVertex *v = fVertexBlock->Vertex(i);
        if ( fVerbose>0 ) cout << "vx=" << v->X() << " vy=" << v->Y() << " vz=" << v->Z() << " sumpt=" << v->SumPt() << endl;
        if ( v->VClass()>=12 ) {
            nVertex++;
            if ( v->SumPt()>sumpt ) {
                sumpt = v->SumPt();
                iVtx = i;
            }
        }
    }
    //float vz;
    if ( iVtx>=0 ) {
        TStnVertex *vtx = fVertexBlock->Vertex(iVtx);
	vz = vtx->Z();
    }    
    else
      vz = 0;
    h.vz->Fill(vz);
    if ( fabs(vz)>60. ) return 1;
// CosmicRay
    if ( fCosmicBlock->HasCosmicRay() ) cosmicflag |= 0x1;
    if ( fCosmicBlock->HasOOTCaloE() ) cosmicflag |= 0x2;
    if ( fCosmicBlock->HasLowQFrac() ) cosmicflag |= 0x4;
    if ( fCosmicBlock->HasNotOutgoingPair() ) cosmicflag |= 0x8;
    if ( getbits(cosmicflag,3,1) ) {
      nEvtCosmic++;
      return 1;
    }
// Trigger Names
    TTl3d* fTl3d = fTriggerBlock->Tl3d();
    int npaths = fTl3d->NPaths();
    for ( int i=0; i<npaths; i++ ) {
        if ( fTl3d->PathPassed(i) ) {
            if ( i < fTriggerTable->NL3Triggers() ) {
                const TStnTrigger* trig = fTriggerTable->GetTrigger(3,i);
                if ( fVerbose>0 ) printf("%10d%10d%s%s \n",fHeaderBlock->RunNumber(),fHeaderBlock->EventNumber()," Trigger Name: ",trig->Name());
                for ( int j=0; j<15; j++ ) {
                    if ( strcmp(L3_name[j],trig->GetName())==0 ) trigger_bit |= 0x1 << j;
                }
            }
        }
    }    
    if ( _Data==1 && getbits(trigger_bit,0,1)!=1 ) return 1;    


//  Missing Et variables:
    missing_et.set_ex(fMetBlock->MetX(3));
    missing_et.set_ey(fMetBlock->MetY(3));
    missing_et.set_et(fMetBlock->Met(3));
    missing_et.set_phi(fMetBlock->MetPhi(3));

//
    nleptons=0; nleptons_temp=0; ngood_leptons=0;
    int nele = 0;
    int ntightele = 0;
//    int track_number[100];
    unsigned quality;
    for ( int i=0; i<fElectronBlock->NElectrons(); i++ ) {
        TStnElectron *el = fElectronBlock->Electron(i);
	if ( ! el->IsCentral() ) continue;
	int it = el->TrackNumber();
	if ( it<0 ) continue;
	TStnTrack *trk = fTrackBlock->Track(it);
	TLorentzVector *mom = new TLorentzVector();
	int status = trk->GetBcMomentum(mom);
//            float corrpt = curvCorr(mom->Pt(),el->Charge(),mom->Phi(),trk->Lam0());
//            float corrp = corrpt*cosh(mom->Eta());
//            float EoP = el->EmE()*(el->Etcor()>0. ? el->Etcor(): 1.)/corrp;
	float d0 = CorrectedD0(trk);
        float eta = el->Momentum()->Eta();
        float phi = TVector2::Phi_0_2pi(el->Momentum()->Phi());
        float et = el->EmEt()*(el->Etcor()>0. ? el->Etcor(): 1.);
        float pt = el->TrackPt();
        if ( et<4. && pt<4. ) continue;
        //if ( ! matched_hepg_ele(el) ) continue;
        quality = electron_quality(el);
	if ( isConversion(el) ) {
	  if ( getbits(quality,23,1) && matched_hepg_ele(el) )
	    h.EleEt_Conv->Fill(el->EmEt()*(el->Etcor()>0. ? el->Etcor(): 1.));
	  continue;
	}
	Fill_eID(el,quality);
	if( (getbits(quality,23,1)==1 || getbits(quality,22,1)==1) && nele<50 ){
	    lep_temp[nele].set_px(pt*cos(phi));
            lep_temp[nele].set_py(pt*sin(phi));
            lep_temp[nele].set_pz(pt*sinh(eta));
            lep_temp[nele].set_pt(pt);
            lep_temp[nele].set_en(el->EmE()*(el->Etcor()>0. ? el->Etcor(): 1.));
            lep_temp[nele].set_et(el->EmEt()*(el->Etcor()>0. ? el->Etcor(): 1.));
            lep_temp[nele].set_eta(eta);
            lep_temp[nele].set_phi(phi);
            lep_temp[nele].set_z0(el->Z0());
            lep_temp[nele].set_id(-11*el->Charge());
            lep_temp[nele].set_qual(quality);
            lep_temp[nele].set_index(i);
            lep_temp[nele].set_mc(matched_hepg_ele(el));
	    //track_number[nele] = it;
            nele++;
	    if(getbits(quality,23,1)==1)
	      ntightele++;
	}

    }
    //cout<<"nele="<<nele<<" , "<<fElectronBlock->NElectrons()<<endl;
    int nmuon = 0;
    for ( int i=0; i<fMuonBlock->NMuons(); i++ ) {
        TStnMuon *mu = fMuonBlock->Muon(i);
        int itrk = mu->TrackNumber();
// Make sure that the muon doesn't coincide with an ele
// that is already in lep_temp array
//	bool already_used_track = false;
//	for(int cnt=0;cnt<nele;cnt++ ) {
//	  if ( track_number[cnt] = itrk ) {
//	    already_used_track = true;
//	    break;
//	  }
//	}
//if ( already_used_track ) continue;
        TStnTrack* trk = fTrackBlock->Track(itrk);
//        TLorentzVector *mom = new TLorentzVector();
//        int status1 = trk->GetBcMomentum(mom);
//        float eta = mom->Eta();
//        float phi = TVector2::Phi_0_2pi(mom->Phi());
//        if ( mom->Pt() < 5. ) continue;
        //if ( ! matched_hepg_mu(mu) ) continue;
        float pt = mu->TrackPt();
        float eta = mu->Momentum()->Eta();
        float phi = TVector2::Phi_0_2pi(mu->Momentum()->Phi());
        if ( pt<4. ) continue;
	if(_OptScheme==3)
	  quality = muon_quality_hipt(mu);
	else
	  quality = muon_quality(mu);
	Fill_mID(mu,quality);
	//	if ( getbits(quality,17,1)!=1 && getbits(quality,18,1)!=1 && getbits(quality,19,1)!=1 && getbits(quality,20,1)!=1 && getbits(quality,21,1)!=1 ) continue;  // requires everything	
	//if ( getbits(quality,17,1)==1 ) continue;
	if ( getbits(quality,17,1)==1 || getbits(quality,18,1)==1 || getbits(quality,19,1)==1 || getbits(quality,20,1)==1 || getbits(quality,21,1)==1 ){
	  if( getbits(quality,23,1)==1 && nele+nmuon<50 ) {
	    lep_temp[nele+nmuon].set_px(pt*cos(phi));
            lep_temp[nele+nmuon].set_py(pt*sin(phi));
            lep_temp[nele+nmuon].set_pz(pt*sinh(eta));
            lep_temp[nele+nmuon].set_pt(pt);
            lep_temp[nele+nmuon].set_en(sqrt(pow(pt*cosh(eta),2)+pow(0.1057,2)));
            lep_temp[nele+nmuon].set_et(pt);
            lep_temp[nele+nmuon].set_eta(eta);
            lep_temp[nele+nmuon].set_phi(phi);
            lep_temp[nele+nmuon].set_z0(mu->Z0());
            lep_temp[nele+nmuon].set_id(-13*mu->Charge());
            lep_temp[nele+nmuon].set_qual(quality);
            lep_temp[nele+nmuon].set_qual2(0x0);
            lep_temp[nele+nmuon].set_index(i);
            lep_temp[nele+nmuon].set_mc(matched_hepg_mu(mu));
            nmuon++;
	  }
	}
    }

//    Order good leptons in et/pt
    nleptons_temp = nele + nmuon;
    if(fVerbose<=-2 && nleptons_temp>1){
      cout<<"nele_temp = "<<nele<<" , nmu_temp = "<<nmuon<<" , nleptons_temp = "<<nleptons_temp<<endl;
    }
    if(nele>1)
      nev_2ele++;
    if(nmuon>1)
      nev_2mu++;
    if(nele>1 && nleptons_temp>2){
      nev_3obj++;
      for(int it=2; it<nleptons_temp; it++){
	if( abs(lep_temp[it].id())==11 ){
	  if( pass_track_cut(lep_temp[it].index(),1) ){
	    nev_3e++;
	    //cout<<"eee evt run/event "<<fHeaderBlock->RunNumber()<<" "<<fHeaderBlock->EventNumber()<<endl;
	    break;
	  }
	}
      }
      for(int it=2; it<nleptons_temp; it++){
	if( abs(lep_temp[it].id())==13 ){
	  if( pass_track_cut(lep_temp[it].index(),2) ){
	    nev_3m++;
	    //cout<<"eem evt run/event "<<fHeaderBlock->RunNumber()<<" "<<fHeaderBlock->EventNumber()<<endl;
	    break;
	  }
	}
      }
    }
    
    Lepton lep_hold;
    for ( int i=0; i<nleptons_temp-1; i++ ) {
        for ( int j=i+1; j<nleptons_temp; j++ ) {
	    if ( lep_temp[i].et() < lep_temp[j].et() ) {
	        lep_hold = lep_temp[i];
		lep_temp[i] = lep_temp[j];
		lep_temp[j] = lep_hold;
	    }
	}
    }
    if(nele>1 && nmuon>0){
      if( abs(lep_temp[0].id())==11 && abs(lep_temp[1].id())==11 ){
	nev_eem++;
      }
    }
    if(nmuon>2){
      if( abs(lep_temp[0].id())==13 && abs(lep_temp[1].id())==13 && abs(lep_temp[2].id())==13 ){
	nev_mmm++;
      }
    }
    
//=================================================================================
    float met_before_corrections = missing_et.et();
    float met_phi_before_corrections = missing_et.phi();
    njet = 0; njet15 = 0; njet_nominal = 0; njet_plus = 0; njet_minus = 0;
    //int lep_subtracted[50]={0};
    int syscode=0;
    int coneSize; // 0=0.4, 1=0.7 and 2=1.0
    if ( _coneSize == 4 )
        coneSize = 0;
    else if ( _coneSize == 7 )
        coneSize = 1;
    JetEnergyCorrections b3 = JetEnergyCorrections("JetCorrections","JetCorrections",3,nVertex,coneSize,_JES_version,syscode,fHeaderBlock->RunNumber(),_JES_IMode);
    JetEnergyCorrections b4 = JetEnergyCorrections("JetCorrections","JetCorrections",4,nVertex,coneSize,_JES_version,syscode,fHeaderBlock->RunNumber(),_JES_IMode);
    JetEnergyCorrections b5 = JetEnergyCorrections("JetCorrections","JetCorrections",5,nVertex,coneSize,_JES_version,syscode,fHeaderBlock->RunNumber(),_JES_IMode);
    JetEnergyCorrections b6 = JetEnergyCorrections("JetCorrections","JetCorrections",6,nVertex,coneSize,_JES_version,syscode,fHeaderBlock->RunNumber(),_JES_IMode);
    JetEnergyCorrections b7 = JetEnergyCorrections("JetCorrections","JetCorrections",7,nVertex,coneSize,_JES_version,syscode,fHeaderBlock->RunNumber(),_JES_IMode);
    b3.setCemCorON(false);
    b4.setCemCorON(false);
    b5.setCemCorON(false);
    b6.setCemCorON(false);
    b7.setCemCorON(false);
    
    int ntemp_jets = 0;
    for ( int i=0; i<fJetBlock->NJets(); i++ ) {
        TStnJet *jet = fJetBlock->Jet(i);
        if ( jet->Et()<5. || fabs(jet->DetEta())>2.5 ) continue;
       HepLorentzVector P4Jet;
        HepLorentzVector P4JetCor;
        P4Jet.setPx(jet->Momentum()->Px());
        P4Jet.setPy(jet->Momentum()->Py());
        P4Jet.setPz(jet->Momentum()->Pz());
        P4Jet.setE(jet->Momentum()->E());

        float emf=jet->Emfr();
        float &emFraction=emf;
        float DetJet=jet->DetEta();

        float scale3=1.; float scale4=1.; float scale5=1.; float scale6=1.; float scale7=1.;
        scale3 = b3.doEnergyCorrections(P4Jet,emFraction,DetJet);
        scale4 = b4.doEnergyCorrections(P4Jet,emFraction,DetJet);
        scale5 = b5.doEnergyCorrections(P4Jet,emFraction,DetJet);
        scale6 = b6.doEnergyCorrections(P4Jet,emFraction,DetJet);
        scale7 = b7.doEnergyCorrections(P4Jet,emFraction,DetJet);
        //h.jet_scaleEta->Fill(DetJet,scale5);
	
	P4JetCor= P4Jet*scale5;
        float corrEt = P4JetCor.e()*sin(P4JetCor.theta());
        float phi = TVector2::Phi_0_2pi(P4JetCor.phi());
        float eta = -log(tan(P4JetCor.theta()/2.));
//     Fill with the uncorrected Jet energy for now.
        if ( pass_jet_cut(jet->Et(),DetJet) && emf<0.9 && ntemp_jets<100 ) {
//        if ( pass_jet_cut(jet->Et(),DetJet) && ntemp_jets<100 ) {
            jt_temp[ntemp_jets].set_ex(jet->Et()*cos(phi));
            jt_temp[ntemp_jets].set_ey(jet->Et()*sin(phi));
            jt_temp[ntemp_jets].set_et(jet->Et());
            jt_temp[ntemp_jets].set_pt(sqrt(pow(P4Jet.px(),2)+pow(P4Jet.py(),2)));
            jt_temp[ntemp_jets].set_eta(eta);
            jt_temp[ntemp_jets].set_phi(phi);
            jt_temp[ntemp_jets].set_en(P4Jet.e());
            jt_temp[ntemp_jets].set_deta(DetJet);
            jt_temp[ntemp_jets].set_scale3(scale3);
            jt_temp[ntemp_jets].set_scale4(scale4);
            jt_temp[ntemp_jets].set_scale5(scale5);
            jt_temp[ntemp_jets].set_scale6(scale6);
            jt_temp[ntemp_jets].set_scale7(scale7);
            jt_temp[ntemp_jets].set_emf(emf);
            float ez = pow(jt_temp[ntemp_jets].en(),2) - pow(jt_temp[ntemp_jets].et(),2);
            if ( ez < 0. ) {
                cout << "Jet En " << jt_temp[i].en() << " < Jet Et " << jt_temp[i].et() << endl;
                ez = 0.;
            }
            jt_temp[ntemp_jets].set_ez(sqrt(ez));
            jt_temp[ntemp_jets].set_index(i);
            ntemp_jets++;
        }
    }
//  remove leptons from a jet if the leptons are inside a jet
    for ( int i=0; i<ntemp_jets; i++ ) {
        TStnJet *jet = fJetBlock->Jet(jt_temp[i].index());
        HepLorentzVector P4Jet;
        HepLorentzVector P4JetCor;
        P4Jet.setPx(jet->Momentum()->Px());
        P4Jet.setPy(jet->Momentum()->Py());
        P4Jet.setPz(jet->Momentum()->Pz());
        P4Jet.setE(jet->Momentum()->E());
        if ( jet->Et()<5. || fabs(jet->DetEta())>2.5 ) continue;
	float scale3=1.;float scale4=1.;float scale5=1.;float scale6=1.;float scale7=1.;
	float new_emf = jet->Emfr();
	float new_DetEta = jet->DetEta();
	float eta = v2eta(P4Jet.px(),P4Jet.py(),P4Jet.pz());
	float phi = v2phi(P4Jet.px(),P4Jet.py());
	float et = P4Jet.e()/cosh(eta);
	if ( et>5. ) {
                    scale3 = b3.doEnergyCorrections(P4Jet,new_emf,new_DetEta);
                    scale4 = b4.doEnergyCorrections(P4Jet,new_emf,new_DetEta);
                    scale5 = b5.doEnergyCorrections(P4Jet,new_emf,new_DetEta);
                    scale6 = b6.doEnergyCorrections(P4Jet,new_emf,new_DetEta);
                    scale7 = b7.doEnergyCorrections(P4Jet,new_emf,new_DetEta);
	}
	if ( pass_jet_cut(et,new_DetEta) ) {
	    jt_temp[i].set_ex(et*cos(phi));
	    jt_temp[i].set_ey(et*sin(phi));
	    jt_temp[i].set_et(et);
	    jt_temp[i].set_pt(sqrt(pow(P4Jet.px(),2)+pow(P4Jet.py(),2)));
	    jt_temp[i].set_eta(eta);
	    jt_temp[i].set_phi(phi);
	    jt_temp[i].set_en(P4Jet.e());
	    jt_temp[i].set_deta(new_DetEta);
	    jt_temp[i].set_scale3(scale3);
	    jt_temp[i].set_scale4(scale4);
	    jt_temp[i].set_scale5(scale5);
	    jt_temp[i].set_scale6(scale6);
	    jt_temp[i].set_scale7(scale7);
	    jt_temp[i].set_emf(new_emf);
	    float ez = pow(jt_temp[i].en(),2) - pow(jt_temp[i].et(),2);
	    if ( ez < 0. ) {
	      cout << "Jet En " << jt_temp[i].en() << " < Jet Et " << jt_temp[i].et() << endl;
	      ez = 0.;
	    }
	    jt_temp[i].set_ez(sqrt(ez));
	    jt_temp[i].set_index(jt_temp[i].index());
	}
	else {
	  jt_temp[i].set_index(1000000+jt_temp[i].index());
	}   
    }
    ngood_jets = 0;
//    for ( int i=0; i<ntemp_jets; i++ ) {
//        jt[ngood_jets] = jt_temp[i];
//        ngood_jets++;
//    }
    for ( int i=0; i<ntemp_jets; i++ ) {
      if  ( jt_temp[i].index()>=1000000 ) continue;
      float eta = jt_temp[i].eta();
      float phi = jt_temp[i].phi();
//      bool jt_lep_overlap = false;
      float jt_en = jt_temp[i].et()*cosh(eta);
      float jt_ex = jt_temp[i].et()*cos(phi);
      float jt_ey = jt_temp[i].et()*sin(phi);
      float jt_ez = jt_temp[i].et()*sinh(eta);
// Since there are no isolated tracks here, I need not check if a trk is 
// contained within a jet. I should go ahead and assume all leptons are
// away from jets and correct the jets --> then later after I have a
// a corrected jet, I apply dR between leptons and jets to ensure
// separation. (=> jt_lep_overlap is always false)
      float jt_et, deta;
      float corrEt, corrEt_plus, corrEt_minus, factor, factor_plus, factor_minus;
//    if ( !jt_lep_overlap ) {
	jt_et = jt_temp[i].et();
	deta = jt_temp[i].deta();
	factor = jt_temp[i].scale5();
//    }
//    else {
//	jt_et = sqrt(pow(jt_ex,2)+pow(jt_ey,2));
//	deta = DetEta(eta,fEvtHdr_fVz);
//	if ( jt_en>0. ) factor = b5.doEnergyCorrections(jt_et,fJets_fEmf[i],deta);
//    }
      if ( jt_en<0. || fabs(deta)>2.5 ) continue;
      b5.setTotalSysUncertainties(1);
      factor_plus = b5.doEnergyCorrections(jt_et,jt_temp[i].emf(),deta);
      b5.setTotalSysUncertainties(-1);
      factor_minus = b5.doEnergyCorrections(jt_et,jt_temp[i].emf(),deta);
      corrEt = jt_et*factor;
      corrEt_plus = jt_et*factor_plus;
      corrEt_minus = jt_et*factor_minus;
          if ( pass_jet_cut(corrEt,deta) && njet_nominal<100 ) {
               jt_nominal[njet_nominal].set_ex(corrEt*cos(phi));
               jt_nominal[njet_nominal].set_ey(corrEt*sin(phi));
               jt_nominal[njet_nominal].set_et(corrEt);
               jt_nominal[njet_nominal].set_pt(jt_et/jt_temp[i].et()*jt_temp[i].pt()*factor);
               jt_nominal[njet_nominal].set_eta(eta);
               jt_nominal[njet_nominal].set_phi(phi);
               jt_nominal[njet_nominal].set_deta(deta);
               jt_nominal[njet_nominal].set_en(corrEt*cosh(eta));
               float ez = pow(jt_nominal[njet_nominal].en(),2) - pow(jt_nominal[njet_nominal].et(),2);
               if ( ez < 0. ) {
                   cout << "Jet En " << jt_nominal[i].en() << " < Jet Et " << jt_nominal[i].et() << endl;
                   ez = 0.;
               }
               jt_nominal[njet_nominal].set_ez(sqrt(ez));
               jt_nominal[njet_nominal].set_emf(jt_temp[i].emf()); // needs fix
               jt_nominal[njet_nominal].set_index(jt_temp[i].index());
               jt_nominal[njet_nominal].set_scale3(jt_temp[i].scale3());
               jt_nominal[njet_nominal].set_scale4(jt_temp[i].scale4());
               jt_nominal[njet_nominal].set_scale5(factor);
               jt_nominal[njet_nominal].set_scale6(jt_temp[i].scale6());
               jt_nominal[njet_nominal].set_scale7(jt_temp[i].scale7());
               njet_nominal++;
           }

           if ( pass_jet_cut(corrEt_plus,deta) && njet_plus<100 ) {
               jt_plus[njet_plus].set_ex(corrEt_plus*cos(phi));
               jt_plus[njet_plus].set_ey(corrEt_plus*sin(phi));
               jt_plus[njet_plus].set_et(corrEt_plus);
               jt_plus[njet_plus].set_pt(jt_et/jt_temp[i].et()*jt_temp[i].pt()*factor_plus);
               jt_plus[njet_plus].set_eta(eta);
               jt_plus[njet_plus].set_phi(phi);
               jt_plus[njet_plus].set_deta(deta);
               jt_plus[njet_plus].set_en(corrEt_plus*cosh(eta));
               float ez = pow(jt_plus[njet_plus].en(),2) - pow(jt_plus[njet_plus].et(),2);
               if ( ez < 0. ) {
                   cout << "Jet En " << jt_plus[i].en() << " < Jet Et " << jt_plus[i].et() << endl;
                   ez = 0.;
               }
               jt_plus[njet_plus].set_ez(sqrt(ez));
               jt_plus[njet_plus].set_emf(jt_temp[i].emf());  // needs fix
               jt_plus[njet_plus].set_scale3(jt_temp[i].scale3());
               jt_plus[njet_plus].set_scale4(jt_temp[i].scale4());
               jt_plus[njet_plus].set_scale5(factor_plus);
               jt_plus[njet_plus].set_scale6(jt_temp[i].scale6());
               jt_plus[njet_plus].set_scale7(jt_temp[i].scale7());
               jt_plus[njet_plus].set_index(i);
               njet_plus++;
           }
           if ( pass_jet_cut(corrEt_minus,deta) && njet_minus<100 ) {
               jt_minus[njet_minus].set_ex(corrEt_minus*cos(phi));
               jt_minus[njet_minus].set_ey(corrEt_minus*sin(phi));
               jt_minus[njet_minus].set_et(corrEt_minus);
               jt_minus[njet_minus].set_pt(jt_et/jt_temp[i].et()*jt_temp[i].pt()*factor_minus);
               jt_minus[njet_minus].set_eta(eta);
               jt_minus[njet_minus].set_phi(phi);
               jt_minus[njet_minus].set_deta(deta);
               jt_minus[njet_minus].set_en(corrEt_minus*cosh(eta));
               float ez = pow(jt_minus[njet_minus].en(),2) - pow(jt_minus[njet_minus].et(),2);
               if ( ez < 0. ) {
                   cout << "Jet En " << jt_minus[i].en() << " < Jet Et " << jt_minus[i].et() << endl;
                   ez = 0.;
               }
               jt_minus[njet_minus].set_ez(sqrt(ez));
               jt_minus[njet_minus].set_emf(jt_temp[i].emf());  // needs fix
               jt_minus[njet_minus].set_scale3(jt_temp[i].scale3());
               jt_minus[njet_minus].set_scale4(jt_temp[i].scale4());
               jt_minus[njet_minus].set_scale5(factor_minus);
               jt_minus[njet_minus].set_scale6(jt_temp[i].scale6());
               jt_minus[njet_minus].set_scale7(jt_temp[i].scale7());
               jt_minus[njet_minus].set_index(i);
               njet_minus++;
           }
    }
     
    if ( _JEScale == 0 ) {
      for ( int i=0; i<njet_nominal; i++ ) jt[i] = jt_nominal[i];
      njet = njet_nominal;
    }
    else if ( _JEScale == 1 ) {
      for ( int i=0; i<njet_plus; i++ ) jt[i] = jt_plus[i];
      njet = njet_plus;
    }
    else if ( _JEScale == -1 ) {
      for ( int i=0; i<njet_minus; i++ ) jt[i] = jt_minus[i];
      njet = njet_minus;
    }
    ht = 0.;
    for ( int i=0; i<njet; i++ ) {
      ht += jt[i].et();
      //if ( jt[i].et()>15. ) {
      if ( jt[i].et()>15. ) {
	jt15[njet15] = jt[i];
	njet15++;
      }
    }
    
    if(njet15>0){
      h.dphi_j1met_nocorr->Fill(delta_phi(jt15[0].phi(),met_phi_before_corrections));
      h.raw_phi_jet1->Fill(jt15[0].phi());
    }
       h.met_noCorr->Fill(met_before_corrections);
       h.met_phi_noCorr->Fill(met_phi_before_corrections);
       correct_met_jets(FINAL_JET_CUT,2.5);
       h.met_JetCorr->Fill(missing_et.et());
       h.met_phi_JetCorr->Fill(missing_et.phi());
       correct_met_muons();
       h.met_MipCorr->Fill(missing_et.et());
       h.met_phi_MipCorr->Fill(missing_et.phi());
// Do the hepg leptons
// The order is: ele (pt sorted) then mu (pt sorted)

//Event Cuts
/*
  HepgEventReset();
  if ( pass_HepgEventReset(_NLepAna) ) {
  Set_HepgEventType(_NLepAna);
  if ( (_EventType_hepg & 0x524) == 0x524 ) {
  Nevents_3goodele_hepg++;
  Nevents_3goodlep_hepg++;
  }
  else if ( (_EventType_hepg & 0x522) == 0x522 ) {
  Nevents_3goodlep_hepg++;
  }
  if ( fVerbose>=10 ) cout << "nlep_hepg=" << nlep_hepg << " ngood_ele_hepg=" << ngood_ele_hepg << " _EventType_hepg=" << hex <<_EventType_hepg << dec << " 3lep=" << Nevents_3goodlep_hepg << " 3ele=" << Nevents_3goodele_hepg << endl;
  }
*/
       // -----------------------------------
       // EventReset() does the following: 
       //   writes lep[] (good leps already sorted in et, not including the track)
       //   counts ngood_leptons (ele for now)
       //   DeltaR's (lep-jet, lep-lep, including the track) are taken care of
       EventReset();
       // Plotting jet_mass/jet_et vs jet_et
       for ( int i=0; i<njet; i++ ) {
	 TStnJet *jet = fJetBlock->Jet(jt[i].index());
	 h.jet_width[0]->Fill(jet->M()/jet->Et(),jet->Et());
       }
       if ( njet>0 ) {
	 TStnJet *jet = fJetBlock->Jet(jt[0].index());
	 h.jet_width[1]->Fill(jet->M()/jet->Et(),jet->Et());
       }

       if ( nleptons<2 )
	 Mlep12 = -1.;
       else {
	 vector<TLorentzVector> v_lep;
	 TLorentzVector vi(lep[0].px(),lep[0].py(),lep[0].pz(),lep[0].en());
	 v_lep.push_back(vi);
	 TLorentzVector vj(lep[1].px(),lep[1].py(),lep[1].pz(),lep[1].en());
	 v_lep.push_back(vj);
	 Mlep12 = inv_mass(v_lep);
       }
       _EventType = 0x0;
       // -----------------------------------
       // This requires:
       //   2 good leptons (ele for now) + 1 track or lepton
       //   satisfying pt requirement (_LEP_PT_CUT3[3])
       //   Lepton charge combination is not required
       //   no pt re-ordering (of good leptons and the 3rd track) is done
       if ( pass_EventReset(_NLepAna) ) {
           // Define the event as ee+e or ee+track
           // flag the correct charge combination
           Set_EventType(_NLepAna);
           // -----------------------------------
           // pass_BaseCuts() = correct charge combination
           // pass_MassCuts(): whether m_os_1 and m_os_2 satisfy requirements
	   // has_radiated_photon = checks is either of the leading two ele
	   // has radiated a photon. This is done only for DY background, for
	   // all other backgrounds this is true. ALthough maybe for WW also it
	   // should be checked.
           //   After Pass_Trilep_BaseCuts, 3 leptons are identified,
           //   charge combination is correct.
           //   We alter the lepton order in pt later in FillHistograms(),
           //   but pass_MassCuts() will not change its decision.
	   Pass_Trilep_BaseCuts = pass_BaseCuts();
	   if ( Pass_Trilep_BaseCuts ) {
	      Pass_Trilep_MassCuts = pass_MassCuts();
	      Pass_RadiatePhoton = has_radiated_photon();
	   }
	   else {
	      Pass_Trilep_MassCuts = false;
	      Pass_RadiatePhoton = false;
	   }
       }
       else {
	 Pass_Trilep_BaseCuts = false;
	 Pass_Trilep_MassCuts = false;
	 Pass_RadiatePhoton = false;
       }
       if ( pass_EventReset(2) ){
	 if ( abs(lep[0].id())==11 && abs(lep[1].id())==11 )
	   Pass_Dilep_BaseCuts = true;
       }
       else
	 Pass_Dilep_BaseCuts = false;
       // -----------------------------------
       // Weight the event with trig/ID effciency
       // and if a dilep event without 3rd hepg_e or hepg_mu matching track
       // weight the event with random track rate
       int return_code = Set_EventWeight();
       // The event weights should be summed here.
       // One problem is that later on we may have cuts which involves
       // the third track. We may apply an event weight on a DY MC event with
       // two tight electrons but with no isotrack. The efficiency of the
       // cut which involves a 3rd track is not applicable for such an event.
       // However, we can count the actual number of dilep+track events
       // which pass base cuts and obtain an efficiency of the cut.
       if ( _Signal || return_code==2 ) {
       // both base and mass cuts have been applied
       // either a signal event or a dilep+track event
           EvtWeightTotal += _EventWeight;
           EvtWeightErrorPlus += _EventWeightErrorPlus;
           EvtWeightErrorMinus += _EventWeightErrorMinus;
       }       
       if ( return_code==-5 ) {
           fakeEvtWeightTotal += _fakeEventWeight;
           fakeEvtWeightErrorPlus += _fakeEventWeightErrorPlus;
           fakeEvtWeightErrorMinus += _fakeEventWeightErrorMinus;
	   fakeEvtWeightOptTotal += _fakeOptEventWeight;
	   //for(int iRegion=0; iRegion<4; i++)
	   //  fakeEvtWeightControlRegion[iRegion] += _fakeEventWeightControlRegion[iRegion];
       }       
       else if ( return_code == 3 ) {  // a dilep-only event
           Nevents_PassDilepOnly++;
           EvtTotalWeight2Lep += _EventWeight;
           EvtTotalWeight2LepErrorPlus += _EventWeightErrorPlus;
           EvtTotalWeight2LepErrorMinus += _EventWeightErrorMinus;
       }
       // -----------------------------------
       // The rest of the code picks the actual dilep+track events.
       // -----------------------------------
       if ( Pass_Trilep_BaseCuts && Pass_RadiatePhoton ) {
           // -----------------------------------
           // The leptons are re-ordered according to their et
           // Decide not to do it on Feb.9,2005.
           // It does NOT alter _EventType
           // Require mass cuts and put weights in histograms, 20050205
           if ( Pass_Trilep_MassCuts ) FillHistograms();
           if ( Pass_Trilep_MassCuts ) {
               Nevents_BaseGoodMass++;
               for ( int iRegion=0; iRegion<4; iRegion++ ) {
                   if ( pass_ControlRegionCuts(iRegion) ) {
                       FillControlRegionHistograms(iRegion);
                       Nevents_ControlRegion[iRegion]++;
                       if ( _Signal || return_code==2 ) {
                           EvtWeightControlRegion[iRegion] += _EventWeight;
                           EvtWeightControlRegionErrorPlus[iRegion] += _EventWeightErrorPlus;
                           EvtWeightControlRegionErrorMinus[iRegion] += _EventWeightErrorMinus;
                       }
                   }
               }
               if ( pass_OptimizedCuts() ) {
                   if ( ! pass_ControlRegionCuts(4) ) {//veto events that pass HiPt cuts
                       FillControlRegionHistograms(4);
                       Nevents_ControlRegion[4]++;
                   }
                   if ( _Signal || return_code==2 ) {
                       EvtWeightOptTotal += _EventWeight;
                       EvtWeightOptErrorPlus += _EventWeightErrorPlus;
                       EvtWeightOptErrorMinus += _EventWeightErrorMinus;
                       if ( ! pass_ControlRegionCuts(4) ) {
                    // We keep events that pass our opt cuts but fail the HiPt cuts
                           EvtWeightControlRegion[4] += _EventWeight;
                           EvtWeightControlRegionErrorPlus[4] += _EventWeightErrorPlus;
                           EvtWeightControlRegionErrorMinus[4] += _EventWeightErrorMinus;
                       }
                   }
	 
                   Nevents_PassAllCuts++;
		   if ( fVerbose<-2 ) {
		     for ( int i=0; i<nleptons; i++ ) {
		       cout << "lep " << i << " " << lep[i].et() << " " << lep[i].eta() << " " << lep[i].phi() << " " << lep[i].id() << endl;
		     }
		   }

               }  // pass all  optimized cuts
               if ( fVerbose>=2 ) cout << " _EventType=" << hex << _EventType << dec << " 3lep=" << Nevents_3goodlep_base << " 3ele=" << Nevents_3goodele_base << " 3match=" << Nevents_3goodlep_matched << endl;
           }
       }  // pass base cuts
       else if ( fVerbose>10 ) {
       // only wrong charge combination or bad opposite-sign masses 
       // ee+e or ee+track events enter here
           cout << "This can only happen for wrong charge combinations or bad opposite-sign masses! 2005/1/23" << endl;
           for ( int i=0; i<nleptons; i++ ) {
               cout << "lep " << i << " " << lep[i].et() << " " << lep[i].eta() << " " << lep[i].phi() << " " << lep[i].id() << endl;
           }
       }



//----------------------------
    _HstFile->cd();
    
    return 0;
}
//______________________________________________________________________________________________

int TSeel::EndJob()
{

    float eff[3] = {0.};
    float eff0, err_eff0;
    if ( Nevents_Total==0 || Nevents_BaseGoodMass==Nevents_Total ) {
        TSeel::Bayesian_Error(Nevents_BaseGoodMass,Nevents_Total,eff);
        eff0 = eff[0];
        err_eff0 = max(fabs(eff[2]),fabs(eff[1]));
    }
    else {
        TSeel::Binomial_Error(Nevents_BaseGoodMass,Nevents_Total,eff);
        eff0 = eff[0];
        err_eff0 = (eff[2]-eff[1])/2.;
    }
//  Fraction of dilep+track events that pass mass cuts
    for ( int i=0; i<3; i++ ) eff[i]=0.;
    float eff1, err_eff1;
    if ( Nevents_BaseGoodMass==0 || Nevents_BaseGoodMass==Nevents_3goodlep_base ) {
        TSeel::Bayesian_Error(Nevents_BaseGoodMass,Nevents_3goodlep_base,eff);
        eff1 = eff[0];
        err_eff1 = max(fabs(eff[2]),fabs(eff[1]));
    }
    else {
        TSeel::Binomial_Error(Nevents_BaseGoodMass,Nevents_3goodlep_base,eff);
        eff1 = eff[0];
        err_eff1 = (eff[2]-eff[1])/2.;
    }
    EvtWeightTotal3LepGoodMass = EvtWeightTotal;
    EvtWeight3LepGoodMassErrorPlus = EvtWeightErrorPlus;
    EvtWeight3LepGoodMassErrorMinus = EvtWeightErrorMinus;
    /*
    if ( _Signal==0 ) {
        EvtWeightTotal = EvtWeightTotal + EvtTotalWeight2Lep*eff1;
        EvtWeightErrorPlus = sqrt(pow(EvtWeightErrorPlus,2)
                                  + pow(EvtTotalWeight2Lep*err_eff1,2)
                                  + pow(EvtTotalWeight2LepErrorPlus*eff1,2));
        EvtWeightErrorMinus = sqrt(pow(EvtWeightErrorMinus,2)
                                  + pow(EvtTotalWeight2Lep*err_eff1,2)
                                  + pow(EvtTotalWeight2LepErrorMinus*eff1,2));
    }
    */
//  Fraction of dilep+track events that pass the optimization cuts
    for ( int i=0; i<3; i++ ) eff[i]=0.;
    float eff2, err_eff2;
    if ( Nevents_3goodlep_base==0 || Nevents_3goodlep_base==Nevents_PassAllCuts ) {
        TSeel::Bayesian_Error(Nevents_PassAllCuts,Nevents_3goodlep_base,eff);
        eff2 = eff[0];
        err_eff2 = max(fabs(eff[2]),fabs(eff[1]));
    }
    else {
        TSeel::Binomial_Error(Nevents_PassAllCuts,Nevents_3goodlep_base,eff);
        eff2 = eff[0];
        err_eff2 = (eff[2]-eff[1])/2.;
    }
    EvtWeightOptTotal3LepGoodMass = EvtWeightOptTotal;
    EvtWeightOpt3LepGoodMassErrorPlus = EvtWeightOptErrorPlus;
    EvtWeightOpt3LepGoodMassErrorMinus = EvtWeightOptErrorMinus;
    /*
    if ( _Signal==0 ) {
        EvtWeightOptTotal = EvtWeightOptTotal + EvtTotalWeight2Lep*eff2;
        EvtWeightOptErrorPlus = sqrt(pow(EvtWeightOptErrorPlus,2)
                                  + pow(EvtTotalWeight2Lep*err_eff2,2)
                                  + pow(EvtTotalWeight2LepErrorPlus*eff2,2));
        EvtWeightOptErrorMinus = sqrt(pow(EvtWeightOptErrorMinus,2)
                                  + pow(EvtTotalWeight2Lep*err_eff2,2)
                                  + pow(EvtTotalWeight2LepErrorMinus*eff2,2));
//        cout << eff2 << "+/-" << err_eff2 << " " << EvtTotalWeight2Lep << "+/-" << EvtTotalWeight2LepErrorPlus << " " << EvtWeightOptTotal3LepGoodMass << "+/-" << EvtWeightOpt3LepGoodMassErrorPlus << " " << EvtWeightOptTotal << "+/-" << EvtWeightOptErrorPlus << " " << Nevents_PassAllCuts << " " << Nevents_3goodlep_base << endl;
    }
    */
//  Fraction of dilep+track events that pass the control region cuts
    for ( int iRegion=0; iRegion<5; iRegion++ ) {
        for ( int i=0; i<3; i++ ) eff[i]=0.;
        float eff3, err_eff3;
        if ( Nevents_3goodlep_base==0 || Nevents_3goodlep_base==Nevents_ControlRegion[iRegion] ) {
            TSeel::Bayesian_Error(Nevents_ControlRegion[iRegion],Nevents_3goodlep_base,eff);
            eff3 = eff[0];
            err_eff3 = max(fabs(eff[2]),fabs(eff[1]));
        }
        else {
            TSeel::Binomial_Error(Nevents_ControlRegion[iRegion],Nevents_3goodlep_base,eff);
            eff3 = eff[0];
            err_eff3 = (eff[2]-eff[1])/2.;
        }
        EvtWeightControlRegion3LepGoodMass[iRegion] = EvtWeightControlRegion[iRegion];
        EvtWeightControlRegion3LepGoodMassErrorPlus[iRegion] = EvtWeightControlRegionErrorPlus[iRegion];
        EvtWeightControlRegion3LepGoodMassErrorMinus[iRegion] = EvtWeightControlRegionErrorMinus[iRegion];
	/*
        if ( _Signal==0 ) {
            EvtWeightControlRegion[iRegion] = EvtWeightControlRegion[iRegion] + EvtTotalWeight2Lep*eff3;
            EvtWeightControlRegionErrorPlus[iRegion] =
                sqrt(pow(EvtWeightControlRegionErrorPlus[iRegion],2)
                     + pow(EvtTotalWeight2Lep*err_eff3,2)
                     + pow(EvtTotalWeight2LepErrorPlus*eff3,2));
            EvtWeightControlRegionErrorMinus[iRegion] =
                sqrt(pow(EvtWeightControlRegionErrorMinus[iRegion],2)
                     + pow(EvtTotalWeight2Lep*err_eff3,2)
                     + pow(EvtTotalWeight2LepErrorMinus*eff3,2));
        }
	*/
    }
    ofstream fout(_SumFileName);
    fout << " Trilepton Analysis eel " << endl;

    if ( _Data == 0 ) {
        AvgEvtTotalWeight2Lep                 = Nevents_PassDilepOnly>0 ? EvtTotalWeight2Lep/((float) Nevents_PassDilepOnly) : 0.;
        AvgEvtTotalWeight2LepErrorPlus        = Nevents_PassDilepOnly>0 ? EvtTotalWeight2LepErrorPlus/((float) Nevents_PassDilepOnly) : 0.;
        AvgEvtTotalWeight2LepErrorMinus       = Nevents_PassDilepOnly>0 ? EvtTotalWeight2LepErrorMinus/((float) Nevents_PassDilepOnly) : 0;
        AvgEvtWeightTotal3LepGoodMass         = Nevents_BaseGoodMass>0 ? EvtWeightTotal3LepGoodMass/((float) Nevents_BaseGoodMass) : 0.;
        AvgEvtWeight3LepGoodMassErrorPlus     = Nevents_BaseGoodMass>0 ? EvtWeight3LepGoodMassErrorPlus/((float) Nevents_BaseGoodMass) : 0.;
        AvgEvtWeight3LepGoodMassErrorMinus    = Nevents_BaseGoodMass>0 ? EvtWeight3LepGoodMassErrorMinus/((float) Nevents_BaseGoodMass) : 0.;
        AvgEvtWeightOptTotal3LepGoodMass      = Nevents_PassAllCuts>0 ? EvtWeightOptTotal3LepGoodMass/((float) Nevents_PassAllCuts) : 0.;
        AvgEvtWeightOpt3LepGoodMassErrorPlus  = Nevents_PassAllCuts>0 ? EvtWeightOpt3LepGoodMassErrorPlus/((float) Nevents_PassAllCuts) : 0.;
        AvgEvtWeightOpt3LepGoodMassErrorMinus = Nevents_PassAllCuts>0 ? EvtWeightOpt3LepGoodMassErrorMinus/((float) Nevents_PassAllCuts) : 0.;
        for ( int i=0; i<5; i++ ) {
            AvgEvtWeightControlRegion3LepGoodMass[i]           = Nevents_ControlRegion[i]>0 ? EvtWeightControlRegion3LepGoodMass[i]/((float) Nevents_ControlRegion[i]) : 0.;
            AvgEvtWeightControlRegion3LepGoodMassErrorPlus[i]  = Nevents_ControlRegion[i]>0 ? EvtWeightControlRegion3LepGoodMassErrorPlus[i]/((float) Nevents_ControlRegion[i]) : 0.;
            AvgEvtWeightControlRegion3LepGoodMassErrorMinus[i] = Nevents_ControlRegion[i]>0 ? EvtWeightControlRegion3LepGoodMassErrorMinus[i]/((float) Nevents_ControlRegion[i]) : 0.;
        }

        fout << "Nevents (cosmics)                  = " << nEvtCosmic << endl;
        fout << "Nevents (trilep hepg)              = " << Nevents_3goodlep_hepg << endl;
        fout << "Nevents (triele hepg)              = " << Nevents_3goodele_hepg << endl;
        fout << "Nevents (triele, base)             = " << Nevents_3goodele_base << endl;
        fout << "Nevents (diele+track, matched)     = " << Nevents_3goodlep_matched << endl;
        fout << endl;
        fout << "Nevents (total MC)                 = " << Nevents_Total << endl;
        fout << "Nevents (diele only)               = " << Nevents_PassDilepOnly << endl;
        fout << "Nevents (diele+track, base)        = " << Nevents_3goodlep_base << endl;
        fout << "Nevents (diele+track, base+mass)   = " << Nevents_BaseGoodMass << endl;
        fout << "Nevents (base+mass+opt)            = " << Nevents_PassAllCuts << endl;
        fout << "Nevents (base+mass+c1 region)      = " << Nevents_ControlRegion[0] << endl;
        fout << "Nevents (base+mass+c2 region)      = " << Nevents_ControlRegion[1] << endl;
        fout << "Nevents (base+mass+c3 region)      = " << Nevents_ControlRegion[2] << endl;
        fout << "Nevents (base+mass+c4 region)      = " << Nevents_ControlRegion[3] << endl;
        fout << "Nevents (opt+HiPt_veto region)     = " << Nevents_ControlRegion[4] << endl;
        fout << endl;
        fout << "AvgEvtWgt (diele only)             = " << AvgEvtTotalWeight2Lep << " + " << AvgEvtTotalWeight2LepErrorPlus << " - " << AvgEvtTotalWeight2LepErrorMinus << endl;
        fout << "AvgEvtWgt (diele+track, base+mass) = " << AvgEvtWeightTotal3LepGoodMass << " + " << AvgEvtWeight3LepGoodMassErrorPlus << " - " << AvgEvtWeight3LepGoodMassErrorMinus << endl;
        fout << "AvgEvtWgt (base+mass+opt)          = " << AvgEvtWeightOptTotal3LepGoodMass << " + " << AvgEvtWeightOpt3LepGoodMassErrorPlus << " - " << AvgEvtWeightOpt3LepGoodMassErrorMinus << endl;
        fout << "AvgEvtWgt (base+mass+c1 region)    = " << AvgEvtWeightControlRegion3LepGoodMass[0] << " + " << AvgEvtWeightControlRegion3LepGoodMassErrorPlus[0] << " - " << AvgEvtWeightControlRegion3LepGoodMassErrorMinus[0] << endl;
        fout << "AvgEvtWgt (base+mass+c2 region)    = " << AvgEvtWeightControlRegion3LepGoodMass[1] << " + " << AvgEvtWeightControlRegion3LepGoodMassErrorPlus[1] << " - " << AvgEvtWeightControlRegion3LepGoodMassErrorMinus[1] << endl;
        fout << "AvgEvtWgt (base+mass+c3 region)    = " << AvgEvtWeightControlRegion3LepGoodMass[2] << " + " << AvgEvtWeightControlRegion3LepGoodMassErrorPlus[2] << " - " << AvgEvtWeightControlRegion3LepGoodMassErrorMinus[2] << endl;
        fout << "AvgEvtWgt (base+mass+c4 region)    = " << AvgEvtWeightControlRegion3LepGoodMass[3] << " + " << AvgEvtWeightControlRegion3LepGoodMassErrorPlus[3] << " - " << AvgEvtWeightControlRegion3LepGoodMassErrorMinus[3] << endl;
        fout << "AvgEvtWgt (opt+HiPt_veto region)   = " << AvgEvtWeightControlRegion3LepGoodMass[4] << " + " << AvgEvtWeightControlRegion3LepGoodMassErrorPlus[4] << " - " << AvgEvtWeightControlRegion3LepGoodMassErrorMinus[4] << endl;
        fout << endl;
        fout << "TotEvtWgt (dilep only or dilep+track):" << endl;
        fout << "TotEvtWgt (dilep only)             = " << EvtTotalWeight2Lep << " + " << EvtTotalWeight2LepErrorPlus << " - " << EvtTotalWeight2LepErrorMinus << endl;
        fout << "TotEvtWgt (diele+track, base+mass) = " << EvtWeightTotal3LepGoodMass << " + " << EvtWeight3LepGoodMassErrorPlus << " - " << EvtWeight3LepGoodMassErrorMinus << endl;
        fout << "TotEvtWgt (base+mass+opt)          = " << EvtWeightOptTotal3LepGoodMass << " + " << EvtWeightOpt3LepGoodMassErrorPlus << " - " << EvtWeightOpt3LepGoodMassErrorMinus << endl;
        fout << "TotEvtWgt (base+mass+c1 region)    = " << EvtWeightControlRegion3LepGoodMass[0] << " + " << EvtWeightControlRegion3LepGoodMassErrorPlus[0] << " - " << EvtWeightControlRegion3LepGoodMassErrorMinus[0] << endl;
        fout << "TotEvtWgt (base+mass+c2 region)    = " << EvtWeightControlRegion3LepGoodMass[1] << " + " << EvtWeightControlRegion3LepGoodMassErrorPlus[1] << " - " << EvtWeightControlRegion3LepGoodMassErrorMinus[1] << endl;
        fout << "TotEvtWgt (base+mass+c3 region)    = " << EvtWeightControlRegion3LepGoodMass[2] << " + " << EvtWeightControlRegion3LepGoodMassErrorPlus[2] << " - " << EvtWeightControlRegion3LepGoodMassErrorMinus[2] << endl;
        fout << "TotEvtWgt (base+mass+c4 region)    = " << EvtWeightControlRegion3LepGoodMass[3] << " + " << EvtWeightControlRegion3LepGoodMassErrorPlus[3] << " - " << EvtWeightControlRegion3LepGoodMassErrorMinus[3] << endl;
        fout << "TotEvtWgt (opt+HiPt_veto region)   = " << EvtWeightControlRegion3LepGoodMass[4] << " + " << EvtWeightControlRegion3LepGoodMassErrorPlus[4] << " - " << EvtWeightControlRegion3LepGoodMassErrorMinus[4] << endl;
//        fout << "Frac. of ctrl (ctrl/base+mass)    = " << eff3 << " +/- " << err_eff3 << endl;
//        fout << "Total Wgt Rel. Err. (base+mass)   = " << (EvtWeightTotal>0. ? (EvtWeightErrorPlus+EvtWeightErrorMinus)/2./EvtWeightTotal : 0.) << endl;
//        fout << "Total Wgt Rel. Err. (opt)         = " << (EvtWeightOptTotal>0. ? (EvtWeightOptErrorPlus+EvtWeightOptErrorMinus)/2./EvtWeightOptTotal : 0.) << endl;
//        fout << "Total Wgt Rel. Err. (ctrl region) = " << (EvtWeightControlRegion>0. ? (EvtWeightControlRegionErrorPlus+EvtWeightControlRegionErrorMinus)/2./EvtWeightControlRegion : 0.) << endl;
        fout << endl;
        fout << "Total Event Weights:" << endl;
        fout << "Total Event Weight (base+mass)     = " << EvtWeightTotal << " + " << EvtWeightErrorPlus << " - " << EvtWeightErrorMinus << endl;
        fout << "Total Event Weight (opt)           = " << EvtWeightOptTotal << " + " << EvtWeightOptErrorPlus << " - " << EvtWeightOptErrorMinus << endl;
        fout << "Total Event Weight (c1 region)     = " << EvtWeightControlRegion[0] << " + " << EvtWeightControlRegionErrorPlus[0] << " - " << EvtWeightControlRegionErrorMinus[0] << endl;
        fout << "Total Event Weight (c2 region)     = " << EvtWeightControlRegion[1] << " + " << EvtWeightControlRegionErrorPlus[1] << " - " << EvtWeightControlRegionErrorMinus[1] << endl;
        fout << "Total Event Weight (c3 region)     = " << EvtWeightControlRegion[2] << " + " << EvtWeightControlRegionErrorPlus[2] << " - " << EvtWeightControlRegionErrorMinus[2] << endl;
        fout << "Total Event Weight (c4 region)     = " << EvtWeightControlRegion[3] << " + " << EvtWeightControlRegionErrorPlus[3] << " - " << EvtWeightControlRegionErrorMinus[3] << endl;
        fout << "Total Event Weight (opt+HiPt_Veto) = " << EvtWeightControlRegion[4] << " + " << EvtWeightControlRegionErrorPlus[4] << " - " << EvtWeightControlRegionErrorMinus[4] << endl;
        fout << endl;
        fout << "Frac. of base (base+mass/total)    = " << eff0 << " +/- " << err_eff0 << endl;
        fout << "Frac. of opt (opt/base+mass)       = " << eff2 << " +/- " << err_eff2 << endl;
        fout << "Nevents rejected due to 3rd ele    = " << NeventsGoodEleReject << endl;
	fout << "Optimization Scheme used is        = " << _OptScheme << endl;
    }
    else {
        fout << "Nevents (cosmics)                 = " << nEvtCosmic << endl;
        fout << "Nevents (total data)              = " << Nevents_Total << endl;
        fout << "Nevents (base+mass+c1 region)     = " << Nevents_ControlRegion[0] << endl;
        fout << "Nevents (base+mass+c2 region)     = " << Nevents_ControlRegion[1] << endl;
        fout << "Nevents (base+mass+c3 region)     = " << Nevents_ControlRegion[2] << endl;
        fout << "Nevents (base+mass+c4 region)     = " << Nevents_ControlRegion[3] << endl;
	fout << "Nevents (fake events)             = " << Nevents_FakeIncremented << endl;
	fout << "TotalEventWeight Fakes            = " << fakeEvtWeightTotal << " + " << fakeEvtWeightErrorPlus << " + " << fakeEvtWeightErrorMinus << endl;
	fout << "Nevents (fake opt events)         = " << Nevents_FakeIncrementedOpt << endl;
	fout << "TotalEventWeight Fakes (mass+opt) = " << fakeEvtWeightOptTotal << endl;
        fout << "Nevents (base+mass+opt)           = " << Nevents_PassAllCuts << endl;
        fout << "Nevents (opt+HiPt_veto)           = " << Nevents_ControlRegion[4] << endl;
        fout << endl;
        fout << "Nevents rejected due to 3rd ele    = " << NeventsGoodEleReject << endl;
	fout << "Optimization Scheme used is        = " << _OptScheme << endl;
    }

     GetFolder()->Write();
    _HstFile->Write();
    _HstFile->Close();

    cout<<"Nevents with 2 or more mu = "<<nev_2mu<<endl;
    cout<<"Nevents with 2 or more ele = "<<nev_2ele<<endl;
    cout<<"Nevents with 3 or more obj = "<<nev_3obj<<endl;
    cout<<"e  +  m   "<<nev_3e<<"  "<<nev_3m<<endl;
    cout<<"NeV sametrack = "<<nev_sametrack<<endl;
    cout<<"nev_eem="<<nev_eem<<" nev_mmm="<<nev_mmm<<endl;
    printf("----- end job: ---- %s\n",GetName());
    return 0;
}
//______________________________________________________________________________________________

bool TSeel::good_run(int run)
{
        // Check cached results for speed
    if (run==_prevRun) {
        return _prevGood;
    }
  
        //Binary search on the good run STL vector
    if (_initialized) {
        _prevRun = run;
        _prevGood = std::binary_search(_good.begin(), _good.end(), run);
        return _prevGood;
    }
    else {
        std::cout << "There were problems making good run list.  Good run set false " << std::endl;
        return false;
    }
}

unsigned TSeel::getbits(unsigned x, int p, int n)
{
    return (x >> (p+1-n)) & ~(~0 << n);
}

float TSeel::v2phi(float x, float y)
{
    if ( y>=0. )
        return atan2(y,x);
    else
        return 2.*M_PI + atan2(y,x);
}

//______________________________________________________________________________
float TSeel::v2eta(float px, float py, float pz)
{
    float pt = sqrt(pow(px,2)+pow(py,2));
    float p = sqrt(pow(pt,2)+pow(pz,2));
    if ( pz>0 )
        return log((p+pz)/pt);
    else
        return log(pt/(p-pz));
}

float TSeel::delta_phi(float phi1, float phi2)
{
    phi1 = TVector2::Phi_0_2pi(phi1);
    phi2 = TVector2::Phi_0_2pi(phi2);
    float dphi=fabs(phi1-phi2);
    if ( dphi>M_PI ) dphi = 2.*M_PI - dphi;
    return dphi;
}

float TSeel::delta_eta(float eta1, float eta2)
{
    return fabs(eta1-eta2);
}

float TSeel::delta_r(float eta1, float phi1, float eta2, float phi2)
{
    float dphi=delta_phi(phi1,phi2);
    return sqrt(pow((eta1-eta2),2)+pow(dphi,2));
}

bool TSeel::isConversion(TStnElectron *ele)
{
  float sepMax = 0.5;//0.2;
    //deltaLam is actually delta(cot[theta])
  float deltaLamMax = 0.1;//0.04;
    double r1, r2, x1, x2, y1, y2;
    double ff, dd, sep, deltaLam, convCharge;
    double xconv, yconv, rconv, dx, dy, sqr;
    TStnTrack* t1 = fTrackBlock->Track(ele->TrackNumber());

    bool susy_ele = false;
    float dr_min = 999.;
    TGenParticle *hepg_e;
    for (int i=0; i<fGenpBlock->NParticles(); i++) {
        TGenParticle *p = fGenpBlock->Particle(i);
        if ( p->GetStatusCode()==1 && p->IsElectron() ) {
            float phi1 = TVector2::Phi_0_2pi(p->Phi());
            float phi2 = TVector2::Phi_0_2pi(t1->Phi0());       // not beam constrained
            float dr = delta_r(p->Eta(),phi1,t1->Eta(),phi2);  // not beam constrained
            if ( dr<dr_min ) {
                dr_min = dr;
                hepg_e = p;
            }
        }
    }
    if ( dr_min<0.05 && fabs(hepg_e->P()-ele->EmE()*(ele->Etcor()>0. ? ele->Etcor(): 1.))/hepg_e->P()<0.1 ) {
        int im = hepg_e->GetFirstMother();
        TGenParticle *pm = fGenpBlock->Particle(im);
        if ( abs(pm->GetPdgCode())==1000023 || abs(pm->GetPdgCode())==1000024 )
            susy_ele = true;
    }
    if ( susy_ele ) h.EleEt_susy_ele->Fill(ele->EmEt()*(ele->Etcor()>0. ? ele->Etcor(): 1.));

    bool found_pair = false;
//    if ( ! pass_track_cut(t1) ) return false;
    for ( int j=0; j<fTrackBlock->NTracks(); j++ ) {
        if( j == t1->Number() ) continue;
        TStnTrack* t2 = fTrackBlock->Track(j);
//        if ( ! pass_track_cut(t2) ) continue;
		//Lam0 is defined as cot(theta) of track - check line 121 and 39 of TStnTrack.hh
        deltaLam = t1->Lam0() - t2->Lam0();

		r1 = 1.0/TMath::Abs(2*t1->C0());
		ff = t1->Charge()*r1+t1->D0();
		x1 = -ff*sin(t1->Phi0());
		y1 =  ff*cos(t1->Phi0());
		r2 = 1.0/TMath::Abs(2*t2->C0());
		ff = t2->Charge()*r2+t2->D0();
		x2 = -ff*sin(t2->Phi0());
		y2 =  ff*cos(t2->Phi0());
        dd  = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
	sep = dd-r1-r2;
        dx  = x2-x1;
        dy  = y2-y1;
        sqr = sqrt(dx*dx+dy*dy);

        xconv = x1+dx/sqr*(r1+fabs(sep)/2);
        yconv = y1+dy/sqr*(r1+fabs(sep)/2);
        rconv = sqrt(xconv*xconv+yconv*yconv);

        convCharge = t1->Charge()*t2->Charge();

        float dkln = 1000.;
//        dkln = isvivekFit(t1,t2);

        h.conv_R->Fill(rconv);
        h.conv_Sxy_dCotTheta->Fill(sep,deltaLam);
        if ( susy_ele && ele->EmEt()*(ele->Etcor()>0. ? ele->Etcor(): 1.)>4. ) {
            h.conv_R_susy_ele->Fill(rconv);
            h.conv_Sxy_dCotTheta_susy_ele->Fill(sep,deltaLam);
        }
        h.conv_Lxy->Fill(dkln);
        if ( fabs(deltaLam)<deltaLamMax
             && fabs(sep)<sepMax
             && convCharge<0.
//             && rconv>2. && rconv<42.
             && dkln>0. ) {
            found_pair = true;
            h.conv_trk2_pt->Fill(t2->Pt());
            if ( susy_ele ) {
                h.EleEt_Conv_susy_ele->Fill(ele->EmEt()*(ele->Etcor()>0. ? ele->Etcor(): 1.));
                h.conv_trk2_pt_susy_ele->Fill(t2->Pt());
            }
            break;
        }
    }

    if ( found_pair )
        return true;
    else
        return false;
}

double TSeel::CorrectedD0(TStnTrack *trk)
{
    double corrected_d0;
    corrected_d0 = trk->D0();
    if ( fBeamPos->Status() >= 0 )
        corrected_d0 = corrected_d0-fBeamPos->D0()*sin(fBeamPos->Phi0()-trk->Phi0());

    return corrected_d0;
}

bool TSeel::matched_hepg_ele(TStnElectron *el)
{
    if ( _JES_IMode == 1 ) return false;
    double dr_min = 999.;
    TGenParticle *hepg_e;
    TStnTrack *trk = fTrackBlock->Track(el->TrackNumber());
    for (int i=0; i<fGenpBlock->NParticles(); i++) {
        TGenParticle* p = fGenpBlock->Particle(i);
        if ( p->GetStatusCode()==1 && p->IsElectron() ) {
            float phi1 = TVector2::Phi_0_2pi(p->Phi());
            float phi2 = TVector2::Phi_0_2pi(trk->Phi0());       // not beam constrained
            double dr = delta_r(p->Eta(),phi1,trk->Eta(),phi2);  // not beam constrained
            if ( dr<dr_min ) {
                dr_min = dr;
                hepg_e = p;
            }
        }
    }
    if ( dr_min < 999. ) {
        h.hepg_ele_matching_dr->Fill(dr_min);
        h.hepg_ele_matching_dp->Fill((hepg_e->P()-el->EmE()*(el->Etcor()>0. ? el->Etcor(): 1.))/hepg_e->P());
        if ( dr_min<0.05 && fabs(hepg_e->P()-el->EmE()*(el->Etcor()>0. ? el->Etcor(): 1.))/hepg_e->P()<0.3 )
            return true;
        else
            return false;
    }
    else
        return false;
}

bool TSeel::matched_hepg_mu(TStnMuon *mu)
{
    if ( _JES_IMode == 1 ) return false;
    int it = mu->TrackNumber();
    TStnTrack* trk = fTrackBlock->Track(it);
//    TLorentzVector *mom = new TLorentzVector();
//    int status = trk->GetBcMomentum(mom);
    double dr_min = 999.;
    TGenParticle *hepg_mu;
    for (int i=0; i<fGenpBlock->NParticles(); i++) {
        TGenParticle* p = fGenpBlock->Particle(i);
        if ( p->GetStatusCode()==1 && p->IsMuon() ) {
            float phi1 = TVector2::Phi_0_2pi(p->Phi());
            float phi2 = TVector2::Phi_0_2pi(trk->Phi0());
            double dr = delta_r(p->Eta(),phi1,trk->Eta(),phi2);
            if ( dr<dr_min ) {
                dr_min = dr;
                hepg_mu = p;
            }
        }
    }
    bool matched;
    if ( dr_min < 999. ) {
        h.hepg_mu_matching_dr->Fill(dr_min);
        h.hepg_mu_matching_dp->Fill((hepg_mu->P()-trk->Pt()*cosh(trk->Eta()))/hepg_mu->P());
        if ( dr_min<0.05 && fabs(hepg_mu->P()-trk->Pt()*cosh(trk->Eta()))/hepg_mu->P()<0.3 )
            matched = true;
        else
            matched = false;
    }
    else
        matched = false;
//    delete mom;
    return matched;
}

unsigned TSeel::electron_quality(TStnElectron *el)
{
    unsigned quality = 0x000000;
    int it = el->TrackNumber();
    TStnTrack *trk;
//    TLorentzVector *mom = new TLorentzVector();
    if ( it >= 0 ) {
        trk = fTrackBlock->Track(it);
//        int status = trk->GetBcMomentum(mom);
    }
    float eChi2Strip = 0.1792*pow(2.11,double(log(el->EmE()*(el->Etcor()>0. ? el->Etcor(): 1.))))*el->Chi2Strip();
//  bit 20 - plug
//  bit 21 - central
//  bit 22 - loose
//  bit 23 - tight
    if ( el->IsCentral() ) {  // Central
//  geometry - fiducial, |Z0|<60. cm, pt>2. GeV
//  bit 0 - HadEm
//  bit 1 - EoverP (loose)
//  bit 2 - Lshr (loose)
//  bit 3 - chi2strip (loose)
//  bit 4 - qDx (loose)
//  bit 5 - Dz (loose)
//  bit 6 - NCotAxSeg
//  bit 7 - NCotStSeg
//  bit 8 - COT chi2/ndof
//  bit 9 - fiso
        quality = quality | 0x200000;
        if ( el->FidEleSmx()!=1
             || trk->Pt()<2.
             || fabs(el->Z0())>60. )
            return quality;
        if ( el->HadEm()<0.055+0.00045*el->EmE()*(el->Etcor()>0. ? el->Etcor(): 1.) ) quality = quality | 0x000001;
        if ( (el->EOverP()<2. && trk->Pt()<50.) || trk->Pt()>50. ) quality = quality | 0x000002;
        if ( el->Lshr()<0.2 ) quality = quality | 0x000004;
        if ( eChi2Strip<10. ) quality = quality | 0x000008;
        if ( el->DelX()*el->Charge()>-3. && el->DelX()*el->Charge()<1.5 ) quality = quality | 0x000010;
        if ( fabs(el->DelZ())<3. ) quality = quality | 0x000020;
        if ( it>=0 && trk->NCotAxSeg(5)>=3 ) quality = quality | 0x000040;
        if ( it>=0 && trk->NCotStSeg(5)>=2 ) quality = quality | 0x000080;
        if ( it>=0 && trk->Chi2Cot()/(trk->NCotHitsTot()-5.)<4. ) quality = quality | 0x000100;
        //if ( el->Iso1Corr()<0.1 ) quality = quality | 0x000200;
	float et = el->EmEt()*(el->Etcor()>0. ? el->Etcor(): 1.);
	float fiso = el->Iso1Corr();
	if ( (et>20 && fiso<0.1)||(et<20 && fiso*et<2) ) quality = quality | 0x000200;
        //if ( (quality & 0x00023F)==0x00023F ) quality = quality | 0x400000; // no COT
        if ( (quality & 0x0002C1)==0x0002C1 && el->Iso1Corr()<0.1 ) quality = quality | 0x400000; // loose hi-pt
        if ( (quality & 0x0002FF)==0x0002FF ) {
//            quality = quality & 0xBFFFFF;
            quality = quality | 0x800000;
        }
    }
/*
    else if ( el->IsPlug() ) {
//  bit 0 - HadEm
//  bit 1 - fiso
//  bit 2 - 3x3FitTower: 0 or 1
//  bit 3 - 3x3Chi2
//  bit 4 - 5x9u
//  bit 5 - 5x9v
//  bit 6 - nSiHits (loose)
//  bit 7 - Dr_PesPem (loose)
//  bit 8 - Z0 (loose)
        if ( fabs(el->PesEta())<1.2 || fabs(el->PesEta())>2.8 ) return quality;
        quality = quality | 0x100000;
        if ( (el->EmEt()*(el->Etcor()>0. ? el->Etcor(): 1.)<100. && el->HadEm()<0.05) ||
             (el->EmEt()*(el->Etcor()>0. ? el->Etcor(): 1.)>=100. && el->HadEm()<0.05+0.026*log(el->EmE()*(el->Etcor()>0. ? el->Etcor(): 1.)/100)) )
            quality = quality | 0x000001;
        if ( el->Iso1Corr()<0.1 ) quality = quality | 0x000002;
        if ( el->Pem3x3FitTower()!=0 ) quality = quality | 0x000004;
        if ( el->Chi2Three()<10. ) quality = quality | 0x000008;
        if ( (el->EmEt()*(el->Etcor()>0. ? el->Etcor(): 1.)<100. && el->Pes5x9(0)>0.65)
             || el->EmEt()*(el->Etcor()>0. ? el->Etcor(): 1.)>100. )
            quality = quality | 0x000010;
        if ( (el->EmEt()*(el->Etcor()>0. ? el->Etcor(): 1.)<100. && el->Pes5x9(1)>0.65)
             || el->EmEt()*(el->Etcor()>0. ? el->Etcor(): 1.)>100. )
            quality = quality | 0x000020;
        if ( it>=0 && trk->NSvxHits()>=3 ) quality = quality | 0x000040;
        if ( el->PesTrkDeltaR()<3. ) quality = quality | 0x000080;
        if ( fabs(el->Z0())<60. ) quality = quality | 0x000100; 
        if ( (quality & 0x00003F)==0x00003F ) quality = quality | 0x400000;
        if ( (quality & 0x0001FF)==0x0001FF ) {
//            quality = quality & 0xBFFFFF;
            quality = quality | 0x800000;
        }
    }
*/
//    delete mom;
    return quality;
}

unsigned TSeel::muon_quality_hipt(TStnMuon *mu)
{
    unsigned quality = 0x000000;
    int it = mu->TrackNumber();
    if ( it<0 ) return quality;
    TStnTrack *trk = fTrackBlock->Track(it);
//    TLorentzVector *mom = new TLorentzVector();
//    int status = trk->GetBcMomentum(mom);
    float d0_corr = CorrectedD0(trk);
    float rho = track_exit_radius(trk->Eta(), trk->Z0());
    float chi2 = trk->Chi2Cot()/(trk->NCotHitsTot()-4.999999);
    if ( mu->TrackPt()<2. || fabs(trk->Z0())>60.) return quality;
    //if ( mu->HasCmxStub() && rho<140. ) return quality;
    //if ( mu->TrackPt()<2. ) return quality;
    float mup = mu->TrackPt()*trk->Eta();
    //replacing "mu->Momentum()->P()" with "mup"
//  bit 0 - Em energy
//  bit 1 - Had energy
//  bit 2 - EM+Had energy
//  bit 3 - Dx
//  bit 4 - d0
//  bit 5 - nAxSegs
//  bit 6 - nStSegs
//  bit 7 - chisq/ndof
//  bit 8 - chisq match
//  bit 9 - fiso
//  bit 17 - CMIO
//  bit 18 - CMU
//  bit 19 - CMP
//  bit 20 - CMUP
//  bit 21 - CMX
//  bit 22 - loose
//  bit 23 - tight

// 0:CMU; 1:CMP; 2:CMX; 3:BMU; 4:CSP; 5:CSX; 6:BSU; 7:TSU

    if ( mu->HasCmuStub() && !mu->HasCmpStub() ) {
        if ( mu->EmEnergy()<2.+max(0.,0.0115*(mup-100.)) )
            quality = quality | 0x000001;
        if ( mu->TrackPt()>20. ) {
            if ( mu->HadEnergy()<6.+max(0.,0.028*(mup-100.)) )
                quality = quality | 0x000002;
        }
        else {
            if ( mu->HadEnergy()<3.5+mu->TrackPt()/8. )
                quality = quality | 0x000002;
        }
        if ( mu->EmEnergy()+mu->HadEnergy()>0.1 )
            quality = quality | 0x000004;
        //if ( fabs(mu->CmuDelX())<3. ) 
	quality = quality | 0x000008;
        //if ( mu->TrackPt()>20. ) {
        //    if ( (it>=0 && trk->NSvxHits()>0 && fabs(d0_corr)<0.02) || (it>=0 && trk->NSvxHits()<=0 && fabs(d0_corr)<0.2) ) quality = quality | 0x000010;
        //}
        //else
	if ( fabs(d0_corr)<0.2 ) quality = quality | 0x000010;
        if ( it>=0 && trk->NCotAxSeg(5)>=3 ) quality = quality | 0x000020;
        if ( it>=0 && trk->NCotStSeg(5)>=2 ) quality = quality | 0x000040;
        //if ( chi2<4. ) 
	quality = quality | 0x000080;
        //if ( mu->CmuChi2Link()<9. ) 
	quality = quality | 0x000100;
        if ( (mu->Iso()/mu->TrackPt()<0.1) )
	  quality = quality | 0x000200;
        quality = quality | 0x040000;
    }
    else if ( mu->HasCmpStub() && !mu->HasCmuStub() ) {
        float CmpChi2 = mu->CmpChi2Link();
        if ( CmpChi2<-998. ) CmpChi2 = cmpChi2(mu);
        if ( mu->EmEnergy()<2.+max(0.,0.0115*(mup-100.)) )
            quality = quality | 0x000001;
        if ( mu->TrackPt()>20. ) {
            if ( mu->HadEnergy()<6.+max(0.,0.028*(mup-100.)) )
                quality = quality | 0x000002;
        }
        else {
            if ( mu->HadEnergy()<3.5+mu->TrackPt()/8. )
                quality = quality | 0x000002;
        }
        if ( mu->EmEnergy()+mu->HadEnergy()>0.1 )
            quality = quality | 0x000004;
        //if ( fabs(mu->CmpDelX())<5. ) 
	quality = quality | 0x000008;
	if ( fabs(d0_corr)<0.2 ) quality = quality | 0x000010;
        if ( it>=0 && trk->NCotAxSeg(5)>=3 ) quality = quality | 0x000020;
        if ( it>=0 && trk->NCotStSeg(5)>=2 ) quality = quality | 0x000040;
        //if ( chi2<4. ) 
	quality = quality | 0x000080;
        //if ( CmpChi2<9. ) 
	quality = quality | 0x000100;
        if ( mu->Iso()/mu->TrackPt()<0.1 ) quality = quality | 0x000200;
        quality = quality | 0x080000;
    }
    else if ( mu->HasCmuStub() && mu->HasCmpStub() ) {
        if ( mu->EmEnergy()<2.+max(0.,0.0115*(mup-100.)) )
            quality = quality | 0x000001;
        if ( mu->TrackPt()>20. ) {
            if ( mu->HadEnergy()<6.+max(0.,0.028*(mup-100.)) )
                quality = quality | 0x000002;
        }
        else {
            if ( mu->HadEnergy()<3.5+mu->TrackPt()/8. )
                quality = quality | 0x000002;
        }
        if ( mu->EmEnergy()+mu->HadEnergy()>0.1 )
            quality = quality | 0x000004;
        //if ( fabs(mu->CmuDelX())<3. && fabs(mu->CmpDelX())<5. ) 
	  quality = quality | 0x000008;
	if ( fabs(d0_corr)<0.2 ) quality = quality | 0x000010;
        if ( it>=0 && trk->NCotAxSeg(5)>=3 ) quality = quality | 0x000020;
        if ( it>=0 && trk->NCotStSeg(5)>=2 ) quality = quality | 0x000040;
        //if ( chi2<4. ) 
	quality = quality | 0x000080;
        //if ( mu->CmuChi2Link()<9. ) 
	quality = quality | 0x000100;
        if ( mu->Iso()/mu->TrackPt()<0.1 ) quality = quality | 0x000200;
        quality = quality | 0x100000;
    }
    else if ( mu->HasCmxStub() ) {
        float CmxChi2 = mu->CmxChi2Link();
        if ( CmxChi2<-998. ) CmxChi2 = cmxChi2(mu);
        if ( mu->EmEnergy()<2.+max(0.,0.0115*(mup-100.)) )
            quality = quality | 0x000001;
        if ( mu->TrackPt()>20. ) {
            if ( mu->HadEnergy()<6.+max(0.,0.028*(mup-100.)) )
                quality = quality | 0x000002;
        }
        else {
            if ( mu->HadEnergy()<3.5+mu->TrackPt()/8. )
                quality = quality | 0x000002;
        }
        if ( mu->EmEnergy()+mu->HadEnergy()>0.1 )
            quality = quality | 0x000004;
        //if ( fabs(mu->CmxDelX())<6. ) 
	quality = quality | 0x000008;
	if ( fabs(d0_corr)<0.2 ) quality = quality | 0x000010;
        if ( it>=0 && trk->NCotAxSeg(5)>=3 ) quality = quality | 0x000020;
        if ( it>=0 && trk->NCotStSeg(5)>=2 ) quality = quality | 0x000040;
        //if ( chi2<4. ) 
	quality = quality | 0x000080;
        //if ( CmxChi2<9. ) 
	quality = quality | 0x000100;
        if ( mu->Iso()/mu->TrackPt() <0.1 ) quality = quality | 0x000200;  // fiso
        quality = quality | 0x200000;
    }
    else if ( mu->IsStubless() ) {
        if ( mu->EmEnergy()<2.+max(0.,0.0115*(mup-100.)) )
              quality = quality | 0x000001;
        if ( mu->TrackPt()>20. ) {
            if ( mu->HadEnergy()<6.+max(0.,0.028*(mup-100.)) )
                quality = quality | 0x000002;
        }
        else {
            if ( mu->HadEnergy()<3.5+mu->TrackPt()/8. )
                quality = quality | 0x000002;
        }
        if ( mu->EmEnergy()+mu->HadEnergy()>0.1 )
            quality = quality | 0x000004;
        quality = quality | 0x000008;  // no dx
	if ( fabs(d0_corr)<0.2 ) quality = quality | 0x000010;
        if ( it>=0 && trk->NCotAxSeg(5)>=3 ) quality = quality | 0x000020;
        if ( it>=0 && trk->NCotStSeg(5)>=3 ) quality = quality | 0x000040;
        if ( 1/*chi2<4.*/ ) quality = quality | 0x000080;
        quality = quality | 0x000100;  // no chi2 link
        if ( mu->Iso()/mu->TrackPt()<0.1 ) quality = quality | 0x000200;  // fiso
        quality = quality | 0x020000;
    }

    if ( (quality & 0x00031F)==0x00031F ) quality = quality | 0x400000;  // no COT
    if ( (quality & 0x00037F)==0x00037F ) quality = quality | 0x800000;

//    delete mom;
    return quality;
}


unsigned TSeel::muon_quality(TStnMuon *mu)
{
    unsigned quality = 0x000000;
    int it = mu->TrackNumber();
    if ( it<0 ) return quality;
    TStnTrack *trk = fTrackBlock->Track(it);
//    TLorentzVector *mom = new TLorentzVector();
//    int status = trk->GetBcMomentum(mom);
    float d0_corr = CorrectedD0(trk);
    float rho = track_exit_radius(trk->Eta(), trk->Z0());
    float chi2 = trk->Chi2Cot()/(trk->NCotHitsTot()-4.999999);
    if ( mu->TrackPt()<2. || fabs(trk->Z0())>60.) return quality;
    if ( mu->HasCmxStub() && rho<140. ) return quality;
    //if ( mu->TrackPt()<2. ) return quality;
    float mup = mu->TrackPt()*trk->Eta();
    //replacing "mu->Momentum()->P()" with "mup"
//  bit 0 - Em energy
//  bit 1 - Had energy
//  bit 2 - EM+Had energy
//  bit 3 - Dx
//  bit 4 - d0
//  bit 5 - nAxSegs
//  bit 6 - nStSegs
//  bit 7 - chisq/ndof
//  bit 8 - chisq match
//  bit 9 - fiso
//  bit 17 - CMIO
//  bit 18 - CMU
//  bit 19 - CMP
//  bit 20 - CMUP
//  bit 21 - CMX
//  bit 22 - loose
//  bit 23 - tight

// 0:CMU; 1:CMP; 2:CMX; 3:BMU; 4:CSP; 5:CSX; 6:BSU; 7:TSU

    if ( mu->HasCmuStub() && !mu->HasCmpStub() ) {
        if ( mu->EmEnergy()<2.+max(0.,0.0115*(mup-100.)) )
            quality = quality | 0x000001;
        if ( mu->TrackPt()>20. ) {
            if ( mu->HadEnergy()<6.+max(0.,0.028*(mup-100.)) )
                quality = quality | 0x000002;
        }
        else {
            if ( mu->HadEnergy()<3.5+mu->TrackPt()/8. )
                quality = quality | 0x000002;
        }
        if ( mu->EmEnergy()+mu->HadEnergy()>0.1 )
            quality = quality | 0x000004;
        if ( fabs(mu->CmuDelX())<3. ) quality = quality | 0x000008;
        //if ( mu->TrackPt()>20. ) {
        //    if ( (it>=0 && trk->NSvxHits()>0 && fabs(d0_corr)<0.02) || (it>=0 && trk->NSvxHits()<=0 && fabs(d0_corr)<0.2) ) quality = quality | 0x000010;
        //}
        //else
	if ( fabs(d0_corr)<0.2 ) quality = quality | 0x000010;
        if ( it>=0 && trk->NCotAxSeg(5)>=3 ) quality = quality | 0x000020;
        if ( it>=0 && trk->NCotStSeg(5)>=2 ) quality = quality | 0x000040;
        if ( chi2<4. ) quality = quality | 0x000080;
        if ( mu->CmuChi2Link()<9. ) quality = quality | 0x000100;
        if ( (mu->TrackPt()>20 && mu->Iso()/mu->TrackPt()<0.1) ||
	     (mu->TrackPt()<20 && mu->Iso()<2) ) quality = quality | 0x000200;
        quality = quality | 0x040000;
    }
    else if ( mu->HasCmpStub() && !mu->HasCmuStub() ) {
        float CmpChi2 = mu->CmpChi2Link();
        if ( CmpChi2<-998. ) CmpChi2 = cmpChi2(mu);
        if ( mu->EmEnergy()<2.+max(0.,0.0115*(mup-100.)) )
            quality = quality | 0x000001;
        if ( mu->TrackPt()>20. ) {
            if ( mu->HadEnergy()<6.+max(0.,0.028*(mup-100.)) )
                quality = quality | 0x000002;
        }
        else {
            if ( mu->HadEnergy()<3.5+mu->TrackPt()/8. )
                quality = quality | 0x000002;
        }
        if ( mu->EmEnergy()+mu->HadEnergy()>0.1 )
            quality = quality | 0x000004;
        if ( fabs(mu->CmpDelX())<5. ) quality = quality | 0x000008;
        //if ( mu->TrackPt()>20. ) {
        //    if ( (it>=0 && trk->NSvxHits()>0 && fabs(d0_corr)<0.02) || (it>=0 && trk->NSvxHits()<=0 && fabs(d0_corr)<0.2) ) quality = quality | 0x000010;
        //}
        //else
	if ( fabs(d0_corr)<0.2 ) quality = quality | 0x000010;
        if ( it>=0 && trk->NCotAxSeg(5)>=3 ) quality = quality | 0x000020;
        if ( it>=0 && trk->NCotStSeg(5)>=2 ) quality = quality | 0x000040;
        if ( chi2<4. ) quality = quality | 0x000080;
        if ( CmpChi2<9. ) quality = quality | 0x000100;
        //if ( mu->Iso()/mu->TrackPt()<0.1 ) quality = quality | 0x000200;
        if ( (mu->TrackPt()>20 && mu->Iso()/mu->TrackPt()<0.1) ||
	     (mu->TrackPt()<20 && mu->Iso()<2) ) quality = quality | 0x000200;
        quality = quality | 0x080000;
    }
    else if ( mu->HasCmuStub() && mu->HasCmpStub() ) {
        if ( mu->EmEnergy()<2.+max(0.,0.0115*(mup-100.)) )
            quality = quality | 0x000001;
        if ( mu->TrackPt()>20. ) {
            if ( mu->HadEnergy()<6.+max(0.,0.028*(mup-100.)) )
                quality = quality | 0x000002;
        }
        else {
            if ( mu->HadEnergy()<3.5+mu->TrackPt()/8. )
                quality = quality | 0x000002;
        }
        if ( mu->EmEnergy()+mu->HadEnergy()>0.1 )
            quality = quality | 0x000004;
        if ( fabs(mu->CmuDelX())<3. /*&& fabs(mu->CmpDelX())<5.*/ ) quality = quality | 0x000008;
        //if ( mu->TrackPt()>20. ) {
        //    if ( (it>=0 && trk->NSvxHits()>0 && fabs(d0_corr)<0.02) || (it>=0 && trk->NSvxHits()<=0 && fabs(d0_corr)<0.2) ) quality = quality | 0x000010;
        //}
        //else
	if ( fabs(d0_corr)<0.2 ) quality = quality | 0x000010;
        if ( it>=0 && trk->NCotAxSeg(5)>=3 ) quality = quality | 0x000020;
        if ( it>=0 && trk->NCotStSeg(5)>=2 ) quality = quality | 0x000040;
        if ( chi2<4. ) quality = quality | 0x000080;
        if ( mu->CmuChi2Link()<9. ) quality = quality | 0x000100;
        //if ( mu->Iso()/mu->TrackPt()<0.1 ) quality = quality | 0x000200;
        if ( (mu->TrackPt()>20 && mu->Iso()/mu->TrackPt()<0.1) ||
	     (mu->TrackPt()<20 && mu->Iso()<2) ) quality = quality | 0x000200;
        quality = quality | 0x100000;
    }
    else if ( mu->HasCmxStub() ) {
        float CmxChi2 = mu->CmxChi2Link();
        if ( CmxChi2<-998. ) CmxChi2 = cmxChi2(mu);
        if ( mu->EmEnergy()<2.+max(0.,0.0115*(mup-100.)) )
            quality = quality | 0x000001;
        if ( mu->TrackPt()>20. ) {
            if ( mu->HadEnergy()<6.+max(0.,0.028*(mup-100.)) )
                quality = quality | 0x000002;
        }
        else {
            if ( mu->HadEnergy()<3.5+mu->TrackPt()/8. )
                quality = quality | 0x000002;
        }
        if ( mu->EmEnergy()+mu->HadEnergy()>0.1 )
            quality = quality | 0x000004;
        if ( fabs(mu->CmxDelX())<6. ) quality = quality | 0x000008;
        //if ( mu->TrackPt()>20. ) {
        //    if ( (it>=0 && trk->NSvxHits()>0 && fabs(d0_corr)<0.02) || (it>=0 && trk->NSvxHits()<=0 && fabs(d0_corr)<0.2) ) quality = quality | 0x000010;
        //}
        //else
	if ( fabs(d0_corr)<0.2 ) quality = quality | 0x000010;
        if ( it>=0 && trk->NCotAxSeg(5)>=3 ) quality = quality | 0x000020;
        if ( it>=0 && trk->NCotStSeg(5)>=2 ) quality = quality | 0x000040;
        if ( chi2<4. ) quality = quality | 0x000080;
        if ( CmxChi2<9. ) quality = quality | 0x000100;
        //if ( mu->Iso()/mu->TrackPt() <0.1 ) quality = quality | 0x000200;  // fiso
        if ( (mu->TrackPt()>20 && mu->Iso()/mu->TrackPt()<0.1) ||
	     (mu->TrackPt()<20 && mu->Iso()<2) ) quality = quality | 0x000200;
        quality = quality | 0x200000;
    }
    else if ( mu->IsStubless() ) {
        if ( mu->EmEnergy()<2.+max(0.,0.0115*(mup-100.)) )
              quality = quality | 0x000001;
        if ( mu->TrackPt()>20. ) {
            if ( mu->HadEnergy()<6.+max(0.,0.028*(mup-100.)) )
                quality = quality | 0x000002;
        }
        else {
            if ( mu->HadEnergy()<3.5+mu->TrackPt()/8. )
                quality = quality | 0x000002;
        }
        if ( mu->EmEnergy()+mu->HadEnergy()>0.1 )
            quality = quality | 0x000004;
        quality = quality | 0x000008;  // no dx
        //if ( mu->TrackPt()>20. ) {
	//  if ( (it>=0 && trk->NSvxHits()>0 && fabs(d0_corr)<0.02) || (it>=0 && trk->NSvxHits()<=0 && fabs(d0_corr)<0.2) ) quality = quality | 0x000010;
        //}
        //else
	if ( fabs(d0_corr)<0.2 ) quality = quality | 0x000010;
        if ( it>=0 && trk->NCotAxSeg(5)>=3 ) quality = quality | 0x000020;
        if ( it>=0 && trk->NCotStSeg(5)>=3 ) quality = quality | 0x000040;
        if ( 1/*chi2<4.*/ ) quality = quality | 0x000080;
        quality = quality | 0x000100;  // no chi2 link
        if ( mu->Iso()/mu->TrackPt()<0.1 ) quality = quality | 0x000200;  // fiso
        quality = quality | 0x020000;
    }

    if ( (quality & 0x00031F)==0x00031F ) quality = quality | 0x400000;  // no COT
    if ( (quality & 0x00037F)==0x00037F ) quality = quality | 0x800000;

//    delete mom;
    return quality;
}


void TSeel::Fill_eID(TStnElectron *el, unsigned quality)
{
    bool pass_all;
    int it = el->TrackNumber();
    TStnTrack *trk;
//    TLorentzVector *mom = new TLorentzVector();
    if ( it >= 0 ) {
       trk = fTrackBlock->Track(it);
//        int status = trk->GetBcMomentum(mom);
    }
    float eChi2Strip = 0.1792*pow(2.11,double(log(el->EmE()*(el->Etcor()>0. ? el->Etcor(): 1.))))*el->Chi2Strip();
    float phi = TVector2::Phi_0_2pi(el->Momentum()->Phi());
    if ( el->IsCentral() ) {
        h.ce_et_raw->Fill(el->EmEt()*(el->Etcor()>0. ? el->Etcor(): 1.));
        h.ce_eta_raw->Fill(el->Momentum()->Eta());
        h.ce_phi_raw->Fill(phi);
        h.ce_id_raw[0]->Fill(el->FidEleSmx());
        h.ce_id_raw[1]->Fill(el->HadEm());
        h.ce_id_raw[2]->Fill(el->Iso1Corr());
        h.ce_id_raw[3]->Fill(el->TrackPt());
        h.ce_id_raw[4]->Fill(el->EOverP());
        if ( it>0 ) h.ce_id_raw[5]->Fill(trk->NCotAxSeg(5));
        if ( it>0 ) h.ce_id_raw[6]->Fill(trk->NCotStSeg(5));
        h.ce_id_raw[7]->Fill(el->Z0());
        h.ce_id_raw[8]->Fill(el->Lshr());
        h.ce_id_raw[9]->Fill(eChi2Strip);
        h.ce_id_raw[10]->Fill(el->DelX());
        h.ce_id_raw[11]->Fill(el->DelX()*el->Charge());
        h.ce_id_raw[12]->Fill(el->DelZ());

        pass_all = (quality & 0x0002FF)==0x0002FF;
        if ( pass_all ) {
            h.ce_et->Fill(el->EmEt()*(el->Etcor()>0. ? el->Etcor(): 1.));
            h.ce_eta->Fill(el->Momentum()->Eta());
            h.ce_phi->Fill(phi);
        }
        if (1/* pass_all || (quality & 0x001FFE)==0x001FFE*/ )
	    h.ce_id[0]->Fill(el->FidEleSmx());
        if ( pass_all || (quality & 0x0002FE)==0x0002FE )
            h.ce_id[1]->Fill(el->HadEm());
        if ( pass_all || (quality & 0x0001FF)==0x0001FF )
            h.ce_id[2]->Fill(el->Iso1Corr());
        if ( 1/*pass_all || (quality & 0x001FF7)==0x001FF7*/ )
            h.ce_id[3]->Fill(el->TrackPt());
        if ( pass_all || (quality & 0x0002FD)==0x0002FD )
            h.ce_id[4]->Fill(el->EOverP());
        if ( pass_all || (quality & 0x0002BF)==0x0002BF )
            h.ce_id[5]->Fill(trk->NCotAxSeg(5));
        if ( pass_all || (quality & 0x00027F)==0x00027F )
            h.ce_id[6]->Fill(trk->NCotStSeg(5));
        if ( 1/*pass_all || (quality & 0x001F7F)==0x001F7F*/ )
	    h.ce_id[7]->Fill(el->Z0());
        if ( pass_all || (quality & 0x0002FB)==0x0002FB )
            h.ce_id[8]->Fill(el->Lshr());
        if ( pass_all || (quality & 0x0002F7)==0x0002F7 )
            h.ce_id[9]->Fill(eChi2Strip);
        if (1/* pass_all || (quality & 0x001BFF)==0x001BFF*/ )
            h.ce_id[10]->Fill(el->DelX());
        if ( pass_all || (quality & 0x0002EF)==0x0002EF )
            h.ce_id[11]->Fill(el->DelX()*el->Charge());
        if ( pass_all || (quality & 0x0002DF)==0x0002DF )
            h.ce_id[12]->Fill(el->DelZ());
    }
    /*
    else if ( el->IsPlug() ) {
        h.pe_et_raw->Fill(el->EmEt()*(el->Etcor()>0. ? el->Etcor(): 1.));
        h.pe_eta_raw->Fill(el->Momentum()->Eta());
        h.pe_phi_raw->Fill(phi);
        h.pe_id_raw[0]->Fill(el->HadEm());
        h.pe_id_raw[1]->Fill(el->Iso1Corr());
        h.pe_id_raw[2]->Fill(el->Pem3x3FitTower());
        h.pe_id_raw[3]->Fill(el->Chi2Three());
        h.pe_id_raw[4]->Fill(el->Pes5x9(0));
        h.pe_id_raw[5]->Fill(el->Pes5x9(1));
        if ( it >= 0 ) h.pe_id_raw[6]->Fill(trk->NSvxHits());
        h.pe_id_raw[7]->Fill(el->PesTrkDeltaR());
        h.pe_id_raw[7]->Fill(el->Z0());

        pass_all = (quality & 0x0000FF)==0x0000FF;
        if ( pass_all ) {
            h.pe_et->Fill(el->EmEt()*(el->Etcor()>0. ? el->Etcor(): 1.));
            h.pe_eta->Fill(el->Momentum()->Eta());
            h.pe_phi->Fill(phi);
        }
        if ( pass_all || (quality & 0x0000FE)==0x0000FE )
            h.pe_id[0]->Fill(el->HadEm());
        if ( pass_all || (quality & 0x0000FD)==0x0000FD )
            h.pe_id[1]->Fill(el->Iso1Corr());
        if ( pass_all || (quality & 0x0000FB)==0x0000FB )
            h.pe_id[2]->Fill(el->Pem3x3FitTower());
        if ( pass_all || (quality & 0x0000F7)==0x0000F7 )
            h.pe_id[3]->Fill(el->Chi2Three());
        if ( pass_all || (quality & 0x0000EF)==0x0000EF )
            h.pe_id[4]->Fill(el->Pes5x9(0));
        if ( pass_all || (quality & 0x0000DF)==0x0000DF )
            h.pe_id[5]->Fill(el->Pes5x9(1));
        if ( pass_all || (quality & 0x0000BF)==0x0000BF )
            h.pe_id[6]->Fill(trk->NSvxHits());
        if ( pass_all || (quality & 0x00007F)==0x00007F )
            h.pe_id[7]->Fill(el->Z0());
    }
    */
//    delete mom;
}

void TSeel::Fill_mID(TStnMuon *mu, unsigned quality)
{
    int it = mu->TrackNumber();
    TStnTrack* trk = fTrackBlock->Track(it);
//    TLorentzVector *mom = new TLorentzVector();
//    int status = trk->GetBcMomentum(mom);
    float d0_corr = CorrectedD0(trk);
//    float rho = track_exit_radius(mom->Eta(), trk->Z0());
    float chi2 = trk->Chi2Cot()/(trk->NCotHitsTot()-4.9999);
    float phi = TVector2::Phi_0_2pi(mu->Momentum()->Phi());

    bool pass_all;
    if ( getbits(quality,18,1)==1 && getbits(quality,19,1)==0 ) {  // CMU
        h.cmu_pt_raw->Fill(mu->Momentum()->Pt());
        h.cmu_eta_raw->Fill(mu->Momentum()->Eta());
        h.cmu_phi_raw->Fill(phi);
        h.cmu_id_raw[0]->Fill(mu->EmEnergy());
        h.cmu_id_raw[1]->Fill(mu->HadEnergy());
        h.cmu_id_raw[2]->Fill(mu->CmuDelX());
        h.cmu_id_raw[3]->Fill(mu->Iso()/mu->TrackPt());
        h.cmu_id_raw[4]->Fill(d0_corr);
        if ( it >= 0 ) h.cmu_id_raw[5]->Fill(trk->NCotAxSeg(5));
        if ( it >= 0 ) h.cmu_id_raw[6]->Fill(trk->NCotStSeg(5));
        h.cmu_id_raw[7]->Fill(chi2);
        h.cmu_id_raw[8]->Fill(mu->CmuChi2Link());

        pass_all = getbits(quality,23,1);
        if ( pass_all ) {
            h.cmu_pt->Fill(mu->Momentum()->Pt());
            h.cmu_eta->Fill(mu->Momentum()->Eta());
            h.cmu_phi->Fill(phi);
        }
        if ( pass_all || (quality & 0x0001FE)==0x0001FE )
            h.cmu_id[0]->Fill(mu->EmEnergy());
        if ( pass_all || (quality & 0x0001FD)==0x0001FD )
            h.cmu_id[1]->Fill(mu->HadEnergy());
        if ( pass_all || (quality & 0x0001FB)==0x0001FB )
            h.cmu_id[2]->Fill(mu->CmuDelX());
        if ( pass_all || (quality & 0x0001F7)==0x0001F7 )
            h.cmu_id[3]->Fill(mu->Iso()/mu->TrackPt());
        if ( pass_all || (quality & 0x0001EF)==0x0001EF )
            h.cmu_id[4]->Fill(d0_corr);
        if ( pass_all || (quality & 0x0001DF)==0x0001DF )
            h.cmu_id[5]->Fill(trk->NCotAxSeg(5));
        if ( pass_all || (quality & 0x0001BF)==0x0001BF )
            h.cmu_id[6]->Fill(trk->NCotStSeg(5));
        if ( pass_all || (quality & 0x00017F)==0x0017F )
            h.cmu_id[7]->Fill(chi2);
        if ( pass_all || (quality & 0x0000FF)==0x0000FF )
            h.cmu_id[8]->Fill(mu->CmuChi2Link());
    }
    else if ( getbits(quality,18,1)==0 && getbits(quality,19,1)==1 ) {  // CMP
        float CmpChi2 = mu->CmpChi2Link();
        if ( CmpChi2<-998. ) CmpChi2 = cmpChi2(mu);
        h.cmp_pt_raw->Fill(mu->Momentum()->Pt());
        h.cmp_eta_raw->Fill(mu->Momentum()->Eta());
        h.cmp_phi_raw->Fill(phi);
        h.cmp_id_raw[0]->Fill(mu->EmEnergy());
        h.cmp_id_raw[1]->Fill(mu->HadEnergy());
        h.cmp_id_raw[2]->Fill(mu->CmpDelX());
        h.cmp_id_raw[3]->Fill(mu->Iso()/mu->TrackPt());
        h.cmp_id_raw[4]->Fill(d0_corr);
        if ( it >= 0 ) h.cmp_id_raw[5]->Fill(trk->NCotAxSeg(5));
        if ( it >= 0 ) h.cmp_id_raw[6]->Fill(trk->NCotStSeg(5));
        h.cmp_id_raw[7]->Fill(chi2);
        h.cmp_id_raw[8]->Fill(CmpChi2);

        pass_all = getbits(quality,23,1);
        if ( pass_all ) {
            h.cmp_pt->Fill(mu->Momentum()->Pt());
            h.cmp_eta->Fill(mu->Momentum()->Eta());
            h.cmp_phi->Fill(phi);
        }
        if ( pass_all || (quality & 0x0001FE)==0x0001FE )
            h.cmp_id[0]->Fill(mu->EmEnergy());
        if ( pass_all || (quality & 0x0001FD)==0x0001FD )
            h.cmp_id[1]->Fill(mu->HadEnergy());
        if ( pass_all || (quality & 0x0001FB)==0x0001FB )
            h.cmp_id[2]->Fill(mu->CmpDelX());
        if ( pass_all || (quality & 0x0001F7)==0x0001F7 )
            h.cmp_id[3]->Fill(mu->Iso()/mu->TrackPt());
        if ( pass_all || (quality & 0x0001EF)==0x0001EF )
            h.cmp_id[4]->Fill(d0_corr);
        if ( pass_all || (quality & 0x0001DF)==0x0001DF )
            h.cmp_id[5]->Fill(trk->NCotAxSeg(5));
        if ( pass_all || (quality & 0x0001BF)==0x0001BF )
            h.cmp_id[6]->Fill(trk->NCotStSeg(5));
        if ( pass_all || (quality & 0x00017F)==0x0017F )
            h.cmp_id[7]->Fill(chi2);
        if ( pass_all || (quality & 0x0000FF)==0x0000FF )
            h.cmp_id[8]->Fill(CmpChi2);
    }
    else if ( getbits(quality,20,1)==1 ) {  // CMUP
        h.cmup_pt_raw->Fill(mu->Momentum()->Pt());
        h.cmup_eta_raw->Fill(mu->Momentum()->Eta());
        h.cmup_phi_raw->Fill(phi);
        h.cmup_id_raw[0]->Fill(mu->EmEnergy());
        h.cmup_id_raw[1]->Fill(mu->HadEnergy());
        h.cmup_id_raw[2]->Fill(mu->CmuDelX());
        h.cmup_id_raw[3]->Fill(mu->Iso()/mu->TrackPt());
        h.cmup_id_raw[4]->Fill(d0_corr);
        if ( it >= 0 ) h.cmup_id_raw[5]->Fill(trk->NCotAxSeg(5));
        if ( it >= 0 ) h.cmup_id_raw[6]->Fill(trk->NCotStSeg(5));
        h.cmup_id_raw[7]->Fill(chi2);
        h.cmup_id_raw[8]->Fill(mu->CmuChi2Link());

        pass_all = getbits(quality,23,1);
        if ( pass_all ) {
            h.cmup_pt->Fill(mu->Momentum()->Pt());
            h.cmup_eta->Fill(mu->Momentum()->Eta());
            h.cmup_phi->Fill(phi);
        }
        if ( pass_all || (quality & 0x0001FE)==0x0001FE )
            h.cmup_id[0]->Fill(mu->EmEnergy());
        if ( pass_all || (quality & 0x0001FD)==0x0001FD )
            h.cmup_id[1]->Fill(mu->HadEnergy());
        if ( pass_all || (quality & 0x0001FB)==0x0001FB )
            h.cmup_id[2]->Fill(mu->CmuDelX());
        if ( pass_all || (quality & 0x0001F7)==0x0001F7 )
            h.cmup_id[3]->Fill(mu->Iso()/mu->TrackPt());
        if ( pass_all || (quality & 0x0001EF)==0x0001EF )
            h.cmup_id[4]->Fill(d0_corr);
        if ( pass_all || (quality & 0x0001DF)==0x0001DF )
            h.cmup_id[5]->Fill(trk->NCotAxSeg(5));
        if ( pass_all || (quality & 0x0001BF)==0x0001BF )
            h.cmup_id[6]->Fill(trk->NCotStSeg(5));
        if ( pass_all || (quality & 0x00017F)==0x0017F )
            h.cmup_id[7]->Fill(chi2);
        if ( pass_all || (quality & 0x0000FF)==0x0000FF )
            h.cmup_id[8]->Fill(mu->CmuChi2Link());
    }
    else if ( getbits(quality,21,1)==1 ) {  // CMX
        float CmxChi2 = mu->CmxChi2Link();
        if ( CmxChi2<-998. ) CmxChi2 = cmxChi2(mu);
        h.cmx_pt_raw->Fill(mu->Momentum()->Pt());
        h.cmx_eta_raw->Fill(mu->Momentum()->Eta());
        h.cmx_phi_raw->Fill(phi);
        h.cmx_id_raw[0]->Fill(mu->EmEnergy());
        h.cmx_id_raw[1]->Fill(mu->HadEnergy());
        h.cmx_id_raw[2]->Fill(mu->CmxDelX());
        h.cmx_id_raw[3]->Fill(mu->Iso()/mu->TrackPt());
        h.cmx_id_raw[4]->Fill(d0_corr);
        if ( it >= 0 ) h.cmx_id_raw[5]->Fill(trk->NCotAxSeg(5));
        if ( it >= 0 ) h.cmx_id_raw[6]->Fill(trk->NCotStSeg(5));
        h.cmx_id_raw[7]->Fill(chi2);
        h.cmx_id_raw[8]->Fill(CmxChi2);

        pass_all = getbits(quality,23,1);
        if ( pass_all ) {
            h.cmx_pt->Fill(mu->Momentum()->Pt());
            h.cmx_eta->Fill(mu->Momentum()->Eta());
            h.cmx_phi->Fill(phi);
        }
        if ( pass_all || (quality & 0x0001FE)==0x0001FE )
            h.cmx_id[0]->Fill(mu->EmEnergy());
        if ( pass_all || (quality & 0x0001FD)==0x0001FD )
            h.cmx_id[1]->Fill(mu->HadEnergy());
        if ( pass_all || (quality & 0x0001FB)==0x0001FB )
            h.cmx_id[2]->Fill(mu->CmxDelX());
        if ( pass_all || (quality & 0x0001F7)==0x0001F7 )
            h.cmx_id[3]->Fill(mu->Iso()/mu->TrackPt());
        if ( pass_all || (quality & 0x0001EF)==0x0001EF )
            h.cmx_id[4]->Fill(d0_corr);
        if ( pass_all || (quality & 0x0001DF)==0x0001DF )
            h.cmx_id[5]->Fill(trk->NCotAxSeg(5));
        if ( pass_all || (quality & 0x0001BF)==0x0001BF )
            h.cmx_id[6]->Fill(trk->NCotStSeg(5));
        if ( pass_all || (quality & 0x00017F)==0x0017F )
            h.cmx_id[7]->Fill(chi2);
        if ( pass_all || (quality & 0x0000FF)==0x0000FF )
            h.cmx_id[8]->Fill(CmxChi2);
    }

//    delete mom;
}

bool TSeel::pass_jet_cut(float Et, float DetEta)
{
    return Et>_JET_ET_CUT && fabs(DetEta)<2.5;
}

void TSeel::correct_met_jets(float et_threshold, float deta_threshold)
{
    float met_ex = missing_et.ex();
    float met_ey = missing_et.ey();
//    for ( int i=0; i<fNJet; i++ ) {
//        if ( fJets_fEt[i]*fJets_fJEScale5[i]>et_threshold
//             && fabs(fJets_fDetEta[i])<deta_threshold
//             && fJets_fIndex[i]<1000000 ) {
//            float phi = fJets_fPhi[i];
//            met_ex = met_ex - fJets_fEt[i]*cos(phi)*(fJets_fJEScale5[i]-1);
//            met_ey = met_ey - fJets_fEt[i]*sin(phi)*(fJets_fJEScale5[i]-1);
//        }
//    }
//  jzz20050205
    for ( int i=0; i<njet; i++ ) {
//  "bug", jzz20050409
//        if ( jt[i].et()*jt[i].scale5()>et_threshold  // wrong but has no impact
        if ( jt[i].et()>et_threshold   // correct
             && fabs(jt[i].deta())<deta_threshold ) {
            float phi = jt[i].phi();
            met_ex = met_ex - jt[i].et()*cos(phi)*(jt[i].scale5()-1.);
            met_ey = met_ey - jt[i].et()*sin(phi)*(jt[i].scale5()-1.);
        }
    }
    float met_phi = atan2(met_ey, met_ex);
    if ( met_phi < 0 ) met_phi += 2*M_PI;
    missing_et.set_ex(met_ex);
    missing_et.set_ey(met_ey);
    missing_et.set_et(sqrt(pow(met_ex,2)+pow(met_ey,2)));
    missing_et.set_phi(met_phi);
}

void TSeel::correct_met_muons()
{
    float met_ex = missing_et.ex();
    float met_ey = missing_et.ey();

    for (int i=0; i<nleptons_temp; i++ ) {
      if ( abs(lep_temp[i].id()) != 13 ) continue;
      TStnMuon *mu = fMuonBlock->Muon(lep_temp[i].index());
      //int it = mu->TrackNumber();
      //TStnTrack *trk = fTrackBlock->Track(it);
      //TLorentzVector *mom = new TLorentzVector();
      //int status = trk->GetBcMomentum(mom);
      float eta = mu->Momentum()->Eta();
      float phi = mu->Momentum()->Phi();
      //float eta = mom->Eta();
      //float phi = mom->Phi();
      met_ex = met_ex - (lep_temp[i].et() - (mu->EmEnergy() + mu->HadEnergy())/cosh(eta))*cos(phi);
      met_ey = met_ey - (lep_temp[i].et() - (mu->EmEnergy() + mu->HadEnergy())/cosh(eta))*sin(phi);
    }
    float met_phi = atan2(met_ey, met_ex);
    if ( met_phi < 0 ) met_phi += 2*M_PI;
    missing_et.set_ex(met_ex);
    missing_et.set_ey(met_ey);
    missing_et.set_et(sqrt(pow(met_ex,2)+pow(met_ey,2)));
    missing_et.set_phi(met_phi);
}

void TSeel::EventReset()
{
// Make raw histograms
  if ( nleptons_temp>=3 )
     if ( lep_temp[0].et()>5. && lep_temp[1].et()>5. && lep_temp[2].et()>5. ) {
       h.raw_dr_lep12->Fill(delta_r(lep_temp[0].eta(),lep_temp[0].phi(),lep_temp[1].eta(),lep_temp[1].phi()));
       h.raw_dr_lep13->Fill(delta_r(lep_temp[0].eta(),lep_temp[0].phi(),lep_temp[2].eta(),lep_temp[2].phi()));
       h.raw_dr_lep23->Fill(delta_r(lep_temp[1].eta(),lep_temp[1].phi(),lep_temp[2].eta(),lep_temp[2].phi()));
     }
  
//================================================================================
//  Make deltaR cuts
    int lep_dr_ok[50];
    for ( int i=0; i<50; i++ ) lep_dr_ok[i] = 1;
    for ( int i=0; i<nleptons_temp-1; i++ ) {  // lepton-lepton dr
        for ( int j=i+1; j<nleptons_temp; j++ ) {
            if ( lep_temp[i].et()>4. && lep_temp[j].et()>4. && delta_r(lep_temp[i].eta(),lep_temp[i].phi(),lep_temp[j].eta(),lep_temp[j].phi())<0.4 ) {
                lep_dr_ok[i] = 0;
                lep_dr_ok[j] = 0;
            }
        }
    }

    for ( int i=0; i<nleptons_temp; i++ ) {   // lepton-jet dr
        float drmin = 999.;
        if ( lep_temp[i].et()>4. ) {
	  int jt_index = -1;
	  for ( int j=0; j<njet; j++ ) {
	    float dr=delta_r(lep_temp[i].eta(),lep_temp[i].phi(),jt[j].eta(),jt[j].phi());
	    if ( jt[j].et()>FINAL_JET_CUT && dr<drmin ) {
	      jt_index = j;
	      drmin = dr;
	    }
	  }
	  if ( drmin<0.4 ) {
	    lep_dr_ok[i]=0;
	  }
	}
    }

//  Write the much anticipated lep[] array
//  An account for the lep_temp[] are in order here:
//  It's first defined as tight leptons(electrons and muons).
//  Then it's sorted in Et.
//  All these done in the main Loop
    nleptons=0; ngood_leptons=0; n_ele=0; n_mu=0;
    
    for ( int i=0; i<nleptons_temp; i++ ) {  // lep_temp[] contains ele and isotrk
        if ( lep_dr_ok[i] ) {
	  lep[nleptons] = lep_temp[i];
	  if ( abs(lep_temp[i].id())==11 ) {
	    lep_ele[n_ele] = lep_temp[i];
	    n_ele++;
	  }
	  if ( abs(lep_temp[i].id())==13 ) {
	    lep_mu[n_mu] = lep_temp[i];
	    n_mu++;
	  }
	  nleptons++;
	  if ( getbits(lep_temp[i].qual(),23,1)==1 )
	    ngood_leptons++;
	}
    }
    if(fVerbose<=-2 && nleptons>1){
      cout<<"n_ele = "<<n_ele<<" ,  n_mu = "<<n_mu<<"  , nleptons = "<<nleptons<<endl;
      cout<<"---------------------------------------------"<<endl;
    }

}

float TSeel::transverse_mass_pt_phi(float pt1,float phi1,float pt2,float phi2)
{
    float dphi = delta_phi(phi1,phi2);
    return sqrt(2.*pt1*pt2*dphi);
}

float TSeel::transverse_mass_px_py(float px1, float py1, float px2, float py2)
{
    float px = px1 + px2;
    float py = py1 + py2;
    float pt = sqrt(pow(px1,2)+pow(py1,2)) + sqrt(pow(px2,2)+pow(py2,2));
    float mt = pow(pt,2) - pow(px,2) - pow(py,2);
    return mt>0. ? sqrt(mt) : 0;
}

float TSeel::transverse_mass(const vector<TLorentzVector> &v)
{
    float sum_pt = 0.;
    TLorentzVector vsum;
    for ( vector<TLorentzVector>::const_iterator i=v.begin(); i!=v.end(); i++ ) {
        vsum += *i;
        sum_pt += i->Pt();
    }
    float mt2 = pow(sum_pt,2)-pow(vsum.Pt(),2);
    float mt = 0.;
    if ( mt2>0. ) mt = sqrt(mt2);
    return mt;
}
float TSeel::Thrust(const vector<TLorentzVector> &vin)
{
    vector<TVector3> p;

    TVector3 qtbo;
    TVector3 zero(0.,0.,0.);
    float vnew;
    float thrust = 0.;

  // select the three-momenta
    for ( vector<TLorentzVector>::const_iterator i=vin.begin(); i!=vin.end(); i++ ) 
        p.push_back(i->Vect());
    int Np = p.size();
    
    if ( Np < 0 ) {
        qtbo = zero;
        return thrust;
    }
    else if ( Np==1 )
        qtbo = p[0];
    else if ( Np==2 ) {
        if (p[0].Dot(p[1]) >= 0.) 
            qtbo = p[0] + p[1];
        else
            qtbo = p[0] - p[1];
    }
    else if ( Np>2 ) { // for more than 2 tracks
        float vmax  = 0.;
        TVector3 vn, vm, vc, vl;
        for(int i=0; i<Np-1; i++) 
            for(int j=i+1; j < Np; j++) {
                vc = p[i].Cross(p[j]);
                vl = zero; 
                for(int k=0; k<Np; k++)
                    if ((k != i) && (k != j))
                        if (p[k].Dot(vc) >= 0.) vl = vl + p[k];
                        else vl = vl - p[k];
            // make all four sign-combinations for i,j
                vn = vl + p[j] + p[i];
                vnew = vn.Mag2();
                if (vnew > vmax) {  
                    vmax = vnew;
                    vm = vn;
                }
                vn = vl + p[j] - p[i];
                vnew = vn.Mag2();
                if (vnew > vmax) {  
                    vmax = vnew;
                    vm = vn;
                }
                vn = vl - p[j] + p[i];
                vnew = vn.Mag2();
                if (vnew > vmax) {
                    vmax = vnew;
                    vm = vn;
                }
                vn = vl - p[j] - p[i];
                vnew = vn.Mag2();
                if (vnew > vmax) {
                    vmax = vnew;
                    vm = vn;
                }
            }
            // sum momenta of all particles and iterate
      
        for(int iter=1; iter<=4; iter++)
        {  
            qtbo = zero;
            for(int i=0; i< Np; i++)
                if (vm.Dot(p[i]) >= 0.) 
                    qtbo = qtbo + p[i];
                else 
                    qtbo = qtbo - p[i];
            vnew = qtbo.Mag2();
            if (vnew  == vmax) break;
            cout << endl << " I have to iterate again " <<  endl << endl;
            vmax = vnew;
            vm = qtbo;
        }
    }  // of if Np > 2
  // normalize thrust -division by total momentum-
    float vsum = 0.;
    for(int i=0; i < Np; i++) vsum = vsum + p[i].Mag();
    vnew  = qtbo.Mag();
    float v = vnew/vsum;
    float x = qtbo.X()/vnew;
    float y = qtbo.Y()/vnew;
    float z = qtbo.Z()/vnew;
    thrust = v;
    return thrust;
}

bool TSeel::Sphericity_Aplanarity(const vector<TLorentzVector> &v, float &sphe, float &aplan)
{
    if ( v.size()<2 ) {
        cout << "ERROR: cannot compute Sphericity on an empty or one vector." << endl;
        sphe = -1.;
        aplan = -1.;
        return false;
    }

/* Comment out jzz20050623
    vector<TVector3> p;
  // select the three-momenta
    for ( vector<TLorentzVector>::const_iterator i=v.begin(); i!=v.end(); i++ ) 
        p.push_back(i->Vect());
*/

    double x[3][3];
    for ( int i=0; i<3; i++ )
        for ( int j=0; j<3; j++ )
            x[i][j] = 0.;
//    for ( vector<TVector3>::const_iterator ip=p.begin(); ip!=p.end(); ip++ )
    for ( vector<TLorentzVector>::const_iterator ip=v.begin(); ip!=v.end(); ip++ )
        for ( int i=0; i<3; i++ )
            for ( int j=0; j<3; j++ ) {
                x[i][j] += (*ip)(i) * (*ip)(j);
//                cout << "[" << i << "][" << j << "]: " << (*ip)(i) << " " << (*ip)(j) << " " << (*ip)(i) * (*ip)(j) << " " << x[i][j] << endl;
            }
    float sum_diagonal = x[0][0]+x[1][1]+x[2][2];
    for ( int i=0; i<3; i++ ) {
        for ( int j=0; j<3; j++ ) {
            x[i][j] /= sum_diagonal;
//            cout << " x[" << i << "][" << j << "]=" << x[i][j];
        }
//        cout << endl;
    }

    TMatrixD Mab(3,3);
    for (int i = Mab.GetRowLwb(); i <= Mab.GetRowUpb(); i++)
        for (int j = Mab.GetColLwb(); j <= Mab.GetColUpb(); j++)
            Mab(i,j) = x[i][j];
    double z[] = {0., 0., 0.};
    TVectorD eigen_values(3,z);
    TMatrixD M_eigen = Mab.EigenVectors(eigen_values);
    vector<double> ez;
    // sort eigen_values in non-decreasing order
    ez.push_back(eigen_values(0));
    ez.push_back(eigen_values(1));
    ez.push_back(eigen_values(2));
    sort(ez.begin(), ez.end());
    sphe = 1.5*(1.-ez[2]);
    aplan = 1.5*ez[0];

    return true;
}
//****************************************************************
//   Baysian_error
// 
//   Purpose: It returns the turn-on error according CDF-5894
//   Input: num, den
//   Output: myeff[3]: [0]=eff, [1]=low error (neg), [2]=high error (pos)
//   REMARKS: by Tony Munor
//
//***************************************************************
void TSeel::Bayesian_Error(Float_t num,Float_t den,Float_t* myeff)
{
    static Int_t numID=0;
    Char_t tit[200];

    Int_t numI;
    Int_t denI;

    numI=(Int_t) num;
    denI=(Int_t) den;

    if(den==0) {
        myeff[0]=0;
        myeff[1]=0;
        myeff[2]=0;
    }
// For very little efficiencies I found the optimal bin size
// to be this one (the area of the distribution is quite close to 1)
    else {
        Double_t bin_width=0.;
        Double_t nbin=0;

        bin_width=1./(Float_t)denI/100.;  
        nbin=1./bin_width; 

        if(nbin<10000) {
            nbin=10000;
            bin_width=0.0001;
        }

//    printf("########################################################\n");
//    printf("## Number of observed: %d Number of trials: %d Eff: %f\n",
//           numI,denI,num/den);
//    printf("## Number of bins for Efficiency Histogram: %d\n",nbin);

// Book the histogram for the efficiencies

        sprintf(tit,"Efficiency Probability Density. Trials: %d Observed: %d Nbins: %d",
                denI,numI,nbin);
        TH1F *hef=new TH1F("ef",tit,nbin,0.,1.);
        hef->SetXTitle("#epsilon");

//Compute the distribution
        Double_t k=0.;
  //lgamma
        k=TMath::LnGamma(denI+1)-TMath::LnGamma(numI+1)-TMath::LnGamma(denI-numI+1);
        Float_t ef=bin_width/2.;
        Float_t d=0.;
        for(Int_t ii=1;ii<=nbin;ii++) {
            d=k+numI*TMath::Log(ef)+(denI-numI)*TMath::Log(1.-ef)+TMath::Log(denI+1);
            d=TMath::Exp(d);
            hef->SetBinContent(ii,d);
            ef=bin_width*ii+bin_width/2.;
        }
//    printf("## Distribution Area: %f\n",hef->Integral(1,nbin,"width"));  

  //printf("k: %f\n",k);
        ef=bin_width/2.;
        d=0.;
        for(Int_t ii=1;ii<=nbin;ii++) {
            d=k+numI*TMath::Log(ef)+(denI-numI)*TMath::Log(1.-ef)+TMath::Log(denI+1);
            d=TMath::Exp(d);
            hef->SetBinContent(ii,d);
            ef=bin_width*ii+bin_width/2.;
        }

//________________________________________________________________________________
  // Compute the errors
  // The scheme is as follows
  // compute the bin number for the maximum and
  // the probability for such bin
        Int_t max_bin=hef->GetMaximumBin();
        Float_t max_bin_p=hef->GetMaximum()*bin_width;
        Int_t nBin=hef->GetNbinsX();
//    printf("-->Number of bins: %d\n",nBin);
//    printf("--> bin the maximum probability: %d Max prob: %f\n",max_bin,max_bin_p);
        sprintf(tit,"h%d",numID);numID++;
        hef->SetName(tit);

//        TFile *fout=new TFile("fdump.root","UPDATE");
//        fout->cd();
//        hef->Write();
//        fout->Close();

        Float_t ls=0.;
        Int_t ls_bin=1;
 
  //Begin with bin 1, continue up to the bin where the probability is 34%
  //Begin with bin of maximum probability. Continue downwards up to 34%
        for(Int_t ii=max_bin-1;ii>0;ii--) {
            ls_bin=ii;
            ls=hef->Integral(ii,max_bin-1,"width")+max_bin_p/2.;
            if(ls>=0.34) break;
        }
  //printf("Best estimation NOBS/NTRIALS: %f bin: %d bin width: %f\n",
  //     hef->GetBinLowEdge(max_bin),max_bin,bin_width);
  //printf("Left Side: at bin: %d we have %f of the events\n",ls_bin,ls);
 
        Float_t rs=0.;
        Int_t rs_bin=nbin;
 
  //Begin with bin max+1, continue up to probability is 34%
        for(Int_t ii=max_bin+1;ii<=nbin;ii++) {
            rs_bin=ii;
            rs=hef->Integral(max_bin+1,ii,"width")+max_bin_p/2.;
            if(rs>=0.34) break;
        }
  //printf("Right Side: at bin: %d we have %f of the events\n",rs_bin,rs);
 
  //Now, if in one side the probability is less than 34%, make the other
  //side bigger
        Float_t prob=0;
        if(rs<0.34) {
            prob=ls+0.34-rs;
            for(Int_t ii=max_bin-1;ii>0;ii--) {
                ls_bin=ii;
                ls=hef->Integral(ii,max_bin-1,"width")+max_bin_p/2.;
                if(ls>=prob) break;
            }
    //printf("Second Pass:\n\t");
    //printf("Left Side: at bin: %d we have %f of the events\n",ls_bin,ls);  
        }
 
        if(ls<0.34) {
            prob=rs+0.34-ls;
            for(Int_t ii=max_bin+1;ii<=nbin;ii++) {
                rs_bin=ii;
                rs=hef->Integral(max_bin+1,ii,"width")+max_bin_p/2.;
                if(rs>=prob) break;
            }
    //printf("Second Pass:\n\t");
    //printf("Right Side: at bin: %d we have %f of the events\n",rs_bin,rs);
        }
 
//        printf("num: %f, den: %f, eff=%f, el_eff=%f, eh_eff=%f\n",
//               num, den,
//               hef->GetBinLowEdge(max_bin),
//               hef->GetBinLowEdge(ls_bin)-hef->GetBinLowEdge(max_bin),
//               hef->GetBinLowEdge(rs_bin)-hef->GetBinLowEdge(max_bin));
//        printf("--------------------------------------------------------------------------\n");
        myeff[0] =  hef->GetBinLowEdge(max_bin);
        myeff[1] =  hef->GetBinLowEdge(ls_bin)-hef->GetBinLowEdge(max_bin);
        myeff[2] =  hef->GetBinLowEdge(rs_bin)-hef->GetBinLowEdge(max_bin);

        hef->Delete();
    }
}

//****************************************************************
//   Binomial_error
// 
//   Purpose: It returns the turn-on error according CDF-5894
//   Input: num, den
//   Output: myeff[3]: [0]=eff, [1]=low error (neg), [2]=high error (pos)
//
//***************************************************************
void TSeel::Binomial_Error(Float_t num,Float_t den,Float_t* myeff)
{
    if ( den==0 ) {
        myeff[0]=0;
        myeff[1]=0;
        myeff[2]=0;
    }
    else if ( num==den ) {
        if ( den==1 ) {
            if ( num==1 ) {
                myeff[0]=1.;
            }
            else if ( num==0 ) {
                myeff[0]=0.;
            }
            myeff[1]=-1.;
            myeff[2]=1.;
        }
    }
    else {
        if ( num==0 ) {
            myeff[0]=0.;
            myeff[1]=-sqrt(1.-1./den)/den;
            myeff[2]=-myeff[1];
        }
        else {
            myeff[0]=num/den;
            myeff[1]=-sqrt((1.-myeff[0])*myeff[0]/den);
            myeff[2]=-myeff[1];
        }
    }
}
float TSeel::inv_mass(const float p1[4], const float p2[4])
{
    float en = p1[0] + p2[0];
    float px = p1[1] + p2[1];
    float py = p1[2] + p2[2];
    float pz = p1[3] + p2[3];

    float mass2 = pow(en,2)-pow(px,2)-pow(py,2)-pow(pz,2);
    if ( mass2>=0. )
        return sqrt(mass2);
    else
        return -sqrt(-mass2);
}
float TSeel::inv_mass(const vector<TLorentzVector> &v)
{
    TLorentzVector vsum;
    for ( vector<TLorentzVector>::const_iterator i=v.begin(); i!=v.end(); i++ )
        vsum += *i;
    return vsum.M();
}

bool TSeel::pass_EventReset(int nlep)
{
//  Require that the leading two leptons in lep[] are electrons or muons
//  and all three leading lep[] pass the Et cut of tri-lep analysis.
    if ( nleptons<nlep ) return false;
    if ( ngood_leptons<2 ) return false;

    bool pass = true;
    if ( nlep==2 ) {
      for ( int i=0; i<2; i++ ) {
	pass = pass && lep[i].et()>_LEP_PT_CUT2[i];
      }
    }
    else if ( nlep==3 ) {
      for ( int i=0; i<3; i++ ) { // First two leptons being "true" and good leptons
	pass = pass && lep[i].et()>_LEP_PT_CUT3[i];
      }
    }
    bool leading_lep_are_tight_ele = false;

    if(_OptScheme==3){
      if( abs(lep[0].id())==11 && abs(lep[1].id())==11
	  && getbits(lep[0].qual(),23,1)==1 )
	leading_lep_are_tight_ele = true;
    }
    else {
      if( abs(lep[0].id())==11 && abs(lep[1].id())==11
	  && getbits(lep[0].qual(),23,1)==1 && getbits(lep[1].qual(),23,1)==1 )
	leading_lep_are_tight_ele = true;
    }

    pass = pass && leading_lep_are_tight_ele;

    
    return pass;
}

void TSeel::Set_EventType(int nlep)
{
//  See JZT.h for the bit definition.
//  For events with good triele: _EventType=0x524, diele+track: _EventType=0x521
//  Events with good eem:_EventType=0x522
  _EventType = 0x0;
  if ( nleptons>=2 ) {
        if ( abs(lep[0].id())==11 )     _EventType |= 0x1 << 8;
        else if ( abs(lep[0].id())==13 ) _EventType |= 0x1 << 7;
        if ( abs(lep[1].id())==11 )     _EventType |= 0x1 << 5;
        else if ( abs(lep[1].id())==13 ) _EventType |= 0x1 << 4;
        if ( nleptons>=3 ) {
            if ( abs(lep[2].id())==11 )     _EventType |= 0x1 << 2;
            else if ( abs(lep[2].id())==13 ) _EventType |= 0x1 << 1;
	    if( (_EventType & 0x124) == 0x124 ){
	      int sum_charge = 0;
	      for ( int i=0; i<3; i++ ) sum_charge += lep[i].id()/abs(lep[i].id());
	      if ( abs(sum_charge)==1 ) _EventType |= 0x1 << 10;
	    }
	    else if ( (_EventType & 0x122) == 0x122 ){
	      int sum_charge = 0;
	      for ( int i=0; i<2; i++ ) sum_charge += lep[i].id()/abs(lep[i].id());
	      if ( abs(sum_charge)==0 ) _EventType |= 0x1 << 10;
	    }
	}
  }
}     

bool TSeel::pass_BaseCuts()
{
    bool pass_charge_combination = false;
    if ( ((_EventType & 0x524) == 0x524) || ((_EventType & 0x522) == 0x522) ) {
        pass_charge_combination = true;
    }

    return pass_charge_combination;
}

bool TSeel::pass_MassCuts()
{
  bool result = false;

    vector<TLorentzVector> v_lep;
    for ( int i=0; i<3; i++ ) {
        TLorentzVector v(lep[i].px(),lep[i].py(),lep[i].pz(),lep[i].en());
        v_lep.push_back(v);
    }
    vector<TLorentzVector> v_os_1,v_os_2;
    if ( lep[0].id()*lep[1].id()<0 ) {
        v_os_1.push_back(v_lep[0]);
        v_os_1.push_back(v_lep[1]);
        if ( lep[0].id()*lep[2].id()<0 ) {
            v_os_2.push_back(v_lep[0]);
            v_os_2.push_back(v_lep[2]);
        }
        else if ( lep[1].id()*lep[2].id()<0 ) {
            v_os_2.push_back(v_lep[1]);
            v_os_2.push_back(v_lep[2]);
        }
    }
    else if ( lep[0].id()*lep[2].id()<0 ) {
        v_os_1.push_back(v_lep[0]);
        v_os_1.push_back(v_lep[2]);
        v_os_2.push_back(v_lep[1]);
        v_os_2.push_back(v_lep[2]);
    }
    if ( v_os_1.empty() )
        m_os_1 = 0.;
    else
        m_os_1 = inv_mass(v_os_1);
    if ( v_os_2.empty() )
        m_os_2 = 0.;
    else
        m_os_2 = inv_mass(v_os_2);

    result = m_os_1>15. && m_os_2>15.;
    
    if( (_EventType & 0x524) == 0x524){ 
      if ( m_os_1 < m_os_2 ) {
	float temp = m_os_1;
	m_os_1 = m_os_2;
	m_os_2 = temp;
      }
    }
    return result;
}

void TSeel::get_CEM_ID_scale_factor(float et, float &x, float &ex, int &iBin)
{
    if ( et>4. && et<6. ) {
//        x = 0.80; ex = 0.0860; iBin=0;
//        x = 0.80; ex = 0.20; iBin=0;
        x = 1.0; ex = 0.20; iBin = 0;
    }
    else if ( et>6. && et<8. ) {
//        x = 0.88; ex = 0.0781; iBin=1;
//        x = 0.88; ex = 0.20; iBin=1;
        x = 1.0; ex = 0.20; iBin = 1;
    }
    else if ( et>8. && et<20. ) {
        x = 1.03; ex = 0.02; iBin=2;
    }
    else if ( et>20. ) {
        x = 0.999; ex = 0.0072; iBin=3;
    }
}

void TSeel::get_CEM4_trig_eff(float et, float &x, float &err0, float &err1, float &err2)
{
// SOURABH changes on 20JUL2005
// These are the parameters from initial study. CDF7095, jet20->conversions
//    float p0 = 4.902;  float ep0 = 0.1639;
//    float p1 = 0.4741; float ep1 = 0.1902;
//    float p2 = 0.956; float ep2 = 0.02826;
// These are parameters from new trigger study, blpc0d->electrons
// IF run<181801 use set=1 CES2 at Level2
// IF run>181839 use set=2 CES3 at Level2
  float p0 = 4.799;  float ep0 = 0.01576;
  float p1 = 0.5025; float ep1 = 0.02029;
  float p2 = 0.9633; float ep2 = 0.001475;
// Following is set2
// p0 = 5.004  +/- 0.02074
// p1 = 0.6731 +/- 0.02666
// p2 = 0.9621 +/- 0.001544

//    float err_sys = 0.03;
    float arg = 0.5*(et-p0)/p1;
    x = 0.5*p2*(1.+TMath::Erf(arg));
    err0 = p2/sqrt(M_PI)*exp(-pow(arg,2))*0.5*ep0/p1;
    err1 = p2/sqrt(M_PI)*exp(-pow(arg,2))*0.5*(et-p0)*ep1/pow(p1,2);
    err2 = 0.5*ep2*(1.+TMath::Erf(arg));
//    ex = sqrt(pow(err0,2)+pow(err1,2)+pow(err2,2));
//  No systematic uncertainty convoluted in the total uncertainty
//    ex = sqrt(pow(err0,2)+pow(err1,2)+pow(err2,2)+pow(err_sys,2));
//  New finding
//    x = 0.9795; ex = 0.0018;

//    if ( fVerbose>10 ) cout << "arg=" << arg << " x=" << x << " ex=" << ex << " err0=" << err0 << " err1=" << err1 << " err2=" << err2 << endl;
}

void TSeel::get_CMU_ID_scale_factor(float pt, float &x, float &ex)
{
    x = 1.0;
    ex = 0.02;
}

void TSeel::get_CMX_ID_scale_factor(float pt, float &x, float &ex)
{
    x = 1.0;
    ex = 0.02;
}

void TSeel::get_CMP_ID_scale_factor(float pt, float &x, float &ex)
{
    x = 1.0;
    ex = 0.02;
}

void TSeel::get_CMUP_ID_scale_factor(float pt, float &x, float &ex)
{
    x = 1.0;
    ex = 0.02;
}

float TSeel::track_exit_radius(float eta, float z0)
{
    float rho;
    const float zcot = 155.;
    const float theta = 2.*atan(exp(-eta));
    const float lambda = TMath::PiOver2()-theta;

    if ( eta>0 )
        rho = (zcot-z0)/tan(lambda);
    else
        rho = -(zcot+z0)/tan(lambda);

    return rho;
}

float TSeel::cmpChi2(TStnMuon *mu)  // from MuonMatchCalculator.cc
{
    float chiXPos = -999.;
    if ( mu->HasCmpStub() ) {
        const double P0 = -49.33;
        const double P1 =  4.621;
        const double P2 = -0.02634;
        const double P3 =  1.526;
        const double P4 =  2.8;
        const double P5 = -0.27;

        double ptMuon = mu->TrackPt();
        double phi0Muon = mu->Momentum()->Phi();

       //pT dependent term that modulates the amplitude of the phi0 dependence
        const double ptDepCoeff = max(P4+P5*ptMuon,0.);
 
       //drphi parameterisation
        const double cmpDrphiParm = (P0+exp(P1+P2*ptMuon))/ptMuon + P3 + 
            ptDepCoeff*(1.-sin(4.*phi0Muon+M_PI/2));
       
        chiXPos = float( fabs(mu->CmpDelX()) / cmpDrphiParm );
        chiXPos *= chiXPos;
    }

    return chiXPos;
}

float TSeel::cmxChi2(TStnMuon *mu) // from MuonMatchCalculator.cc
{
    float chiXPos = -999.;
    if ( mu->HasCmxStub() ) {
        const double P0 = 0.320;
        const double P1 = 34.4;
        const double P2 = 35.81;

        double ptMuon = mu->TrackPt();
        //drphi parameterisation
        const double cmxDrphiParm = P0 + P1/ptMuon + P2/pow(ptMuon,2);  
       
        chiXPos = float( fabs(mu->CmxDelX()) / cmxDrphiParm );
        chiXPos *= chiXPos;
    }

    return chiXPos;
}


int TSeel::Set_EventWeight()
{
    int return_code = 0;
    if ( ngood_leptons<2 ) {
        _EventWeight=0.;
        _EventWeightErrorPlus=0.;
        _EventWeightErrorMinus=0.;
        return -1;
    }
    bool pass0 = true;
    for ( int i=0; i<2; i++ ) { // First two leptons being "true" and good leptons
        pass0 = pass0 && lep[i].et()>_LEP_PT_CUT3[i];
    }
    //bool pass3 = pass0 && lep[2].et()>_LEP_PT_CUT3[2];
    if ( ! pass0 ) {
        _EventWeight = 0.;
        _EventWeightErrorPlus = 0.;
        _EventWeightErrorMinus = 0.;
        return -2;
    }

    _EventWeight = 0.;
    _EventWeightErrorPlus = 0.;
    _EventWeightErrorMinus = 0.;

    float p1[4], p2[4];
    p1[0]=lep[0].en(); p2[0]=lep[1].en();
    p1[1]=lep[0].px(); p2[1]=lep[1].px();
    p1[2]=lep[0].py(); p2[2]=lep[1].py();
    p1[3]=lep[0].pz(); p2[3]=lep[1].pz();

    if ( lep[0].id()*lep[1].id()==-121 )
        h.dilep_mass_ee_os->Fill(inv_mass(p1,p2));
    else if ( lep[0].id()*lep[1].id()==121 )
        h.dilep_mass_ee_ss->Fill(inv_mass(p1,p2));
    else if ( lep[0].id()*lep[1].id()==-169 )
        h.dilep_mass_mm_os->Fill(inv_mass(p1,p2));
    else if ( lep[0].id()*lep[1].id()==169 )
        h.dilep_mass_mm_ss->Fill(inv_mass(p1,p2));
    else if ( lep[0].id()*lep[1].id()==-143 )
        h.dilep_mass_em_os->Fill(inv_mass(p1,p2));
    else if ( lep[0].id()*lep[1].id()==143 )
        h.dilep_mass_em_ss->Fill(inv_mass(p1,p2));

    if ( _Data==1 ) {
// If this is a data event we have two things to do -
// 1) We are running the analysis and wish to see how many data events we see, as in the final step.
// 2) We are running over data and we wish to estimate fakes.
// This is step 1)
        _EventWeight = 1.;
        _EventWeightErrorPlus = 0.;
        _EventWeightErrorMinus = 0.;
	//return 0;
    }
    if ( _Data==1 ) {
// This is step 2) Here this is what we do:
// For every dilepton event, we count the number of fakeable objects. We get the fakerate for each
// object and we sum over it. This gives us the event weight, (fake event's weight.)
// The sum of _fakeEventWeight at the end of whole data is the number of fakes we expect.
      if ( Pass_Dilep_BaseCuts ) {
	float fake_evtwgt=0;
	float fake_optevtwgt=0;
	//fake_evtwgt=fake_optevtwgt=0;
	fake_optevtwgt = get_Fake_Weight(fake_evtwgt);
	float err_fake_evtwgt = 0;
	_fakeEventWeight = fake_evtwgt;
	_fakeEventWeightErrorPlus  = err_fake_evtwgt; 
	_fakeEventWeightErrorMinus = err_fake_evtwgt;
	_fakeOptEventWeight = fake_optevtwgt;
        return -5;
      }
    }

    float trig_eff[10], err_trig_eff[10], ep0[10], ep1[10], ep2[10];
    float id_eff[10],   err_id_eff[10];
    for ( int i=0; i<10; i++ ) {
        trig_eff[i]=0.;
        id_eff[i]=0.;
        ep0[i] = ep1[i] = ep2[i] = 0.;
    }



//  Specific for diele+track channel
// For ID, there are four bins:
//   4-6, 6-8, 8-20, >20. Check if any of the two eles fall into the same bin.
//  ID weight, jzz20050313.
//  Assume track ID scale factor = 1 and ignore its uncertainty.
    int eID_BinContent[4]={0};
    for ( int i=0; i<3; i++ ) {  // We require only two good electrons.
        if ( abs(lep[i].id())==11 ) {
            int iBin = -1;
            get_CEM_ID_scale_factor(lep[i].et(),id_eff[i],err_id_eff[i],iBin);
            eID_BinContent[iBin]++;
        }
	if ( abs(lep[i].id())==13 ) {
	  //int iBin = -1;
	  get_CMU_ID_scale_factor(lep[i].et(),id_eff[i],err_id_eff[i]);
	}
    }
    float id_weight = 0.; float err_id_weight = 0.;
    id_weight = id_eff[0]*id_eff[1]*id_eff[2];
    err_id_weight = id_weight*sqrt(pow(err_id_eff[0]/id_eff[0],2)+pow(err_id_eff[1]/id_eff[1],2)+pow(err_id_eff[2]/id_eff[2],2));
    for ( int i=0; i<4; i++ ) {
        if ( eID_BinContent[i]==2 ) {
            err_id_weight = 2.*id_eff[0]*err_id_eff[0]; // same eff in the same bin
            break;
        }
    }

// Trigger weight, jzz20050313
// Assign trigger efficiency and uncertainty on one CEM4 term by each electron
    for ( int i=0; i<min(nleptons,10); i++ ) {
        if ( abs(lep[i].id())==11 ) {
            get_CEM4_trig_eff(lep[i].et(),trig_eff[i],ep0[i],ep1[i],ep2[i]);
        }
	/*
        else if ( abs(lep[i].id())==1 ) {
            if ( getbits(lep[i].qual(),12,1) ) { // This is an emobject track
                if ( EmTrack_pass_L3EM(lep[i].index()) ) {
                    get_CEM4_trig_eff(fElectrons_fEt[lep[i].index()],trig_eff[i],ep0[i],ep1[i],ep2[i]);
        // add 10% error since we use offline info
                    ep0[i] *= 1.1;
                    ep1[i] *= 1.1;
                    ep2[i] *= 1.1;
                }
            }
        }
	*/
    }
// Calculate the AllFail and AllFailButOne efficiency.
    float trig_weight = 0.; float err_trig_weight = 0.;
    float ep0_total = 0.; float ep1_total = 0.; float ep2_total = 0.;
    float AllFailButOne = 0.; float AllFail = 1.;
    for ( int i=0; i<min(nleptons,10); i++ ) {
        AllFail *= 1.-trig_eff[i];
        float AllFailButOneI = 1.;
        for ( int j=0; j<min(nleptons,10); j++ ) {
            if ( j==i ) continue;
            AllFailButOneI *= 1.-trig_eff[j];
        }
        AllFailButOne += AllFailButOneI;
    }
// Now calculate the uncertainties
    for ( int i=0; i<min(nleptons,10); i++ ) {
        float AllFailButOneI = 1.;
        for ( int j=0; j<min(nleptons,10); j++ ) {
            if ( j==i ) continue;
            AllFailButOneI *= 1.-trig_eff[j];
        }
    // (AllFailButOne - (min(nleptons,10)-1)*AllFail - AllFailButOneI)
    //  gives the terms that has (1-trig_eff[i]) factor in the efficiecy formula.
    //  Take the derivative of these terms is equivalent to dividing
    //  (1-trig_eff[i]) and times -1, which is omitted in quadratic sum
    //  of errors.
        ep0_total += pow((AllFailButOne - (min(nleptons,10)-1)*AllFail - AllFailButOneI)/(1.-trig_eff[i])*ep0[i],2);
        ep1_total += pow((AllFailButOne - (min(nleptons,10)-1)*AllFail - AllFailButOneI)/(1.-trig_eff[i])*ep1[i],2);
        ep2_total += pow((AllFailButOne - (min(nleptons,10)-1)*AllFail - AllFailButOneI)/(1.-trig_eff[i])*ep2[i],2);
    }
    trig_weight = 1. - AllFailButOne + (min(nleptons,10)-1)*AllFail;
    err_trig_weight = sqrt(ep0_total + ep1_total + ep2_total);


//  EventWeight depends on the sample

    if ( _Signal ) {
        if ( Pass_Trilep_BaseCuts && Pass_Trilep_MassCuts ) {
            _EventWeight = trig_weight*id_weight;
            if ( trig_weight>0 && id_weight>0. ) {
                _EventWeightErrorPlus = _EventWeight*sqrt(pow(err_trig_weight/trig_weight,2)+pow(err_id_weight/id_weight,2));
                _EventWeightErrorMinus = _EventWeight*sqrt(pow(err_trig_weight/trig_weight,2)+pow(err_id_weight/id_weight,2));
            }
            else {
                _EventWeightErrorPlus = 0.;
                _EventWeightErrorMinus = 0.;
            }
        }
        else {
            _EventWeight = 0.;
            _EventWeightErrorPlus = 0.;
            _EventWeightErrorMinus = 0.;
        }
        return_code = 1;
    }
    
    else if ( Pass_Trilep_BaseCuts ) {
      if ( Pass_Trilep_MassCuts ) {
                _EventWeight = trig_weight*id_weight;
                if ( trig_weight>0 && id_weight>0. ) {
                    _EventWeightErrorPlus = _EventWeight*sqrt(pow(err_trig_weight/trig_weight,2)+pow(err_id_weight/id_weight,2));
                    _EventWeightErrorMinus = _EventWeight*sqrt(pow(err_trig_weight/trig_weight,2)+pow(err_id_weight/id_weight,2));
                }
                else {
                    _EventWeightErrorPlus = 0.;
                    _EventWeightErrorMinus = 0.;
		}
		return_code = 2;
      }	
    }
    else {
        _EventWeight = trig_weight*id_weight;
        if ( trig_weight>0 && id_weight>0. ) {
            _EventWeightErrorPlus = _EventWeight*sqrt(pow(err_trig_weight/trig_weight,2)+pow(err_id_weight/id_weight,2));
            _EventWeightErrorMinus = _EventWeight*sqrt(pow(err_trig_weight/trig_weight,2)+pow(err_id_weight/id_weight,2));
        }
        else {
            _EventWeightErrorPlus = 0;
            _EventWeightErrorMinus = 0;
        }
        if ( lep[0].id()*lep[1].id()>0 ) {
        // take care of prob of failing the cut on charge combination
            _EventWeight /=2.;
            _EventWeightErrorPlus /=2.;
            _EventWeightErrorMinus /=2.;
        }
        return_code = 3;
    }      
    /*
    if ( nleptons>=3 ) {
        h.trilep_trig_weight->Fill(trig_weight);
        h.trilep_lepID_weight->Fill(id_weight);
        h.trilep_event_weight->Fill(_EventWeight);
        h.trilep_trig_weight_err->Fill(err_trig_weight);
        h.trilep_lepID_weight_err->Fill(err_id_weight);
        h.trilep_event_weight_err->Fill((_EventWeightErrorPlus+_EventWeightErrorMinus)/2.);
    }
    else if ( nleptons>=2 ) {
        h.dilep_trig_weight->Fill(trig_weight);
        h.dilep_lepID_weight->Fill(id_weight);
        h.dilep_event_weight->Fill(_EventWeight);
        h.dilep_trig_weight_err->Fill(err_trig_weight);
        h.dilep_lepID_weight_err->Fill(err_id_weight);
        h.dilep_event_weight_err->Fill((_EventWeightErrorPlus+_EventWeightErrorMinus)/2.);
    }
    */
    if ( fVerbose>=2 ) {
        for ( int i=0; i<nleptons; i++ ) {
            cout << i << " id=" << lep[i].id() << " et=" << lep[i].et() << " trig_eff=" << trig_eff[i] << "+/-" << err_trig_eff[i] << " id_eff=" << id_eff[i] << "+/-" << err_id_eff[i] << endl;
        }
        cout << " nleptons=" << nleptons << " ngood_leptons=" << ngood_leptons << endl;
        cout << " trig_weight=" << trig_weight << "+/-" << err_trig_weight << endl;
        cout << " id_weight=" << id_weight << "+/-" << err_id_weight << endl;
        cout << " EventWeight=" << _EventWeight << " + " << _EventWeightErrorPlus << " - " << _EventWeightErrorMinus << endl;
    }
    return return_code;
}

float TSeel::get_Fake_Weight(float &x)
{


  float y=0;
  // loop over jets to find fakeable objects
  for(int ijet=0; ijet<fJetBlock->NJets(); ijet++) {
    TStnJet* jet = fJetBlock->Jet(ijet);
    
    //stop primary leptons in event being double-counted as fakeable objects
    float deltar=999;
    for( int ilep=0; ilep<nleptons; ilep++ ) {
      deltar = delta_r(lep[ilep].eta(),lep[ilep].phi(),jet->Eta(),jet->Phi());
      if(deltar<0.4)
	break;
    }

    //et of fake electrons (different from jet Et)
    float ettlep=-999;
    
    //fake rate for electrons
    float frt=0;
    
    if(deltar>0.4){
      //get tight elec fakes
      //frt = GetCentralElecTightFake(jet,ettlep);
      //get loose elec fakes
      frt = GetCentralElecLooseFake(jet,ettlep);
      if(ettlep>_LEP_PT_CUT3[2]){
	x+= frt;
	Nevents_FakeIncremented++;
      }
      if(frt>0){
	// Now apply the mass and optimization cuts and sum in y 
	// (this is fake event optimization weight)
	// Transverse Masses
	vector<float> mt_lep;
	mt_lep.push_back(transverse_mass_px_py(lep[0].px(),lep[0].py(),missing_et.ex(),missing_et.ey()));
	mt_lep.push_back(transverse_mass_px_py(lep[1].px(),lep[1].py(),missing_et.ex(),missing_et.ey()));
	mt_lep.push_back(transverse_mass_px_py(jet->Momentum()->Px(),jet->Momentum()->Py(),missing_et.ex(),missing_et.ey()));
	sort(mt_lep.begin(),mt_lep.end());
	mtmin = mt_lep[0];
	// Invariant Masses
	float wgt_jet_is_positive,wgt_jet_is_negative;
	wgt_jet_is_positive=wgt_jet_is_negative=0;
	vector<TLorentzVector> v_lep;
	for ( int i=0; i<2; i++ ) {
	  TLorentzVector v(lep[i].px(),lep[i].py(),lep[i].pz(),lep[i].en());
	  v_lep.push_back(v);
	}
	TLorentzVector v(jet->Momentum()->Px(),jet->Momentum()->Py(),jet->Momentum()->Pz(),jet->Momentum()->E());
	v_lep.push_back(v);
	vector<TLorentzVector> v_os_1,v_os_2,v_os_3;
	if ( lep[0].id()*lep[1].id()<0 ) {
	  v_os_1.push_back(v_lep[0]);
	  v_os_1.push_back(v_lep[1]);
	  v_os_2.push_back(v_lep[0]);
	  v_os_2.push_back(v_lep[2]);
	  v_os_3.push_back(v_lep[1]);
	  v_os_3.push_back(v_lep[2]);
	}	
	else {
	  v_os_2.push_back(v_lep[0]);
	  v_os_2.push_back(v_lep[2]);
	  v_os_3.push_back(v_lep[1]);
	  v_os_3.push_back(v_lep[2]);
	}	  
	if ( v_os_1.empty() ) mos1_fake=0;
	else mos1_fake = inv_mass(v_os_1); 
	
	if ( v_os_2.empty() ) mos2_fake=0;
	else mos2_fake = inv_mass(v_os_2);
	if(mos1_fake<mos2_fake){
	  float temp = mos1_fake; mos1_fake=mos2_fake; mos2_fake=temp;
	}
	//cout<<missing_et.et()<<" "<<
	//  mos1_fake<<" "<<mos2_fake<<" "<<mtmin<<" "<<
	//  delta_phi(lep[0].phi(),lep[1].phi())<<" "<<lep[0].et()<<
	//  " "<<lep[1].et()<<" "<<ettlep<<endl;
	if ( pass_opt_fake(mos1_fake,mos2_fake,mtmin,ettlep) )
	  wgt_jet_is_positive = frt;
	
	if ( v_os_3.empty() ) mos2_fake=0;
	else mos2_fake = inv_mass(v_os_3);
	if(mos1_fake<mos2_fake){
	  float temp = mos1_fake; mos1_fake=mos2_fake; mos2_fake=temp;
	}
	if ( pass_opt_fake(mos1_fake,mos2_fake,mtmin,ettlep) )
	  wgt_jet_is_negative = frt;
	
	y += ((wgt_jet_is_positive+wgt_jet_is_negative)/2);
	if( (wgt_jet_is_positive+wgt_jet_is_negative)/2 > 0 )
	  Nevents_FakeIncrementedOpt++;
	// Also apply the various control region cuts and sum directly in
	// _fakeEventWeightControlRegion[iRegion] (this is the control region
	// fake contribution

      }// if positive fake rate
      
    }
  }
  
  //loop over tracks to find fakeable objects
  for(int itrk=0; itrk<fTrackBlock->NTracks(); itrk++) {
    TStnTrack* trk = fTrackBlock->Track(itrk);
    TLorentzVector* mom = new TLorentzVector();
    int status = trk->GetBcMomentum(mom);
    //stop primary leptons in event being double-counted as fakeable objects
    float deltar=999;
    for( int ilep=0; ilep<nleptons; ilep++ ) {
      deltar = delta_r(lep[ilep].eta(),lep[ilep].phi(),mom->Eta(),mom->Phi());
      if(deltar<0.4)
	break;
    }

    //pt of fake muons
    float ptcmuplep=-999;float ptcmxlep=-999;float ptcmiolep=-999;
    
    //fake rate for muons
    float frcmup=0; float frcmx=0; float frcmio=0;
    if(deltar>0.4 && CorrectedD0(trk)<0.02){
      //get CMUP fakes
      frcmup = GetCMUPFake(trk,fCalDataBlock,ptcmuplep);
      if(ptcmuplep>_LEP_PT_CUT3[2]){
	x+= frcmup;
	Nevents_FakeIncremented++;
      }
      //get CMX fakes
      frcmx = GetCMXFake(trk,fCalDataBlock,ptcmxlep);
      if(ptcmxlep>_LEP_PT_CUT3[2]){
	x+= frcmx;
	Nevents_FakeIncremented++;
      }
      //get CMIO fakes
      frcmio = GetCMIOFake(trk,fCalDataBlock,ptcmiolep);
      if(ptcmiolep>_LEP_PT_CUT3[2]){
	x+= frcmio;
	Nevents_FakeIncremented++;
      }
      if(frcmup>0 || frcmx>0 || frcmio>0){
	// Now apply the mass and optimization cuts and sum in y 
	// (this is fake event optimization weight)
	// Transverse Masses
	vector<float> mt_lep;
	mt_lep.push_back(transverse_mass_px_py(lep[0].px(),lep[0].py(),missing_et.ex(),missing_et.ey()));
	mt_lep.push_back(transverse_mass_px_py(lep[1].px(),lep[1].py(),missing_et.ex(),missing_et.ey()));
	mt_lep.push_back(transverse_mass_px_py(trk->Momentum()->Px(),trk->Momentum()->Py(),missing_et.ex(),missing_et.ey()));
	sort(mt_lep.begin(),mt_lep.end());
	mtmin = mt_lep[0];
	// Invariant Masses
	vector<TLorentzVector> v_lep;
	for ( int i=0; i<2; i++ ) {
	  TLorentzVector v(lep[i].px(),lep[i].py(),lep[i].pz(),lep[i].en());
	  v_lep.push_back(v);
	}
	TLorentzVector v(trk->Momentum()->Px(),trk->Momentum()->Py(),trk->Momentum()->Pz(),trk->Momentum()->E());
	v_lep.push_back(v);
	vector<TLorentzVector> v_os_1,v_os_2;     
	if ( lep[0].id()*lep[1].id()<0 ) {
	  v_os_1.push_back(v_lep[0]);
	  v_os_1.push_back(v_lep[1]);
	  if ( (lep[0].id()/-11)*trk->Charge()<0 ){
	    v_os_2.push_back(v_lep[0]);
	    v_os_2.push_back(v_lep[2]);
	  }
	  else if ( (lep[1].id()/-11)*trk->Charge()<0 ){
	    v_os_2.push_back(v_lep[1]);
	    v_os_2.push_back(v_lep[2]);
	  }
	}
	if ( v_os_1.empty() ) mos1_fake=0;
	else mos1_fake = inv_mass(v_os_1);
	if ( v_os_2.empty() ) mos2_fake=0;
	else mos2_fake = inv_mass(v_os_2);

	if ( pass_opt_fake(mos1_fake,mos2_fake,mtmin,ptcmuplep) ){
	  y += frcmup;
	  Nevents_FakeIncrementedOpt++;
	}
	if ( pass_opt_fake(mos1_fake,mos2_fake,mtmin,ptcmxlep) ){
	  y += frcmx;
	  Nevents_FakeIncrementedOpt++;
	}
	if ( pass_opt_fake(mos1_fake,mos2_fake,mtmin,ptcmiolep) ){
	  y += frcmio;	
	  Nevents_FakeIncrementedOpt++;
	}
      }// if positive fake rate
    }
  }

  return y;

}
//convert jet Et -> electron Et
float TSeel::Smear(float mean,float sigma){
//-------------------------------------------------------
//
// this function does a random gaussian smearing around the "mean"
// with a sigma of "sigma" and returns the smeared value
//
//--------------------------------------------------------

  float smear=1.0;
  float t1=0.0;
  float s1=0.0;
  float t2=0.0;
  float rmax=RAND_MAX;
  for(int i=0;i<100 ;i++){
    t1 = 0.0;
    t1 = 3.*sigma*(2.0*rand()/rmax -1.) + mean;
    
    s1 = exp(-pow((t1-mean),2)/(2.*pow(sigma,2)));
    t2 = 0.0;
    t2 = rand()/rmax;
    if(t2<s1){
      smear=t1;
      return smear;
    }
  }
   return mean;
}
bool TSeel::isFakeableCentralJet(TStnJet* jet) {
  return (fabs(jet->Momentum()->Eta())<1.1 && jet->Et()>4.);
}
bool TSeel::isFakeableCMUPTrack(TStnTrack* trk, TCalDataBlock* fCalDataBlock, float eoverpcut){
  return (trk->IsCMUFid() && trk->IsCMPFid() && isFakeableTrack(trk, fCalDataBlock, eoverpcut));
}

bool TSeel::isFakeableCMXTrack(TStnTrack* trk, TCalDataBlock* fCalDataBlock, float eoverpcut){
  return (trk->IsCMXFid() && isFakeableTrack(trk, fCalDataBlock, eoverpcut));
}
bool TSeel::isFakeableCMIOTrack(TStnTrack* trk, TCalDataBlock* fCalDataBlock, float eoverpcut){
  bool fid = (!(trk->IsCMUFid() && trk->IsCMPFid()) && !trk->IsCMXFid());
  int nCotSt = 0;
  for(int i=0; i<8; i+=2) {
    if(trk->NCotHits(i)>=5) nCotSt++;
  }
  bool stsegcut = nCotSt>2;
  bool ptcut = trk->Pt()>10;
  return (fid && stsegcut && ptcut && isFakeableTrack(trk, fCalDataBlock, eoverpcut));
}
float TSeel::GetCentralElecTightFake(TStnJet* jet, float &etlep){
  if(isFakeableCentralJet(jet)) {      
    //printf("hi\n");
    float fakerate=_centralelectightfakerate->Eval(jet->Et());
    //printf("evaluated fakerate\n");
    //map jet Et to fake lepton Et
    etlep = Smear(jet->Et()*_smearMean,jet->Et()*_smearRMS);
    return fakerate;
  }
  return 0;
}
//main function to return central loose electron fakes
float TSeel::GetCentralElecLooseFake(TStnJet* jet, float &etlep){
  if(isFakeableCentralJet(jet)) {      
    float fakerate=_centralelecloosefakerate->Eval(jet->Et());
    //cout<<"fakerate = "<<fakerate<<endl;
    etlep = Smear(jet->Et()*_smearMean,jet->Et()*_smearRMS);
    return fakerate;
  }
  return 0;
}
//main functions to return muon fakes in each category
float TSeel::GetCMUPFake(TStnTrack* trk, TCalDataBlock* fCalDataBlock, float &etlep){
  if(isFakeableCMUPTrack(trk, fCalDataBlock)){
    float fakerate=_CMUPfakerate->Eval(trk->Pt());
    etlep = trk->Pt();
    return fakerate;
  }
  return 0;
}

float TSeel::GetCMXFake(TStnTrack* trk, TCalDataBlock* fCalDataBlock, float &etlep){
  if(isFakeableCMXTrack(trk, fCalDataBlock)){
    float fakerate=_CMXfakerate->Eval(trk->Pt());
    etlep = trk->Pt();
    return fakerate;
  }
  return 0;
}
float TSeel::GetCMIOFake(TStnTrack* trk, TCalDataBlock* fCalDataBlock, float &etlep){
  if(isFakeableCMIOTrack(trk, fCalDataBlock)){
    float fakerate=_CMIOfakerate->Eval(trk->Pt());
    etlep = trk->Pt();
    return fakerate;
  }
  return 0;
}
//check track cuts common to all categories 
bool TSeel::isFakeableTrack(TStnTrack* trk, TCalDataBlock* fCalDataBlock, float eoverpcut){

  //track quality cuts
  bool d0cut = true; //have to apply d0cut outside this function at present because function doesn't know about the correct beam position
  
  bool axsegcut = trk->NCotAxSeg(5)>2;
  int nCotSt = 0;
  for(int i=0; i<8; i+=2) {
    if(trk->NCotHits(i)>=5) nCotSt++;
  }
  bool stsegcut = nCotSt>1;
  bool z0cut = fabs(trk->Z0())<60;
  


  if (trk->Pt()>4 && axsegcut && stsegcut && z0cut && d0cut){


    
    float iso4 = TStntuple::CalIso(fCalDataBlock,trk->Momentum()->Eta(),trk->Phi0(),trk->Z0(),0.4);

    float sec_xCes = -999.;
    float sec_zCes = -999.;
    int   ieta = -999;
    int   iphi = -999;
      
    Double_t sec_gxyz[3];
    Double_t sec_mom[3];
    Double_t sec_mtot = -999.;

    int sec_status = GetCesInfo(trk,sec_xCes,sec_zCes,ieta,iphi,sec_gxyz,sec_mom,sec_mtot);
    if(sec_status!=0) return false;



    // reject electrons by cutting on E/p in that tower
    float eoverp=0;float ettow=0;
    if (fCalDataBlock->Tower(ieta, iphi)) {
      TCalTower* to = fCalDataBlock->Tower(ieta,iphi);
      ettow = to->Et();
      eoverp= to->Energy()/trk->Momentum()->P();
    } 
    else {
      eoverp=0;
    }
    iso4-=ettow;
    // cut at isolation cut of <4 GeV
    if(eoverp<eoverpcut && (iso4<4 || iso4/trk->Pt()<0.1)){

      return true;
    }
  }
  return false;
}
int TSeel::GetCesInfo(TStnTrack* track, float &xCes, float &zCes, int &ieta, int &iphi, 
			      Double_t *gxyz, Double_t *mom, Double_t &mtot) {

  TSimpleExtrapolator* fExtrapolator = new TSimpleExtrapolator();
  TTrajectoryPoint p0;
  double xyz[8], xw, zw;
  int side, wedge;
  
  xyz[0] = -track->D0() * TMath::Sin(track->Phi0());
  xyz[1] =  track->D0() * TMath::Cos(track->Phi0());
  xyz[2] =  track->Z0();
  xyz[3] =  track->Momentum()->Px()/track->Momentum()->P();
  xyz[4] =  track->Momentum()->Py()/track->Momentum()->P();
  xyz[5] =  track->Momentum()->Pz()/track->Momentum()->P();
  xyz[6] =  0;
  xyz[7] =  track->Momentum()->P();
  
  p0.SetPoint(xyz);
      
  int status = fExtrapolator->SwimToCes(&p0,track->Charge(),side,wedge,xw,zw);
  //  if(status!=0) return status; //you don't delete fExtrapolator

  fExtrapolator->GetFinalResults(gxyz,mom,mtot);
  fExtrapolator->GetIEtaIPhi(side,wedge,zw,ieta,iphi);

  if(status==0) {
    if(side == 0) zCes = -zw;
    else          zCes =  zw;
    xCes = xw;
  }

  delete fExtrapolator;

  return status;
}

void TSeel::FillHistograms()
{
// Examination Variables
//h.beam_xy->Fill(fEvtHdr_fVx,fEvtHdr_fVy);
    for ( int i=0; i<njet; i++ ) {
        h.jet_scaleEta->Fill(jt[i].eta(),jt[i].scale5());
        h.jet20_et->Fill(jt[i].et());
        h.jet20_eta->Fill(jt[i].eta());
        h.jet20_phi->Fill(jt[i].phi());
    }

    if ( nleptons>=2 ) {
        for ( int i=0; i<1; i++ ) {
            h.dilep_pt[i]->Fill(lep[i].et());
            h.dilep_eta[i]->Fill(lep[i].eta());
            h.dilep_phi[i]->Fill(lep[i].phi());
        }
    }
    if ( nleptons>=3 ) {
//  Trilepton kinematics
        for ( int i=0; i<3; i++ ) {
            h.trilep_pt[i]->Fill(lep[i].et());
            h.trilep_eta[i]->Fill(lep[i].eta());
            h.trilep_phi[i]->Fill(lep[i].phi());
        }
        h.vz_lep1->Fill(vz-lep[0].z0());
        h.vz_lep2->Fill(vz-lep[1].z0());
        h.vz_lep3->Fill(vz-lep[2].z0());
//  st
//        float st = ht;
        float st = 0.;
        for ( int i=0; i<3; i++ ) st += lep[i].et();
        float ht2 = 0.; if ( njet>0 ) ht2 = ht - jt[0].et();
        float st2 = 0.; if ( njet>0 ) st2 = st - jt[0].et();
    
        if ( getbits(_EventType,10,1) ) { // correct charge combination
            h.jet20_n->Fill(njet);

            h.met->Fill(missing_et.et());
            h.met_phi->Fill(missing_et.phi());
            h.st->Fill(st);
            h.st2->Fill(st2);
            h.st2_over_st->Fill((st>0. ? st2/st : 0.) );
            h.ht->Fill(ht);
            h.ht2->Fill(ht2);
            h.ht2_over_ht->Fill((ht>0. ? ht2/ht : 0.) );
            h.trk_met->Fill(lep[2].pt()*missing_et.et());
            h.st_met->Fill(st+missing_et.et());
        }


//  invariant masses
        vector<TLorentzVector> v_lep;
        for ( int i=0; i<3; i++ ) {
            TLorentzVector v(lep[i].px(),lep[i].py(),lep[i].pz(),lep[i].en());
            v_lep.push_back(v);
        }
        float m_leps = inv_mass(v_lep);
        vector<TLorentzVector> v_leps_jets = v_lep;
        for ( int i=0; i<njet; i++ ) {
            TLorentzVector v(jt[i].ex(),jt[i].ey(),jt[i].ez(),jt[i].en());
            v_leps_jets.push_back(v);
        }
        float m_leps_jets = inv_mass(v_leps_jets);
        vector<TLorentzVector> v_all = v_leps_jets;
        TLorentzVector v(missing_et.ex(),missing_et.ey(),0.,missing_et.et());
        v_all.push_back(v);
        float m_all = inv_mass(v_all);
        float mt_all = transverse_mass(v_all);

        vector<TLorentzVector> v_os_1,v_os_2,v_ss,v_ss_1,v_ss_2,v_ss_3;
        if ( lep[0].id()*lep[1].id()<0 ) {
            v_os_1.push_back(v_lep[0]);
            v_os_1.push_back(v_lep[1]);
            if ( lep[0].id()*lep[2].id()<0 ) {
                v_os_2.push_back(v_lep[0]);
                v_os_2.push_back(v_lep[2]);
            }
            else if ( lep[1].id()*lep[2].id()<0 ) {
                v_os_2.push_back(v_lep[1]);
                v_os_2.push_back(v_lep[2]);
            }
            if ( lep[0].id()*lep[2].id()>0 ) {
                v_ss.push_back(v_lep[0]);
                v_ss.push_back(v_lep[2]);
            }
            else if ( lep[1].id()*lep[2].id()>0 ) {
                v_ss.push_back(v_lep[1]);
                v_ss.push_back(v_lep[2]);
            }
        }
        else if ( lep[0].id()*lep[2].id()<0 ) {
            v_os_1.push_back(v_lep[0]);
            v_os_1.push_back(v_lep[2]);
            v_os_2.push_back(v_lep[1]);
            v_os_2.push_back(v_lep[2]);
            v_ss.push_back(v_lep[0]);
            v_ss.push_back(v_lep[1]);
        }
        else { // all three leptons same sign
            v_ss_1.push_back(v_lep[0]); v_ss_1.push_back(v_lep[1]);
            v_ss_2.push_back(v_lep[0]); v_ss_2.push_back(v_lep[2]);
            v_ss_3.push_back(v_lep[1]); v_ss_3.push_back(v_lep[2]);
        }
//        m_os_1 = inv_mass(v_os_1);
//        m_os_2 = inv_mass(v_os_2);
//        if ( m_os_1 < m_os_2 ) {
//            float temp = m_os_1;
//            m_os_1 = m_os_2;
//            m_os_2 = temp;
//        }

        vector<float> m_ss;
        m_ss.push_back(inv_mass(v_ss_1));
        m_ss.push_back(inv_mass(v_ss_2));
        m_ss.push_back(inv_mass(v_ss_3));
        sort(m_ss.begin(), m_ss.end());
        if ( getbits(_EventType,10,1) ) { // correct charge combination
            h.m_leps->Fill(m_leps);
            h.m_leps_jets->Fill(m_leps_jets);
            h.m_all->Fill(m_all);
            h.m_os_1->Fill(m_os_1);
            h.m_os_2->Fill(m_os_2);
            h.m_ss->Fill(inv_mass(v_ss));
            h.mt_all->Fill(mt_all);
            h.met_m_os_1->Fill(m_os_1,missing_et.et());
        }
        else {
            h.m_ss_1->Fill(m_ss[2]);
            h.m_ss_2->Fill(m_ss[1]);
            h.m_ss_3->Fill(m_ss[0]);
        }

//  Transverse masses
        vector<float> mt_lep;
        mt_lep.push_back(transverse_mass_px_py(lep[0].px(),lep[0].py(),missing_et.ex(),missing_et.ey()));
        mt_lep.push_back(transverse_mass_px_py(lep[1].px(),lep[1].py(),missing_et.ex(),missing_et.ey()));
        mt_lep.push_back(transverse_mass_px_py(lep[2].px(),lep[2].py(),missing_et.ex(),missing_et.ey()));
        if ( getbits(_EventType,10,1) ) { // correct charge combination
            h.mt_lep1_met->Fill(mt_lep[0]);
            h.mt_lep2_met->Fill(mt_lep[1]);
            h.mt_lep3_met->Fill(mt_lep[2]);
            sort(mt_lep.begin(),mt_lep.end());
            h.mt_lep_met_max->Fill(mt_lep[2]);
            h.mt_lep_met_min->Fill(mt_lep[0]);
        }

//  Sphericity, Alpanarity, and Thrust
        float sphe=0.; float apla=0.;
        if ( getbits(_EventType,10,1) ) { // correct charge combination
            float thrust = Thrust(v_all); h.thrust->Fill(thrust);
            if ( Sphericity_Aplanarity(v_all, sphe, apla) ) {
                h.sphericity->Fill(sphe);
                h.aplanarity->Fill(apla);
            }
            else {
                h.sphericity->Fill(-999.);
                h.aplanarity->Fill(-999.);
            }
        }
        vector<TLorentzVector> vtrans_all = v_all;
        for ( vector<TLorentzVector>::iterator i=vtrans_all.begin(); i!=vtrans_all.end(); i++ ) i->SetPz(0.);
        float sphe_2d=0.; float apla_2d=0.;
        if ( getbits(_EventType,10,1) ) { // correct charge combination
            float thrust_2d = Thrust(vtrans_all); h.thrust_2d->Fill(thrust_2d);
            if ( Sphericity_Aplanarity(vtrans_all, sphe_2d, apla_2d) ) {
                h.sphericity_2d->Fill(sphe_2d);
                h.aplanarity_2d->Fill(apla_2d);
            }
            else {
                h.sphericity->Fill(-999.);
                h.aplanarity->Fill(-999.);
            }
        }
//  Object Separations
        h.dphi_lep12->Fill(delta_phi(lep[0].phi(),lep[1].phi()));
        h.dphi_lep13->Fill(delta_phi(lep[0].phi(),lep[2].phi()));
        h.dphi_lep23->Fill(delta_phi(lep[1].phi(),lep[2].phi()));
        h.dphi_lep1met->Fill(delta_phi(lep[0].phi(),missing_et.phi()));
        h.dphi_lep2met->Fill(delta_phi(lep[1].phi(),missing_et.phi()));
        h.dphi_lep3met->Fill(delta_phi(lep[2].phi(),missing_et.phi()));
        h.dr_lep12->Fill(delta_r(lep[0].eta(),lep[0].phi(),lep[1].eta(),lep[1].phi()));
        h.dr_lep13->Fill(delta_r(lep[0].eta(),lep[0].phi(),lep[2].eta(),lep[2].phi()));
        h.dr_lep23->Fill(delta_r(lep[1].eta(),lep[1].phi(),lep[2].eta(),lep[2].phi()));
	if ( njet15>0 ) {
	    h.dphi_j1met_j15->Fill(delta_phi(jt15[0].phi(),missing_et.phi()));
	}
        if ( njet>0 ) {
            h.dphi_j1met->Fill(delta_phi(jt[0].phi(),missing_et.phi()));
            float drmin[3], dr[3]; for (int i=0; i<3; i++) drmin[i]=999.;
            for ( int i=0; i<3; i++ ) {
                for ( int j=0; j<njet; j++ ) {
                    dr[i] = delta_r(jt[j].eta(),jt[j].phi(),lep[i].eta(),lep[i].phi());
                    if ( dr[i]<drmin[i] ) drmin[i]=dr[i];
                }
            }
            h.dr_jt_lep1->Fill(drmin[0]);
            h.dr_jt_lep2->Fill(drmin[1]);
            h.dr_jt_lep3->Fill(drmin[2]);
            if ( njet>1 ) h.dphi_j2met->Fill(delta_phi(jt[1].phi(),missing_et.phi()));
        }
        h.met_trkpt->Fill(lep[2].pt(),missing_et.et());
        opt_met = missing_et.et();
        opt_st = st;
        opt_m_os_1 = m_os_1;
        opt_m_os_2 = m_os_2;
        opt_dphi_lep12 = delta_phi(lep[0].phi(),lep[1].phi());
        opt_ht = ht;
        opt_trk_met = lep[2].pt()*missing_et.et();
        opt_mt_min = mt_lep[0];

	opt_ele1_et = lep[0].et();
	opt_ele2_et = lep[1].et();
	opt_trk_pt  = lep[2].et();
    } // nleptons>=3
}

bool TSeel::pass_ControlRegionCuts(int iRegion)
{
// iRegion=0: c1 - MET<10 & Diele + track  --- This is implemented already.
// iRegion=1: c2 - MET>10 & Diele + track & 76<Mee<106
// iRegion=2: c3 - MET<10 & Diele + track & 2 jets
// iRegion=3: c4 - MET>10 & Diele + track & 2 jets & 76<Mee<106
// iRegion=4: veto events HiPt analysis would have selected
    vector<TLorentzVector> v_lep;
    for ( int i=0; i<2; i++ ) {
        TLorentzVector v(lep[i].px(),lep[i].py(),lep[i].pz(),lep[i].en());
        v_lep.push_back(v);
    }
    float m_leps = inv_mass(v_lep);
    if ( iRegion==0 )
        return ( Pass_Trilep_BaseCuts && Pass_Trilep_MassCuts && missing_et.et()<10. );
    else if ( iRegion==1 )
        return ( Pass_Trilep_BaseCuts && Pass_Trilep_MassCuts && missing_et.et()>10. && m_leps>76. && m_leps<106. );
    else if ( iRegion==2 )
        return ( Pass_Trilep_BaseCuts && Pass_Trilep_MassCuts && missing_et.et()<10. && njet15>=2 );
    else if ( iRegion==3 )
        return ( Pass_Trilep_BaseCuts && Pass_Trilep_MassCuts && missing_et.et()>10. && m_leps>76. && m_leps<106. && njet15>=2 );
    else if ( iRegion==4 ) {
//  HiPt cuts that are not part of my cuts:
//  Tight 20 GeV electron (mine already tight, just need to look at Et)
//  Loose 8 GeV electron (mine already tight, just need to look at Et)
//  Third is in Table 3, on p9 of CDF7499.
//  reject events with two 20 GeV jets.
//  To make our analysis orthogonal, we veto the event if any of the above is
//  true, while
        bool good_3rd_lepton = false;
        for ( int i=2; i<nleptons; i++ ) {
            if ( lep[i].pt()<5. ) continue;
            if ( abs(lep[i].id())==11 && lep[i].et()>5. )
                good_3rd_lepton = true;
            else if ( getbits(lep[i].qual(),13,1) ) {  // an CDF Muon
                good_3rd_lepton = true;
	    }
	}
        
        int njet20 = 0;
        for ( int i=0; i<njet; i++ ) {
            if ( jt[i].et()>20. )
                njet20++;
        }

        bool pass_HiPt_cuts = false;
        if ( lep[0].et()>20.
             && lep[1].et()>8.
             && delta_phi(lep[0].phi(),lep[1].phi())*180/M_PI<160.
             && good_3rd_lepton
             && njet20<2 )
            pass_HiPt_cuts = true;
        return pass_HiPt_cuts;
    }
    else // no such control region
        return false;
}

void TSeel::FillControlRegionHistograms(int iRegion)
{
// ====================================================================
//  ControlRegion variables
//  Control Regions:
// c1 - MET<10 & Diele + track  --- This is implemented already.
// c2 - MET>10 & Diele + track & 76<Mee<106
// c3 - MET<10 & Diele + track & 2 jets
// c4 - MET>10 & Diele + track & 2 jets & 76<Mee<106
// ---------------------------------------
    if ( Pass_Trilep_BaseCuts && Pass_Trilep_MassCuts ) {
//  Trilepton kinematics
        for ( int i=0; i<3; i++ ) {
            h.c_trilep_pt[iRegion][i]->Fill(lep[i].et());
            h.c_trilep_eta[iRegion][i]->Fill(lep[i].eta());
            h.c_trilep_phi[iRegion][i]->Fill(lep[i].phi());
        }
//  st
//        float st = ht;
        float st = 0.;
        for ( int i=0; i<3; i++ ) st += lep[i].et();
        float ht2 = 0.; if ( njet>0 ) ht2 = ht - jt[0].et();
        float st2 = 0.; if ( njet>0 ) st2 = st - jt[0].et();
    
        if ( getbits(_EventType,10,1) ) { // correct charge combination
            h.c_met[iRegion]->Fill(missing_et.et());
            h.c_st[iRegion]->Fill(st);
            h.c_st2[iRegion]->Fill(st2);
            h.c_st2_over_st[iRegion]->Fill((st>0. ? st2/st : 0.) );
            h.c_ht[iRegion]->Fill(ht);
            h.c_ht2[iRegion]->Fill(ht2);
            h.c_ht2_over_ht[iRegion]->Fill((ht>0. ? ht2/ht : 0.) );
            h.c_trk_met[iRegion]->Fill(lep[2].pt()*missing_et.et());
            h.c_st_met[iRegion]->Fill(st+missing_et.et());
            h.c_njet15[iRegion]->Fill(njet15);
        }

//  invariant masses
        vector<TLorentzVector> v_lep;
        for ( int i=0; i<3; i++ ) {
            TLorentzVector v(lep[i].px(),lep[i].py(),lep[i].pz(),lep[i].en());
            v_lep.push_back(v);
        }
        float m_leps = inv_mass(v_lep);
        vector<TLorentzVector> v_leps_jets = v_lep;
        for ( int i=0; i<njet; i++ ) {
            TLorentzVector v(jt[i].ex(),jt[i].ey(),jt[i].ez(),jt[i].en());
            v_leps_jets.push_back(v);
        }
        float m_leps_jets = inv_mass(v_leps_jets);
        vector<TLorentzVector> v_all = v_leps_jets;
        TLorentzVector v(missing_et.ex(),missing_et.ey(),0.,missing_et.et());
        v_all.push_back(v);
        float m_all = inv_mass(v_all);
        float mt_all = transverse_mass(v_all);

        vector<TLorentzVector> v_os_1,v_os_2,v_ss,v_ss_1,v_ss_2,v_ss_3;
        if ( lep[0].id()*lep[1].id()<0 ) {
            v_os_1.push_back(v_lep[0]);
            v_os_1.push_back(v_lep[1]);
            if ( lep[0].id()*lep[2].id()<0 ) {
                v_os_2.push_back(v_lep[0]);
                v_os_2.push_back(v_lep[2]);
            }
            else if ( lep[1].id()*lep[2].id()<0 ) {
                v_os_2.push_back(v_lep[1]);
                v_os_2.push_back(v_lep[2]);
            }
            if ( lep[0].id()*lep[2].id()>0 ) {
                v_ss.push_back(v_lep[0]);
                v_ss.push_back(v_lep[2]);
            }
            else if ( lep[1].id()*lep[2].id()>0 ) {
                v_ss.push_back(v_lep[1]);
                v_ss.push_back(v_lep[2]);
            }
        }
        else if ( lep[0].id()*lep[2].id()<0 ) {
            v_os_1.push_back(v_lep[0]);
            v_os_1.push_back(v_lep[2]);
            v_os_2.push_back(v_lep[1]);
            v_os_2.push_back(v_lep[2]);
            v_ss.push_back(v_lep[0]);
            v_ss.push_back(v_lep[1]);
        }
        else { // all three leptons same sign
            v_ss_1.push_back(v_lep[0]); v_ss_1.push_back(v_lep[1]);
            v_ss_2.push_back(v_lep[0]); v_ss_2.push_back(v_lep[2]);
            v_ss_3.push_back(v_lep[1]); v_ss_3.push_back(v_lep[2]);
        }
//        m_os_1 = inv_mass(v_os_1);
//        m_os_2 = inv_mass(v_os_2);
//        if ( m_os_1 < m_os_2 ) {
//            float temp = m_os_1;
//            m_os_1 = m_os_2;
//            m_os_2 = temp;
//        }
        vector<float> m_ss;
        m_ss.push_back(inv_mass(v_ss_1));
        m_ss.push_back(inv_mass(v_ss_2));
        m_ss.push_back(inv_mass(v_ss_3));
        sort(m_ss.begin(), m_ss.end());
        if ( getbits(_EventType,10,1) ) { // correct charge combination
            h.c_m_leps[iRegion]->Fill(m_leps);
            h.c_m_leps_jets[iRegion]->Fill(m_leps_jets);
            h.c_m_all[iRegion]->Fill(m_all);
            h.c_m_os_1[iRegion]->Fill(m_os_1);
            h.c_m_os_2[iRegion]->Fill(m_os_2);
            h.c_m_ss[iRegion]->Fill(inv_mass(v_ss));
            h.c_mt_all[iRegion]->Fill(mt_all);
            h.c_met_m_os_1[iRegion]->Fill(m_os_1,missing_et.et());
        }
        else {
            h.c_m_ss_1[iRegion]->Fill(m_ss[2]);
            h.c_m_ss_2[iRegion]->Fill(m_ss[1]);
            h.c_m_ss_3[iRegion]->Fill(m_ss[0]);
        }

//  Transverse masses
        vector<float> mt_lep;
        mt_lep.push_back(transverse_mass_px_py(lep[0].px(),lep[0].py(),missing_et.ex(),missing_et.ey()));
        mt_lep.push_back(transverse_mass_px_py(lep[1].px(),lep[1].py(),missing_et.ex(),missing_et.ey()));
        mt_lep.push_back(transverse_mass_px_py(lep[2].px(),lep[2].py(),missing_et.ex(),missing_et.ey()));
        if ( getbits(_EventType,10,1) ) { // correct charge combination
            h.c_mt_lep1_met[iRegion]->Fill(mt_lep[0]);
            h.c_mt_lep2_met[iRegion]->Fill(mt_lep[1]);
            h.c_mt_lep3_met[iRegion]->Fill(mt_lep[2]);
            sort(mt_lep.begin(),mt_lep.end());
            h.c_mt_lep_met_max[iRegion]->Fill(mt_lep[2]);
            h.c_mt_lep_met_min[iRegion]->Fill(mt_lep[0]);
        }

//  Sphericity, Alpanarity, and Thrust
        float sphe=0.; float apla=0.;
        if ( getbits(_EventType,10,1) ) { // correct charge combination
            float thrust = Thrust(v_all); h.c_thrust[iRegion]->Fill(thrust);
            if ( Sphericity_Aplanarity(v_all, sphe, apla) ) {
                h.c_sphericity[iRegion]->Fill(sphe);
                h.c_aplanarity[iRegion]->Fill(apla);
            }
            else {
                h.c_sphericity[iRegion]->Fill(-999.);
                h.c_aplanarity[iRegion]->Fill(-999.);
            }
        }
        vector<TLorentzVector> vtrans_all = v_all;
        for ( vector<TLorentzVector>::iterator i=vtrans_all.begin(); i!=vtrans_all.end(); i++ ) i->SetPz(0.);
        float sphe_2d=0.; float apla_2d=0.;
        if ( getbits(_EventType,10,1) ) { // correct charge combination
            float thrust_2d = Thrust(vtrans_all); h.c_thrust_2d[iRegion]->Fill(thrust_2d);
            if ( Sphericity_Aplanarity(vtrans_all, sphe_2d, apla_2d) ) {
                h.c_sphericity_2d[iRegion]->Fill(sphe_2d);
                h.c_aplanarity_2d[iRegion]->Fill(apla_2d);
            }
            else {
                h.c_sphericity[iRegion]->Fill(-999.);
                h.c_aplanarity[iRegion]->Fill(-999.);
            }
        }
//  Object Separations
        h.c_dphi_lep12[iRegion]->Fill(delta_phi(lep[0].phi(),lep[1].phi()));
        h.c_dphi_lep13[iRegion]->Fill(delta_phi(lep[0].phi(),lep[2].phi()));
        h.c_dphi_lep23[iRegion]->Fill(delta_phi(lep[1].phi(),lep[2].phi()));
        h.c_dphi_lep1met[iRegion]->Fill(delta_phi(lep[0].phi(),missing_et.phi()));
        h.c_dphi_lep2met[iRegion]->Fill(delta_phi(lep[1].phi(),missing_et.phi()));
        h.c_dphi_lep3met[iRegion]->Fill(delta_phi(lep[2].phi(),missing_et.phi()));
        h.c_dr_lep12[iRegion]->Fill(delta_r(lep[0].eta(),lep[0].phi(),lep[1].eta(),lep[1].phi()));
        h.c_dr_lep13[iRegion]->Fill(delta_r(lep[0].eta(),lep[0].phi(),lep[2].eta(),lep[2].phi()));
        h.c_dr_lep23[iRegion]->Fill(delta_r(lep[1].eta(),lep[1].phi(),lep[2].eta(),lep[2].phi()));
        if ( njet15>0 ) {
            h.c_dphi_j1met[iRegion]->Fill(delta_phi(jt15[0].phi(),missing_et.phi()));
            float drmin[3], dr[3]; for (int i=0; i<3; i++) drmin[i]=999.;
            for ( int i=0; i<3; i++ ) {
                for ( int j=0; j<njet15; j++ ) {
                    dr[i] = delta_r(jt15[j].eta(),jt15[j].phi(),lep[i].eta(),lep[i].phi());
                    if ( dr[i]<drmin[i] ) drmin[i]=dr[i];
                }
            }
            h.c_dr_jt_lep1[iRegion]->Fill(drmin[0]);
            h.c_dr_jt_lep2[iRegion]->Fill(drmin[1]);
            h.c_dr_jt_lep3[iRegion]->Fill(drmin[2]);
            if ( njet15>1 ) h.c_dphi_j2met[iRegion]->Fill(delta_phi(jt15[1].phi(),missing_et.phi()));
        }
    } // Pass_Trilep_BaseCuts && Pass_Trilep_MassCuts
}

bool TSeel::pass_opt_fake(float mos1_fake,float mos2_fake,float mtmin,float et3)
{
//  met, m_os_1, m_os_2, ht, dphi_lep12
//  et1, et2,  et3

  // Raw distributions
  h.fopt_met->Fill(missing_et.et());
  h.fopt_mos1->Fill(mos1_fake);
  h.fopt_mos2->Fill(mos2_fake);
  h.fopt_ht->Fill(ht);
  h.fopt_dphi->Fill(delta_phi(lep[0].phi(),lep[1].phi()));
  h.fopt_mt->Fill(mtmin);
  h.fopt_et[0]->Fill(lep[0].et());
  h.fopt_et[1]->Fill(lep[1].et());
  h.fopt_et[2]->Fill(et3);

    int njet20 = 0;
    for ( int i=0; i<njet; i++ ) {
      if ( jt[i].et()>20. ) njet20++;
    }
    bool cut2[10];
    if ( _OptScheme==1 ) {
      cut2[0] = missing_et.et()>15.;
      cut2[1] = mos1_fake<76. || mos1_fake>106.;
      cut2[2] = mos2_fake<60.;
      cut2[3] = ht<80;
      cut2[4] = delta_phi(lep[0].phi(),lep[1].phi())<2.9;
      cut2[5] = true;
      cut2[6] = mtmin>10;
      cut2[7] = lep[0].et()>10;
      cut2[8] = lep[1].et()>5;
      cut2[9] = et3>4;
    }
    else if ( _OptScheme==3 ){
      cut2[0] = missing_et.et()>15.;
      cut2[1] = mos1_fake<76. || mos1_fake>106.;
      cut2[2] = true;//mos2_fake<60.;
      cut2[3] = njet20<2;//ht<80;
      cut2[4] = delta_phi(lep[0].phi(),lep[1].phi())<2.792527;
      cut2[5] = true;
      cut2[6] = true;//mtmin>10;
      cut2[7] = lep[0].et()>20 && lep[0].pt()>10;
      cut2[8] = lep[1].et()>8 && lep[1].pt()>8;
      cut2[9] = et3>5;
    }
    else if ( _OptScheme==9 ) {
      cut2[0] = missing_et.et()>20.;
      cut2[1] = mos1_fake<76. || mos1_fake>106.;
      cut2[2] = mos2_fake<60.;
      cut2[3] = ht<100;
      cut2[4] = delta_phi(lep[0].phi(),lep[1].phi())<2.8;
      cut2[5] = true;
      cut2[6] = mtmin>10;
      cut2[7] = lep[0].et()>10;
      cut2[8] = lep[1].et()>5;
      cut2[9] = et3>4;
    }

    bool pass_all = true;

    for ( int i=0; i<10; i++ ){

      pass_all = pass_all && cut2[i];

      bool pass = true;
      for ( int j=0;j<10;j++ ) {
	if(j!=i) {
	  pass = pass && cut2[j];
	}
      }
      if ( pass ) {
	if (i==0)
	  h.fopt_n1_met->Fill(missing_et.et());
	else if (i==1) 
	  h.fopt_n1_mos1->Fill(mos1_fake);
	else if (i==2)
	  h.fopt_n1_mos2->Fill(mos2_fake);
	else if (i==3)
	  h.fopt_n1_ht->Fill(ht);
	else if (i==4)
	  h.fopt_dphi->Fill(delta_phi(lep[0].phi(),lep[1].phi()));
	else if (i==6)
	  h.fopt_n1_mt->Fill(mtmin);
	else if (i==7)
	  h.fopt_n1_et[0]->Fill(lep[0].et());
	else if (i==8)
	  h.fopt_n1_et[1]->Fill(lep[1].et());
	else if (i==9)
	  h.fopt_n1_et[2]->Fill(et3);
      }
    

    }
    pass_all = pass_all && mos1_fake>15 && mos2_fake>15;

    return pass_all;

}

bool TSeel::pass_OptimizedCuts()
{
//  met, m_os_1, m_os_2, ht, dphi_lep12
//  et1, et2,  et3

    int njet20 = 0;
    for ( int i=0; i<njet; i++ ) {
      if ( jt[i].et()>20. ) njet20++;
    }
    bool cut2[10];

    if ( _OptScheme==1 ) { // sourabh scheme Original
        cut2[0] = opt_met>15.;
        cut2[1] = opt_m_os_1<76. || opt_m_os_1>106.;
        cut2[2] = opt_m_os_2<60.;
        cut2[3] = opt_ht<80.;
        cut2[4] = opt_dphi_lep12<2.9;
        cut2[5] = true;  // opt_trk_met>0.;
        cut2[6] = opt_mt_min>10.;
	cut2[7] = opt_ele1_et>10;
	cut2[8] = opt_ele2_et>5;
	cut2[9] = opt_trk_pt>4;
    }
    else if ( _OptScheme==3 ) { //sourabh scheme HiPt analysis
        cut2[0] = opt_met>15.;
	cut2[1] = opt_m_os_1<76. || opt_m_os_1>106.;
        cut2[2] = true;//opt_m_os_2<60.;
        cut2[3] = njet20<2;// opt_ht<100.; 
        cut2[4] = opt_dphi_lep12<2.792527;//160.000011 degrees = 2.792527 rad
        cut2[5] = true;  // opt_trk_met>0.;
        cut2[6] = true;//opt_mt_min>10.;
	cut2[7] = opt_ele1_et>20 && lep[0].pt()>10;
	cut2[8] = opt_ele2_et>8 && lep[1].pt()>8;
	cut2[9] = opt_trk_pt>5 && lep[2].pt()>5;
    }
    else if ( _OptScheme==5 ) { //sourabh scheme 5
        cut2[0] = opt_met>20.;
	cut2[1] = opt_m_os_1<76. || opt_m_os_1>106.;
        cut2[2] = opt_m_os_2<60.;
        cut2[3] = opt_ht<100.; 
        cut2[4] = opt_dphi_lep12<2.8;
        cut2[5] = true;  // opt_trk_met>0.;
        cut2[6] = true;//opt_mt_min>10.;
	cut2[7] = opt_ele1_et>10;
	cut2[8] = opt_ele2_et>8;
	cut2[9] = opt_trk_pt>4;
    }
    else if ( _OptScheme==7 ) { //sourabh scheme 7
        cut2[0] = opt_met>20.;
	cut2[1] = opt_m_os_1<76. || opt_m_os_1>106.;
        cut2[2] = opt_m_os_2<60.;
        cut2[3] = opt_ht<100.; 
        cut2[4] = opt_dphi_lep12<2.8;
        cut2[5] = true;  // opt_trk_met>0.;
        cut2[6] = true;//opt_mt_min>10.;
	cut2[7] = opt_ele1_et>18;
	cut2[8] = opt_ele2_et>5;
	cut2[9] = opt_trk_pt>4;
    }
    else if ( _OptScheme==9 ) { //sourabh scheme 9
        cut2[0] = opt_met>20.;
	cut2[1] = opt_m_os_1<76. || opt_m_os_1>106.;
        cut2[2] = opt_m_os_2<60.;
        cut2[3] = opt_ht<100.; 
        cut2[4] = opt_dphi_lep12<2.8;
        cut2[5] = true;  // opt_trk_met>0.;
        cut2[6] = opt_mt_min>10.;
	cut2[7] = opt_ele1_et>10;
	cut2[8] = opt_ele2_et>5;
	cut2[9] = opt_trk_pt>4;
    }
	
    bool pass_all = true;

    for ( int i=0; i<10; i++ ) {
        bool pass = true;
        for ( int j=0; j<10; j++ ) {
            if ( j!=i ) {
                pass = pass && cut2[j];
            }
        }
        if ( pass ) {
            if ( i==0 )
                h.n1_opt_met->Fill(opt_met);
            else if ( i==1 )
                h.n1_opt_m_os_1->Fill(opt_m_os_1);
            else if ( i==2 )
                h.n1_opt_m_os_2->Fill(opt_m_os_2);
            else if ( i==3 )
                h.n1_opt_ht->Fill(opt_ht);
            else if ( i==4 )
                h.n1_opt_dphi_lep12->Fill(opt_dphi_lep12);
            else if ( i==5 )
                h.n1_opt_trk_met->Fill(opt_trk_met);
            else if ( i==6 )
                h.n1_opt_mt_min->Fill(opt_mt_min);
	    else if ( i==7 )
	        h.n1_opt_ele1_et->Fill(opt_ele1_et);
	    else if ( i==8 )
	        h.n1_opt_ele2_et->Fill(opt_ele2_et);
	    else if ( i==9 )
	        h.n1_opt_trk_pt->Fill(opt_trk_pt);
	    
        }
        pass_all = pass_all && cut2[i];
    }
    
    if ( pass_all ) {
      //if ( _OptScheme==5 ) 
      h.n1_opt_met_trkpt->Fill(lep[2].pt(),opt_met);
      //FillOptRegionHistograms();
      if(abs(lep[0].id())==11 && abs(lep[1].id())==11){
	int trk1 = fElectronBlock->Electron(lep[0].index())->TrackNumber();
	int trk2 = fElectronBlock->Electron(lep[1].index())->TrackNumber();
	int trk3 = -1;
	if(abs(lep[2].id())==13)
	  trk3 = fMuonBlock->Muon(lep[2].index())->TrackNumber();
	if(trk1==trk3)
	  nev_sametrack++;
	if(trk2==trk3)
	  nev_sametrack++;
      }
    }


    return pass_all;
}

void TSeel::FillOptRegionHistograms()
{
  
  if ( nleptons>=3 ) {
    //cout<<"Nleptons = "<<nleptons<<endl;
//  Trilepton kinematics
    for ( int i=0; i<3; i++ ) {
      h.opt_trilep_pt[i]->Fill(lep[i].et());
      h.opt_trilep_eta[i]->Fill(lep[i].eta());
      h.opt_trilep_phi[i]->Fill(lep[i].phi());
    }
//  st
//        float st = ht;
    float st = 0.;
    for ( int i=0; i<3; i++ ) st += lep[i].et();
    float ht2 = 0.; if ( njet>0 ) ht2 = ht - jt[0].et();
    float st2 = 0.; if ( njet>0 ) st2 = st - jt[0].et();
    if ( getbits(_EventType,10,1) ) { // correct charge combination
            h.opt_met->Fill(missing_et.et());
            h.opt_st->Fill(st);
            h.opt_st2->Fill(st2);
            h.opt_st2_over_st->Fill((st>0. ? st2/st : 0.) );
            h.opt_ht->Fill(ht);
            h.opt_ht2->Fill(ht2);
            h.opt_ht2_over_ht->Fill((ht>0. ? ht2/ht : 0.) );
            h.opt_trk_met->Fill(lep[2].pt()*missing_et.et());
            h.opt_st_met->Fill(st+missing_et.et());
            h.opt_njet15->Fill(njet15);
    }

//  invariant masses
    vector<TLorentzVector> v_lep;
    for ( int i=0; i<3; i++ ) {
      TLorentzVector v(lep[i].px(),lep[i].py(),lep[i].pz(),lep[i].en());
      v_lep.push_back(v);
    }
    float m_leps = inv_mass(v_lep);
    vector<TLorentzVector> v_leps_jets = v_lep;
    for ( int i=0; i<njet; i++ ) {
      TLorentzVector v(jt[i].ex(),jt[i].ey(),jt[i].ez(),jt[i].en());
      v_leps_jets.push_back(v);
    }
    float m_leps_jets = inv_mass(v_leps_jets);
    vector<TLorentzVector> v_all = v_leps_jets;
    TLorentzVector v(missing_et.ex(),missing_et.ey(),0.,missing_et.et());
    v_all.push_back(v);
    float m_all = inv_mass(v_all);
    float mt_all = transverse_mass(v_all);

    vector<TLorentzVector> v_os_1,v_os_2,v_ss,v_ss_1,v_ss_2,v_ss_3;
    if ( lep[0].id()*lep[1].id()<0 ) {
            v_os_1.push_back(v_lep[0]);
            v_os_1.push_back(v_lep[1]);
            if ( lep[0].id()*lep[2].id()<0 ) {
                v_os_2.push_back(v_lep[0]);
                v_os_2.push_back(v_lep[2]);
            }
            else if ( lep[1].id()*lep[2].id()<0 ) {
	        v_os_2.push_back(v_lep[1]);
                v_os_2.push_back(v_lep[2]);
            }
            if ( lep[0].id()*lep[2].id()>0 ) {
                v_ss.push_back(v_lep[0]);
                v_ss.push_back(v_lep[2]);
            }
            else if ( lep[1].id()*lep[2].id()>0 ) {
                v_ss.push_back(v_lep[1]);
                v_ss.push_back(v_lep[2]);
            }
    }
    else if ( lep[0].id()*lep[2].id()<0 ) {
            v_os_1.push_back(v_lep[0]);
            v_os_1.push_back(v_lep[2]);
            v_os_2.push_back(v_lep[1]);
            v_os_2.push_back(v_lep[2]);
            v_ss.push_back(v_lep[0]);
            v_ss.push_back(v_lep[1]);
    }
    else { // all three leptons same sign
            v_ss_1.push_back(v_lep[0]); v_ss_1.push_back(v_lep[1]);
            v_ss_2.push_back(v_lep[0]); v_ss_2.push_back(v_lep[2]);
            v_ss_3.push_back(v_lep[1]); v_ss_3.push_back(v_lep[2]);
    }
//        m_os_1 = inv_mass(v_os_1);
//        m_os_2 = inv_mass(v_os_2);
//        if ( m_os_1 < m_os_2 ) {
//            float temp = m_os_1;
//            m_os_1 = m_os_2;
//            m_os_2 = temp;
//        }
    vector<float> m_ss;
    m_ss.push_back(inv_mass(v_ss_1));
    m_ss.push_back(inv_mass(v_ss_2));
    m_ss.push_back(inv_mass(v_ss_3));
    sort(m_ss.begin(), m_ss.end());
    if ( getbits(_EventType,10,1) ) { // correct charge combination
            h.opt_m_leps->Fill(m_leps);
            h.opt_m_leps_jets->Fill(m_leps_jets);
            h.opt_m_all->Fill(m_all);
            h.opt_m_os_1->Fill(m_os_1);
            h.opt_m_os_2->Fill(m_os_2);
            h.opt_m_ss->Fill(inv_mass(v_ss));
            h.opt_mt_all->Fill(mt_all);
            h.opt_met_m_os_1->Fill(m_os_1,missing_et.et());
    }
    else {
            h.opt_m_ss_1->Fill(m_ss[2]);
            h.opt_m_ss_2->Fill(m_ss[1]);
            h.opt_m_ss_3->Fill(m_ss[0]);
    }

//  Transverse masses
    vector<float> mt_lep;
    mt_lep.push_back(transverse_mass_px_py(lep[0].px(),lep[0].py(),missing_et.ex(),missing_et.ey()));
    mt_lep.push_back(transverse_mass_px_py(lep[1].px(),lep[1].py(),missing_et.ex(),missing_et.ey()));
    mt_lep.push_back(transverse_mass_px_py(lep[2].px(),lep[2].py(),missing_et.ex(),missing_et.ey()));
    if ( getbits(_EventType,10,1) ) { // correct charge combination
            h.opt_mt_lep1_met->Fill(mt_lep[0]);
            h.opt_mt_lep2_met->Fill(mt_lep[1]);
            h.opt_mt_lep3_met->Fill(mt_lep[2]);
            sort(mt_lep.begin(),mt_lep.end());
            h.opt_mt_lep_met_max->Fill(mt_lep[2]);
            h.opt_mt_lep_met_min->Fill(mt_lep[0]);
    }

//  Sphericity, Alpanarity, and Thrust
    float sphe=0.; float apla=0.;
    if ( getbits(_EventType,10,1) ) { // correct charge combination
            float thrust = Thrust(v_all); h.opt_thrust->Fill(thrust);
            if ( Sphericity_Aplanarity(v_all, sphe, apla) ) {
                h.opt_sphericity->Fill(sphe);
                h.opt_aplanarity->Fill(apla);
            }
            else {
                h.opt_sphericity->Fill(-999.);
                h.opt_aplanarity->Fill(-999.);
            }
    }
    vector<TLorentzVector> vtrans_all = v_all;
    for ( vector<TLorentzVector>::iterator i=vtrans_all.begin(); i!=vtrans_all.end(); i++ ) i->SetPz(0.);
    float sphe_2d=0.; float apla_2d=0.;
    if ( getbits(_EventType,10,1) ) { // correct charge combination
      float thrust_2d = Thrust(vtrans_all); h.opt_thrust_2d->Fill(thrust_2d);
      if ( Sphericity_Aplanarity(vtrans_all, sphe_2d, apla_2d) ) {
	h.opt_sphericity_2d->Fill(sphe_2d);
	h.opt_aplanarity_2d->Fill(apla_2d);
      }
      else {
	h.opt_sphericity->Fill(-999.);
	h.opt_aplanarity->Fill(-999.);
      }
    }
    
  } // nleptons>=3


}



void TSeel::BookHistograms()
{
  DeleteHistograms();



  //Examine
    h.jet_scaleEta = new TH2F("jet_scaleEta","jet_scaleEta",50,-2.5,2.5,20,1.0,2.0);
    AddHistogram(h.jet_scaleEta);
  h.beam_xy = new TH2F("Beam_XY","Beam_XY",100,-5.,5.,100,-5.,5.);AddHistogram(h.beam_xy); 

  h.EleEt_Conv = new TH1F("ele_et_conv","Conversion Et",200,0.,50.); AddHistogram(h.EleEt_Conv);
  h.EleEt_susy_ele = new TH1F("ele_et_susy_ele","Et (SUSY ele before conv)",200,0.,50.); AddHistogram(h.EleEt_susy_ele);
  h.conv_R = new TH1F("conv_R","conv_R",50,0.,50.);AddHistogram(h.conv_R);
  h.conv_Sxy_dCotTheta = new TH2F("conv_Sxy_dCotTheta","conv_Sxy_dCotTheta",200,-0.2,0.2,100,-0.1,0.1);AddHistogram(h.conv_Sxy_dCotTheta);
  h.conv_R_susy_ele = new TH1F("conv_R_susy_ele","conv_R_susy_ele",50,0.,50.);AddHistogram(h.conv_R_susy_ele);
  h.conv_Sxy_dCotTheta_susy_ele = new TH2F("conv_Sxy_dCotTheta_susy_ele","conv_Sxy_dCotTheta_susy_ele",200,-0.2,0.2,100,-0.1,0.1);AddHistogram(h.conv_Sxy_dCotTheta_susy_ele);
  h.conv_trk2_pt = new TH1F("conv_trk2_pt","conv_trk2_pt",100,0.,50.);AddHistogram(h.conv_trk2_pt);
  h.conv_trk2_pt_susy_ele = new TH1F("conv_trk2_pt_susy_ele","conv_trk2_pt_susy_ele",100,0.,50.);AddHistogram(h.conv_trk2_pt_susy_ele);
  h.conv_Lxy = new TH1F("conv_Lxy","conv_Lxy",60,-10.,50.);AddHistogram(h.conv_Lxy);


    h.hepg_ele_matching_dr = new TH1F("hepg_ele_matching_dr","hepg_ele_matching_dr",100,0.,0.2); AddHistogram(h.hepg_ele_matching_dr);
    h.hepg_mu_matching_dr = new TH1F("hepg_mu_matching_dr","hepg_mu_matching_dr",100,0.,0.2); AddHistogram(h.hepg_mu_matching_dr);
    h.hepg_ele_matching_dp = new TH1F("hepg_ele_matching_dp","hepg_ele_matching_dp",200,-0.2,0.2); AddHistogram(h.hepg_ele_matching_dp);
    h.hepg_mu_matching_dp = new TH1F("hepg_mu_matching_dp","hepg_mu_matching_dp",200,-0.2,0.2); AddHistogram(h.hepg_mu_matching_dp);


    h.dphi_j1met_nocorr = new TH1F("dphi_j1met_nocorr","dphi_j1met_nocorr",64,0.,3.2);AddHistogram(h.dphi_j1met_nocorr);
    h.raw_phi_jet1= new TH1F("raw_phi_jet1","raw_phi_jet1",140,0,7);AddHistogram(h.raw_phi_jet1);
    h.met_noCorr = new TH1F("met_noCorr","met_noCorr",200,0.,200.);AddHistogram(h.met_noCorr);
    h.met_phi_noCorr = new TH1F("met_phi_noCorr","met_phi_noCorr",128,0.,6.4);AddHistogram(h.met_phi_noCorr);
    h.met_JetCorr = new TH1F("met_JetCorr","met_JetCorr",200,0.,200.);AddHistogram(h.met_JetCorr);
    h.met_phi_JetCorr = new TH1F("met_phi_JetCorr","met_phi_JetCorr",128,0.,6.4);AddHistogram(h.met_phi_JetCorr);
    h.met_MipCorr = new TH1F("met_MipCorr","met_MipCorr",200,0.,200.);AddHistogram(h.met_MipCorr);
    h.met_phi_MipCorr = new TH1F("met_phi_MipCorr","met_phi_MipCorr",128,0.,6.4);AddHistogram(h.met_phi_MipCorr);


    h.raw_dr_lep12 = new TH1F("raw_dr_lep12","raw_dr_lep12",100,0.,5.);AddHistogram(h.raw_dr_lep12);
    h.raw_dr_lep13 = new TH1F("raw_dr_lep13","raw_dr_lep13",100,0.,5.);AddHistogram(h.raw_dr_lep13);
    h.raw_dr_lep23 = new TH1F("raw_dr_lep23","raw_dr_lep23",100,0.,5.);AddHistogram(h.raw_dr_lep23);

    h.jet_width[0] = new TH2F("jet_width_all","M/Et vs Et ,all jets",200,0,2,200,0,100);AddHistogram(h.jet_width[0]);
    h.jet_width[1] = new TH2F("jet_width_1","M/Et vs Et ,lead jet",200,0,2,200,0,100);AddHistogram(h.jet_width[1]);

//  Central Electron
    h.ce_id_raw[0] = new TH1F("raw_ceFid","raw_ceFid",10,0.,10.);
    h.ce_id_raw[1] = new TH1F("raw_ceHadEm","raw_ceHadEm",50,0.,0.5);
    h.ce_id_raw[2] = new TH1F("raw_ceIsoFr","raw_ceIsoFr",50,0.,0.5);
    h.ce_id_raw[3] = new TH1F("raw_cePt","raw_cePt",200,0.,200.);
    h.ce_id_raw[4] = new TH1F("raw_ceEoP","raw_ceEoP",100,0.,5.);
    h.ce_id_raw[5] = new TH1F("raw_ceNAxSeg","raw_ceNAxSeg",5,0.,5.);
    h.ce_id_raw[6] = new TH1F("raw_ceNStSeg","raw_ceNStSeg",5,0.,5.);
    h.ce_id_raw[7] = new TH1F("raw_ceZ0","raw_ceZ0",100,-100.,100.);
    h.ce_id_raw[8] = new TH1F("raw_ceLshr","raw_ceLshr",100,-0.5,0.5);
    h.ce_id_raw[9] = new TH1F("raw_ceChi2Strip","raw_ceChi2Strip",100,0.,50.);
    h.ce_id_raw[10] = new TH1F("raw_ceDx","raw_ceDx",100,-5.,5.);
    h.ce_id_raw[11] = new TH1F("raw_ceQDx","raw_cdQDx",100,-5.,5.);
    h.ce_id_raw[12] = new TH1F("raw_ceDz","raw_ceDz",100,-10.,10.);
    h.ce_et_raw = new TH1F("raw_ce_et","raw_ce_et",200,0.,200.);
    h.ce_eta_raw = new TH1F("raw_ce_eta","raw_ce_eta",200,-4.,4.);
    h.ce_phi_raw = new TH1F("raw_ce_phi","raw_ce_phi",160,0.,6.4);
    for ( int i=0; i<13; i++ ) AddHistogram(h.ce_id_raw[i]);
    AddHistogram(h.ce_et_raw); AddHistogram(h.ce_eta_raw); AddHistogram(h.ce_phi_raw);
    h.ce_id[0] = new TH1F("ceFid","ceFid (N-1)",10,0.,10.);
    h.ce_id[1] = new TH1F("ceHadEm","ceHadEm (N-1)",50,0.,0.5);
    h.ce_id[2] = new TH1F("ceIsoFr","ceIsoFr (N-1)",50,0.,0.5);
    h.ce_id[3] = new TH1F("cePt","cePt (N-1)",200,0.,200.);
    h.ce_id[4] = new TH1F("ceEoP","ceEoP (N-1)",100,0.,5.);
    h.ce_id[5] = new TH1F("ceNAxSeg","ceNAxSeg (N-1)",5,0.,5.);
    h.ce_id[6] = new TH1F("ceNStSeg","ceNStSeg (N-1)",5,0.,5.);
    h.ce_id[7] = new TH1F("ceZ0","ceZ0 (N-1)",100,-100.,100.);
    h.ce_id[8] = new TH1F("ceLshr","ceLshr (N-1)",100,-0.5,0.5);
    h.ce_id[9] = new TH1F("ceChi2Strip","ceChi2Strip (N-1)",100,0.,50.);
    h.ce_id[10] = new TH1F("ceDx","ceDx (N-1)",100,-5.,5.);
    h.ce_id[11] = new TH1F("ceQDx","cdQDx (N-1)",100,-5.,5.);
    h.ce_id[12] = new TH1F("ceDz","ceDz (N-1)",100,-10.,10.);
    h.ce_et = new TH1F("ce_et","ce_et (pass all)",200,0.,200.);
    h.ce_eta = new TH1F("ce_eta","ce_eta (pass all)",200,-4.,4.);
    h.ce_phi = new TH1F("ce_phi","ce_phi (pass all)",160,0.,6.4);
    for ( int i=0; i<13; i++ ) AddHistogram(h.ce_id[i]);
    AddHistogram(h.ce_et); AddHistogram(h.ce_eta); AddHistogram(h.ce_phi);
//    Central muon
// Raw distributions
//_____________________
//  CMU
    h.cmu_id_raw[0] = new TH1F("raw_cmu_Em_en","raw_cmu_Em_en",200,0.,20.);
    h.cmu_id_raw[1] = new TH1F("raw_cmu_Had_en","raw_cmu_Had_en",200,0.,20.);
    h.cmu_id_raw[2] = new TH1F("raw_cmu_Dx","raw_cmu_Dx",500,-50.,50.);
    h.cmu_id_raw[3] = new TH1F("raw_cmu_fiso","raw_cmu_fiso",100,0.,10.);
    h.cmu_id_raw[4] = new TH1F("raw_cmu_d0","raw_cmu_d0",100,-0.5,0.5);
    h.cmu_id_raw[5] = new TH1F("raw_cmu_nAxSegs","raw_cmu_nAxSegs",5,0.,5.);
    h.cmu_id_raw[6] = new TH1F("raw_cmu_nStSegs","raw_cmu_nStSegs",5,0.,5.);
    h.cmu_id_raw[7] = new TH1F("raw_cmu_trk_chi2ndof","raw_cmu_trk_chi2ndof",100,0.,5.);
    h.cmu_id_raw[8] = new TH1F("raw_cmu_chi2match","raw_cmu_chi2match",250,0.,50.);
    h.cmu_pt_raw = new TH1F("raw_cmu_pt","raw_cmu_pt",200,0.,200.);
    h.cmu_eta_raw = new TH1F("raw_cmu_eta","raw_cmu_eta",100,-2.,2.);
    h.cmu_phi_raw = new TH1F("raw_cmu_phi","raw_cmu_phi",128,0.,6.4);
    for ( int i=0; i<9; i++ ) AddHistogram(h.cmu_id_raw[i]);
    AddHistogram(h.cmu_pt_raw); AddHistogram(h.cmu_eta_raw); AddHistogram(h.cmu_phi_raw);

//  CMUP
    h.cmup_id_raw[0] = new TH1F("raw_cmup_Em_en","raw_cmup_Em_en",200,0.,20.);
    h.cmup_id_raw[1] = new TH1F("raw_cmup_Had_en","raw_cmup_Had_en",200,0.,20.);
    h.cmup_id_raw[2] = new TH1F("raw_cmup_Dx","raw_cmup_Dx",500,-50.,50.);
    h.cmup_id_raw[3] = new TH1F("raw_cmup_fiso","raw_cmup_fiso",100,0.,10.);
    h.cmup_id_raw[4] = new TH1F("raw_cmup_d0","raw_cmup_d0",100,-0.5,0.5);
    h.cmup_id_raw[5] = new TH1F("raw_cmup_nAxSegs","raw_cmup_nAxSegs",5,0.,5.);
    h.cmup_id_raw[6] = new TH1F("raw_cmup_nStSegs","raw_cmup_nStSegs",5,0.,5.);
    h.cmup_id_raw[7] = new TH1F("raw_cmup_trk_chi2ndof","raw_cmup_trk_chi2ndof",100,0.,5.);
    h.cmup_id_raw[8] = new TH1F("raw_cmup_chi2match","raw_cmup_chi2match",250,0.,50.);
    h.cmup_pt_raw = new TH1F("raw_cmup_pt","raw_cmup_pt",200,0.,200.);
    h.cmup_eta_raw = new TH1F("raw_cmup_eta","raw_cmup_eta",100,-2.,2.);
    h.cmup_phi_raw = new TH1F("raw_cmup_phi","raw_cmup_phi",128,0.,6.4);
    for ( int i=0; i<9; i++ ) AddHistogram(h.cmup_id_raw[i]);
    AddHistogram(h.cmup_pt_raw); AddHistogram(h.cmup_eta_raw); AddHistogram(h.cmup_phi_raw);

//  CMP
    h.cmp_id_raw[0] = new TH1F("raw_cmp_Em_en","raw_cmp_Em_en",200,0.,20.);
    h.cmp_id_raw[1] = new TH1F("raw_cmp_Had_en","raw_cmp_Had_en",200,0.,20.);
    h.cmp_id_raw[2] = new TH1F("raw_cmp_Dx","raw_cmp_Dx",500,-50.,50.);
    h.cmp_id_raw[3] = new TH1F("raw_cmp_fiso","raw_cmp_fiso",100,0.,10.);
    h.cmp_id_raw[4] = new TH1F("raw_cmp_d0","raw_cmp_d0",100,-0.5,0.5);
    h.cmp_id_raw[5] = new TH1F("raw_cmp_nAxSegs","raw_cmp_nAxSegs",5,0.,5.);
    h.cmp_id_raw[6] = new TH1F("raw_cmp_nStSegs","raw_cmp_nStSegs",5,0.,5.);
    h.cmp_id_raw[7] = new TH1F("raw_cmp_trk_chi2ndof","raw_cmp_trk_chi2ndof",100,0.,5.);
    h.cmp_id_raw[8] = new TH1F("raw_cmp_chi2match","raw_cmp_chi2match",250,0.,50.);
    h.cmp_pt_raw = new TH1F("raw_cmp_pt","raw_cmp_pt",200,0.,200.);
    h.cmp_eta_raw = new TH1F("raw_cmp_eta","raw_cmp_eta",100,-2.,2.);
    h.cmp_phi_raw = new TH1F("raw_cmp_phi","raw_cmp_phi",128,0.,6.4);
    for ( int i=0; i<9; i++ ) AddHistogram(h.cmp_id_raw[i]);
    AddHistogram(h.cmp_pt_raw); AddHistogram(h.cmp_eta_raw); AddHistogram(h.cmp_phi_raw);

//  CMX
    h.cmx_id_raw[0] = new TH1F("raw_cmx_Em_en","raw_cmx_Em_en",200,0.,20.);
    h.cmx_id_raw[1] = new TH1F("raw_cmx_Had_en","raw_cmx_Had_en",200,0.,20.);
    h.cmx_id_raw[2] = new TH1F("raw_cmx_Dx","raw_cmx_Dx",500,-50.,50.);
    h.cmx_id_raw[3] = new TH1F("raw_cmx_fiso","raw_cmx_fiso",100,0.,10.);
    h.cmx_id_raw[4] = new TH1F("raw_cmx_d0","raw_cmx_d0",100,-0.5,0.5);
    h.cmx_id_raw[5] = new TH1F("raw_cmx_nAxSegs","raw_cmx_nAxSegs",5,0.,5.);
    h.cmx_id_raw[6] = new TH1F("raw_cmx_nStSegs","raw_cmx_nStSegs",5,0.,5.);
    h.cmx_id_raw[7] = new TH1F("raw_cmx_trk_chi2ndof","raw_cmx_trk_chi2ndof",100,0.,5.);
    h.cmx_id_raw[8] = new TH1F("raw_cmx_chi2match","raw_cmx_chi2match",250,0.,50.);
    h.cmx_pt_raw = new TH1F("raw_cmx_pt","raw_cmx_pt",200,0.,200.);
    h.cmx_eta_raw = new TH1F("raw_cmx_eta","raw_cmx_eta",100,-2.,2.);
    h.cmx_phi_raw = new TH1F("raw_cmx_phi","raw_cmx_phi",128,0.,6.4);
    for ( int i=0; i<9; i++ ) AddHistogram(h.cmx_id_raw[i]);
    AddHistogram(h.cmx_pt_raw); AddHistogram(h.cmx_eta_raw); AddHistogram(h.cmx_phi_raw);
//  CMU
    h.cmu_id[0] = new TH1F("cmu_Em_en","cmu_Em_en (N-1)",200,0.,20.);
    h.cmu_id[1] = new TH1F("cmu_Had_en","cmu_Had_en (N-1)",200,0.,20.);
    h.cmu_id[2] = new TH1F("cmu_Dx","cmu_Dx (N-1)",500,-50.,50.);
    h.cmu_id[3] = new TH1F("cmu_fiso","cmu_fiso (N-1)",100,0.,10.);
    h.cmu_id[4] = new TH1F("cmu_d0","cmu_d0 (N-1)",500,-0.5,0.5);
    h.cmu_id[5] = new TH1F("cmu_nAxSegs","cmu_nAxSegs (N-1)",5,0.,5.);
    h.cmu_id[6] = new TH1F("cmu_nStSegs","cmu_nStSegs (N-1)",5,0.,5.);
    h.cmu_id[7] = new TH1F("cmu_trk_chi2ndof","cmu_trk_chi2ndof (N-1)",100,0.,5.);
    h.cmu_id[8] = new TH1F("cmu_chi2match","cmu_chi2match (N-1)",250,0.,50.);
    h.cmu_pt = new TH1F("cmu_pt","cmu_pt (pass all)",200,0.,200.);
    h.cmu_eta = new TH1F("cmu_eta","cmu_eta (pass all)",100,-2.,2.);
    h.cmu_phi = new TH1F("cmu_phi","cmu_phi (pass all)",128,0.,6.4);
    for ( int i=0; i<9; i++ ) AddHistogram(h.cmu_id[i]);
    AddHistogram(h.cmu_pt); AddHistogram(h.cmu_eta); AddHistogram(h.cmu_phi);

//  CMUP
    h.cmup_id[0] = new TH1F("cmup_Em_en","cmup_Em_en (N-1)",200,0.,20.);
    h.cmup_id[1] = new TH1F("cmup_Had_en","cmup_Had_en (N-1)",200,0.,20.);
    h.cmup_id[2] = new TH1F("cmup_Dx","cmup_Dx (N-1)",500,-50.,50.);
    h.cmup_id[3] = new TH1F("cmup_fiso","cmup_fiso (N-1)",100,0.,10.);
    h.cmup_id[4] = new TH1F("cmup_d0","cmup_d0 (N-1)",500,-0.5,0.5);
    h.cmup_id[5] = new TH1F("cmup_nAxSegs","cmup_nAxSegs (N-1)",5,0.,5.);
    h.cmup_id[6] = new TH1F("cmup_nStSegs","cmup_nStSegs (N-1)",5,0.,5.);
    h.cmup_id[7] = new TH1F("cmup_trk_chi2ndof","cmup_trk_chi2ndof (N-1)",100,0.,5.);
    h.cmup_id[8] = new TH1F("cmup_chi2match","cmup_chi2match (N-1)",250,0.,50.);
    h.cmup_pt = new TH1F("cmup_pt","cmup_pt (pass all)",200,0.,200.);
    h.cmup_eta = new TH1F("cmup_eta","cmup_eta (pass all)",100,-2.,2.);
    h.cmup_phi = new TH1F("cmup_phi","cmup_phi (pass all)",128,0.,6.4);
    for ( int i=0; i<9; i++ ) AddHistogram(h.cmup_id[i]);
    AddHistogram(h.cmup_pt); AddHistogram(h.cmup_eta); AddHistogram(h.cmup_phi);

//  CMP
    h.cmp_id[0] = new TH1F("cmp_Em_en","cmp_Em_en (N-1)",200,0.,20.);
    h.cmp_id[1] = new TH1F("cmp_Had_en","cmp_Had_en (N-1)",200,0.,20.);
    h.cmp_id[2] = new TH1F("cmp_Dx","cmp_Dx (N-1)",500,-50.,50.);
    h.cmp_id[3] = new TH1F("cmp_fiso","cmp_fiso (N-1)",100,0.,10.);
    h.cmp_id[4] = new TH1F("cmp_d0","cmp_d0 (N-1)",500,-0.5,0.5);
    h.cmp_id[5] = new TH1F("cmp_nAxSegs","cmp_nAxSegs (N-1)",5,0.,5.);
    h.cmp_id[6] = new TH1F("cmp_nStSegs","cmp_nStSegs (N-1)",5,0.,5.);
    h.cmp_id[7] = new TH1F("cmp_trk_chi2ndof","cmp_trk_chi2ndof (N-1)",100,0.,5.);
    h.cmp_id[8] = new TH1F("cmp_chi2match","cmp_chi2match (N-1)",250,0.,50.);
    h.cmp_pt = new TH1F("cmp_pt","cmp_pt (pass all)",200,0.,200.);
    h.cmp_eta = new TH1F("cmp_eta","cmp_eta (pass all)",100,-2.,2.);
    h.cmp_phi = new TH1F("cmp_phi","cmp_phi (pass all)",128,0.,6.4);
    for ( int i=0; i<9; i++ ) AddHistogram(h.cmp_id[i]);
    AddHistogram(h.cmp_pt); AddHistogram(h.cmp_eta); AddHistogram(h.cmp_phi);

//  CMX
    h.cmx_id[0] = new TH1F("cmx_Em_en","cmx_Em_en (N-1)",200,0.,20.);
    h.cmx_id[1] = new TH1F("cmx_Had_en","cmx_Had_en (N-1)",200,0.,20.);
    h.cmx_id[2] = new TH1F("cmx_Dx","cmx_Dx (N-1)",500,-50.,50.);
    h.cmx_id[3] = new TH1F("cmx_fiso","cmx_fiso (N-1)",100,0.,10.);
    h.cmx_id[4] = new TH1F("cmx_d0","cmx_d0 (N-1)",500,-0.5,0.5);
    h.cmx_id[5] = new TH1F("cmx_nAxSegs","cmx_nAxSegs (N-1)",5,0.,5.);
    h.cmx_id[6] = new TH1F("cmx_nStSegs","cmx_nStSegs (N-1)",5,0.,5.);
    h.cmx_id[7] = new TH1F("cmx_trk_chi2ndof","cmx_trk_chi2ndof (N-1)",100,0.,5.);
    h.cmx_id[8] = new TH1F("cmx_chi2match","cmx_chi2match (N-1)",250,0.,50.);
    h.cmx_pt = new TH1F("cmx_pt","cmx_pt (pass all)",200,0.,200.);
    h.cmx_eta = new TH1F("cmx_eta","cmx_eta (pass all)",100,-2.,2.);
    h.cmx_phi = new TH1F("cmx_phi","cmx_phi (pass all)",128,0.,6.4);
    for ( int i=0; i<9; i++ ) AddHistogram(h.cmx_id[i]);
    AddHistogram(h.cmx_pt); AddHistogram(h.cmx_eta); AddHistogram(h.cmx_phi);

//=================================
// Full Selection

    h.dilep_mass_ee_os = new TH1F("dilep_mass_ee_os","dilep_mass_ee_os",200,0.,200.);AddHistogram(h.dilep_mass_ee_os);
    h.dilep_mass_ee_ss = new TH1F("dilep_mass_ee_ss","dilep_mass_ee_ss",200,0.,200.);AddHistogram(h.dilep_mass_ee_ss);
    h.dilep_mass_mm_os = new TH1F("dilep_mass_mm_os","dilep_mass_mm_os",200,0.,200.);AddHistogram(h.dilep_mass_mm_os);
    h.dilep_mass_mm_ss = new TH1F("dilep_mass_mm_ss","dilep_mass_mm_ss",200,0.,200.);AddHistogram(h.dilep_mass_mm_ss);
    h.dilep_mass_em_os = new TH1F("dilep_mass_em_os","dilep_mass_em_os",200,0.,200.);AddHistogram(h.dilep_mass_em_os);
    h.dilep_mass_em_ss = new TH1F("dilep_mass_em_ss","dilep_mass_em_ss",200,0.,200.);AddHistogram(h.dilep_mass_em_ss);
    

    h.trilep_trig_weight = new TH1F("trilep_trig_weight","trilep_trig_weight",100,0.,1.);AddHistogram(h.trilep_trig_weight);
    h.trilep_lepID_weight = new TH1F("trilep_lepID_weight","trilep_lepID_weight",100,0.,1.);AddHistogram(h.trilep_lepID_weight);
    h.trilep_event_weight = new TH1F("trilep_event_weight","trilep_event_weight",100,0.,1.);AddHistogram(h.trilep_event_weight);
    h.trilep_trig_weight_err = new TH1F("trilep_trig_weight_err","trilep_trig_weight_err",100,0.,1.);AddHistogram(h.trilep_trig_weight_err);
    h.trilep_lepID_weight_err = new TH1F("trilep_lepID_weight_err","trilep_lepID_weight_err",100,0.,1.);AddHistogram(h.trilep_lepID_weight_err);
    h.trilep_event_weight_err = new TH1F("trilep_event_weight_err","trilep_event_weight_err",100,0.,1.);AddHistogram(h.trilep_event_weight_err);

//=====================================================================
//  ControlRegion variables
//  Control Regions:
// c1 - MET<10 & Diele + track  --- This is implemented already.
// c2 - MET>10 & Diele + track & 76<Mee<106
// c3 - MET<10 & Diele + track & 2 jets
// c4 - MET>10 & Diele + track & 2 jets & 76<Mee<106
// ---------------------------------------
// c1
// ---------------------------------------
    h.c_trilep_pt[0][0] = new TH1F("c1_trilep_pt1","c1_trilep pt 1",100,0.,100.);
    h.c_trilep_pt[0][1] = new TH1F("c1_trilep_pt2","c1_trilep pt 2",100,0.,100.);
    h.c_trilep_pt[0][2] = new TH1F("c1_trilep_pt3","c1_trilep pt 3",100,0.,100.);
    h.c_trilep_eta[0][0] = new TH1F("c1_trilep_eta1","c1_trilep eta 1",120,-3.,3.);
    h.c_trilep_eta[0][1] = new TH1F("c1_trilep_eta2","c1_trilep eta 2",120,-3.,3.);
    h.c_trilep_eta[0][2] = new TH1F("c1_trilep_eta3","c1_trilep eta 3",120,-3.,3.);
    h.c_trilep_phi[0][0] = new TH1F("c1_trilep_phi1","c1_trilep phi 1",128,0.,6.4);
    h.c_trilep_phi[0][1] = new TH1F("c1_trilep_phi2","c1_trilep phi 2",128,0.,6.4);
    h.c_trilep_phi[0][2] = new TH1F("c1_trilep_phi3","c1_trilep phi 3",128,0.,6.4);

    h.c_met[0] = new TH1F("c1_met","c1_met",400,0.,400.);
    h.c_st[0] = new TH1F("c1_st","c1_st",400,0.,400.);
    h.c_st2[0] = new TH1F("c1_st2","c1_st2",400,0.,400.);
    h.c_st2_over_st[0] = new TH1F("c1_st2_over_st","c1_st2_over_st",100,0.,1.);
    h.c_ht[0] = new TH1F("c1_ht","c1_ht",200,0.,200.);
    h.c_ht2[0] = new TH1F("c1_ht2","c1_ht2",200,0.,200.);
    h.c_ht2_over_ht[0] = new TH1F("c1_ht2_over_ht","c1_ht2_over_ht",100,0.,1.);

    h.c_m_leps[0] = new TH1F("c1_m_leps","c1_m_leps",400,0.,400.);
    h.c_m_leps_jets[0] = new TH1F("c1_m_leps_jets","c1_m_leps_jets",400,0.,400.);
    h.c_m_all[0] = new TH1F("c1_m_all","c1_m_all",400,0.,400.);
    h.c_m_os_1[0] = new TH1F("c1_m_os_1","c1_m_os_1",200,0.,200.); // highest opposite sign mass
    h.c_m_os_2[0] = new TH1F("c1_m_os_2","c1_m_os_2",200,0.,200.);
    h.c_m_ss[0] = new TH1F("c1_m_ss","c1_m_ss",200,0.,200.);
    h.c_m_ss_1[0] = new TH1F("c1_m_ss_1","c1_m_ss_1",200,0.,200.); // highest same sign mass
    h.c_m_ss_2[0] = new TH1F("c1_m_ss_2","c1_m_ss_2",200,0.,200.);
    h.c_m_ss_3[0] = new TH1F("c1_m_ss_3","c1_m_ss_3",200,0.,200.);

    h.c_mt_all[0] = new TH1F("c1_mt_all","c1_mt_all",500,0.,500.);
    h.c_mt_lep1_met[0] = new TH1F("c1_mt_lep1_met","c1_mt_lep1_met",200,0.,200.);
    h.c_mt_lep2_met[0] = new TH1F("c1_mt_lep2_met","c1_mt_lep2_met",200,0.,200.);
    h.c_mt_lep3_met[0] = new TH1F("c1_mt_lep3_met","c1_mt_lep3_met",200,0.,200.);
    h.c_mt_lep_met_max[0] = new TH1F("c1_mt_lep_met_max","c1_mt_lep_met_max",200,0.,200.);
    h.c_mt_lep_met_min[0] = new TH1F("c1_mt_lep_met_min","c1_mt_lep_met_min",200,0.,200.);

    h.c_sphericity[0] = new TH1F("c1_sphericity","c1_sphericity",40,0.,2.);
    h.c_aplanarity[0] = new TH1F("c1_aplanarity","c1_aplanarity",40,0.,2.);
    h.c_thrust[0] = new TH1F("c1_thrust","c1_thrust",50,0.,1.);
    h.c_sphericity_2d[0] = new TH1F("c1_sphericity_2d","c1_sphericity_2d",40,0.,2.);
    h.c_aplanarity_2d[0] = new TH1F("c1_aplanarity_2d","c1_aplanarity_2d",40,0.,2.);
    h.c_thrust_2d[0] = new TH1F("c1_thrust_2d","c1_thrust_2d",50,0.,1.);

    h.c_dphi_lep12[0] = new TH1F("c1_dphi_lep12","c1_dphi_lep12",64,0.,3.2);
    h.c_dphi_lep13[0] = new TH1F("c1_dphi_lep13","c1_dphi_lep13",64,0.,3.2);
    h.c_dphi_lep23[0] = new TH1F("c1_dphi_lep23","c1_dphi_lep23",64,0.,3.2);
    h.c_dphi_lep1met[0] = new TH1F("c1_dphi_lep1met","c1_dphi_lep1met",64,0.,3.2);
    h.c_dphi_lep2met[0] = new TH1F("c1_dphi_lep2met","c1_dphi_lep2met",64,0.,3.2);
    h.c_dphi_lep3met[0] = new TH1F("c1_dphi_lep3met","c1_dphi_lep3met",64,0.,3.2);
    h.c_dr_lep12[0] = new TH1F("c1_dr_lep12","c1_dr_lep12",100,0.,5.);
    h.c_dr_lep13[0] = new TH1F("c1_dr_lep13","c1_dr_lep13",100,0.,5.);
    h.c_dr_lep23[0] = new TH1F("c1_dr_lep23","c1_dr_lep23",100,0.,5.);

    h.c_dphi_j1met[0] = new TH1F("c1_dphi_j1met","c1_dphi_j1met",64,0.,3.2);
    h.c_dphi_j2met[0] = new TH1F("c1_dphi_j2met","c1_dphi_j2met",64,0.,3.2);
    h.c_dr_jt_lep1[0] = new TH1F("c1_dr_jt_lep1","c1_dr_jt_lep1",100,0.,5.);
    h.c_dr_jt_lep2[0] = new TH1F("c1_dr_jt_lep2","c1_dr_jt_lep2",100,0.,5.);
    h.c_dr_jt_lep3[0] = new TH1F("c1_dr_jt_lep3","c1_dr_jt_lep3",100,0.,5.);

    h.c_met_m_os_1[0] = new TH2F("c1_met_m_os_1","c1_met_m_os_1",50,0.,200.,50,0.,200.);
    h.c_trk_met[0] = new TH1F("c1_trk_met","c1_trk_met",1000,0.,1000.);
    h.c_st_met[0] = new TH1F("c1_st_met","c1_st_met",500,0.,500.);
    h.c_njet15[0] = new TH1F("c1_njet15","c1_njet15",5,0.,5.);
// ---------------------------------------
// c2
// ---------------------------------------
    h.c_trilep_pt[1][0] = new TH1F("c2_trilep_pt1","c2_trilep pt 1",100,0.,100.);
    h.c_trilep_pt[1][1] = new TH1F("c2_trilep_pt2","c2_trilep pt 2",100,0.,100.);
    h.c_trilep_pt[1][2] = new TH1F("c2_trilep_pt3","c2_trilep pt 3",100,0.,100.);
    h.c_trilep_eta[1][0] = new TH1F("c2_trilep_eta1","c2_trilep eta 1",120,-3.,3.);
    h.c_trilep_eta[1][1] = new TH1F("c2_trilep_eta2","c2_trilep eta 2",120,-3.,3.);
    h.c_trilep_eta[1][2] = new TH1F("c2_trilep_eta3","c2_trilep eta 3",120,-3.,3.);
    h.c_trilep_phi[1][0] = new TH1F("c2_trilep_phi1","c2_trilep phi 1",128,0.,6.4);
    h.c_trilep_phi[1][1] = new TH1F("c2_trilep_phi2","c2_trilep phi 2",128,0.,6.4);
    h.c_trilep_phi[1][2] = new TH1F("c2_trilep_phi3","c2_trilep phi 3",128,0.,6.4);

    h.c_met[1] = new TH1F("c2_met","c2_met",400,0.,400.);
    h.c_st[1] = new TH1F("c2_st","c2_st",400,0.,400.);
    h.c_st2[1] = new TH1F("c2_st2","c2_st2",400,0.,400.);
    h.c_st2_over_st[1] = new TH1F("c2_st2_over_st","c2_st2_over_st",100,0.,1.);
    h.c_ht[1] = new TH1F("c2_ht","c2_ht",200,0.,200.);
    h.c_ht2[1] = new TH1F("c2_ht2","c2_ht2",200,0.,200.);
    h.c_ht2_over_ht[1] = new TH1F("c2_ht2_over_ht","c2_ht2_over_ht",100,0.,1.);

    h.c_m_leps[1] = new TH1F("c2_m_leps","c2_m_leps",400,0.,400.);
    h.c_m_leps_jets[1] = new TH1F("c2_m_leps_jets","c2_m_leps_jets",400,0.,400.);
    h.c_m_all[1] = new TH1F("c2_m_all","c2_m_all",400,0.,400.);
    h.c_m_os_1[1] = new TH1F("c2_m_os_1","c2_m_os_1",200,0.,200.); // highest opposite sign mass
    h.c_m_os_2[1] = new TH1F("c2_m_os_2","c2_m_os_2",200,0.,200.);
    h.c_m_ss[1] = new TH1F("c2_m_ss","c2_m_ss",200,0.,200.);
    h.c_m_ss_1[1] = new TH1F("c2_m_ss_1","c2_m_ss_1",200,0.,200.); // highest same sign mass
    h.c_m_ss_2[1] = new TH1F("c2_m_ss_2","c2_m_ss_2",200,0.,200.);
    h.c_m_ss_3[1] = new TH1F("c2_m_ss_3","c2_m_ss_3",200,0.,200.);

    h.c_mt_all[1] = new TH1F("c2_mt_all","c2_mt_all",500,0.,500.);
    h.c_mt_lep1_met[1] = new TH1F("c2_mt_lep1_met","c2_mt_lep1_met",200,0.,200.);
    h.c_mt_lep2_met[1] = new TH1F("c2_mt_lep2_met","c2_mt_lep2_met",200,0.,200.);
    h.c_mt_lep3_met[1] = new TH1F("c2_mt_lep3_met","c2_mt_lep3_met",200,0.,200.);
    h.c_mt_lep_met_max[1] = new TH1F("c2_mt_lep_met_max","c2_mt_lep_met_max",200,0.,200.);
    h.c_mt_lep_met_min[1] = new TH1F("c2_mt_lep_met_min","c2_mt_lep_met_min",200,0.,200.);

    h.c_sphericity[1] = new TH1F("c2_sphericity","c2_sphericity",40,0.,2.);
    h.c_aplanarity[1] = new TH1F("c2_aplanarity","c2_aplanarity",40,0.,2.);
    h.c_thrust[1] = new TH1F("c2_thrust","c2_thrust",50,0.,1.);
    h.c_sphericity_2d[1] = new TH1F("c2_sphericity_2d","c2_sphericity_2d",40,0.,2.);
    h.c_aplanarity_2d[1] = new TH1F("c2_aplanarity_2d","c2_aplanarity_2d",40,0.,2.);
    h.c_thrust_2d[1] = new TH1F("c2_thrust_2d","c2_thrust_2d",50,0.,1.);

    h.c_dphi_lep12[1] = new TH1F("c2_dphi_lep12","c2_dphi_lep12",64,0.,3.2);
    h.c_dphi_lep13[1] = new TH1F("c2_dphi_lep13","c2_dphi_lep13",64,0.,3.2);
    h.c_dphi_lep23[1] = new TH1F("c2_dphi_lep23","c2_dphi_lep23",64,0.,3.2);
    h.c_dphi_lep1met[1] = new TH1F("c2_dphi_lep1met","c2_dphi_lep1met",64,0.,3.2);
    h.c_dphi_lep2met[1] = new TH1F("c2_dphi_lep2met","c2_dphi_lep2met",64,0.,3.2);
    h.c_dphi_lep3met[1] = new TH1F("c2_dphi_lep3met","c2_dphi_lep3met",64,0.,3.2);
    h.c_dr_lep12[1] = new TH1F("c2_dr_lep12","c2_dr_lep12",100,0.,5.);
    h.c_dr_lep13[1] = new TH1F("c2_dr_lep13","c2_dr_lep13",100,0.,5.);
    h.c_dr_lep23[1] = new TH1F("c2_dr_lep23","c2_dr_lep23",100,0.,5.);

    h.c_dphi_j1met[1] = new TH1F("c2_dphi_j1met","c2_dphi_j1met",64,0.,3.2);
    h.c_dphi_j2met[1] = new TH1F("c2_dphi_j2met","c2_dphi_j2met",64,0.,3.2);
    h.c_dr_jt_lep1[1] = new TH1F("c2_dr_jt_lep1","c2_dr_jt_lep1",100,0.,5.);
    h.c_dr_jt_lep2[1] = new TH1F("c2_dr_jt_lep2","c2_dr_jt_lep2",100,0.,5.);
    h.c_dr_jt_lep3[1] = new TH1F("c2_dr_jt_lep3","c2_dr_jt_lep3",100,0.,5.);

    h.c_met_m_os_1[1] = new TH2F("c2_met_m_os_1","c2_met_m_os_1",50,0.,200.,50,0.,200.);
    h.c_trk_met[1] = new TH1F("c2_trk_met","c2_trk_met",1000,0.,1000.);
    h.c_st_met[1] = new TH1F("c2_st_met","c2_st_met",500,0.,500.);
    h.c_njet15[1] = new TH1F("c2_njet15","c2_njet15",5,0.,5.);
// ---------------------------------------
// c3
// ---------------------------------------
    h.c_trilep_pt[2][0] = new TH1F("c3_trilep_pt1","c3_trilep pt 1",100,0.,100.);
    h.c_trilep_pt[2][1] = new TH1F("c3_trilep_pt2","c3_trilep pt 2",100,0.,100.);
    h.c_trilep_pt[2][2] = new TH1F("c3_trilep_pt3","c3_trilep pt 3",100,0.,100.);
    h.c_trilep_eta[2][0] = new TH1F("c3_trilep_eta1","c3_trilep eta 1",120,-3.,3.);
    h.c_trilep_eta[2][1] = new TH1F("c3_trilep_eta2","c3_trilep eta 2",120,-3.,3.);
    h.c_trilep_eta[2][2] = new TH1F("c3_trilep_eta3","c3_trilep eta 3",120,-3.,3.);
    h.c_trilep_phi[2][0] = new TH1F("c3_trilep_phi1","c3_trilep phi 1",128,0.,6.4);
    h.c_trilep_phi[2][1] = new TH1F("c3_trilep_phi2","c3_trilep phi 2",128,0.,6.4);
    h.c_trilep_phi[2][2] = new TH1F("c3_trilep_phi3","c3_trilep phi 3",128,0.,6.4);

    h.c_met[2] = new TH1F("c3_met","c3_met",400,0.,400.);
    h.c_st[2] = new TH1F("c3_st","c3_st",400,0.,400.);
    h.c_st2[2] = new TH1F("c3_st2","c3_st2",400,0.,400.);
    h.c_st2_over_st[2] = new TH1F("c3_st2_over_st","c3_st2_over_st",100,0.,1.);
    h.c_ht[2] = new TH1F("c3_ht","c3_ht",200,0.,200.);
    h.c_ht2[2] = new TH1F("c3_ht2","c3_ht2",200,0.,200.);
    h.c_ht2_over_ht[2] = new TH1F("c3_ht2_over_ht","c3_ht2_over_ht",100,0.,1.);

    h.c_m_leps[2] = new TH1F("c3_m_leps","c3_m_leps",400,0.,400.);
    h.c_m_leps_jets[2] = new TH1F("c3_m_leps_jets","c3_m_leps_jets",400,0.,400.);
    h.c_m_all[2] = new TH1F("c3_m_all","c3_m_all",400,0.,400.);
    h.c_m_os_1[2] = new TH1F("c3_m_os_1","c3_m_os_1",200,0.,200.); // highest opposite sign mass
    h.c_m_os_2[2] = new TH1F("c3_m_os_2","c3_m_os_2",200,0.,200.);
    h.c_m_ss[2] = new TH1F("c3_m_ss","c3_m_ss",200,0.,200.);
    h.c_m_ss_1[2] = new TH1F("c3_m_ss_1","c3_m_ss_1",200,0.,200.); // highest same sign mass
    h.c_m_ss_2[2] = new TH1F("c3_m_ss_2","c3_m_ss_2",200,0.,200.);
    h.c_m_ss_3[2] = new TH1F("c3_m_ss_3","c3_m_ss_3",200,0.,200.);

    h.c_mt_all[2] = new TH1F("c3_mt_all","c3_mt_all",500,0.,500.);
    h.c_mt_lep1_met[2] = new TH1F("c3_mt_lep1_met","c3_mt_lep1_met",200,0.,200.);
    h.c_mt_lep2_met[2] = new TH1F("c3_mt_lep2_met","c3_mt_lep2_met",200,0.,200.);
    h.c_mt_lep3_met[2] = new TH1F("c3_mt_lep3_met","c3_mt_lep3_met",200,0.,200.);
    h.c_mt_lep_met_max[2] = new TH1F("c3_mt_lep_met_max","c3_mt_lep_met_max",200,0.,200.);
    h.c_mt_lep_met_min[2] = new TH1F("c3_mt_lep_met_min","c3_mt_lep_met_min",200,0.,200.);

    h.c_sphericity[2] = new TH1F("c3_sphericity","c3_sphericity",40,0.,2.);
    h.c_aplanarity[2] = new TH1F("c3_aplanarity","c3_aplanarity",40,0.,2.);
    h.c_thrust[2] = new TH1F("c3_thrust","c3_thrust",50,0.,1.);
    h.c_sphericity_2d[2] = new TH1F("c3_sphericity_2d","c3_sphericity_2d",40,0.,2.);
    h.c_aplanarity_2d[2] = new TH1F("c3_aplanarity_2d","c3_aplanarity_2d",40,0.,2.);
    h.c_thrust_2d[2] = new TH1F("c3_thrust_2d","c3_thrust_2d",50,0.,1.);

    h.c_dphi_lep12[2] = new TH1F("c3_dphi_lep12","c3_dphi_lep12",64,0.,3.2);
    h.c_dphi_lep13[2] = new TH1F("c3_dphi_lep13","c3_dphi_lep13",64,0.,3.2);
    h.c_dphi_lep23[2] = new TH1F("c3_dphi_lep23","c3_dphi_lep23",64,0.,3.2);
    h.c_dphi_lep1met[2] = new TH1F("c3_dphi_lep1met","c3_dphi_lep1met",64,0.,3.2);
    h.c_dphi_lep2met[2] = new TH1F("c3_dphi_lep2met","c3_dphi_lep2met",64,0.,3.2);
    h.c_dphi_lep3met[2] = new TH1F("c3_dphi_lep3met","c3_dphi_lep3met",64,0.,3.2);
    h.c_dr_lep12[2] = new TH1F("c3_dr_lep12","c3_dr_lep12",100,0.,5.);
    h.c_dr_lep13[2] = new TH1F("c3_dr_lep13","c3_dr_lep13",100,0.,5.);
    h.c_dr_lep23[2] = new TH1F("c3_dr_lep23","c3_dr_lep23",100,0.,5.);

    h.c_dphi_j1met[2] = new TH1F("c3_dphi_j1met","c3_dphi_j1met",64,0.,3.2);
    h.c_dphi_j2met[2] = new TH1F("c3_dphi_j2met","c3_dphi_j2met",64,0.,3.2);
    h.c_dr_jt_lep1[2] = new TH1F("c3_dr_jt_lep1","c3_dr_jt_lep1",100,0.,5.);
    h.c_dr_jt_lep2[2] = new TH1F("c3_dr_jt_lep2","c3_dr_jt_lep2",100,0.,5.);
    h.c_dr_jt_lep3[2] = new TH1F("c3_dr_jt_lep3","c3_dr_jt_lep3",100,0.,5.);

    h.c_met_m_os_1[2] = new TH2F("c3_met_m_os_1","c3_met_m_os_1",50,0.,200.,50,0.,200.);
    h.c_trk_met[2] = new TH1F("c3_trk_met","c3_trk_met",1000,0.,1000.);
    h.c_st_met[2] = new TH1F("c3_st_met","c3_st_met",500,0.,500.);
    h.c_njet15[2] = new TH1F("c3_njet15","c3_njet15",5,0.,5.);
// ---------------------------------------
// c4
// ---------------------------------------
    h.c_trilep_pt[3][0] = new TH1F("c4_trilep_pt1","c4_trilep pt 1",100,0.,100.);
    h.c_trilep_pt[3][1] = new TH1F("c4_trilep_pt2","c4_trilep pt 2",100,0.,100.);
    h.c_trilep_pt[3][2] = new TH1F("c4_trilep_pt3","c4_trilep pt 3",100,0.,100.);
    h.c_trilep_eta[3][0] = new TH1F("c4_trilep_eta1","c4_trilep eta 1",120,-3.,3.);
    h.c_trilep_eta[3][1] = new TH1F("c4_trilep_eta2","c4_trilep eta 2",120,-3.,3.);
    h.c_trilep_eta[3][2] = new TH1F("c4_trilep_eta3","c4_trilep eta 3",120,-3.,3.);
    h.c_trilep_phi[3][0] = new TH1F("c4_trilep_phi1","c4_trilep phi 1",128,0.,6.4);
    h.c_trilep_phi[3][1] = new TH1F("c4_trilep_phi2","c4_trilep phi 2",128,0.,6.4);
    h.c_trilep_phi[3][2] = new TH1F("c4_trilep_phi3","c4_trilep phi 3",128,0.,6.4);

    h.c_met[3] = new TH1F("c4_met","c4_met",400,0.,400.);
    h.c_st[3] = new TH1F("c4_st","c4_st",400,0.,400.);
    h.c_st2[3] = new TH1F("c4_st2","c4_st2",400,0.,400.);
    h.c_st2_over_st[3] = new TH1F("c4_st2_over_st","c4_st2_over_st",100,0.,1.);
    h.c_ht[3] = new TH1F("c4_ht","c4_ht",200,0.,200.);
    h.c_ht2[3] = new TH1F("c4_ht2","c4_ht2",200,0.,200.);
    h.c_ht2_over_ht[3] = new TH1F("c4_ht2_over_ht","c4_ht2_over_ht",100,0.,1.);

    h.c_m_leps[3] = new TH1F("c4_m_leps","c4_m_leps",400,0.,400.);
    h.c_m_leps_jets[3] = new TH1F("c4_m_leps_jets","c4_m_leps_jets",400,0.,400.);
    h.c_m_all[3] = new TH1F("c4_m_all","c4_m_all",400,0.,400.);
    h.c_m_os_1[3] = new TH1F("c4_m_os_1","c4_m_os_1",200,0.,200.); // highest opposite sign mass
    h.c_m_os_2[3] = new TH1F("c4_m_os_2","c4_m_os_2",200,0.,200.);
    h.c_m_ss[3] = new TH1F("c4_m_ss","c4_m_ss",200,0.,200.);
    h.c_m_ss_1[3] = new TH1F("c4_m_ss_1","c4_m_ss_1",200,0.,200.); // highest same sign mass
    h.c_m_ss_2[3] = new TH1F("c4_m_ss_2","c4_m_ss_2",200,0.,200.);
    h.c_m_ss_3[3] = new TH1F("c4_m_ss_3","c4_m_ss_3",200,0.,200.);

    h.c_mt_all[3] = new TH1F("c4_mt_all","c4_mt_all",500,0.,500.);
    h.c_mt_lep1_met[3] = new TH1F("c4_mt_lep1_met","c4_mt_lep1_met",200,0.,200.);
    h.c_mt_lep2_met[3] = new TH1F("c4_mt_lep2_met","c4_mt_lep2_met",200,0.,200.);
    h.c_mt_lep3_met[3] = new TH1F("c4_mt_lep3_met","c4_mt_lep3_met",200,0.,200.);
    h.c_mt_lep_met_max[3] = new TH1F("c4_mt_lep_met_max","c4_mt_lep_met_max",200,0.,200.);
    h.c_mt_lep_met_min[3] = new TH1F("c4_mt_lep_met_min","c4_mt_lep_met_min",200,0.,200.);

    h.c_sphericity[3] = new TH1F("c4_sphericity","c4_sphericity",40,0.,2.);
    h.c_aplanarity[3] = new TH1F("c4_aplanarity","c4_aplanarity",40,0.,2.);
    h.c_thrust[3] = new TH1F("c4_thrust","c4_thrust",50,0.,1.);
    h.c_sphericity_2d[3] = new TH1F("c4_sphericity_2d","c4_sphericity_2d",40,0.,2.);
    h.c_aplanarity_2d[3] = new TH1F("c4_aplanarity_2d","c4_aplanarity_2d",40,0.,2.);
    h.c_thrust_2d[3] = new TH1F("c4_thrust_2d","c4_thrust_2d",50,0.,1.);

    h.c_dphi_lep12[3] = new TH1F("c4_dphi_lep12","c4_dphi_lep12",64,0.,3.2);
    h.c_dphi_lep13[3] = new TH1F("c4_dphi_lep13","c4_dphi_lep13",64,0.,3.2);
    h.c_dphi_lep23[3] = new TH1F("c4_dphi_lep23","c4_dphi_lep23",64,0.,3.2);
    h.c_dphi_lep1met[3] = new TH1F("c4_dphi_lep1met","c4_dphi_lep1met",64,0.,3.2);
    h.c_dphi_lep2met[3] = new TH1F("c4_dphi_lep2met","c4_dphi_lep2met",64,0.,3.2);
    h.c_dphi_lep3met[3] = new TH1F("c4_dphi_lep3met","c4_dphi_lep3met",64,0.,3.2);
    h.c_dr_lep12[3] = new TH1F("c4_dr_lep12","c4_dr_lep12",100,0.,5.);
    h.c_dr_lep13[3] = new TH1F("c4_dr_lep13","c4_dr_lep13",100,0.,5.);
    h.c_dr_lep23[3] = new TH1F("c4_dr_lep23","c4_dr_lep23",100,0.,5.);

    h.c_dphi_j1met[3] = new TH1F("c4_dphi_j1met","c4_dphi_j1met",64,0.,3.2);
    h.c_dphi_j2met[3] = new TH1F("c4_dphi_j2met","c4_dphi_j2met",64,0.,3.2);
    h.c_dr_jt_lep1[3] = new TH1F("c4_dr_jt_lep1","c4_dr_jt_lep1",100,0.,5.);
    h.c_dr_jt_lep2[3] = new TH1F("c4_dr_jt_lep2","c4_dr_jt_lep2",100,0.,5.);
    h.c_dr_jt_lep3[3] = new TH1F("c4_dr_jt_lep3","c4_dr_jt_lep3",100,0.,5.);

    h.c_met_m_os_1[3] = new TH2F("c4_met_m_os_1","c4_met_m_os_1",50,0.,200.,50,0.,200.);
    h.c_trk_met[3] = new TH1F("c4_trk_met","c4_trk_met",1000,0.,1000.);
    h.c_st_met[3] = new TH1F("c4_st_met","c4_st_met",500,0.,500.);
    h.c_njet15[3] = new TH1F("c4_njet15","c4_njet15",5,0.,5.);
// ---------------------------------------
// Veto HiPt selection
// ---------------------------------------
    h.c_trilep_pt[4][0] = new TH1F("veto_trilep_pt1","veto_trilep pt 1",100,0.,100.);
    h.c_trilep_pt[4][1] = new TH1F("veto_trilep_pt2","veto_trilep pt 2",100,0.,100.);
    h.c_trilep_pt[4][2] = new TH1F("veto_trilep_pt3","veto_trilep pt 3",100,0.,100.);
    h.c_trilep_eta[4][0] = new TH1F("veto_trilep_eta1","veto_trilep eta 1",120,-3.,3.);
    h.c_trilep_eta[4][1] = new TH1F("veto_trilep_eta2","veto_trilep eta 2",120,-3.,3.);
    h.c_trilep_eta[4][2] = new TH1F("veto_trilep_eta3","veto_trilep eta 3",120,-3.,3.);
    h.c_trilep_phi[4][0] = new TH1F("veto_trilep_phi1","veto_trilep phi 1",128,0.,6.4);
    h.c_trilep_phi[4][1] = new TH1F("veto_trilep_phi2","veto_trilep phi 2",128,0.,6.4);
    h.c_trilep_phi[4][2] = new TH1F("veto_trilep_phi3","veto_trilep phi 3",128,0.,6.4);

    h.c_met[4] = new TH1F("veto_met","veto_met",400,0.,400.);
    h.c_st[4] = new TH1F("veto_st","veto_st",400,0.,400.);
    h.c_st2[4] = new TH1F("veto_st2","veto_st2",400,0.,400.);
    h.c_st2_over_st[4] = new TH1F("veto_st2_over_st","veto_st2_over_st",100,0.,1.);
    h.c_ht[4] = new TH1F("veto_ht","veto_ht",200,0.,200.);
    h.c_ht2[4] = new TH1F("veto_ht2","veto_ht2",200,0.,200.);
    h.c_ht2_over_ht[4] = new TH1F("veto_ht2_over_ht","veto_ht2_over_ht",100,0.,1.);

    h.c_m_leps[4] = new TH1F("veto_m_leps","veto_m_leps",400,0.,400.);
    h.c_m_leps_jets[4] = new TH1F("veto_m_leps_jets","veto_m_leps_jets",400,0.,400.);
    h.c_m_all[4] = new TH1F("veto_m_all","veto_m_all",400,0.,400.);
    h.c_m_os_1[4] = new TH1F("veto_m_os_1","veto_m_os_1",200,0.,200.); // highest opposite sign mass
    h.c_m_os_2[4] = new TH1F("veto_m_os_2","veto_m_os_2",200,0.,200.);
    h.c_m_ss[4] = new TH1F("veto_m_ss","veto_m_ss",200,0.,200.);
    h.c_m_ss_1[4] = new TH1F("veto_m_ss_1","veto_m_ss_1",200,0.,200.); // highest same sign mass
    h.c_m_ss_2[4] = new TH1F("veto_m_ss_2","veto_m_ss_2",200,0.,200.);
    h.c_m_ss_3[4] = new TH1F("veto_m_ss_3","veto_m_ss_3",200,0.,200.);

    h.c_mt_all[4] = new TH1F("veto_mt_all","veto_mt_all",500,0.,500.);
    h.c_mt_lep1_met[4] = new TH1F("veto_mt_lep1_met","veto_mt_lep1_met",200,0.,200.);
    h.c_mt_lep2_met[4] = new TH1F("veto_mt_lep2_met","veto_mt_lep2_met",200,0.,200.);
    h.c_mt_lep3_met[4] = new TH1F("veto_mt_lep3_met","veto_mt_lep3_met",200,0.,200.);
    h.c_mt_lep_met_max[4] = new TH1F("veto_mt_lep_met_max","veto_mt_lep_met_max",200,0.,200.);
    h.c_mt_lep_met_min[4] = new TH1F("veto_mt_lep_met_min","veto_mt_lep_met_min",200,0.,200.);

    h.c_sphericity[4] = new TH1F("veto_sphericity","veto_sphericity",40,0.,2.);
    h.c_aplanarity[4] = new TH1F("veto_aplanarity","veto_aplanarity",40,0.,2.);
    h.c_thrust[4] = new TH1F("veto_thrust","veto_thrust",50,0.,1.);
    h.c_sphericity_2d[4] = new TH1F("veto_sphericity_2d","veto_sphericity_2d",40,0.,2.);
    h.c_aplanarity_2d[4] = new TH1F("veto_aplanarity_2d","veto_aplanarity_2d",40,0.,2.);
    h.c_thrust_2d[4] = new TH1F("veto_thrust_2d","veto_thrust_2d",50,0.,1.);

    h.c_dphi_lep12[4] = new TH1F("veto_dphi_lep12","veto_dphi_lep12",64,0.,3.2);
    h.c_dphi_lep13[4] = new TH1F("veto_dphi_lep13","veto_dphi_lep13",64,0.,3.2);
    h.c_dphi_lep23[4] = new TH1F("veto_dphi_lep23","veto_dphi_lep23",64,0.,3.2);
    h.c_dphi_lep1met[4] = new TH1F("veto_dphi_lep1met","veto_dphi_lep1met",64,0.,3.2);
    h.c_dphi_lep2met[4] = new TH1F("veto_dphi_lep2met","veto_dphi_lep2met",64,0.,3.2);
    h.c_dphi_lep3met[4] = new TH1F("veto_dphi_lep3met","veto_dphi_lep3met",64,0.,3.2);
    h.c_dr_lep12[4] = new TH1F("veto_dr_lep12","veto_dr_lep12",100,0.,5.);
    h.c_dr_lep13[4] = new TH1F("veto_dr_lep13","veto_dr_lep13",100,0.,5.);
    h.c_dr_lep23[4] = new TH1F("veto_dr_lep23","veto_dr_lep23",100,0.,5.);

    h.c_dphi_j1met[4] = new TH1F("veto_dphi_j1met","veto_dphi_j1met",64,0.,3.2);
    h.c_dphi_j2met[4] = new TH1F("veto_dphi_j2met","veto_dphi_j2met",64,0.,3.2);
    h.c_dr_jt_lep1[4] = new TH1F("veto_dr_jt_lep1","veto_dr_jt_lep1",100,0.,5.);
    h.c_dr_jt_lep2[4] = new TH1F("veto_dr_jt_lep2","veto_dr_jt_lep2",100,0.,5.);
    h.c_dr_jt_lep3[4] = new TH1F("veto_dr_jt_lep3","veto_dr_jt_lep3",100,0.,5.);

    h.c_met_m_os_1[4] = new TH2F("veto_met_m_os_1","veto_met_m_os_1",50,0.,200.,50,0.,200.);
    h.c_trk_met[4] = new TH1F("veto_trk_met","veto_trk_met",1000,0.,1000.);
    h.c_st_met[4] = new TH1F("veto_st_met","veto_st_met",500,0.,500.);
    h.c_njet15[4] = new TH1F("veto_njet15","veto_njet15",5,0.,5.);


    h.jet20_et = new TH1F("jet20_et","jet20_et",100,0.,100.);
    h.jet20_eta = new TH1F("jet20_eta","jet20_eta",100,-2.5,2.5);
    h.jet20_phi = new TH1F("jet20_phi","jet20_phi",160,0.,6.4);
    h.jet20_n = new TH1F("jet20_n","jet20_n",5,0.,5.);
    h.jet_scaleEta = new TH2F("jet_scaleEta","jet_scaleEta",50,-2.5,2.5,75,0.5,2.0);
    h.beam_xy = new TH2F("Beam_XY","Beam_XY",100,-5.,5.,100,-5.,5.);
    h.vz = new TH1F("ZVertex","ZVertex",200,-100.,100.);

    h.dilep_pt[0] = new TH1F("dilep_pt1","dilep pt 1",100,0.,100.);
    h.dilep_pt[1] = new TH1F("dilep_pt2","dilep pt 2",100,0.,100.);
    h.dilep_eta[0] = new TH1F("dilep_eta1","dilep eta 1",100,-2.,2.);
    h.dilep_eta[1] = new TH1F("dilep_eta2","dilep eta 2",100,-2.,2.);
    h.dilep_phi[0] = new TH1F("dilep_phi1","dilep phi 1",100,0.,100.);
    h.dilep_phi[1] = new TH1F("dilep_phi2","dilep phi 2",100,0.,100.);
    h.trilep_pt[0] = new TH1F("trilep_pt1","trilep pt 1",100,0.,100.);
    h.trilep_pt[1] = new TH1F("trilep_pt2","trilep pt 2",100,0.,100.);
    h.trilep_pt[2] = new TH1F("trilep_pt3","trilep pt 3",100,0.,100.);
    h.trilep_eta[0] = new TH1F("trilep_eta1","trilep eta 1",120,-3.,3.);
    h.trilep_eta[1] = new TH1F("trilep_eta2","trilep eta 2",120,-3.,3.);
    h.trilep_eta[2] = new TH1F("trilep_eta3","trilep eta 3",120,-3.,3.);
    h.trilep_phi[0] = new TH1F("trilep_phi1","trilep phi 1",128,0.,6.4);
    h.trilep_phi[1] = new TH1F("trilep_phi2","trilep phi 2",128,0.,6.4);
    h.trilep_phi[2] = new TH1F("trilep_phi3","trilep phi 3",128,0.,6.4);

    h.vz_lep1 = new TH1F("vz_lep1","vz_lep1",200,-10.,10.);
    h.vz_lep2 = new TH1F("vz_lep2","vz_lep2",200,-10.,10.);
    h.vz_lep3 = new TH1F("vz_lep3","vz_lep3",200,-10.,10.);

//  Advanced variables
    h.st = new TH1F("st","st",400,0.,400.);
    h.st2 = new TH1F("st2","st2",400,0.,400.);
    h.st2_over_st = new TH1F("st2_over_st","st2_over_st",100,0.,1.);
    h.ht = new TH1F("ht","ht",200,0.,200.);
    h.ht2 = new TH1F("ht2","ht2",200,0.,200.);
    h.ht2_over_ht = new TH1F("ht2_over_ht","ht2_over_ht",100,0.,1.);
    h.met = new TH1F("met","met",200,0.,200.);
    h.met_phi = new TH1F("met_phi","met_phi",160,0.,6.4);
    h.trk_met = new TH1F("trk_met","trk_met",1000,0.,1000.);
    h.st_met = new TH1F("st_met","st_met",500,0.,500.);
    h.met_trkpt = new TH2F("met_trkpt","met_trkpt",100,0.,100.,100,0.,100.);
    h.m_leps = new TH1F("m_leps","m_leps",400,0.,400.);
    h.m_leps_jets = new TH1F("m_leps_jets","m_leps_jets",400,0.,400.);
    h.m_all = new TH1F("m_all","m_all",400,0.,400.);
    h.m_os_1 = new TH1F("m_os_1","m_os_1",200,0.,200.); // highest opposite sign mass
    h.m_os_2 = new TH1F("m_os_2","m_os_2",200,0.,200.);
    h.m_ss   = new TH1F("m_ss","m_ss",200,0.,200.);
    h.m_ss_1 = new TH1F("m_ss_1","m_ss_1",200,0.,200.); // highest same sign mass
    h.m_ss_2 = new TH1F("m_ss_2","m_ss_2",200,0.,200.);
    h.m_ss_3 = new TH1F("m_ss_3","m_ss_3",200,0.,200.);
    h.mt_all = new TH1F("mt_all","mt_all",500,0.,500.);
    h.mt_lep1_met = new TH1F("mt_lep1_met","mt_lep1_met",200,0.,200.);
    h.mt_lep2_met = new TH1F("mt_lep2_met","mt_lep2_met",200,0.,200.);
    h.mt_lep3_met = new TH1F("mt_lep3_met","mt_lep3_met",200,0.,200.);
    h.mt_lep_met_max = new TH1F("mt_lep_met_max","mt_lep_met_max",200,0.,200.);
    h.mt_lep_met_min = new TH1F("mt_lep_met_min","mt_lep_met_min",200,0.,200.);
    h.met_m_os_1 = new TH2F("met_m_os_1","met_m_os_1",50,0.,200.,50,0.,200.);
    
    h.sphericity = new TH1F("sphericity","sphericity",40,0.,2.);
    h.aplanarity = new TH1F("aplanarity","aplanarity",40,0.,2.);
    h.thrust = new TH1F("thrust","thrust",50,0.,1.);
    h.sphericity_2d = new TH1F("sphericity_2d","sphericity_2d",40,0.,2.);
    h.aplanarity_2d = new TH1F("aplanarity_2d","aplanarity_2d",40,0.,2.);
    h.thrust_2d = new TH1F("thrust_2d","thrust_2d",50,0.,1.);
 
    h.dphi_lep12 = new TH1F("dphi_lep12","dphi_lep12",64,0.,3.2);
    h.dphi_lep13 = new TH1F("dphi_lep13","dphi_lep13",64,0.,3.2);
    h.dphi_lep23 = new TH1F("dphi_lep23","dphi_lep23",64,0.,3.2);
    h.dphi_lep1met = new TH1F("dphi_lep1met","dphi_lep1met",64,0.,3.2);
    h.dphi_lep2met = new TH1F("dphi_lep2met","dphi_lep2met",64,0.,3.2);
    h.dphi_lep3met = new TH1F("dphi_lep3met","dphi_lep3met",64,0.,3.2);
    h.dr_lep12 = new TH1F("dr_lep12","dr_lep12",100,0.,5.);
    h.dr_lep13 = new TH1F("dr_lep13","dr_lep13",100,0.,5.);
    h.dr_lep23 = new TH1F("dr_lep23","dr_lep23",100,0.,5.);
    h.dphi_j1met_j15 = new TH1F("dphi_j1met_j15","dphi_j1met_j15",64,0.,3.2);
    h.dphi_j1met = new TH1F("dphi_j1met","dphi_j1met",64,0.,3.2);
    h.dphi_j2met = new TH1F("dphi_j2met","dphi_j2met",64,0.,3.2);
    h.dr_jt_lep1 = new TH1F("dr_jt_lep1","dr_jt_lep1",100,0.,5.);
    h.dr_jt_lep2 = new TH1F("dr_jt_lep2","dr_jt_lep2",100,0.,5.);
    h.dr_jt_lep3 = new TH1F("dr_jt_lep3","dr_jt_lep3",100,0.,5.);

//=====================================================================
//  OptRegion variables
    h.n1_opt_dphi_lep12 = new TH1F("n1_opt_dphi_lep12","n1_opt_dphi_lep12",64,0.,3.2);
    h.n1_opt_met = new TH1F("n1_opt_met","n1_opt_met",200,0.,200.);
    h.n1_opt_m_os_1 = new TH1F("n1_opt_m_os_1","n1_opt_m_os_1",200,0.,200.);
    h.n1_opt_m_os_2 = new TH1F("n1_opt_m_os_2","n1_opt_m_os_2",200,0.,200.);
    h.n1_opt_ht = new TH1F("n1_opt_ht","n1_opt_ht",200,0.,200.);
    h.n1_opt_trk_met = new TH1F("n1_opt_trk_met","n1_opt_trk_met",1000,0.,1000.);
    h.n1_opt_mt_min = new TH1F("n1_opt_mt_min","n1_opt_mt_min",150,0.,150.);
    h.n1_opt_met_trkpt = new TH2F("n1_opt_met_trkpt","n1_opt_met_trkpt",100,0.,100.,100,0.,100.);
    h.n1_opt_ele1_et = new TH1F("n1_opt_ele1_et","n1_opt_ele1_et",200,0,100);
    h.n1_opt_ele2_et = new TH1F("n1_opt_ele2_et","n1_opt_ele2_et",200,0,100);
    h.n1_opt_trk_pt = new TH1F("n1_opt_trk_pt","n1_opt_trk_pt",200,0,100);


    h.opt_met = new TH1F("opt_met","opt_met",400,0.,400.);
    h.opt_st = new TH1F("opt_st","opt_st",400,0.,400.);
    h.opt_st2 = new TH1F("opt_st2","opt_st2",400,0.,400.);
    h.opt_st2_over_st = new TH1F("opt_st2_over_st","opt_st2_over_st",100,0.,1.);
    h.opt_ht = new TH1F("opt_ht","opt_ht",200,0.,200.);
    h.opt_ht2 = new TH1F("opt_ht2","opt_ht2",200,0.,200.);
    h.opt_ht2_over_ht = new TH1F("opt_ht2_over_ht","opt_ht2_over_ht",100,0.,1.);

    h.opt_m_leps = new TH1F("opt_m_leps","opt_m_leps",400,0.,400.);
    h.opt_m_leps_jets = new TH1F("opt_m_leps_jets","opt_m_leps_jets",400,0.,400.);
    h.opt_m_all = new TH1F("opt_m_all","opt_m_all",400,0.,400.);
    h.opt_m_os_1 = new TH1F("opt_m_os_1","opt_m_os_1",200,0.,200.); // highest opposite sign mass
    h.opt_m_os_2 = new TH1F("opt_m_os_2","opt_m_os_2",200,0.,200.);
    h.opt_m_ss   = new TH1F("opt_m_ss","opt_m_ss",200,0.,200.);
    h.opt_m_ss_1 = new TH1F("opt_m_ss_1","opt_m_ss_1",200,0.,200.); // highest same sign mass
    h.opt_m_ss_2 = new TH1F("opt_m_ss_2","opt_m_ss_2",200,0.,200.);
    h.opt_m_ss_3 = new TH1F("opt_m_ss_3","opt_m_ss_3",200,0.,200.);

    h.opt_mt_all = new TH1F("opt_mt_all","opt_mt_all",500,0.,500.);
    h.opt_mt_lep1_met = new TH1F("opt_mt_lep1_met","opt_mt_lep1_met",200,0.,200.);
    h.opt_mt_lep2_met = new TH1F("opt_mt_lep2_met","opt_mt_lep2_met",200,0.,200.);
    h.opt_mt_lep3_met = new TH1F("opt_mt_lep3_met","opt_mt_lep3_met",200,0.,200.);
    h.opt_mt_lep_met_max = new TH1F("opt_mt_lep_met_max","opt_mt_lep_met_max",200,0.,200.);
    h.opt_mt_lep_met_min = new TH1F("opt_mt_lep_met_min","opt_mt_lep_met_min",200,0.,200.);

    h.opt_sphericity = new TH1F("opt_sphericity","opt_sphericity",40,0.,2.);
    h.opt_aplanarity = new TH1F("opt_aplanarity","opt_aplanarity",40,0.,2.);
    h.opt_thrust = new TH1F("opt_thrust","opt_thrust",50,0.,1.);
    h.opt_sphericity_2d = new TH1F("opt_sphericity_2d","opt_sphericity_2d",40,0.,2.);
    h.opt_aplanarity_2d = new TH1F("opt_aplanarity_2d","opt_aplanarity_2d",40,0.,2.);
    h.opt_thrust_2d = new TH1F("opt_thrust_2d","opt_thrust_2d",50,0.,1.);

    h.opt_met_m_os_1 = new TH2F("opt_met_m_os_1","opt_met_m_os_1",50,0.,200.,50,0.,200.);
    h.opt_trk_met = new TH1F("opt_trk_met","opt_trk_met",1000,0.,1000.);
    h.opt_st_met = new TH1F("opt_st_met","opt_st_met",500,0.,500.);
    h.opt_njet15 = new TH1F("opt_njet15","opt_njet15",5,0.,5.);

    h.fopt_met = new TH1F("fopt_met","fopt_met",200,0,100);
    h.fopt_mos1 = new TH1F("fopt_mos1","fopt_mos1",300,0,150);
    h.fopt_mos2 = new TH1F("fopt_mos2","fopt_mos2",300,0,150);
    h.fopt_ht = new TH1F("fopt_ht","fopt_ht",200,0,200);
    h.fopt_dphi = new TH1F("fopt_dphi","fopt_dphi",64,0,3.2);
    h.fopt_mt = new TH1F("fopt_mt","fopt_mt",150,0,150);
    h.fopt_et[0] = new TH1F("fopt_et0","fopt_et0",200,0,100);
    h.fopt_et[1] = new TH1F("fopt_et1","fopt_et1",200,0,100);
    h.fopt_et[2] = new TH1F("fopt_et2","fopt_et2",200,0,100);
    h.fopt_n1_met = new TH1F("fopt_n1_met","fopt_n1_met",200,0,100);
    h.fopt_n1_mos1 = new TH1F("fopt_n1_mos1","fopt_n1_mos1",300,0,150);
    h.fopt_n1_mos2 = new TH1F("fopt_n1_mos2","fopt_n1_mos2",300,0,150);
    h.fopt_n1_ht = new TH1F("fopt_n1_ht","fopt_n1_ht",200,0,200);
    h.fopt_n1_dphi = new TH1F("fopt_n1_dphi","fopt_n1_dphi",64,0,3.2);
    h.fopt_n1_mt = new TH1F("fopt_n1_mt","fopt_n1_mt",150,0,150);
    h.fopt_n1_et[0] = new TH1F("fopt_n1_et0","fopt_n1_et0",200,0,100);
    h.fopt_n1_et[1] = new TH1F("fopt_n1_et1","fopt_n1_et1",200,0,100);
    h.fopt_n1_et[2] = new TH1F("fopt_n1_et2","fopt_n1_et2",200,0,100);


    AddHistogram(h.c_trilep_pt[0][0]);
    AddHistogram(h.c_trilep_pt[0][1]);
    AddHistogram(h.c_trilep_pt[0][2]);
    AddHistogram(h.c_trilep_eta[0][0]);
    AddHistogram(h.c_trilep_eta[0][1]);
AddHistogram(h.c_trilep_eta[0][2]);
AddHistogram(h.c_trilep_phi[0][0]);
AddHistogram(h.c_trilep_phi[0][1]);
AddHistogram(h.c_trilep_phi[0][2]);
AddHistogram(h.c_met[0]);
AddHistogram(h.c_st[0]);
AddHistogram(h.c_st2[0]);
AddHistogram(h.c_st2_over_st[0]);
AddHistogram(h.c_ht[0]);
AddHistogram(h.c_ht2[0]);
AddHistogram(h.c_ht2_over_ht[0]);
AddHistogram(h.c_m_leps[0]);
AddHistogram(h.c_m_leps_jets[0]);
AddHistogram(h.c_m_all[0]);
AddHistogram(h.c_m_os_1[0]);
AddHistogram(h.c_m_os_2[0]);
AddHistogram(h.c_m_ss[0]);
AddHistogram(h.c_m_ss_1[0]);
AddHistogram(h.c_m_ss_2[0]);
AddHistogram(h.c_m_ss_3[0]);
AddHistogram(h.c_mt_all[0]);
AddHistogram(h.c_mt_lep1_met[0]);
AddHistogram(h.c_mt_lep2_met[0]);
AddHistogram(h.c_mt_lep3_met[0]);
AddHistogram(h.c_mt_lep_met_max[0]);
AddHistogram(h.c_mt_lep_met_min[0]);
AddHistogram(h.c_sphericity[0]);
AddHistogram(h.c_aplanarity[0]);
AddHistogram(h.c_thrust[0]);
AddHistogram(h.c_sphericity_2d[0]);
AddHistogram(h.c_aplanarity_2d[0]);
AddHistogram(h.c_thrust_2d[0]);
AddHistogram(h.c_dphi_lep12[0]);
AddHistogram(h.c_dphi_lep13[0]);
AddHistogram(h.c_dphi_lep23[0]);
AddHistogram(h.c_dphi_lep1met[0]);
AddHistogram(h.c_dphi_lep2met[0]);
AddHistogram(h.c_dphi_lep3met[0]);
AddHistogram(h.c_dr_lep12[0]);
AddHistogram(h.c_dr_lep13[0]);
AddHistogram(h.c_dr_lep23[0]);
AddHistogram(h.c_dphi_j1met[0]);
AddHistogram(h.c_dphi_j2met[0]);
AddHistogram(h.c_dr_jt_lep1[0]);
AddHistogram(h.c_dr_jt_lep2[0]);
AddHistogram(h.c_dr_jt_lep3[0]);
AddHistogram(h.c_met_m_os_1[0]);
AddHistogram(h.c_trk_met[0]);
AddHistogram(h.c_st_met[0]);
AddHistogram(h.c_njet15[0]);
AddHistogram(h.c_trilep_pt[1][0]);
AddHistogram(h.c_trilep_pt[1][1]);
AddHistogram(h.c_trilep_pt[1][2]);
AddHistogram(h.c_trilep_eta[1][0]);
AddHistogram(h.c_trilep_eta[1][1]);
AddHistogram(h.c_trilep_eta[1][2]);
AddHistogram(h.c_trilep_phi[1][0]);
AddHistogram(h.c_trilep_phi[1][1]);
AddHistogram(h.c_trilep_phi[1][2]);
AddHistogram(h.c_met[1]);
AddHistogram(h.c_st[1]);
AddHistogram(h.c_st2[1]);
AddHistogram(h.c_st2_over_st[1]);
AddHistogram(h.c_ht[1]);
AddHistogram(h.c_ht2[1]);
AddHistogram(h.c_ht2_over_ht[1]);
AddHistogram(h.c_m_leps[1]);
AddHistogram(h.c_m_leps_jets[1]);
AddHistogram(h.c_m_all[1]);
AddHistogram(h.c_m_os_1[1]);
AddHistogram(h.c_m_os_2[1]);
AddHistogram(h.c_m_ss[1]);
AddHistogram(h.c_m_ss_1[1]);
AddHistogram(h.c_m_ss_2[1]);
AddHistogram(h.c_m_ss_3[1]);
AddHistogram(h.c_mt_all[1]);
AddHistogram(h.c_mt_lep1_met[1]);
AddHistogram(h.c_mt_lep2_met[1]);
AddHistogram(h.c_mt_lep3_met[1]);
AddHistogram(h.c_mt_lep_met_max[1]);
AddHistogram(h.c_mt_lep_met_min[1]);
AddHistogram(h.c_sphericity[1]);
AddHistogram(h.c_aplanarity[1]);
AddHistogram(h.c_thrust[1]);
AddHistogram(h.c_sphericity_2d[1]);
AddHistogram(h.c_aplanarity_2d[1]);
AddHistogram(h.c_thrust_2d[1]);
AddHistogram(h.c_dphi_lep12[1]);
AddHistogram(h.c_dphi_lep13[1]);
AddHistogram(h.c_dphi_lep23[1]);
AddHistogram(h.c_dphi_lep1met[1]);
AddHistogram(h.c_dphi_lep2met[1]);
AddHistogram(h.c_dphi_lep3met[1]);
AddHistogram(h.c_dr_lep12[1]);
AddHistogram(h.c_dr_lep13[1]);
AddHistogram(h.c_dr_lep23[1]);
AddHistogram(h.c_dphi_j1met[1]);
AddHistogram(h.c_dphi_j2met[1]);
AddHistogram(h.c_dr_jt_lep1[1]);
AddHistogram(h.c_dr_jt_lep2[1]);
AddHistogram(h.c_dr_jt_lep3[1]);
AddHistogram(h.c_met_m_os_1[1]);
AddHistogram(h.c_trk_met[1]);
AddHistogram(h.c_st_met[1]);
AddHistogram(h.c_njet15[1]);
AddHistogram(h.c_trilep_pt[2][0]);
AddHistogram(h.c_trilep_pt[2][1]);
AddHistogram(h.c_trilep_pt[2][2]);
AddHistogram(h.c_trilep_eta[2][0]);
AddHistogram(h.c_trilep_eta[2][1]);
AddHistogram(h.c_trilep_eta[2][2]);
AddHistogram(h.c_trilep_phi[2][0]);
AddHistogram(h.c_trilep_phi[2][1]);
AddHistogram(h.c_trilep_phi[2][2]);
AddHistogram(h.c_met[2]);
AddHistogram(h.c_st[2]);
AddHistogram(h.c_st2[2]);
AddHistogram(h.c_st2_over_st[2]);
AddHistogram(h.c_ht[2]);
AddHistogram(h.c_ht2[2]);
AddHistogram(h.c_ht2_over_ht[2]);
AddHistogram(h.c_m_leps[2]);
AddHistogram(h.c_m_leps_jets[2]);
AddHistogram(h.c_m_all[2]);
AddHistogram(h.c_m_os_1[2]);
AddHistogram(h.c_m_os_2[2]);
AddHistogram(h.c_m_ss[2]);
AddHistogram(h.c_m_ss_1[2]);
AddHistogram(h.c_m_ss_2[2]);
AddHistogram(h.c_m_ss_3[2]);
AddHistogram(h.c_mt_all[2]);
AddHistogram(h.c_mt_lep1_met[2]);
AddHistogram(h.c_mt_lep2_met[2]);
AddHistogram(h.c_mt_lep3_met[2]);
AddHistogram(h.c_mt_lep_met_max[2]);
AddHistogram(h.c_mt_lep_met_min[2]);
AddHistogram(h.c_sphericity[2]);
AddHistogram(h.c_aplanarity[2]);
AddHistogram(h.c_thrust[2]);
AddHistogram(h.c_sphericity_2d[2]);
AddHistogram(h.c_aplanarity_2d[2]);
AddHistogram(h.c_thrust_2d[2]);
AddHistogram(h.c_dphi_lep12[2]);
AddHistogram(h.c_dphi_lep13[2]);
AddHistogram(h.c_dphi_lep23[2]);
AddHistogram(h.c_dphi_lep1met[2]);
AddHistogram(h.c_dphi_lep2met[2]);
AddHistogram(h.c_dphi_lep3met[2]);
AddHistogram(h.c_dr_lep12[2]);
AddHistogram(h.c_dr_lep13[2]);
AddHistogram(h.c_dr_lep23[2]);
AddHistogram(h.c_dphi_j1met[2]);
AddHistogram(h.c_dphi_j2met[2]);
AddHistogram(h.c_dr_jt_lep1[2]);
AddHistogram(h.c_dr_jt_lep2[2]);
AddHistogram(h.c_dr_jt_lep3[2]);
AddHistogram(h.c_met_m_os_1[2]);
AddHistogram(h.c_trk_met[2]);
AddHistogram(h.c_st_met[2]);
AddHistogram(h.c_njet15[2]);
AddHistogram(h.c_trilep_pt[3][0]);
AddHistogram(h.c_trilep_pt[3][1]);
AddHistogram(h.c_trilep_pt[3][2]);
AddHistogram(h.c_trilep_eta[3][0]);
AddHistogram(h.c_trilep_eta[3][1]);
AddHistogram(h.c_trilep_eta[3][2]);
AddHistogram(h.c_trilep_phi[3][0]);
AddHistogram(h.c_trilep_phi[3][1]);
AddHistogram(h.c_trilep_phi[3][2]);
AddHistogram(h.c_met[3]);
AddHistogram(h.c_st[3]);
AddHistogram(h.c_st2[3]);
AddHistogram(h.c_st2_over_st[3]);
AddHistogram(h.c_ht[3]);
AddHistogram(h.c_ht2[3]);
AddHistogram(h.c_ht2_over_ht[3]);
AddHistogram(h.c_m_leps[3]);
AddHistogram(h.c_m_leps_jets[3]);
AddHistogram(h.c_m_all[3]);
AddHistogram(h.c_m_os_1[3]);
AddHistogram(h.c_m_os_2[3]);
AddHistogram(h.c_m_ss[3]);
AddHistogram(h.c_m_ss_1[3]);
AddHistogram(h.c_m_ss_2[3]);
AddHistogram(h.c_m_ss_3[3]);
AddHistogram(h.c_mt_all[3]);
AddHistogram(h.c_mt_lep1_met[3]);
AddHistogram(h.c_mt_lep2_met[3]);
AddHistogram(h.c_mt_lep3_met[3]);
AddHistogram(h.c_mt_lep_met_max[3]);
AddHistogram(h.c_mt_lep_met_min[3]);
AddHistogram(h.c_sphericity[3]);
AddHistogram(h.c_aplanarity[3]);
AddHistogram(h.c_thrust[3]);
AddHistogram(h.c_sphericity_2d[3]);
AddHistogram(h.c_aplanarity_2d[3]);
AddHistogram(h.c_thrust_2d[3]);
AddHistogram(h.c_dphi_lep12[3]);
AddHistogram(h.c_dphi_lep13[3]);
AddHistogram(h.c_dphi_lep23[3]);
AddHistogram(h.c_dphi_lep1met[3]);
AddHistogram(h.c_dphi_lep2met[3]);
AddHistogram(h.c_dphi_lep3met[3]);
AddHistogram(h.c_dr_lep12[3]);
AddHistogram(h.c_dr_lep13[3]);
AddHistogram(h.c_dr_lep23[3]);
AddHistogram(h.c_dphi_j1met[3]);
AddHistogram(h.c_dphi_j2met[3]);
AddHistogram(h.c_dr_jt_lep1[3]);
AddHistogram(h.c_dr_jt_lep2[3]);
AddHistogram(h.c_dr_jt_lep3[3]);
AddHistogram(h.c_met_m_os_1[3]);
AddHistogram(h.c_trk_met[3]);
AddHistogram(h.c_st_met[3]);
AddHistogram(h.c_njet15[3]);
AddHistogram(h.c_trilep_pt[4][0]);
AddHistogram(h.c_trilep_pt[4][1]);
AddHistogram(h.c_trilep_pt[4][2]);
AddHistogram(h.c_trilep_eta[4][0]);
AddHistogram(h.c_trilep_eta[4][1]);
AddHistogram(h.c_trilep_eta[4][2]);
AddHistogram(h.c_trilep_phi[4][0]);
AddHistogram(h.c_trilep_phi[4][1]);
AddHistogram(h.c_trilep_phi[4][2]);
AddHistogram(h.c_met[4]);
AddHistogram(h.c_st[4]);
AddHistogram(h.c_st2[4]);
AddHistogram(h.c_st2_over_st[4]);
AddHistogram(h.c_ht[4]);
AddHistogram(h.c_ht2[4]);
AddHistogram(h.c_ht2_over_ht[4]);
AddHistogram(h.c_m_leps[4]);
AddHistogram(h.c_m_leps_jets[4]);
AddHistogram(h.c_m_all[4]);
AddHistogram(h.c_m_os_1[4]);
AddHistogram(h.c_m_os_2[4]);
AddHistogram(h.c_m_ss[4]);
AddHistogram(h.c_m_ss_1[4]);
AddHistogram(h.c_m_ss_2[4]);
AddHistogram(h.c_m_ss_3[4]);
AddHistogram(h.c_mt_all[4]);
AddHistogram(h.c_mt_lep1_met[4]);
AddHistogram(h.c_mt_lep2_met[4]);
AddHistogram(h.c_mt_lep3_met[4]);
AddHistogram(h.c_mt_lep_met_max[4]);
AddHistogram(h.c_mt_lep_met_min[4]);
AddHistogram(h.c_sphericity[4]);
AddHistogram(h.c_aplanarity[4]);
AddHistogram(h.c_thrust[4]);
AddHistogram(h.c_sphericity_2d[4]);
AddHistogram(h.c_aplanarity_2d[4]);
AddHistogram(h.c_thrust_2d[4]);
AddHistogram(h.c_dphi_lep12[4]);
AddHistogram(h.c_dphi_lep13[4]);
AddHistogram(h.c_dphi_lep23[4]);
AddHistogram(h.c_dphi_lep1met[4]);
AddHistogram(h.c_dphi_lep2met[4]);
AddHistogram(h.c_dphi_lep3met[4]);
AddHistogram(h.c_dr_lep12[4]);
AddHistogram(h.c_dr_lep13[4]);
AddHistogram(h.c_dr_lep23[4]);
AddHistogram(h.c_dphi_j1met[4]);
AddHistogram(h.c_dphi_j2met[4]);
AddHistogram(h.c_dr_jt_lep1[4]);
AddHistogram(h.c_dr_jt_lep2[4]);
AddHistogram(h.c_dr_jt_lep3[4]);
AddHistogram(h.c_met_m_os_1[4]);
AddHistogram(h.c_trk_met[4]);
AddHistogram(h.c_st_met[4]);
AddHistogram(h.c_njet15[4]);
AddHistogram(h.jet20_et);
AddHistogram(h.jet20_eta);
AddHistogram(h.jet20_phi);
AddHistogram(h.jet20_n);
AddHistogram(h.jet_scaleEta);
AddHistogram(h.beam_xy);
AddHistogram(h.vz);
AddHistogram(h.dilep_pt[0]);
AddHistogram(h.dilep_pt[1]);
AddHistogram(h.dilep_eta[0]);
AddHistogram(h.dilep_eta[1]);
AddHistogram(h.dilep_phi[0]);
AddHistogram(h.dilep_phi[1]);
AddHistogram(h.trilep_pt[0]);
AddHistogram(h.trilep_pt[1]);
AddHistogram(h.trilep_pt[2]);
AddHistogram(h.trilep_eta[0]);
AddHistogram(h.trilep_eta[1]);
AddHistogram(h.trilep_eta[2]);
AddHistogram(h.trilep_phi[0]);
AddHistogram(h.trilep_phi[1]);
AddHistogram(h.trilep_phi[2]);
AddHistogram(h.vz_lep1);
AddHistogram(h.vz_lep2);
AddHistogram(h.vz_lep3);
AddHistogram(h.st);
AddHistogram(h.st2);
AddHistogram(h.st2_over_st);
AddHistogram(h.ht);
AddHistogram(h.ht2);
AddHistogram(h.ht2_over_ht);
AddHistogram(h.met);
AddHistogram(h.met_phi);
AddHistogram(h.trk_met);
AddHistogram(h.st_met);
AddHistogram(h.met_trkpt);
AddHistogram(h.m_leps);
AddHistogram(h.m_leps_jets);
AddHistogram(h.m_all);
AddHistogram(h.m_os_1);
AddHistogram(h.m_os_2);
AddHistogram(h.m_ss);
AddHistogram(h.m_ss_1);
AddHistogram(h.m_ss_2);
AddHistogram(h.m_ss_3);
AddHistogram(h.mt_all);
AddHistogram(h.mt_lep1_met);
AddHistogram(h.mt_lep2_met);
AddHistogram(h.mt_lep3_met);
AddHistogram(h.mt_lep_met_max);
AddHistogram(h.mt_lep_met_min);
AddHistogram(h.met_m_os_1);
AddHistogram(h.sphericity);
AddHistogram(h.aplanarity);
AddHistogram(h.thrust);
AddHistogram(h.sphericity_2d);
AddHistogram(h.aplanarity_2d);
AddHistogram(h.thrust_2d);
AddHistogram(h.dphi_lep12);
AddHistogram(h.dphi_lep13);
AddHistogram(h.dphi_lep23);
AddHistogram(h.dphi_lep1met);
AddHistogram(h.dphi_lep2met);
AddHistogram(h.dphi_lep3met);
AddHistogram(h.dr_lep12);
AddHistogram(h.dr_lep13);
AddHistogram(h.dr_lep23);
AddHistogram(h.dphi_j1met_j15);
AddHistogram(h.dphi_j1met);
AddHistogram(h.dphi_j2met);
AddHistogram(h.dr_jt_lep1);
AddHistogram(h.dr_jt_lep2);
AddHistogram(h.dr_jt_lep3);
AddHistogram(h.n1_opt_dphi_lep12);
AddHistogram(h.n1_opt_met);
AddHistogram(h.n1_opt_m_os_1);
AddHistogram(h.n1_opt_m_os_2);
AddHistogram(h.n1_opt_ht);
AddHistogram(h.n1_opt_trk_met);
AddHistogram(h.n1_opt_mt_min);
AddHistogram(h.n1_opt_met_trkpt);
AddHistogram(h.n1_opt_ele1_et);
AddHistogram(h.n1_opt_ele2_et);
AddHistogram(h.n1_opt_trk_pt);
AddHistogram(h.opt_met);
AddHistogram(h.opt_st);
AddHistogram(h.opt_st2);
AddHistogram(h.opt_st2_over_st);
AddHistogram(h.opt_ht);
AddHistogram(h.opt_ht2);
AddHistogram(h.opt_ht2_over_ht);
AddHistogram(h.opt_m_leps);
AddHistogram(h.opt_m_leps_jets);
AddHistogram(h.opt_m_all);
AddHistogram(h.opt_m_os_1);
AddHistogram(h.opt_m_os_2);
AddHistogram(h.opt_m_ss);
AddHistogram(h.opt_m_ss_1);
AddHistogram(h.opt_m_ss_2);
AddHistogram(h.opt_m_ss_3);
AddHistogram(h.opt_mt_all);
AddHistogram(h.opt_mt_lep1_met);
AddHistogram(h.opt_mt_lep2_met);
AddHistogram(h.opt_mt_lep3_met);
AddHistogram(h.opt_mt_lep_met_max);
AddHistogram(h.opt_mt_lep_met_min);
AddHistogram(h.opt_sphericity);
AddHistogram(h.opt_aplanarity);
AddHistogram(h.opt_thrust);
AddHistogram(h.opt_sphericity_2d);
AddHistogram(h.opt_aplanarity_2d);
AddHistogram(h.opt_thrust_2d);
AddHistogram(h.opt_met_m_os_1);
AddHistogram(h.opt_trk_met);
AddHistogram(h.opt_st_met);
AddHistogram(h.opt_njet15);

AddHistogram(h.fopt_met);
AddHistogram(h.fopt_mos1);
AddHistogram(h.fopt_mos2);
AddHistogram(h.fopt_ht);
AddHistogram(h.fopt_dphi);
AddHistogram(h.fopt_mt);
AddHistogram(h.fopt_et[0]);
AddHistogram(h.fopt_et[1]);
AddHistogram(h.fopt_et[2]);
AddHistogram(h.fopt_n1_met);
AddHistogram(h.fopt_n1_mos1);
AddHistogram(h.fopt_n1_mos2);
AddHistogram(h.fopt_n1_ht);
AddHistogram(h.fopt_n1_dphi);
AddHistogram(h.fopt_n1_mt);
AddHistogram(h.fopt_n1_et[0]);
AddHistogram(h.fopt_n1_et[1]);
AddHistogram(h.fopt_n1_et[2]);

 
}
// This gets the cos_theta_star of a in the rest frame of a and b.
float TSeel::cos_theta_star(const TLorentzVector &a, const TLorentzVector &b)
{
    TLorentzVector c = -1.* (a + b);
    c.SetE(-c.E());
    TLorentzVector d = a;
    d.Boost(c.BoostVector());
    
    return (float)d.CosTheta();
}


bool TSeel::pass_track_cut(int itrk, int category)
{
    bool pass = false;
    if ( category == 1 ) { // Electron
      TStnElectron *el = fElectronBlock->Electron(itrk);
      int it = el->TrackNumber();
      if(it>=0){
	TStnTrack *trk = fTrackBlock->Track(it);
	if ( trk->NCotAxSeg(5)>=3 && trk->NCotStSeg(5)>=3 )
	  pass = true;
	//             && fElectrons_fCOTChi2[itrk]/(fElectrons_fnCOTHits[itrk]-5)<3.
      }
    }
    else if ( category == 2 ) { // Muon
      TStnMuon *mu = fMuonBlock->Muon(itrk);
      int it = mu->TrackNumber();
      if(it>=0){
	TStnTrack *trk = fTrackBlock->Track(it);
	if ( trk->NCotAxSeg(5)>=3 && trk->NCotStSeg(5)>=3 && TrkIso(trk,0.4,4,0.4)<0.1 )
	  pass = true;
	//             && fMuons_fCOTChi2[itrk]/(fMuons_fnCOTHits[itrk]-5)<3.
      }
    }

    return pass;

}

float TSeel::TrkIso(TStnTrack *track, float cone, float dz, float ptMin) 
{
  // calculate `track' isolation: (sumPt of other tracks)/trackPt
  // from `trackList' within `cone' in eta-phi from the `track'. 
  // Tracks with Pt > `ptMin' are required to have |Ztrk-Z0|<`dz'
  // where Z0 is z-coordinate of the `track' in the point of closest
  // approach to the origin, thus cutting against multiple interactions

  // Modified on 9/30/2004
  // using all tracks around THE TRACK in a cone of 0.4 and in |dz|<5 cm.

    float  pt0, phi0, eta0, z0;
    float  pt, phi, eta, z, dr, trk_iso;
    int    ntrk;
    float sum_p,sum_pt;
    sum_p=sum_pt=0;

//    TLorentzVector *mom = new TLorentzVector();
//    int status = track->GetBcMomentum(mom);
    pt0 = track->Pt();
    phi0 = TVector2::Phi_0_2pi(track->Phi0());
    eta0 = track->Eta();
    z0   = track->Z0();
    ntrk  = fTrackBlock->NTracks();

    float px=0.; float py=0.; float pz=0.;
    for ( int i=0; i<ntrk; i++ ) {
        TStnTrack* ti = (TStnTrack*) fTrackBlock->Track(i);
        if ( ti!=track && pass_track_cut(ti) ) {
//            TLorentzVector* momi = new TLorentzVector();
//            int status = ti->GetBcMomentum(momi);
            pt  = ti->Pt();
            phi = TVector2::Phi_0_2pi(ti->Phi0());
            eta = ti->Eta();
            z   = ti->Z0();
            if ( pt>ptMin && fabs(z-z0)<dz ) {
                dr   = delta_r(eta,phi,eta0,phi0);
                if ( dr < cone ) {
                    sum_pt += pt;
                    px+=pt*cos(phi); py+=pt*sin(phi); pz+=pt0*sinh(eta);
                }
            }
//            delete momi;
        }
    }
    sum_p = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
  // finally calculate isolation
    trk_iso = sum_pt/track->Pt();
//    delete mom;
    return trk_iso;
}

bool TSeel::pass_track_cut(TStnTrack *trk)
{
//    TLorentzVector *mom = new TLorentzVector();
//    int status = trk->GetBcMomentum(mom);
//    delete mom;
    if ( trk->NCotAxSeg(5)>=2
         && trk->NCotStSeg(5)>=2
         && fabs(trk->Z0())<60. )
        return true;
    else
        return false;
}

bool TSeel::has_radiated_photon()
{
  if ( !_DY ) return true; // if Data or Signal return true

  // For each of the two main leptons in lep array
  // check if it has radiated a photon at the hepg
  // level. If so, return true.
  // Ensure before this that lep array has two good ele

  bool has_radiated_photon = false;

  TStnElectron *ele1 = fElectronBlock->Electron(lep[0].index());
  TStnTrack *trk1 = fTrackBlock->Track(ele1->TrackNumber());
  TStnElectron *ele2 = fElectronBlock->Electron(lep[1].index());
  TStnTrack *trk2 = fTrackBlock->Track(ele2->TrackNumber());

  int hepg_ele[100],nhepg_ele=0;
  for (int i=0; i<fGenpBlock->NParticles(); i++) {
      TGenParticle* p = fGenpBlock->Particle(i);  
      if ( p->GetStatusCode()==1 && p->IsPhoton() ) {
	 int im = p->GetFirstMother();
	 TGenParticle *pm = fGenpBlock->Particle(im);
	 if ( pm->IsElectron() ) {
	    hepg_ele[nhepg_ele] = im;
	    nhepg_ele++;
	 }
      }
  }
  //h.nhepg_ele->Fill(nhepg_ele);
  
  // Now hepg_ele[] has all the hepg electrons which radiated a photon
  // we match these to lep[0] and lep[1] and return true is one is matched

  for( int i=0; i<nhepg_ele; i++){

      TGenParticle* p = fGenpBlock->Particle(hepg_ele[i]);  
      float phi1 = TVector2::Phi_0_2pi(p->Phi());
      // check is this hepg ele matched lep[0]
      float phi2 = TVector2::Phi_0_2pi(trk1->Phi0());       // not beam constrained
      double dr = delta_r(p->Eta(),phi1,trk1->Eta(),phi2);  // not beam constrained    
      if( dr<0.05 && 
	  fabs(p->P()-ele1->EmE()*(ele1->Etcor()>0. ? ele1->Etcor(): 1.))/p->P()<0.3 ){
	 has_radiated_photon = true;
	 break;
      }
      // else check is this hepg ele matched lep[1]
      phi2 = TVector2::Phi_0_2pi(trk2->Phi0());       // not beam constrained
      dr = delta_r(p->Eta(),phi1,trk2->Eta(),phi2);  // not beam constrained    
      if( dr<0.05 && 
	  fabs(p->P()-ele2->EmE()*(ele2->Etcor()>0. ? ele2->Etcor(): 1.))/p->P()<0.3 ){
	 has_radiated_photon = true;
	 break;
      }

  }

  return has_radiated_photon;

}
/*
    nleptons_temp = nele + nmuon;
    if(fVerbose<=-2 && nleptons_temp>1){
      cout<<"nele_temp = "<<nele<<" , nmu_temp = "<<nmuon<<" , nleptons_temp = "<<nleptons_temp<<endl;
    }
    if(nele>1)
      nev_2ele++;
    if(nmuon>1)
      nev_2mu++;
    if(nele>1 && nleptons_temp>2){
      nev_3obj++;
      for(int it=2; it<nleptons_temp; it++){
	if( abs(lep_temp[it].id())==11 ){
	  if( pass_track_cut(lep_temp[it].index(),1) ){
	    nev_3e++;
	    //cout<<"eee evt run/event "<<fHeaderBlock->RunNumber()<<" "<<fHeaderBlock->EventNumber()<<endl;
	    break;
	  }
	}
      }
      for(int it=2; it<nleptons_temp; it++){
	if( abs(lep_temp[it].id())==13 ){
	  if( pass_track_cut(lep_temp[it].index(),2) ){
	    nev_3m++;
	    //cout<<"eem evt run/event "<<fHeaderBlock->RunNumber()<<" "<<fHeaderBlock->EventNumber()<<endl;
	    break;
	  }
	}
      }

      if(nele>2)
      if( pass_track_cut(lep_temp[2].index(),1) ){
      nev_3e++;
      cout<<"eee evt run/event "<<fHeaderBlock->RunNumber()<<" "<<fHeaderBlock->EventNumber()<<endl;
      }
      if(nmuon>0){
      if( pass_track_cut(lep_temp[nele].index(),2) ){
      nev_3m++;
      cout<<"eem evt run/event "<<fHeaderBlock->RunNumber()<<" "<<fHeaderBlock->EventNumber()<<endl;
      }
      }
      if(nmuon>1){
      if( !pass_track_cut(lep_temp[nele].index(),2) && pass_track_cut(lep_temp[nele+1].index(),2) ){
      nev_3m++;
      cout<<"eem evt run/event "<<fHeaderBlock->RunNumber()<<" "<<fHeaderBlock->EventNumber()<<endl;
      }
      }

    }
*/
