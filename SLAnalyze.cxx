#include "SLAnalyze.h"
#include <TMath.h>
#include <iostream>
#include <vector>


// This macro reads in a starlight output file and creates histograms of 
// the p_T and rapidity of the daughters, as well as the p_T, rapidity and 
// mass of the parent.  The macro assumes there are only two daughter 
// tracks.  Currently, the macro only accomodates electrons, muons or 
// pions as daughter particles.  The histograms for the daughter particles 
// are called fPt2, fPt2, fRap1, and fRap2.  Parent histograms are created 
// for each possible daughter species (e.g., parent p_T histograms are 
// created with the names fPtEl, fPtMu, and fPtPi), but only the ones 
// corresponding to  the actual daughter particle are filled. The 
// histograms are saved in a  file called histograms.root.
//
// To use this macro, modify the file ana.C to call your input file 
// (as downloaded, it calls slight.out) and the number of events you wish to 
// process (as downloaded, it processes 20 events).  Then open root and type 
// ".x Ana.C" .   

using namespace std;
Analyze::Analyze() :
  fInfile("/Users/Nick/STAR/starlight/utils/slight.out"),

  //running for gold
  //fInfile("/Users/Nick/STAR/docker_mount/JobFilesAu1M/slight.out"),
  fNEvents(1)
{
  //Constructor
  
  //Creating histograms
  mPt = new TH1F("mPt", "Pair Transverse Momentum; Pair P_{T} (GeV/c); Counts", 50, 0, 0.4);
  mPtEl = new TH1F("mPtEl", "Electron pT; pT (GeV/c); Counts", 50, 0, .5);
  mMass = new TH1F("mMass", "Pair Mass; M_{ee} (GeV/c^2); Counts", 50, 0, 3);
  mPhi = new TH1F("mPhi", "Pair #Phi; #Phi (rad); counts", 100, -3.141592, 3.141592);
  mEta = new TH1F("mEta", "Pair Eta; #eta; counts", 100, -1, 1);

  mSLCos1phivsPTvsMass = new TH3F("mSLCos1phivsPTvsMass", "Pair P_{T} (GeV/c); pair Mass (GeV); 2<cos#phi>; counts", 50, -2, 2, 10, 0, 0.4, 20, 0, 1.5);
  mSLCos2phivsPTvsMass = new TH3F("mSLCos2phivsPTvsMass", "Pair P_{T} (GeV/c); pair Mass (GeV); 2<cos2#phi>; counts", 50, -2, 2, 10, 0, 0.4, 20, 0, 1.5);
  mSLCos3phivsPTvsMass = new TH3F("mSLCos3phivsPTvsMass", "Pair P_{T} (GeV/c); pair Mass (GeV); 2<cos3#phi>; counts", 50, -2, 2, 10, 0, 0.4, 20, 0, 1.5);
  mSLCos4phivsPTvsMass = new TH3F("mSLCos4phivsPTvsMass", "Pair P_{T} (GeV/c); pair Mass (GeV); 2<cos4#phi>; counts", 50, -2, 2, 10, 0, 0.4, 20, 0, 1.5);


  /*const Int_t ptbins = 10;
  Double_t edges[ptbins + 1] = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.14, .2};*/

  mcos2phivspt = new TH2F("mCos2phivsPT", "pair pT; 2<cos2#phi>", 50, -2, 2, 7, 0, 0.1);
  mcos4phivspt = new TH2F("mCos4phivsPT", "pair pT; 2<cos4#phi>", 50, -2, 2, 7, 0, 0.1);

  SLphivPt = new TH2F("SLphivPt", "SLphivPt", 20, -3.14, 3.14, 7, 0, 0.1);
}

Analyze::Analyze(TString infile, Int_t nEvents) :
  fInfile(infile),
  fNEvents(nEvents)
{
  //Special constructor

  //Creating histograms
  mPt = new TH1F("mPt", "Pair Transverse Momentum; Pair P_{T} (GeV/c); counts", 100, 0, 0.2);
  mPtEl = new TH1F("mPtEl", "Electron pT; pT (GeV/c); Counts", 50, 0, .5);
  mMass = new TH1F("mMass", "Pair Mass; M_{ee} (GeV); counts", 100, 0, .5);
  mPhi = new TH1F("mPhi", "Pair #Phi; #Phi (rad); counts", 100, -3.141592, 3.141592);
  mEta = new TH1F("mEta", "Pair Eta; #eta; counts", 100, -1, 1);
  //mPionPt = new TH1F("mPionPt", "#pi Transverse Momentum; #pi P_{T} (GeV/c); counts", 50, 0, 0.75);

  mSLCos1phivsPTvsMass = new TH3F("mSLCos1phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos#phi>; counts", 50, -2, 2, 10, 0, 0.4, 20, 0, 1.5);
  mSLCos2phivsPTvsMass = new TH3F("mSLCos2phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos2#phi>; counts", 50, -2, 2, 10, 0, 0.4, 20, 0, 1.5);
  mSLCos3phivsPTvsMass = new TH3F("mSLCos3phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos3#phi>; counts", 50, -2, 2, 10, 0, 0.4, 20, 0, 1.5);
  mSLCos4phivsPTvsMass = new TH3F("mSLCos4phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos4#phi>; counts", 50, -2, 2, 10, 0, 0.4, 20, 0, 1.5);

  const Int_t ptbins = 10;
  Double_t edges[ptbins + 1] = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.14, .2};

  mcos2phivspt = new TH2F("mCos2phivsPT", "pair pT; 2<cos2#phi>", 50, -2, 2, 10, edges);
  mcos4phivspt = new TH2F("mCos4phivsPT", "pair pT; 2<cos4#phi>", 50, -2, 2, 10, edges);
  SLphivPt = new TH2F("SLphivPt", "SLphivPt", 20, -3.14, 3.14, 10, edges);

}

Analyze::~Analyze()
{
  //Destructor
  delete mPt;
  delete mPtEl;
  delete mMass;
  delete mPhi;
  delete mEta;
  delete mPionPt;

  delete mSLCos1phivsPTvsMass;
  delete mSLCos2phivsPTvsMass;
  delete mSLCos3phivsPTvsMass;
  delete mSLCos4phivsPTvsMass;

  delete mcos2phivspt;
  delete mcos4phivspt;
  delete SLphivPt;
}

Int_t Analyze::Init()
{
  
  
  cout << "Opening textfile " << fInfile << endl;
  if( !(filelist=fopen(fInfile.Data(),"r")) ){
    cout<<"Couldn't open input file: "<<fInfile<<endl;
    return -1;
  }
  cout << "Done opening textfile" << endl;
  return 0;
}

Int_t Analyze::NextEvent()
{
  char linelabel[20];
  int i1=0;
  int i2=0;
  int i3=0;
  double x1=0.0;
  double x2=0.0;
  double x3=0.0;
  double x4=0.0;
  int ntrk=0;
  int nvtx=0;
  // Event line 
  fscanf(filelist,"%s %d %d %d ",linelabel,&i1,&ntrk,&i2);
  fNParticles = ntrk;
  cout<<linelabel<<"  "<<i1<<"  "<<ntrk<<"  "<<i2<<"   "<<fNParticles<<endl;
  // Vertex line 
  fscanf(filelist,"%s %lf %lf %lf %lf %d %d %d %d",
         linelabel,&x1,&x2,&x3,&x4,&i1,&i2,&i3,&nvtx);
  cout<<linelabel<<"  "<<x1<<"  "<<x2<<"  "<<x3<<"  "<<x4<<"  "<<i1<<"  "<<i2<<"  "<<i3<<"  "<<nvtx<<endl;
  if(ntrk != nvtx)cout<<"ERROR: ntrk = "<<ntrk<<"  nvtx = "<<nvtx<<endl;
  //
  return fNParticles;
}

TParticle* Analyze::NextParticle()
{
  char tracklabel[20];
  int i1=0;
  int i2=0;
  int i3=0;
  int i4=0;
  Int_t idpart = 0;
  Double_t px = 0.0;
  Double_t py = 0.0;
  Double_t pz = 0.0;
  Double_t ep = 0.0;
  
  cout<<"In NextParticle: fNparticles = "<<fNParticles<<endl;

  // for ( int itk=0; itk < fNParticles; itk++){


    fscanf(filelist,"%s %d %le %le %le %d %d %d %d",
         tracklabel,&i1,&px,&py,&pz,&i2,&i3,&i4,&idpart);
    cout<<"   "<<tracklabel<<"  "<<i1<<"  "<<px<<"  "<<py<<"  "<<pz<<"  "<<i2<<"  "<<i3<<"  "<<i4<<"  "<<idpart<<endl;
  
    TParticle *particle = 
      new TParticle(idpart, 0, -1, -1, -1, -1, px, py, pz, ep, 0., 0., 0., 0.);
      if(idpart == 11 || idpart == -11){particle->SetCalcMass(0.00051099907);}
      else if(idpart == 13 || idpart == -13){particle->SetCalcMass(0.105658389);}
      else if(idpart == 211 || idpart == -211){particle->SetCalcMass(0.13956995);}
      else {cout << "unknown daughter!  please edit the code to accomodate"<< endl;
	      exit(0);
	   }

  return particle;
}

double calc_Phi( TLorentzVector lv1, TLorentzVector lv2) {
    TLorentzVector lvPlus = lv1 + lv2;
    lv1.Boost(-lvPlus.BoostVector());
    lv2.Boost(-lvPlus.BoostVector());
    TLorentzVector lvMinus = lv1 - lv2;
    double Px = lvPlus.Px();
    double Py = lvPlus.Py();
    double Qx = lvMinus.Px();
    double Qy = lvMinus.Py();
    double PcrossQ = (Px*Qy) - (Py*Qx);
    double cosphi = (Px*Qx + Py*Qy) / (lvPlus.Pt()*lvMinus.Pt());
    double PairPhi = acos(cosphi);
    if ( PcrossQ > 0 ){
        PairPhi = PairPhi - 3.141592;
    } else {
        PairPhi = 3.141592 - PairPhi;
    }
    return PairPhi;
}

void Analyze::doAnalysis()
{

  Int_t check = Init();
  if(check < 0) return;
  Double_t mass;
  //Doing the analysis:  loop over events
  for(Int_t ev = 0; ev < fNEvents; ev++){
    	   
    const Int_t ntracks = NextEvent();
    //Array of TLorentzVectors. One vector for each track
    vector<int> vecidpart; 
    TLorentzVector* vecArr[ntracks];
      Int_t idpart = 0;
    //Looping over the tracks of the event
    for(Int_t tr = 0; tr < ntracks; tr++){
      //Getting the TParticle for the track
      TParticle *part = NextParticle();
      mass = part->GetCalcMass();
      idpart = part->GetPdgCode();
      vecidpart.push_back(idpart);
      Double_t energy = TMath::Sqrt(mass*mass+part->Px()*part->Px()+part->Py()*part->Py()+part->Pz()*part->Pz());
      //Creating a new TLorentzVector and setting px, py, pz and E.
      vecArr[tr] = new TLorentzVector;
      vecArr[tr]->SetPxPyPzE(part->Px(), part->Py(), part->Pz(), energy); 
      cout << "particle " << tr << ": px: " << part->Px() << " py: " << part->Py() << " pz: " << part->Pz() << " Energy: " << energy << endl;
    }
    // spin analysis
    TLorentzVector lvp, lvm, lv;
    double M_e = 0.00051;

    if (vecidpart[0] == 11 || vecidpart[1] == -11){
      lvp.SetPtEtaPhiM(vecArr[0]->Pt(), vecArr[0]->Eta(), vecArr[0]->Phi(), M_e);
      lvm.SetPtEtaPhiM(vecArr[1]->Pt(), vecArr[1]->Eta(), vecArr[1]->Phi(), M_e);
    } else if (vecidpart[0] == -11 || vecidpart[1] == 11){
      lvp.SetPtEtaPhiM(vecArr[1]->Pt(), vecArr[1]->Eta(), vecArr[1]->Phi(), M_e);
      lvm.SetPtEtaPhiM(vecArr[0]->Pt(), vecArr[0]->Eta(), vecArr[0]->Phi(), M_e);
    }

    //if ( abs(lvp.Eta()) < 10 ){ mPionPt->Fill(lvp.Pt()); }
    //if ( abs(lvm.Eta()) < 10 ){ mPionPt->Fill(lvm.Pt()); }

    if ( abs(lvp.Eta()) < 1 && abs(lvm.Eta()) <1) {
      lv = lvp+lvm;

      double PhiVal = calc_Phi(lvp, lvm);
      double cos1phi = 2*cos(PhiVal);
      double cos2phi = 2*cos(2*PhiVal);
      double cos3phi = 2*cos(3*PhiVal);
      double cos4phi = 2*cos(4*PhiVal);

      if(lv.Pt() !=0) mPt->Fill(lv.Pt());
      mPtEl->Fill(lvm.Pt());
      //if(lv.M() != 0) 
      mMass->Fill(lv.M());
      mPhi->Fill(PhiVal);
      if(lv.Eta() != 0) mEta->Fill(lv.Eta());

      mSLCos1phivsPTvsMass->Fill(cos1phi, lv.Pt(), lv.M());
      mSLCos2phivsPTvsMass->Fill(cos2phi, lv.Pt(), lv.M());
      mSLCos3phivsPTvsMass->Fill(cos3phi, lv.Pt(), lv.M());
      mSLCos4phivsPTvsMass->Fill(cos4phi, lv.Pt(), lv.M());


      if(lv.M() > 0.45 && lv.M() < 0.76){
        mcos2phivspt->Fill(cos2phi, lv.Pt());
        mcos4phivspt->Fill(cos4phi, lv.Pt());
        SLphivPt->Fill(PhiVal, lv.Pt());
      }

    }
  }

  //Writing the histograms to file
  TFile file("./output_root_files/SL_plotsBetter.root", "RECREATE");
  mMass->Write();
  mPt->Write();
  mPtEl->Write();
  mPhi->Write();
  mEta->Write();
  //mPionPt->Write();

  mSLCos1phivsPTvsMass->Write();
  mSLCos2phivsPTvsMass->Write();
  mSLCos3phivsPTvsMass->Write();
  mSLCos4phivsPTvsMass->Write();

  mcos2phivspt->Write();
  mcos4phivspt->Write();
  SLphivPt->Write();
  std::cout << "wrote";
}
