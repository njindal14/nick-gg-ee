#ifndef ANALYZE_H
#define ANALYZE_H

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TH1F.h"
#include <string>

class Analyze 
{
 public:
  Analyze(); //Constructor
  Analyze(TString infile, Int_t nEvents); //Special constructor
  ~Analyze(); //Destructor
  Int_t Init();
  Int_t NextEvent();
  TParticle* NextParticle();
  void doAnalysis(); //Function doing the analysis
  
 private:
 
  TH1F* mPt;
  TH1F* mPtEl;
  TH1F* mMass;
  TH1F* mPhi;
  TH1F* mEta;
  TH1F* mPionPt;

  TH3F* mSLCos1phivsPTvsMass;
  TH3F* mSLCos2phivsPTvsMass;
  TH3F* mSLCos3phivsPTvsMass;
  TH3F* mSLCos4phivsPTvsMass;

  TH2F* mcos2phivspt;
  TH2F* mcos4phivspt;


  FILE *filelist;
  TString fInfile;
  Int_t fNParticles;
  Int_t fNEvents;
  
};

#endif
