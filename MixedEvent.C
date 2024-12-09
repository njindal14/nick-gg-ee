#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "FemtoPairFormat.h"
#include "TRandom3.h"
#include <vector>
//#include "../Include/FemtoPairFormat.h"
#include "AnalysisFunctions.h"
#include "PlottingFunctions.h"

std::random_device global_rng;
TRandom3 rng(global_rng());

void MixedEvent() {

    TH1F("h1", "ntuple", 100, -4, 4);
    TFile * fo = new TFile( "output_root_files/MixedEventplots_Au.root", "RECREATE" );

    auto * mAllTrackPt = new TH1F("mAllTrackPt", "e Track P_{T}; P_{T} (GeV/c); counts", 50, 0, 0.8);
    auto * mLSTrackPt = new TH1F("mLSTrackPt", "LS e Track P_{T}; P_{T} (GeV/c); counts", 50, 0, 0.8);
    auto * mULSTrackPt = new TH1F("mULSTrackPt", "ULS e Track P_{T}; P_{T} (GeV/c); counts", 50, 0, 0.8);

    auto * mULSLSphipt = new TH2F("mULSLSphipt", "mULSLSphipt; #phi; pT (GeV/c)", 30, -3.14, 3.14, 30, 0, .3);

    auto * mLSLSPtvsM = new TH2F("mLSLSPtvsM", "LSLS pair P_{T}; P_{T} (GeV/c); counts", 300, 0, 5, 300, 0, 2);
    auto * mLSULSPtvsM = new TH2F("mLSULSPtvsM", "LSULS pair P_{T}; P_{T} (GeV/c); counts", 300, 0, 5, 300, 0, 2);
    auto * mULSLSPtvsM = new TH2F("mULSLSPtvsM", "ULSLS pair P_{T}; P_{T} (GeV/c); counts", 300, 0, 5, 300, 0, 2);
    auto * mULSULSPtvsM = new TH2F("mULSULSPtvsM", "ULSULS pair P_{T}; P_{T} (GeV/c); counts", 300, 0, 5, 300, 0, 2);

    auto * mPosPhi = new TH1F("mPosPhi", "#e^{+} #phi; #phi (rad); counts", 50, -3.141592, 3.141592);
    auto * mNegPhi = new TH1F("mNegPhi", "#e^{-} #phi; #phi (rad); counts", 50, -3.141592, 3.141592);
    auto * acoULSULS = new TH1F("acoULSULS", "acoULSULS", 20, 0, .2);

    auto * mLSLSCos1phivsPTvsMass = new TH3F("mLSLSCos1phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos1#phi>; counts", 50, -2, 2, 30, 0, 0.6, 30, 0, 2);
    auto * mLSLSCos2phivsPTvsMass = new TH3F("mLSLSCos2phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos2#phi>; counts", 50, -2, 2, 30, 0, 0.6, 30, 0, 2);
    auto * mLSLSCos3phivsPTvsMass = new TH3F("mLSLSCos3phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos3#phi>; counts", 50, -2, 2, 30, 0, 0.6, 30, 0, 2);
    auto * mLSLSCos4phivsPTvsMass = new TH3F("mLSLSCos4phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos4#phi>; counts", 50, -2, 2, 30, 0, 0.6, 30, 0, 2);

    auto * mLSULSCos1phivsPTvsMass = new TH3F("mLSULSCos1phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos1#phi>; counts", 50, -2, 2, 30, 0, 0.6, 30, 0, 2);
    auto * mLSULSCos2phivsPTvsMass = new TH3F("mLSULSCos2phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos2#phi>; counts", 50, -2, 2, 30, 0, 0.6, 30, 0, 2);
    auto * mLSULSCos3phivsPTvsMass = new TH3F("mLSULSCos3phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos3#phi>; counts", 50, -2, 2, 30, 0, 0.6, 30, 0, 2);
    auto * mLSULSCos4phivsPTvsMass = new TH3F("mLSULSCos4phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos4#phi>; counts", 50, -2, 2, 30, 0, 0.6, 30, 0, 2);

    auto * mULSLSCos1phivsPTvsMass = new TH3F("mULSLSCos1phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos1#phi>; counts", 50, -2, 2, 30, 0, 0.6, 30, 0, 2);
    auto * mULSLSCos2phivsPTvsMass = new TH3F("mULSLSCos2phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos2#phi>; counts", 50, -2, 2, 30, 0, 0.6, 30, 0, 2);
    auto * mULSLSCos3phivsPTvsMass = new TH3F("mULSLSCos3phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos3#phi>; counts", 50, -2, 2, 30, 0, 0.6, 30, 0, 2);
    auto * mULSLSCos4phivsPTvsMass = new TH3F("mULSLSCos4phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos4#phi>; counts", 50, -2, 2, 30, 0, 0.6, 30, 0, 2);

    auto * mULSULSCos1phivsPTvsMass = new TH3F("mULSULSCos1phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos1#phi>; counts", 50, -2, 2, 30, 0, 0.6, 30, 0, 2);
    auto * mULSULSCos2phivsPTvsMass = new TH3F("mULSULSCos2phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos2#phi>; counts", 50, -2, 2, 30, 0, 0.6, 30, 0, 2);
    auto * mULSULSCos3phivsPTvsMass = new TH3F("mULSULSCos3phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos3#phi>; counts", 50, -2, 2, 30, 0, 0.6, 30, 0, 2);
    auto * mULSULSCos4phivsPTvsMass = new TH3F("mULSULSCos4phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos4#phi>; counts", 50, -2, 2, 30, 0, 0.6, 30, 0, 2);


    TChain * ch = new TChain("PairDst");
    ch->Add("/Users/Nick/STAR/breit-wheeler/rootFiles/slim_pair_dst_Run10AuAu.root");
    ch->Add("/Users/Nick/STAR/breit-wheeler/rootFiles/slim_pair_dst_Run11AuAu.root");
    //ch->Add("/Users/Nick/STAR/breit-wheeler/rootFiles/pair_dst_Run12UU.root");

    TTreeReader myReader(ch);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

    vector<TLorentzVector> LSPosBuffer;
    vector<TLorentzVector> LSNegBuffer;
    vector<TLorentzVector> ULSPosBuffer;
    vector<TLorentzVector> ULSNegBuffer;

    while (myReader.Next()) {
        double chipipi = pow( pair->d1_mNSigmaPion, 2) + pow( pair->d2_mNSigmaPion, 2);
        double chiee = pow(pair->d1_mNSigmaElectron,2) + pow(pair->d2_mNSigmaElectron, 2);
        double dca1 = pair->d1_mDCA;
        double dca2 = pair->d2_mDCA;

        int LS_buffer_size = 100;
        int ULS_buffer_size = 100;

        TLorentzVector lv1, lv2;
        lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.00051 );
        lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.00051 );

        double ddTOFval = ddTOF( 0.00051, lv1.P(), pair->d1_mTof, pair->d1_mLength, lv2.P(), pair->d2_mTof, pair->d2_mLength );
        
        if ( chiee < 10 && 3*chiee < chipipi && dca1 <1 && dca2 <1 && fabs(pair->mVertexZ) < 100 &&  pair->mGRefMult <= 4 && lv1.Pt() > 0.2 && lv2.Pt() > 0.2){

//POSITIVE LIKE SIGN SOURCE==========================================================================================================
            if ( pair->mChargeSum==2 ) {
            mAllTrackPt->Fill(lv1.Pt());
            mAllTrackPt->Fill(lv2.Pt());
            mLSTrackPt->Fill(lv1.Pt());
            mLSTrackPt->Fill(lv2.Pt());

            mPosPhi->Fill(lv1.Phi());
            mPosPhi->Fill(lv2.Phi());
//MAKE LS PAIRS
                for ( TLorentzVector mixed_lv1 : LSPosBuffer ) {
                    TLorentzVector mixed_lv = mixed_lv1 + lv2 ;
                    double mixedPairPhi = calc_DeltaPhi(mixed_lv1, lv2);
                    if ( abs(mixed_lv.Rapidity()) < 1 && pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag !=0 && abs(ddTOFval) < 0.4 ) {
                        mLSLSCos1phivsPTvsMass->Fill( 2*cos(1*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSLSCos2phivsPTvsMass->Fill( 2*cos(2*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSLSCos3phivsPTvsMass->Fill( 2*cos(3*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSLSCos4phivsPTvsMass->Fill( 2*cos(4*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                    }
                }
                for ( TLorentzVector mixed_lv2 : LSPosBuffer ) {
                    TLorentzVector mixed_lv = lv1 + mixed_lv2 ;
                    double mixedPairPhi = calc_DeltaPhi(lv1, mixed_lv2);
                    if ( abs(mixed_lv.Rapidity()) < 1 && pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag !=0 && abs(ddTOFval) < 0.4 ) {
                        mLSLSCos1phivsPTvsMass->Fill( 2*cos(1*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSLSCos2phivsPTvsMass->Fill( 2*cos(2*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSLSCos3phivsPTvsMass->Fill( 2*cos(3*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSLSCos4phivsPTvsMass->Fill( 2*cos(4*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                    }
                }
//MAKE ULS PAIRS
                for ( TLorentzVector mixed_lv1 : LSNegBuffer ) {
                    TLorentzVector mixed_lv = mixed_lv1 + lv2 ;
                    double mixedPairPhi = calc_DeltaPhi(mixed_lv1, lv2);
                    if ( abs(mixed_lv.Rapidity()) < 1 && pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag !=0 && abs(ddTOFval) < 0.4 ) {
                        mLSULSCos1phivsPTvsMass->Fill( 2*cos(1*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSULSCos2phivsPTvsMass->Fill( 2*cos(2*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSULSCos3phivsPTvsMass->Fill( 2*cos(3*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSULSCos4phivsPTvsMass->Fill( 2*cos(4*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                    }
                }
                for ( TLorentzVector mixed_lv2 : LSNegBuffer ) {
                    TLorentzVector mixed_lv = lv1 + mixed_lv2 ;
                    double mixedPairPhi = calc_DeltaPhi(lv1, mixed_lv2);
                    if ( abs(mixed_lv.Rapidity()) < 1 && pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag !=0 && abs(ddTOFval) < 0.4 ) {
                        mLSULSCos1phivsPTvsMass->Fill( 2*cos(1*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSULSCos2phivsPTvsMass->Fill( 2*cos(2*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSULSCos3phivsPTvsMass->Fill( 2*cos(3*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSULSCos4phivsPTvsMass->Fill( 2*cos(4*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                    }
                }
                int buff_index = rng.Integer(LS_buffer_size - 1 );
                LSPosBuffer = update_buffer( LSPosBuffer, lv1, LS_buffer_size, buff_index );
                LSPosBuffer = update_buffer( LSPosBuffer, lv2, LS_buffer_size, buff_index );
            }

//NEGATIVE LIKE SIGN SOURCE ==========================================================================================================
            if ( pair->mChargeSum==-2 ) {
            mAllTrackPt->Fill(lv1.Pt());
            mAllTrackPt->Fill(lv2.Pt());
            mLSTrackPt->Fill(lv1.Pt());
            mLSTrackPt->Fill(lv2.Pt());

            mNegPhi->Fill(lv1.Phi());
            mNegPhi->Fill(lv2.Phi());


//MAKE LS PAIRS
                for ( TLorentzVector mixed_lv1 : LSNegBuffer ) {
                    TLorentzVector mixed_lv = mixed_lv1 + lv2 ;
                    double mixedPairPhi = calc_DeltaPhi(mixed_lv1, lv2);
                    if ( abs(mixed_lv.Rapidity()) < 1 && pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag !=0 && abs(ddTOFval) < 0.4 ) {
                        mLSLSPtvsM->Fill(mixed_lv.Pt(), mixed_lv.M());
                        mLSLSCos1phivsPTvsMass->Fill( 2*cos(1*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSLSCos2phivsPTvsMass->Fill( 2*cos(2*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSLSCos3phivsPTvsMass->Fill( 2*cos(3*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSLSCos4phivsPTvsMass->Fill( 2*cos(4*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                    }
                }
                for ( TLorentzVector mixed_lv2 : LSNegBuffer ) {
                    TLorentzVector mixed_lv = lv1 + mixed_lv2 ;
                    double mixedPairPhi = calc_DeltaPhi(lv1, mixed_lv2);
                    if ( abs(mixed_lv.Rapidity()) < 1 && pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag !=0 && abs(ddTOFval) < 0.4 ) {
                        mLSLSPtvsM->Fill(mixed_lv.Pt(), mixed_lv.M());
                        mLSLSCos1phivsPTvsMass->Fill( 2*cos(1*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSLSCos2phivsPTvsMass->Fill( 2*cos(2*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSLSCos3phivsPTvsMass->Fill( 2*cos(3*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSLSCos4phivsPTvsMass->Fill( 2*cos(4*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                    }
                }
//MAKE ULS PAIRS
                for ( TLorentzVector mixed_lv1 : LSPosBuffer ) {
                    TLorentzVector mixed_lv = mixed_lv1 + lv2 ;
                    double mixedPairPhi = calc_DeltaPhi(mixed_lv1, lv2);
                    if ( abs(mixed_lv.Rapidity()) < 1 && pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag !=0 && abs(ddTOFval) < 0.4 ) {
                        mLSULSPtvsM->Fill(mixed_lv.Pt(), mixed_lv.M());
                        mLSULSCos1phivsPTvsMass->Fill( 2*cos(1*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSULSCos2phivsPTvsMass->Fill( 2*cos(2*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSULSCos3phivsPTvsMass->Fill( 2*cos(3*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSULSCos4phivsPTvsMass->Fill( 2*cos(4*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );

                    }
                }
                for ( TLorentzVector mixed_lv2 : LSPosBuffer ) {
                    TLorentzVector mixed_lv = lv1 + mixed_lv2 ;
                    double mixedPairPhi = calc_DeltaPhi(lv1, mixed_lv2);
                    if ( abs(mixed_lv.Rapidity()) < 1 && pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag !=0 && abs(ddTOFval) < 0.4 ) {
                        mLSULSPtvsM->Fill(mixed_lv.Pt(), mixed_lv.M());
                        mLSULSCos1phivsPTvsMass->Fill( 2*cos(1*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSULSCos2phivsPTvsMass->Fill( 2*cos(2*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSULSCos3phivsPTvsMass->Fill( 2*cos(3*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mLSULSCos4phivsPTvsMass->Fill( 2*cos(4*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );

                    }
                }
                int buff_index = rng.Integer(LS_buffer_size - 1 );
                LSNegBuffer = update_buffer( LSNegBuffer, lv1, LS_buffer_size, buff_index );
                LSNegBuffer = update_buffer( LSNegBuffer, lv2, LS_buffer_size, buff_index );
            }

//UNLIKE SIGN SOURCE==================================================================================================================
            if ( pair->mChargeSum==0 ) {
            mAllTrackPt->Fill(lv1.Pt());
            mAllTrackPt->Fill(lv2.Pt());
            mULSTrackPt->Fill(lv1.Pt());
            mULSTrackPt->Fill(lv2.Pt());

            mPosPhi->Fill(lv1.Phi());
            mNegPhi->Fill(lv2.Phi());
//MAKE LS PAIRS
                for ( TLorentzVector mixed_lv1 : ULSNegBuffer ) {
                    TLorentzVector mixed_lv = mixed_lv1 + lv2 ;
                    double mixedPairPhi = calc_DeltaPhi(mixed_lv1, lv2);
                    if ( abs(mixed_lv.Rapidity()) < 1 && pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag !=0 && abs(ddTOFval) < 0.4 ) {
                        mULSLSPtvsM->Fill(mixed_lv.Pt(), mixed_lv.M());
                        mULSLSCos1phivsPTvsMass->Fill( 2*cos(1*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mULSLSCos2phivsPTvsMass->Fill( 2*cos(2*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mULSLSCos3phivsPTvsMass->Fill( 2*cos(3*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mULSLSCos4phivsPTvsMass->Fill( 2*cos(4*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        if(mixed_lv.M() > 0.4 && mixed_lv.M() < 0.76 && mixed_lv.Pt() > 0.1) {mULSLSphipt->Fill(mixedPairPhi, mixed_lv.Pt());}

                    }
                }
                for ( TLorentzVector mixed_lv2 : ULSPosBuffer ) {
                    TLorentzVector mixed_lv = lv1 + mixed_lv2 ;
                    double mixedPairPhi = calc_DeltaPhi(lv1, mixed_lv2);
                    if ( abs(mixed_lv.Rapidity()) < 1 && pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag !=0 && abs(ddTOFval) < 0.4 ) {
                        mULSLSPtvsM->Fill(mixed_lv.Pt(), mixed_lv.M());
                        mULSLSCos1phivsPTvsMass->Fill( 2*cos(1*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mULSLSCos2phivsPTvsMass->Fill( 2*cos(2*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mULSLSCos3phivsPTvsMass->Fill( 2*cos(3*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mULSLSCos4phivsPTvsMass->Fill( 2*cos(4*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        if(mixed_lv.M() < 0.76) {mULSLSphipt->Fill(mixedPairPhi, mixed_lv.Pt());}

                    }
                }
//MAKE ULS PAIRS
                for ( TLorentzVector mixed_lv1 : ULSPosBuffer ) {
                    TLorentzVector mixed_lv = mixed_lv1 + lv2 ;
                    double mixedPairPhi = calc_DeltaPhi(mixed_lv1, lv2);
                    if ( abs(mixed_lv.Rapidity()) < 1 && pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag !=0 && abs(ddTOFval) < 0.4 ) {
                        mULSULSPtvsM->Fill(mixed_lv.Pt(), mixed_lv.M());
                        mULSULSCos1phivsPTvsMass->Fill( 2*cos(1*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mULSULSCos2phivsPTvsMass->Fill( 2*cos(2*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mULSULSCos3phivsPTvsMass->Fill( 2*cos(3*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mULSULSCos4phivsPTvsMass->Fill( 2*cos(4*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        if(mixed_lv.M() > 0.4 && mixed_lv.M() < 1){
                            double acoval = 1 - abs(mixed_lv1.Phi()-lv2.Phi())/M_PI;
                            acoULSULS->Fill(acoval);
                        }
                    }
                }
                for ( TLorentzVector mixed_lv2 : ULSNegBuffer ) {
                    TLorentzVector mixed_lv = lv1 + mixed_lv2 ;
                    double mixedPairPhi = calc_DeltaPhi(lv1, mixed_lv2);
                    if ( abs(mixed_lv.Rapidity()) < 1 && pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag !=0 && abs(ddTOFval) < 0.4 ) {
                        mULSULSPtvsM->Fill(mixed_lv.Pt(), mixed_lv.M());
                        mULSULSCos1phivsPTvsMass->Fill( 2*cos(1*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mULSULSCos2phivsPTvsMass->Fill( 2*cos(2*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mULSULSCos3phivsPTvsMass->Fill( 2*cos(3*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mULSULSCos4phivsPTvsMass->Fill( 2*cos(4*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                    }
                }
                int buff_index = rng.Integer(ULS_buffer_size - 1 );
                ULSPosBuffer = update_buffer( ULSPosBuffer, lv1, ULS_buffer_size, buff_index );
                ULSNegBuffer = update_buffer( ULSNegBuffer, lv2, ULS_buffer_size, buff_index );
            }
        }
    }

fo->cd();

fo->Write();
}
