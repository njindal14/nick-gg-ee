#include <iostream>
#include "FemtoPairFormat.h"
double const pi = 3.14159265358979323;
double c = 3.0e1; //in cm/ns
double me2 = pow(0.00051,2);
double m2 = .00051;

int ican2 = 0;
void makeCanvas()  {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican2++ ), "", 900, 600);
    can->SetTopMargin(0.08);
    can->SetRightMargin(0.15);
    can->SetBottomMargin(0.15);
}

double calc_Phi( TLorentzVector lv1, TLorentzVector lv2) {
    TLorentzVector lvPlus = lv1 + lv2;
    TLorentzVector lvMinus = lv1 - lv2;
    return lvMinus.DeltaPhi(lvPlus);
}

void backgroundCorrection() {

    TFile *myFile = TFile::Open("/Users/Nick/STAR/breit-wheeler/rootFiles/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");
    TLorentzVector lv1, lv2, lv, lvn;

    auto * mPtUnlike = new TH1F("mPtUnlike", "Pair pTUnlike", 20, 0, 1);
    auto * mPtLike = new TH1D("mPtLike", "mPtLike", 20, 0, 1);
    auto * cos2phivPtLike = new TH2F("Cos2#phivPtLike", "A_{2#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2}", 400, -2, 2, 4, 0, 0.3);
    auto * cos4phivPtLike = new TH2F("Cos4#phivPtLike", "A_{4#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2}", 400, -2, 2, 4, 0, 0.3);
    auto * cos2phivPtUnlike = new TH2F("Cos2#phivPtUnlike", "A_{2#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2}", 400, -2, 2, 30, 0, 0.6);
    auto * cos4phivPtUnlike = new TH2F("Cos4#phivPtUnlike", "A_{4#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2}", 400, -2, 2, 30, 0, 0.6);

    TFile * starlight_out = new TFile("/Users/Nick/STAR/starlight/utils/SL_plotsBetter.root");
    TH1F * negpt = (TH1F*)starlight_out->Get("mNegPt");
    TH1F * pospt = (TH1F*)starlight_out->Get("mPosPt"); 

    TH1F * negeta = (TH1F*)starlight_out->Get("mNegEta");
    TH1F * poseta = (TH1F*)starlight_out->Get("mPosEta"); 

    TH1F * negphi = (TH1F*)starlight_out->Get("mNegPhi");
    TH1F * posphi = (TH1F*)starlight_out->Get("mPosPhi"); 


    while(myReader.Next()){

        double chiee = pow( pair->d1_mNSigmaElectron, 2 ) + pow( pair->d2_mNSigmaElectron, 2 );
        double chipipi = pow( pair -> d1_mNSigmaPion, 2) + pow( pair -> d2_mNSigmaPion, 2);

        lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.00051 );
        lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.00051 ); 

        lv = lv1 + lv2;
        lvn = lv1 - lv2;

        Float_t p1_2 = pow(lv1.P(), 2);
        Float_t p2_2 = pow(lv2.P(), 2);
        Float_t dTofVal = pair->d1_mTof - pair->d2_mTof;
        Float_t dTofexpVal = pair->d1_mLength/c * sqrt(1 + me2/p1_2) - pair->d2_mLength/c * sqrt(1 + me2/p2_2);
        Float_t ddTofVal = dTofVal - dTofexpVal;

        if(lv1.Pt() < 0.2 || lv2.Pt() < 0.2) continue;

        if( fabs(lv1.Eta()) > 1 || fabs(lv2.Eta()) > 1 || fabs(lv.Rapidity()) > 1) continue;

        if( /*fabs(pair->mVertexZ) < 100 &&  pair->mGRefMult <= 4 && pair->d1_mDCA < 1 && pair->d2_mDCA < 1 && 
        pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag!=0 &&*/ fabs(ddTofVal) < 0.4 && ddTofVal !=0 && chiee < 10 && 3*chiee < chipipi) {
            
            if(lv.M() > 0.4 && lv.M() < 0.76){
                double phival = calc_Phi(lv1, lv2);

                if(pair->mChargeSum == 0 && fabs(pair->mVertexZ) < 100 &&  pair->mGRefMult <= 4 && pair->d1_mDCA < 1 && pair->d2_mDCA < 1 && 
        pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag!=0) {
                    mPtUnlike->Fill(lv.Pt());
                    cos2phivPtUnlike->Fill(2*cos(2*phival), lv.Pt());
                    cos4phivPtUnlike->Fill(2*cos(4*phival), lv.Pt());
                }
                else if(pair->mChargeSum != 0 && pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag!=0) {
                    mPtLike->Fill(lv.Pt());
                    cos2phivPtLike->Fill(2*cos(2*phival), lv.Pt());
                    cos4phivPtLike->Fill(2*cos(4*phival), lv.Pt());

                }
            }

        }

    }

    //created resampled histograms
    TH2F * cos2phivPtResample= new TH2F("cos2phivPtResample", "A_{2#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2}", 400, -2, 2, 30, 0, 0.6);
    TH2F * cos4phivPtResample = new TH2F("cos4phivPtResample", "A_{4#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2}", 400, -2, 2, 30, 0, 0.6);
    TH1F * resamplePt = new TH1F("resamplePt", "", 20, 0, 1);

    TH2F * cos2phivPtResampleUniform= new TH2F("cos2phivPtResampleUniform", "A_{2#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2} Uni sampling", 400, -2, 2, 30, 0, 0.6);
    TH2F * cos4phivPtResampleUniform = new TH2F("cos4phivPtResampleUniform", "A_{4#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2} Uni sampling", 400, -2, 2, 30, 0, 0.6);
    TH1F * resamplePtUniform = new TH1F("resamplePtUniform", "", 20, 0, 1);


    pospt->GetXaxis()->SetRangeUser(.2, 1);
    negpt->GetXaxis()->SetRangeUser(.2, 1);

    //now sample from the starlight dists
    for ( int i = 0; i<100000000; i++) {
        TLorentzVector pos, neg, pair;
        TLorentzVector posUniform, negUniform, pairUniform;
        pos.SetPtEtaPhiM( pospt->GetRandom(),  poseta->GetRandom(), posphi->GetRandom(), 0.00051);
        neg.SetPtEtaPhiM( negpt->GetRandom(),  negeta->GetRandom(), negphi->GetRandom(), 0.00051);
        posUniform.SetPtEtaPhiM(gRandom->Uniform(0, 1), poseta->GetRandom() /*gRandom->Uniform(-1,1)*/, gRandom->Uniform(-3.14,3.14), .00051);
        negUniform.SetPtEtaPhiM(gRandom->Uniform(0, 1), negeta->GetRandom() /*gRandom->Uniform(-1,1)*/, gRandom->Uniform(-3.14,3.14), .00051);
        pair = pos+neg;
        pairUniform = posUniform + negUniform;
        if(i%1000000 == 0) std::cout << "Sampled event number: " << i << "\n";
        double phi = calc_Phi(pos, neg);
        double uniformPhi = calc_Phi(posUniform, negUniform);
        if(pair.M() > 0.4 && pair.M() < 0.76){
            cos2phivPtResample->Fill( 2*cos(2*phi), pair.Pt() );
            cos4phivPtResample->Fill( 2*cos(4*phi), pair.Pt() );
            resamplePt->Fill(pair.Pt());
        }
        if(pairUniform.M() > 0.4 && pairUniform.M() < 0.76){
            cos2phivPtResampleUniform->Fill( 2*cos(2*uniformPhi), pairUniform.Pt() );
            cos4phivPtResampleUniform->Fill( 2*cos(4*uniformPhi), pairUniform.Pt() );
            resamplePtUniform->Fill(pairUniform.Pt());
        }

    }

    //normalizing to data at a certain pt
    resamplePt->Scale(0.001);

    //look at signal/background and signal/sqrt(signal+background)
    makeCanvas();
    TH1F * SBR = (TH1F*)mPtUnlike->Clone();
    SBR->Divide(resamplePt);
    SBR->SetTitle("Signal/ Rescaled Background Ratio");
    SBR->GetXaxis()->SetTitle("pT (GeV/c)");
    SBR->Draw("PE");

    makeCanvas();
    TH1F *SNR = (TH1F*)mPtUnlike->Clone();
    TH1F* tot = (TH1F*)mPtUnlike->Clone();
    tot->Add(resamplePt);
    for(int i = 0; i <= tot->GetNbinsX(); ++i){
        tot->SetBinContent(i, sqrt(tot->GetBinContent(i)));
    }
    SNR->Divide(tot);
    SNR->SetTitle("Signal/sqrt(signal+background) (significance)");
    SNR->GetXaxis()->SetTitle("pT (GeV/c)");
    SNR->Draw("PE");

    
    
    makeCanvas();
    mPtUnlike->Draw("PE");
    mPtUnlike->SetTitle("pT signal and background, 0.4 < M < 0.76 Gev");
    mPtUnlike->GetXaxis()->SetTitle("pT (GeV/c)");
    mPtUnlike->GetYaxis()->SetTitle("Counts");
    resamplePt->SetLineColor(kRed);
    //resamplePt->Scale(1/resamplePt->Integral());
    resamplePt->Draw("PE;same");
    TLegend * leg = new TLegend(.8, .7, 1, .9);
    leg->SetHeader("Legend","C"); // option "C" allows to center the header
    leg->AddEntry(mPtUnlike,"Unlike sign");
    leg->AddEntry(resamplePt, "Resampled pair pT");
    leg->Draw("same");
    


    auto *m2Ptcos4phimomentsLikesign = cos4phivPtLike->ProfileY("m2Ptcos4phimomentsLike",1, -1);
    auto *m2Ptcos2phimomentsLikesign = cos2phivPtLike->ProfileY("m2Ptcos2phimomentsLike",1, -1);

    auto *m2Ptcos4phimomentsUnlikesign = cos4phivPtUnlike->ProfileY("m2Ptcos4phimomentsUnlike",1, -1);
    auto *m2Ptcos2phimomentsUnlikesign = cos2phivPtUnlike->ProfileY("m2Ptcos2phimomentsUnlike",1, -1);

    auto * m2phiMomentsResample = cos2phivPtResample->ProfileY("m2phiMomentsResample", 1, -1);
    auto * m4phiMomentsResample = cos4phivPtResample->ProfileY("m4phiMomentsResample", 1, -1);

    auto * m2phiMomentsResampleUniform = cos2phivPtResampleUniform->ProfileY("m2phiMomentsResampleUniform", 1, -1);
    auto * m4phiMomentsResampleUniform = cos4phivPtResampleUniform->ProfileY("m4phiMomentsResampleUniform", 1, -1);




    makeCanvas();
    m2Ptcos2phimomentsLikesign->GetXaxis()->SetTitle("pT (GeV/c)");
    m2Ptcos2phimomentsLikesign->GetYaxis()->SetTitle("A_{2#phi}");
    m2Ptcos2phimomentsUnlikesign->SetLineColor(kGreen);
    m2Ptcos2phimomentsUnlikesign->Draw("PE;same");
    m2phiMomentsResample->SetLineColor(kRed);
    m2phiMomentsResample->Draw("PE;same");
    m2phiMomentsResampleUniform->SetLineColor(kBlue);
    m2phiMomentsResampleUniform->Draw("PE;same");
    TLegend * leg2 = new TLegend(.8, .7, 1, .9);
    leg2->SetHeader("Legend","C"); // option "C" allows to center the header
    leg2->AddEntry(m2Ptcos2phimomentsUnlikesign, "Unlike sign, all cuts");
    leg2->AddEntry(m2phiMomentsResample, "slight random sampled");
    leg2->AddEntry(m2phiMomentsResampleUniform, "slight random sampled uniform");
    leg2->Draw("same");

    makeCanvas();
    m2Ptcos4phimomentsLikesign->GetXaxis()->SetTitle("pT (GeV/c)");
    m2Ptcos4phimomentsLikesign->GetYaxis()->SetTitle("A_{4#phi}");
    m2Ptcos4phimomentsUnlikesign->SetLineColor(kGreen);
    m2Ptcos4phimomentsUnlikesign->Draw("PE;same");
    m4phiMomentsResample->SetLineColor(kRed);
    m4phiMomentsResample->Draw("PE;same");
    m4phiMomentsResampleUniform->SetLineColor(kBlue);
    m4phiMomentsResampleUniform->Draw("PE;same");
    TLegend * leg4 = new TLegend(.8, .7, 1, .9);
    leg4->SetHeader("Legend","C"); // option "C" allows to center the header
    leg4->AddEntry(m2Ptcos4phimomentsUnlikesign, "Unlike sign, all cuts");
    leg4->AddEntry(m4phiMomentsResample, "slight random sampled");
    leg4->AddEntry(m4phiMomentsResampleUniform, "slight random sampled uniform");
    leg4->Draw("same");

    makeCanvas();
    pospt->Draw("PE");
    negpt->SetLineColor(kRed);
    negpt->Draw("PE'same");

    makeCanvas();
    poseta->Draw("PE");
    negeta->SetLineColor(kRed);
    negeta->Draw("PE;same");

    makeCanvas();
    posphi->Draw("PE");
    negphi->SetLineColor(kRed);
    negphi->Draw("PE;same");





}
