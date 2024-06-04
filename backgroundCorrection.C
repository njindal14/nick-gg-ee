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

double phiFit(double *x, double *par){
    double fitval;
    fitval = par[0]*(1+par[1]*cos(x[0]) + par[2]*cos(2*x[0]) + par[3]*cos(3*x[0]) + par[4]*cos(4*x[0]));
    return fitval;
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

    auto * cos2phivMassUnlike = new TH2F("Cos2#phivMassUnlike", "", 400, -2, 2, 50, 0, 2);
    auto * cos4phivMassUnlike = new TH2F("Cos4#phivMassUnlike", "", 400, -2, 2, 50, 0, 2);

    auto * cos1phivPtUnlike = new TH2F("Cos1#phivPtUnlike", "A_{1#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2}", 400, -2, 2, 30, 0, 0.6);
    auto * cos3phivPtUnlike = new TH2F("Cos3#phivPtUnlike", "A_{3#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2}", 400, -2, 2, 30, 0, 0.6);

    auto * resampledDeltaPhi = new TH1F("resampledDeltaPhi", "resampledDeltaPhi, 0.4 < M_{ee} < 0.76 GeV/c^{2}", 50, -3.14, 3.14);


    auto * mMassUnlike = new TH1F("mMassUnlike", "Pair Mass Data", 60, 0, 3);


    TFile * starlight_out = new TFile("/Users/Nick/STAR/starlight/utils/SL_plotsBetter.root");
    TH1F * negpt = (TH1F*)starlight_out->Get("mNegPt");
    TH1F * pospt = (TH1F*)starlight_out->Get("mPosPt"); 

    TH1F * negeta = (TH1F*)starlight_out->Get("mNegEta");
    TH1F * poseta = (TH1F*)starlight_out->Get("mPosEta"); 

    TH1F * negphi = (TH1F*)starlight_out->Get("mNegPhi");
    TH1F * posphi = (TH1F*)starlight_out->Get("mPosPhi"); 


    TFile * mixedeventmoments = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/MixedEventMoments.root");
    TH1F * cos1phimixed = (TH1F*)mixedeventmoments->Get("mULSULSCos1phivsPTvsMass_p0");
    TH1F * cos2phimixed = (TH1F*)mixedeventmoments->Get("mULSULSCos2phivsPTvsMass_p1");
    TH1F * cos3phimixed = (TH1F*)mixedeventmoments->Get("mULSULSCos3phivsPTvsMass_p2");
    TH1F * cos4phimixed = (TH1F*)mixedeventmoments->Get("mULSULSCos4phivsPTvsMass_p3");


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
            
            mMassUnlike->Fill(lv.M());

            if(lv.M() > 0.4 && lv.M() < 0.76){
                double phival = calc_Phi(lv1, lv2);

                if(pair->mChargeSum == 0 && fabs(pair->mVertexZ) < 100 &&  pair->mGRefMult <= 4 && pair->d1_mDCA < 1 && pair->d2_mDCA < 1 && 
        pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag!=0) {
                    mPtUnlike->Fill(lv.Pt());
                    cos2phivPtUnlike->Fill(2*cos(2*phival), lv.Pt());
                    cos4phivPtUnlike->Fill(2*cos(4*phival), lv.Pt());
                    cos1phivPtUnlike->Fill(2*cos(phival), lv.Pt());
                    cos3phivPtUnlike->Fill(2*cos(3*phival), lv.Pt());
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

    TH2F * cos1phivPtResample= new TH2F("cos1phivPtResample", "A_{1#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2}", 400, -2, 2, 30, 0, 0.6);
    TH2F * cos3phivPtResample = new TH2F("cos3phivPtResample", "A_{3#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2}", 400, -2, 2, 30, 0, 0.6);

    TH1F * acoResample = new TH1F("acoResample", "acoResample", 20, 0, .2);

    TH1F * resamplePt = new TH1F("resamplePt", "", 20, 0, 1);

    TH2F * cos2phivPtResampleUniform= new TH2F("cos2phivPtResampleUniform", "A_{2#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2} Uni sampling", 400, -2, 2, 20, 0, 1);
    TH2F * cos4phivPtResampleUniform = new TH2F("cos4phivPtResampleUniform", "A_{4#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2} Uni sampling", 400, -2, 2, 20, 0, 1);

    TH2F * cos1phivPtResampleUniform = new TH2F("cos1phivPtResampleUniform", "A_{1#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2} Uni sampling", 400, -2, 2, 30, 0, 0.6);
    TH2F * cos3phivPtResampleUniform = new TH2F("cos3phivPtResampleUniform", "A_{3#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2} Uni sampling", 400, -2, 2, 30, 0, 0.6);


    auto * cos2phivPtPosLS = new TH2F("Cos2#phivPtPosLS", "A_{2#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2}", 400, -2, 2, 30, 0, 0.6);
    auto * cos4phivPtPosLS = new TH2F("Cos4#phivPtPosLS", "A_{4#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2}", 400, -2, 2, 30, 0, 0.6);

    auto * cos1phivPtPosLS = new TH2F("Cos1#phivPtPosLS", "A_{1#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2}", 400, -2, 2, 30, 0, 0.6);
    auto * cos3phivPtPosLS = new TH2F("Cos3#phivPtPosLS", "A_{3#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2}", 400, -2, 2, 30, 0, 0.6);

    TH1F * resamplePtUniform = new TH1F("resamplePtUniform", "", 20, 0, 1);

    TH1F * resampleMass = new TH1F("resampleMass", "", 60, 0, 3);
    TH1F * resampleMassUniform = new TH1F("resampleMassUniform", "", 60, 0, 3);


    pospt->GetXaxis()->SetRangeUser(.2, 1);
    negpt->GetXaxis()->SetRangeUser(.2, 1);

    //now sample from the starlight dists
    for ( int i = 0; i<10000000; i++) {
        TLorentzVector pos, neg, pair;
        TLorentzVector posUniform, negUniform, pairUniform;

        //like sign sampling from starlight tracks
        TLorentzVector pos1, pos2, pairPosLS;
        TLorentzVector neg1, neg2;

        //sample and make sure they pass kinematic cuts
        double posptval, negptval, pos1ptval, pos2ptval, posUptval, negUptval;

        posptval = pospt->GetRandom();
        while(posptval < 0.2){posptval = pospt->GetRandom();}
        negptval = negpt->GetRandom();
        while(negptval < 0.2){negptval = negpt->GetRandom();}

        pos1ptval = pospt->GetRandom();
        while(pos1ptval < 0.2){pos1ptval = pospt->GetRandom();}
        pos2ptval = negpt->GetRandom();
        while(pos2ptval < 0.2){pos2ptval = negpt->GetRandom();}

        pos.SetPtEtaPhiM( posptval,  poseta->GetRandom(), posphi->GetRandom(), 0.00051);
        neg.SetPtEtaPhiM( negptval,  negeta->GetRandom(), negphi->GetRandom(), 0.00051);

        pos1.SetPtEtaPhiM( pos1ptval,  poseta->GetRandom(), posphi->GetRandom(), 0.00051);
        pos2.SetPtEtaPhiM( pos2ptval,  poseta->GetRandom(), posphi->GetRandom(), 0.00051); 

        posUniform.SetPtEtaPhiM(gRandom->Uniform(0, 1), /*poseta->GetRandom()*/ gRandom->Uniform(-1,1), gRandom->Uniform(-3.14,3.14), .00051);
        negUniform.SetPtEtaPhiM(gRandom->Uniform(0, 1), /*negeta->GetRandom()*/ gRandom->Uniform(-1,1), gRandom->Uniform(-3.14,3.14), .00051);

        pair = pos+neg;
        pairUniform = posUniform + negUniform;
        pairPosLS = pos1 + pos2;

        if(i%1000000 == 0) std::cout << "Sampled event number: " << i << "\n";
        double phi = calc_Phi(pos, neg);
        double uniformPhi = calc_Phi(posUniform, negUniform);
        double posLSphi = calc_Phi(pos1, pos2);
        resampleMassUniform->Fill(pairUniform.M());
        resampleMass->Fill(pair.M());

        cos2phivMassUnlike->Fill(2*cos(2*phi), pair.M());
        cos4phivMassUnlike->Fill(2*cos(4*phi), pair.M());

        if(pair.M() > 0.4 && pair.M() < 0.76){
            cos2phivPtResample->Fill( 2*cos(2*phi), pair.Pt() );
            cos4phivPtResample->Fill( 2*cos(4*phi), pair.Pt() );

            cos1phivPtResample->Fill( 2*cos(1*phi), pair.Pt() );
            cos3phivPtResample->Fill( 2*cos(3*phi), pair.Pt() );
            resamplePt->Fill(pair.Pt());
            resampledDeltaPhi->Fill(phi);
        }
        if(pairUniform.M() > 0.4 && pairUniform.M() < 0.76){
            cos2phivPtResampleUniform->Fill( 2*cos(2*uniformPhi), pairUniform.Pt() );
            cos4phivPtResampleUniform->Fill( 2*cos(4*uniformPhi), pairUniform.Pt() );

            cos1phivPtResampleUniform->Fill( 2*cos(1*uniformPhi), pairUniform.Pt() );
            cos3phivPtResampleUniform->Fill( 2*cos(3*uniformPhi), pairUniform.Pt() );

            resamplePtUniform->Fill(pairUniform.Pt());
        }

        if(pairPosLS.M() > 0.4 && pairPosLS.M() < 0.76){
            cos2phivPtPosLS->Fill( 2*cos(2*posLSphi), pairPosLS.Pt());
            cos4phivPtPosLS->Fill( 2*cos(4*posLSphi), pairPosLS.Pt());

            cos1phivPtPosLS->Fill( 2*cos(1*posLSphi), pairPosLS.Pt());
            cos3phivPtPosLS->Fill( 2*cos(3*posLSphi), pairPosLS.Pt());
        }

        if(pair.M() > 0.4 && pair.M() < 1){
            double acoval = 1 - abs(pos.Phi() - neg.Phi())/M_PI;
            acoResample->Fill(acoval);
        }

    }

    //normalizing to data at a certain pt
    resamplePt->Scale( mPtUnlike->GetBinContent( mPtUnlike->FindBin(0.47) ) / resamplePt->GetBinContent( resamplePt->FindBin(.47) ) );

    //look at signal/background and signal/sqrt(signal+background)
    makeCanvas();
    TH1F * SBR = (TH1F*)mPtUnlike->Clone();
    SBR->SetTitle("Signal/ Background Ratio");
    SBR->GetXaxis()->SetTitle("pT (GeV/c)");
    SBR->Draw("PE");
    SBR->SetStats(false);
    gPad->SetLogy();
    gPad->Print("note_plots/mixed_event_plots/SBR.png");

    makeCanvas();
    TH1F *SNR = (TH1F*)mPtUnlike->Clone();
    TH1F* tot = (TH1F*)mPtUnlike->Clone();
    tot->Add(resamplePt);
    for(int i = 0; i <= tot->GetNbinsX(); ++i){
        tot->SetBinContent(i, sqrt(tot->GetBinContent(i)));
    }
    SNR->Divide(tot);
    SNR->SetTitle("Significance");
    SNR->GetXaxis()->SetTitle("pT (GeV/c)");
    SNR->Draw("PE");
    SNR->SetStats(false);
    gPad->SetLogy();
    gPad->Print("note_plots/mixed_event_plots/significance.png");

    
    
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


    makeCanvas();
    TH1F * purity = (TH1F*)mPtUnlike->Clone();
    TH1F * denom = (TH1F*)resamplePt->Clone();
    denom->Add(mPtUnlike);
    purity->Divide(denom);
    purity->SetTitle("Pair Purity Using Mixed Event");
    purity->Draw("PE");
    gPad->Print("note_plots/mixed_event_plots/purity.png");
    

    auto *m2Ptcos4phimomentsLikesign = cos4phivPtLike->ProfileY("m2Ptcos4phimomentsLike",1, -1);
    auto *m2Ptcos2phimomentsLikesign = cos2phivPtLike->ProfileY("m2Ptcos2phimomentsLike",1, -1);

    auto *m2Ptcos4phimomentsUnlikesign = cos4phivPtUnlike->ProfileY("m2Ptcos4phimomentsUnlike",1, -1);
    auto *m2Ptcos2phimomentsUnlikesign = cos2phivPtUnlike->ProfileY("m2Ptcos2phimomentsUnlike",1, -1);

    auto *m2Ptcos1phimomentsUnlikesign = cos1phivPtUnlike->ProfileY("m2Ptcos1phimomentsUnlike",1, -1);
    auto *m2Ptcos3phimomentsUnlikesign = cos3phivPtUnlike->ProfileY("m2Ptcos3phimomentsUnlike",1, -1);

    auto * m2phiMomentsResample = cos2phivPtResample->ProfileY("m2phiMomentsResample", 1, -1);
    auto * m4phiMomentsResample = cos4phivPtResample->ProfileY("m4phiMomentsResample", 1, -1);
    
    auto * m1phiMomentsResample = cos1phivPtResample->ProfileY("m1phiMomentsResample", 1, -1);
    auto * m3phiMomentsResample = cos3phivPtResample->ProfileY("m3phiMomentsResample", 1, -1);

    auto * m2phiMomentsResampleUniform = cos2phivPtResampleUniform->ProfileY("m2phiMomentsResampleUniform", 1, -1);
    auto * m4phiMomentsResampleUniform = cos4phivPtResampleUniform->ProfileY("m4phiMomentsResampleUniform", 1, -1);

    auto * m1phiMomentsResampleUniform = cos1phivPtResampleUniform->ProfileY("m1phiMomentsResampleUniform", 1, -1);
    auto * m3phiMomentsResampleUniform = cos3phivPtResampleUniform->ProfileY("m3phiMomentsResampleUniform", 1, -1);

    auto * m2phiMomentsPosLS = cos2phivPtPosLS->ProfileY("m2phiMomentsPosLS", 1, -1);
    auto * m4phiMomentsPosLS = cos4phivPtPosLS->ProfileY("m4phiMomentsPosLS", 1, -1);

    auto * m1phiMomentsPosLS = cos1phivPtPosLS->ProfileY("m1phiMomentsPosLS", 1, -1);
    auto * m3phiMomentsPosLS = cos3phivPtPosLS->ProfileY("m3phiMomentsPosLS", 1, -1);



    makeCanvas();
    m2Ptcos1phimomentsUnlikesign->SetLineColor(kGreen);
    m2Ptcos1phimomentsUnlikesign->Draw("PE;same");
    m1phiMomentsResample->SetLineColor(kRed);
    m1phiMomentsResample->Draw("PE;same");
    m1phiMomentsResampleUniform->SetLineColor(kBlue);
    m1phiMomentsResampleUniform->Draw("PE;same");
    m1phiMomentsPosLS->SetLineColor(kBlack);
    m1phiMomentsPosLS->Draw("PE;same");
    cos1phimixed->SetLineColor(kPink);
    cos1phimixed->Draw("PE;same");
    m2Ptcos1phimomentsUnlikesign->GetXaxis()->SetTitle("pT (GeV/c)");
    m2Ptcos1phimomentsUnlikesign->GetYaxis()->SetTitle("A_{1#phi}");

    TLegend * leg1 = new TLegend(.8, .7, 1, .9);
    leg1->SetHeader("Legend","C"); // option "C" allows to center the header
    leg1->AddEntry(m2Ptcos1phimomentsUnlikesign, "Unlike sign, all cuts");
    leg1->AddEntry(m1phiMomentsResample, "slight random sampled");
    leg1->AddEntry(m1phiMomentsResampleUniform, "slight random sampled uniform");
    leg1->AddEntry(m1phiMomentsPosLS, "slight LS pos");
    leg1->AddEntry(cos1phimixed, "Mixed Event");
    leg1->Draw("same");
    m2Ptcos1phimomentsUnlikesign->SetStats(false);
    gPad->Print("note_plots/mixed_event_plots/onephi_with_mixed.png");


    makeCanvas();
    //m2Ptcos2phimomentsLikesign->GetXaxis()->SetTitle("pT (GeV/c)");
    //m2Ptcos2phimomentsLikesign->GetYaxis()->SetTitle("A_{2#phi}");
    m2Ptcos2phimomentsUnlikesign->SetLineColor(kGreen);
    m2Ptcos2phimomentsUnlikesign->Draw("PE;same");
    m2phiMomentsResample->SetLineColor(kRed);
    m2phiMomentsResample->Draw("PE;same");
    m2phiMomentsResampleUniform->SetLineColor(kBlue);
    m2phiMomentsResampleUniform->Draw("PE;same");
    m2phiMomentsPosLS->SetLineColor(kBlack);
    m2phiMomentsPosLS->Draw("PE;same");
    cos2phimixed->SetLineColor(kPink);
    cos2phimixed->Draw("PE;same");
    m2Ptcos2phimomentsUnlikesign->GetXaxis()->SetTitle("pT (GeV/c)");
    m2Ptcos2phimomentsUnlikesign->GetYaxis()->SetTitle("A_{2#phi}");

    TLegend * leg2 = new TLegend(.8, .7, 1, .9);
    leg2->SetHeader("Legend","C"); // option "C" allows to center the header
    leg2->AddEntry(m2Ptcos2phimomentsUnlikesign, "Unlike sign, all cuts");
    leg2->AddEntry(m2phiMomentsResample, "slight random sampled");
    leg2->AddEntry(m2phiMomentsResampleUniform, "slight random sampled uniform");
    leg2->AddEntry(m2phiMomentsPosLS, "slight LS pos");
    leg2->AddEntry(cos2phimixed, "Mixed Event");
    leg2->Draw("same");
    m2Ptcos2phimomentsUnlikesign->SetStats(false);
    gPad->Print("note_plots/mixed_event_plots/twophi_with_mixed.png");

    makeCanvas();
    m2Ptcos3phimomentsUnlikesign->SetLineColor(kGreen);
    m2Ptcos3phimomentsUnlikesign->Draw("PE;same");
    m3phiMomentsResample->SetLineColor(kRed);
    m3phiMomentsResample->Draw("PE;same");
    m3phiMomentsResampleUniform->SetLineColor(kBlue);
    m3phiMomentsResampleUniform->Draw("PE;same");
    m3phiMomentsPosLS->SetLineColor(kBlack);
    m3phiMomentsPosLS->Draw("PE;same");
    cos3phimixed->SetLineColor(kPink);
    cos3phimixed->Draw("PE;same");
    m2Ptcos3phimomentsUnlikesign->GetXaxis()->SetTitle("pT (GeV/c)");
    m2Ptcos3phimomentsUnlikesign->GetYaxis()->SetTitle("A_{3#phi}");

    TLegend * leg3 = new TLegend(.8, .7, 1, .9);
    leg3->SetHeader("Legend","C"); // option "C" allows to center the header
    leg3->AddEntry(m2Ptcos3phimomentsUnlikesign, "Unlike sign, all cuts");
    leg3->AddEntry(m3phiMomentsResample, "slight random sampled");
    leg3->AddEntry(m3phiMomentsResampleUniform, "slight random sampled uniform");
    leg3->AddEntry(m3phiMomentsPosLS, "slight LS pos");
    leg3->AddEntry(cos3phimixed, "Mixed Event");
    leg3->Draw("same");
    m2Ptcos3phimomentsUnlikesign->SetStats(false);
    gPad->Print("note_plots/mixed_event_plots/threephi_with_mixed.png");

    makeCanvas();
    m2Ptcos4phimomentsUnlikesign->SetLineColor(kGreen);
    m2Ptcos4phimomentsUnlikesign->Draw("PE;same");
    m4phiMomentsResample->SetLineColor(kRed);
    m4phiMomentsResample->Draw("PE;same");
    m4phiMomentsResampleUniform->SetLineColor(kBlue);
    m4phiMomentsResampleUniform->Draw("PE;same");
    m4phiMomentsPosLS->SetLineColor(kBlack);
    m4phiMomentsPosLS->Draw("PE;same");
    cos4phimixed->SetLineColor(kPink);
    cos4phimixed->Draw("PE;same");
    m2Ptcos4phimomentsUnlikesign->GetXaxis()->SetTitle("pT (GeV/c)");
    m2Ptcos4phimomentsUnlikesign->GetYaxis()->SetTitle("A_{4#phi}");
    TLegend * leg4 = new TLegend(.8, .7, 1, .9);
    leg4->SetHeader("Legend","C"); // option "C" allows to center the header
    leg4->AddEntry(m2Ptcos4phimomentsUnlikesign, "Unlike sign, all cuts");
    leg4->AddEntry(m4phiMomentsResample, "slight random sampled");
    leg4->AddEntry(m4phiMomentsResampleUniform, "slight random sampled uniform");
    leg4->AddEntry(m4phiMomentsPosLS, "slight LS pos");
    leg4->AddEntry(cos4phimixed, "Mixed Event");
    leg4->Draw("same");
    m2Ptcos4phimomentsUnlikesign->SetStats(false);
    gPad->Print("note_plots/mixed_event_plots/fourphi_with_mixed.png");

    
    makeCanvas();
    pospt->Draw("PE");
    negpt->SetLineColor(kRed);
    negpt->Draw("PE'same");
    pospt->GetXaxis()->SetTitle("pT (GeV/c)");
    pospt->GetYaxis()->SetTitle("Counts");
    pospt->SetTitle("Track pT STARLight Tracks");
    pospt->SetStats(false);
    gPad->Print("note_plots/mixed_event_plots/track_pt.png");
    


    makeCanvas();
    poseta->Draw("PE");
    negeta->SetLineColor(kRed);
    negeta->Draw("PE;same");
    poseta->GetXaxis()->SetTitle("#eta");
    poseta->GetYaxis()->SetTitle("Counts");
    poseta->SetTitle("#eta STARLight Tracks");
    poseta->SetStats(false);
    gPad->Print("note_plots/mixed_event_plots/track_etas.png");
    


    makeCanvas();
    posphi->Draw("PE");
    negphi->SetLineColor(kRed);
    negphi->Draw("PE;same");
    posphi->GetXaxis()->SetTitle("#phi");
    posphi->GetYaxis()->SetTitle("Counts");
    posphi->SetTitle("Azimuthal Phi STARLight Tracks");
    posphi->SetStats(false);
    gPad->Print("note_plots/mixed_event_plots/track_phis.png");

    makeCanvas();
    resampledDeltaPhi->SetLineColor(kBlack);
    resampledDeltaPhi->GetXaxis()->SetTitle("#Delta#phi");
    resampledDeltaPhi->GetYaxis()->SetTitle("Counts");
    resampledDeltaPhi->Draw("PE");
    resampledDeltaPhi->SetTitle("Pair #Delta#phi Mixed Event");
    resampledDeltaPhi->SetStats(false);
    gPad->Print("note_plots/mixed_event_plots/mixed_delta_phi.png");



    makeCanvas();
    resampleMass->Scale(1/resampleMass->Integral(resampleMass->FindBin(.4), resampleMass->FindBin(.76), "width"));
    resampleMassUniform->Scale(1/resampleMassUniform->Integral(resampleMassUniform->FindBin(.4), resampleMassUniform->FindBin(.76), "width"));
    mMassUnlike->Scale(1/mMassUnlike->Integral(mMassUnlike->FindBin(.4), mMassUnlike->FindBin(.76), "width"));

    resampleMass->SetLineColor(kBlack);
    resampleMass->Draw("PE");
    mMassUnlike->SetLineColor(kGreen);
    mMassUnlike->Draw("PE;same");

    resampleMassUniform->SetLineColor(kRed);
    resampleMassUniform->Draw("PE;same");
  
    resampleMass->GetXaxis()->SetTitle("M_{ee} (GeV)");
    resampleMass->GetYaxis()->SetTitle("Normalized Counts");

    TLegend * legmass = new TLegend(.8, .7, 1, .9);
    legmass->SetHeader("Legend","C"); // option "C" allows to center the header
    legmass->AddEntry(resampleMassUniform, "Mixed Event Uniform");
    legmass->AddEntry(resampleMass, "Mixed Event slight tracks");
    legmass->AddEntry(mMassUnlike, "Data");
    legmass->Draw("same");
    gPad->SetLogy();
    resampleMass->SetStats(false);
    gPad->Print("note_plots/mixed_event_plots/mass_dists.png");

    
    /*makeCanvas();
    resamplePt->Scale(1/resamplePt->Integral());
    resamplePtUniform->Scale(1/resamplePtUniform->Integral());
    mPtUnlike->Scale(1/mPtUnlike->Integral());

    resamplePt->SetLineColor(kBlack);
    resamplePt->Draw("PE");
    mPtUnlike->SetLineColor(kGreen);
    mPtUnlike->Draw("PE;same");

    resamplePtUniform->SetLineColor(kRed);
    resamplePtUniform->Draw("PE;same");
  
    resamplePt->GetXaxis()->SetTitle("pT (GeV)");
    resamplePt->GetYaxis()->SetTitle("Normalized Counts");

    TLegend * legpt = new TLegend(.8, .7, 1, .9);
    legpt->SetHeader("Legend","C"); // option "C" allows to center the header
    legpt->AddEntry(resamplePtUniform, "Mixed Event Uniform");
    legpt->AddEntry(resamplePt, "Mixed Event slight tracks");
    legpt->AddEntry(mPtUnlike, "Data");
    legpt->Draw("same");
    resamplePt->SetStats(false);
    gPad->SetLogy();
    gPad->Print("note_plots/mixed_event_plots/pair_pts.png");*/

    makeCanvas();
    cos2phivMassUnlike->Draw("colz");
    cos2phivMassUnlike->GetXaxis()->SetTitle("A_{2#phi}");
    cos2phivMassUnlike->GetYaxis()->SetTitle("M_{ee}");
    cos2phivMassUnlike->SetTitle("A_{2#phi} vs pair mass");

    makeCanvas();
    cos4phivMassUnlike->Draw("colz");
    cos4phivMassUnlike->GetXaxis()->SetTitle("A_{4#phi}");
    cos4phivMassUnlike->GetYaxis()->SetTitle("M_{ee}");
    cos4phivMassUnlike->SetTitle("A_{4#phi} vs pair mass");

    makeCanvas();
    acoResample->SetLineColor(kBlack);
    acoResample->SetTitle("#alpha, 0.4 < M < 1");
    acoResample->Scale(1/acoResample->GetEntries());
    for(int i = 0; i < acoResample->GetNbinsX(); ++i){
        double binwidth = acoResample->GetBinWidth(i);
        std::cout << "bin width " << binwidth; 
        std::cout << "bin content: " << acoResample->GetBinContent(i);
        double content = acoResample->GetBinContent(i)/binwidth;
        acoResample->SetBinContent(i, content);
    }
    acoResample->GetXaxis()->SetTitle("#alpha");
    acoResample->GetYaxis()->SetTitle("1/N * dN/d#alpha");
    acoResample->Draw("PE");
    /*MCaco->SetLineColor(kBlue);
    rcaco->SetLineColor(kRed);
    MCaco->Draw("PE;same");
    rcaco->Draw("PE;same");
    auto acolegend = new TLegend(0.75,0.6,1,0.85);
    acolegend->SetHeader("Legend","C"); // option "C" allows to center the header
    acolegend->AddEntry(aco,"Run 12");
    acolegend->AddEntry(MCaco,"Simulation MC");
    acolegend->AddEntry(rcaco,"Simulation Reco");*/



    TFile file("output_root_files/background_plots.root", "RECREATE");
    purity->Write();
    m2phiMomentsResampleUniform->Write();
    m4phiMomentsResampleUniform->Write();

    m2phiMomentsResample->Write();
    m4phiMomentsResample->Write();
    acoResample->Write();





}
