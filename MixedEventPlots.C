#include <iostream>
#include "PlottingFunctions.h"

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


void MixedEventPlots(){
    

    TFile * mixedeventplots = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/MixedEventplots_Au.root");

    TH3F * mULSULSCos1phivsPTvsMass = (TH3F*)mixedeventplots->Get("mULSULSCos1phivsPTvsMass");
    TH3F * mULSULSCos2phivsPTvsMass = (TH3F*)mixedeventplots->Get("mULSULSCos2phivsPTvsMass");    
    TH3F * mULSULSCos3phivsPTvsMass = (TH3F*)mixedeventplots->Get("mULSULSCos3phivsPTvsMass");    
    TH3F * mULSULSCos4phivsPTvsMass = (TH3F*)mixedeventplots->Get("mULSULSCos4phivsPTvsMass");

    TH3F * mULSLSCos1phivsPTvsMass = (TH3F*)mixedeventplots->Get("mULSLSCos1phivsPTvsMass");
    TH3F * mULSLSCos2phivsPTvsMass = (TH3F*)mixedeventplots->Get("mULSLSCos2phivsPTvsMass");    
    TH3F * mULSLSCos3phivsPTvsMass = (TH3F*)mixedeventplots->Get("mULSLSCos3phivsPTvsMass");    
    TH3F * mULSLSCos4phivsPTvsMass = (TH3F*)mixedeventplots->Get("mULSLSCos4phivsPTvsMass");
    

    TH1F * acoULSULS = (TH1F*)mixedeventplots->Get("acoULSULS");
    
    auto *cos1phimoments = (TH1F*)ProjectTH3ToProfile(mULSULSCos1phivsPTvsMass, "yx", "ProfileY", 0, 0, 0, 0, 0.4, 0.76);
    auto *cos2phimoments = (TH1F*)ProjectTH3ToProfile(mULSULSCos2phivsPTvsMass, "yx", "ProfileY", 0, 0, 0, 0, 0.4, 0.76);
    auto *cos3phimoments = (TH1F*)ProjectTH3ToProfile(mULSULSCos3phivsPTvsMass, "yx", "ProfileY", 0, 0, 0, 0, 0.4, 0.76);
    auto *cos4phimoments = (TH1F*)ProjectTH3ToProfile(mULSULSCos4phivsPTvsMass, "yx", "ProfileY", 0, 0, 0, 0, 0.4, 0.76);

    auto *cos1phimomentsMass = (TH1F*)ProjectTH3ToProfile(mULSULSCos1phivsPTvsMass, "zx", "ProfileY", 0, 0, 0, 0.15, 0, 0);
    auto *cos2phimomentsMass = (TH1F*)ProjectTH3ToProfile(mULSULSCos2phivsPTvsMass, "zx", "ProfileY", 0, 0, 0, 0.15, 0, 0);
    auto *cos3phimomentsMass = (TH1F*)ProjectTH3ToProfile(mULSULSCos3phivsPTvsMass, "zx", "ProfileY", 0, 0, 0, 0.15, 0, 0);
    auto *cos4phimomentsMass = (TH1F*)ProjectTH3ToProfile(mULSULSCos4phivsPTvsMass, "zx", "ProfileY", 0, 0, 0, 0.15, 0, 0);


    auto *cos1phimomentsULSLS = (TH1F*)ProjectTH3ToProfile(mULSLSCos1phivsPTvsMass, "yx", "ProfileY", 0, 0, 0, 0, 0.4, 0.76);
    auto *cos2phimomentsULSLS = (TH1F*)ProjectTH3ToProfile(mULSLSCos2phivsPTvsMass, "yx", "ProfileY", 0, 0, 0, 0, 0.4, 0.76);
    auto *cos3phimomentsULSLS = (TH1F*)ProjectTH3ToProfile(mULSLSCos3phivsPTvsMass, "yx", "ProfileY", 0, 0, 0, 0, 0.4, 0.76);
    auto *cos4phimomentsULSLS = (TH1F*)ProjectTH3ToProfile(mULSLSCos4phivsPTvsMass, "yx", "ProfileY", 0, 0, 0, 0, 0.4, 0.76);


    makeCanvas();
    cos1phimoments->GetXaxis()->SetTitle("pT (GeV/c)");
    cos1phimoments->GetYaxis()->SetTitle("A_{1#phi}");
    cos1phimoments->Draw();
    cos1phimomentsULSLS->SetLineColor(kRed);
    cos1phimomentsULSLS->Draw("PE;same");

    makeCanvas();
    cos2phimoments->GetXaxis()->SetTitle("pT (GeV/c)");
    cos2phimoments->GetYaxis()->SetTitle("A_{2#phi}");
    cos2phimoments->Draw(); 
    cos2phimomentsULSLS->SetLineColor(kRed);
    cos2phimomentsULSLS->Draw("PE;same");
    
    makeCanvas();
    cos3phimoments->GetXaxis()->SetTitle("pT (GeV/c)");
    cos3phimoments->GetYaxis()->SetTitle("A_{3#phi}");
    cos3phimoments->Draw();
    cos3phimomentsULSLS->SetLineColor(kRed);
    cos3phimomentsULSLS->Draw("PE;same");

    makeCanvas();
    cos4phimoments->GetXaxis()->SetTitle("pT (GeV/c)");
    cos4phimoments->GetYaxis()->SetTitle("A_{4#phi}");
    cos4phimoments->Draw();
    cos4phimomentsULSLS->SetLineColor(kRed);
    cos4phimomentsULSLS->Draw("PE;same");



    makeCanvas();
    cos1phimomentsMass->GetXaxis()->SetTitle("M (GeV/c^{2})");
    cos1phimomentsMass->GetYaxis()->SetTitle("A_{1#phi}");
    cos1phimomentsMass->Draw();

    makeCanvas();
    cos2phimomentsMass->GetXaxis()->SetTitle("M (GeV/c^{2})");
    cos2phimomentsMass->GetYaxis()->SetTitle("A_{2#phi}");
    cos2phimomentsMass->Draw(); 

    makeCanvas();
    cos3phimomentsMass->GetXaxis()->SetTitle("M (GeV/c^{2})");
    cos3phimomentsMass->GetYaxis()->SetTitle("A_{3#phi}");
    cos3phimomentsMass->Draw();

    makeCanvas();
    cos4phimomentsMass->GetXaxis()->SetTitle("M (GeV/c^{2})");
    cos4phimomentsMass->GetYaxis()->SetTitle("A_{4#phi}");
    cos4phimomentsMass->Draw();

    makeCanvas();
    acoULSULS->SetTitle("#alpha, 0.4 < M < 1");
    acoULSULS->Scale(1/acoULSULS->GetEntries());
    for(int i = 0; i < acoULSULS->GetNbinsX(); ++i){
        double binwidth = acoULSULS->GetBinWidth(i);
        double content = acoULSULS->GetBinContent(i)/binwidth;
        acoULSULS->SetBinContent(i, content);
    }
    acoULSULS->GetXaxis()->SetTitle("#alpha");
    acoULSULS->GetYaxis()->SetTitle("1/N * dN/d#alpha");
    acoULSULS->Draw("PE");

    TFile * fo = new TFile( "output_root_files/MixedEventMoments_Au.root", "RECREATE" );
    cos1phimoments->Write();
    cos2phimoments->Write();
    cos3phimoments->Write();
    cos4phimoments->Write();

    cos1phimomentsULSLS->Write();
    cos2phimomentsULSLS->Write();
    cos3phimomentsULSLS->Write();
    cos4phimomentsULSLS->Write();



}