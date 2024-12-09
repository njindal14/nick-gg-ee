#include <iostream>
#include "FemtoPairFormat.h"
double const pi = 3.14159265358979323;
double c = 3.0e1; //in cm/ns
double me2 = pow(0.00051,2);

int ican2 = 0;
void makeCanvas()  {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican2++ ), "", 900, 600);
    can->SetTopMargin(0.08);
    can->SetRightMargin(0.15);
    can->SetBottomMargin(0.15);
}

void crossSectionsUU() {

    TFile *myFile = TFile::Open("/Users/Nick/STAR/breit-wheeler/rootFiles/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");
    TLorentzVector lv1, lv2, lv, lvn;
    //gSystem->Load("FemtoPairFormat_h.so");

    auto * mMass = new TH1F("mMass", "Pair Mass", 20, 0.3, 2.7);
    auto * mPt = new TH1F("mPt", "Pair pT", 20, 0, 0.1);
    auto * mPt2 = new TH1F("mPt2", "Pair pT^{2}", 20, 0, 0.008);
    auto * mY = new TH1F("mY", "Pair Rapidity", 20, -1, 1);

    //for looking at uranium cross sections using gold corrections
    auto * mMassAu = new TH1F("", "Pair Mass Au Corrections", 20, 0.3, 2.7);
    auto * mPtAu = new TH1F("", "Pair pT Au Corrections", 20, 0, 0.1);
    auto * mPt2Au = new TH1F("", "Pair pT^{2} Au Corrections", 20, 0, 0.008);
    auto * mYAu = new TH1F("", "Pair Rapidity Au Corrections", 20, -1, 1);

    //tofmatch efficiencies
    TFile * tofmatchEfficiencies = new TFile("/Users/Nick/STAR/breit-wheeler/rootFiles/TOF_efficiency_y2017_TpcRS_NOMINAL.root");
    TH1D * tofMatch_eff_mass = (TH1D*)tofmatchEfficiencies->Get("tof_eff_mass");
    TH1D * tofMatch_eff_pt = (TH1D*)tofmatchEfficiencies->Get("tof_eff_pt");
    TH1D * tofMatch_eff_pt2 = (TH1D*)tofmatchEfficiencies->Get("tof_eff_pt2");
    TH1D * tofMatch_eff_y = (TH1D*)tofmatchEfficiencies->Get("tof_eff_y");


    //embedding reco efficiencies 
    TFile * simulationplots = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/simulation_plots_new.root");
    TEfficiency * mass_eff = (TEfficiency*)simulationplots->Get("mEff");
    TEfficiency * pt_eff = (TEfficiency*)simulationplots->Get("hEff");
    TEfficiency * pt2_eff = (TEfficiency*)simulationplots->Get("h2Eff");
    TEfficiency * y_eff = (TEfficiency*)simulationplots->Get("YEff");

    TEfficiency * ptm_eff = (TEfficiency*)simulationplots->Get("PtMEff");

    

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

        if(pair->d1_mNHitsFit < 10 || pair->d2_mNHitsFit < 10) continue;
        if(pair->d1_mNHitsDedx < 15 || pair->d2_mNHitsDedx < 15) continue;


        if( fabs(pair->mVertexZ) < 100 &&  pair->mGRefMult <= 4 && pair->mChargeSum == 0 && pair->d1_mDCA < 1 && pair->d2_mDCA < 1 && 
        pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag!=0 && fabs(ddTofVal) < 0.4 && ddTofVal !=0 && chiee < 10 && 3*chiee < chipipi && pair->mZDCEast < 500 && pair->mZDCWest < 500) {
            
            if(lv.M() > 0.4 && lv.M() < 0.76){
                mPt->Fill(lv.Pt());
                mPt2->Fill(pow(lv.Pt(),2)); 
                mY->Fill(lv.Rapidity());
                
                mPtAu->Fill(lv.Pt());
                mPt2Au->Fill(pow(lv.Pt(),2)); 
                mYAu->Fill(lv.Rapidity());    

            }   
            if(lv.Pt() < 0.1  && abs(lv.Rapidity()) < 1){ 
                mMass->Fill(lv.M()); mMassAu->Fill(lv.M());
            }

        }   
    }

    //scale the histograms using tofmatch efficiencies and the reco efficiencies

    for(int ix = 0; ix < mMass->GetNbinsX(); ix++){
        double binContent = mMass->GetBinContent(ix);
        double binError = mMass->GetBinError(ix);
        double binCenter = mMass->GetBinCenter(ix);
        double tofEffBin = tofMatch_eff_mass->FindBin(binCenter);
        double tof_eff = tofMatch_eff_mass->GetBinContent(tofEffBin);
        double recoEffBin = mass_eff->FindFixBin(binCenter);
        double reco_eff = mass_eff->GetEfficiency(recoEffBin);
        if(tof_eff == 0) tof_eff = tofMatch_eff_mass->GetBinContent(5);
        if(reco_eff == 0) reco_eff = mass_eff->GetEfficiency(20);
        mMass->SetBinContent(ix, binContent/(tof_eff*reco_eff));
        mMass->SetBinError(ix, binError/(tof_eff*reco_eff));
        mMassAu->SetBinContent(ix, binContent/(tof_eff*reco_eff));
        mMassAu->SetBinError(ix, binError/(tof_eff*reco_eff));
    }

    for(int ix = 0; ix < mPt->GetNbinsX(); ix++){
        double binContent = mPt->GetBinContent(ix);
        double binError = mPt->GetBinError(ix);
        double binCenter = mPt->GetBinCenter(ix);
        double tofEffBin = tofMatch_eff_pt->FindBin(binCenter);
        double tof_eff = tofMatch_eff_pt->GetBinContent(tofEffBin);
        double recoEffBin = pt_eff->FindFixBin(binCenter);
        double reco_eff = pt_eff->GetEfficiency(recoEffBin);
        if(tof_eff == 0) tof_eff = tofMatch_eff_pt->GetBinContent(1);
        if(reco_eff == 0) reco_eff = pt_eff->GetEfficiency(20);
        mPt->SetBinContent(ix, binContent/(tof_eff*reco_eff));
        mPt->SetBinError(ix, binError/(tof_eff*reco_eff));
        mPtAu->SetBinContent(ix, binContent/(tof_eff*reco_eff));
        mPtAu->SetBinError(ix, binError/(tof_eff*reco_eff));
    }

    for(int ix = 0; ix < mPt2->GetNbinsX(); ix++){
        double binContent = mPt2->GetBinContent(ix);
        double binError = mPt2->GetBinError(ix);
        double binCenter = mPt2->GetBinCenter(ix);
        double tofEffBin = tofMatch_eff_pt2->FindBin(binCenter);
        double tof_eff = tofMatch_eff_pt2->GetBinContent(tofEffBin);
        double recoEffBin = pt2_eff->FindFixBin(binCenter);
        double reco_eff = pt2_eff->GetEfficiency(recoEffBin);
        if(tof_eff == 0) tof_eff = tofMatch_eff_pt2->GetBinContent(1);
        if(reco_eff == 0) reco_eff = pt_eff->GetEfficiency(20);

        mPt2->SetBinContent(ix, binContent/(tof_eff*reco_eff));
        mPt2->SetBinError(ix, binError/(tof_eff*reco_eff));
        mPt2Au->SetBinContent(ix, binContent/(tof_eff*reco_eff));
        mPt2Au->SetBinError(ix, binError/(tof_eff*reco_eff));
    }

    for(int ix = 0; ix < mY->GetNbinsX(); ix++){
        double binContent = mY->GetBinContent(ix);
        double binError = mY->GetBinError(ix);
        double binCenter = mY->GetBinCenter(ix);
        double tofEffBin = tofMatch_eff_y->FindBin(binCenter);
        double tof_eff = tofMatch_eff_y->GetBinContent(tofEffBin);
        double recoEffBin = y_eff->FindFixBin(binCenter);
        double reco_eff = y_eff->GetEfficiency(recoEffBin);
        if(tof_eff == 0) tof_eff = tofMatch_eff_y->GetBinContent(5);
        if(reco_eff == 0) reco_eff = y_eff->GetEfficiency(20);
        mY->SetBinContent(ix, binContent/(tof_eff*reco_eff));
        mY->SetBinError(ix, binError/(tof_eff*reco_eff));
        
        mYAu->SetBinContent(ix, binContent/(tof_eff*reco_eff));
        mYAu->SetBinError(ix, binError/(tof_eff*reco_eff));
    }

    
    //global correction factors (not bin-by-bin)
    double luminosity = 274170;
    double lumi_fraction = .793; //confirmed for my analysis

    double lumi_fraction_au = .696;

    double bbc_eff = 0.683; //taken from JDB analysis
    double vertex_eff = 0.68; //taken from JDB analysis


    double purity_corrections = 0.975*0.996; //subject to change
    double tpc_eff = 1*1; //100% for each track -- no missing sectors -- this may change after looking at simulation
    double XnXn_correction = 0.552; //still an estimate, need to recalculate. got it from integral of zdcs less than 500 on both sides over the total integral

    double total_eff = luminosity*lumi_fraction*bbc_eff*purity_corrections*vertex_eff*XnXn_correction;

    double total_eff_with_au = luminosity*lumi_fraction_au*bbc_eff*purity_corrections*vertex_eff*XnXn_correction;

    //scale and draw cross sections
    mMass->Scale(1/(total_eff*mMass->GetBinWidth(1)));
    mPt->Scale(1/(total_eff*mPt->GetBinWidth(1)));
    mPt2->Scale(1/(total_eff*mPt2->GetBinWidth(1)));
    mY->Scale(1/(total_eff*mY->GetBinWidth(1)));

    //scale and draw cross sections
    mMassAu->Scale(1/(total_eff_with_au*mMass->GetBinWidth(1)));
    mPtAu->Scale(1/(total_eff_with_au*mPt->GetBinWidth(1)));
    mPt2Au->Scale(1/(total_eff_with_au*mPt2->GetBinWidth(1)));
    mYAu->Scale(1/(total_eff_with_au*mY->GetBinWidth(1)));

    makeCanvas();
    gPad->SetLogy();
    mMass->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
    mMass->GetYaxis()->SetTitle("#frac{d#sigma (#gamma#gamma --> e^{+}e^{-})}{dM} (mb/(GeV/c^{2}))");
    mMass->Draw("PE");
    //mMassAu->SetLineColor(kRed);
    //mMassAu->Draw("PE;same");
    gPad->Print("note_plots/results_plots/UU_MassXSecBothNEW.png");

    makeCanvas();
    mPt->GetXaxis()->SetTitle("pT_{ee} (GeV/c)");
    mPt->GetYaxis()->SetTitle("#frac{d#sigma (#gamma#gamma --> e^{+}e^{-})}{dpT} (mb/(GeV/c))");
    mPt->Draw("PE");
    //mPtAu->SetLineColor(kRed);
    //mPtAu->Draw("PE;same");
    gPad->Print("note_plots/results_plots/UU_PtXSecBothNEW.png");

    makeCanvas();
    gPad->SetLogy();
    mPt2->GetXaxis()->SetTitle("pT_{ee}^{2} (GeV/c)^{2}");
    mPt2->GetYaxis()->SetTitle("#frac{d#sigma (#gamma#gamma --> e^{+}e^{-})}{dpT^{2}} (mb/(GeV/c)^{2})");
    mPt2->Draw("PE");
   // mPt2Au->SetLineColor(kRed);
    //mPt2Au->Draw("PE;same");
    gPad->Print("note_plots/results_plots/UU_Pt2XSecBothNEW.png");


    makeCanvas();
    mY->GetXaxis()->SetTitle("y_{ee}");
    mY->GetYaxis()->SetTitle("#frac{d#sigma (#gamma#gamma --> e^{+}e^{-})}{dy} (mb)");
    gPad->SetLogy();
    mY->Draw("PE");
    //mYAu->SetLineColor(kRed);
    //mYAu->Draw("PE;same");
    gPad->Print("note_plots/results_plots/UU_YXSecBothNEW.png");

   
    // calculate and plot the absolute cross sections
    
    Double_t error_m;
    Double_t error_pt;
    Double_t error_pt2;
    Double_t error_y;
    std::cout << "mass total section : " << mMass->IntegralAndError(1, 20, error_m, "width") << " +- " << error_m << "\n";
    std::cout << "pT cross section integral: " << mPt->IntegralAndError(1, 20, error_pt, "width") << " +- " << error_pt << "\n";
    std::cout << "pT2 cross section integral: " << mPt2->IntegralAndError(1, 20, error_pt2, "width") << " +- " << error_pt2 << "\n";
    std::cout << "Y cross section integral: " << mY->IntegralAndError(1, 20, error_y, "width") << " +- " << error_y << "\n";
    //Int_t n = 4;
    Double_t x[4]  = {1, 2, 3, 4};
    Double_t y[4]  = {mMass->IntegralAndError(1, 20, error_m, "width"), mPt->IntegralAndError(1, 20, error_pt, "width"), mPt2->IntegralAndError(1, 20, error_pt2, "width"), mY->IntegralAndError(1, 20, error_y, "width")};
    Double_t ex[4] = {0.5, 0.5, 0.5, 0.5};
    Double_t ey[4] = {error_m, error_pt, error_pt2, error_y};
    makeCanvas();
    TGraph * gr = new TGraphErrors(4,x,y,ex,ey);
    gr->SetTitle("Absolute Cross Sections");
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->GetXaxis()->SetBinLabel(20,"M_{ee} (GeV/c^{2})");
    gr->GetXaxis()->SetBinLabel(40, "pT (GeV/c)");
    gr->GetXaxis()->SetBinLabel(62, "pT^{2} (GeV/c)^{2}");
    gr->GetXaxis()->SetBinLabel(84, "y_{ee} [1]");  
    gr->GetYaxis()->SetTitle("Absolute Cross Section (mb)");  
    gr->Draw("AP");
    gPad->Print("note_plots/results_plots/UUAbsXSectionsNEW.png");


    //write to root file
    TFile file("output_root_files/crossSectionsUU_NEW.root", "RECREATE");
    mMass->Write();
    mPt->Write();
    mPt2->Write();
    mY->Write();


}