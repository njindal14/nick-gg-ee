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
}

void crossSectionsUU() {

    TFile *myFile = TFile::Open("/Users/Nick/STAR/breit-wheeler/rootFiles/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");
    TLorentzVector lv1, lv2, lv, lvn;
    //gSystem->Load("FemtoPairFormat_h.so");

    auto * mMass = new TH1F("", "Pair Mass", 20, 0.3, 2.7);
    auto * mPt = new TH1F("", "Pair pT", 20, 0, 0.1);
    auto * mPt2 = new TH1F("", "Pair pT^{2}", 20, 0, 0.01);
    auto * mY = new TH1F("", "Pair Rapidity", 20, -1, 1);

    //tofmatch efficiencies
    TFile * tofmatchEfficiencies = new TFile("/Users/Nick/STAR/breit-wheeler/rootFiles/TOF_efficiency_y2017_TpcRS_NOMINAL.root");
    TH1D * tofMatch_eff_mass = (TH1D*)tofmatchEfficiencies->Get("tof_eff_mass");
    TH1D * tofMatch_eff_pt = (TH1D*)tofmatchEfficiencies->Get("tof_eff_pt");
    TH1D * tofMatch_eff_pt2 = (TH1D*)tofmatchEfficiencies->Get("tof_eff_pt2");
    TH1D * tofMatch_eff_y = (TH1D*)tofmatchEfficiencies->Get("tof_eff_y");

    //embedding reco efficiencies 
    TFile * simulationplots = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/simulation_plots.root");
    TEfficiency * mass_eff = (TEfficiency*)simulationplots->Get("mEff");
    TEfficiency * pt_eff = (TEfficiency*)simulationplots->Get("hEff");
    TEfficiency * pt2_eff = (TEfficiency*)simulationplots->Get("h2Eff");
    TEfficiency * y_eff = (TEfficiency*)simulationplots->Get("YEff");

    while(myReader.Next()){

        double chiee = pow( pair->d1_mNSigmaElectron, 2 ) + pow( pair->d2_mNSigmaElectron, 2 );
        double chipipi = pow( pair -> d1_mNSigmaPion, 2) + pow( pair -> d2_mNSigmaPion, 2);
        int chargesumval = pair->mChargeSum;

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

        if( fabs(pair->mVertexZ) < 100 &&  pair->mGRefMult <= 4 && chargesumval == 0 && pair->d1_mDCA < 1 && pair->d2_mDCA < 1 && 
        pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag!=0 && fabs(ddTofVal < 0.4) && ddTofVal !=0 && chiee < 10 && 3*chiee < chipipi) {
            
            if(lv.M() > 0.4 && lv.M() < 0.8){
                mPt->Fill(lv.Pt());
                mPt2->Fill(pow(lv.Pt(),2)); 
                mY->Fill(lv.Rapidity());
            }   
            if(lv.Pt() < 0.15 && fabs(lv.Eta() < 1 && fabs(lv.Rapidity() < 1))) mMass->Fill(lv.M());

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
        if(reco_eff == 0) reco_eff = pt_eff->GetEfficiency(20);
        mY->SetBinContent(ix, binContent/(tof_eff*reco_eff));
        mY->SetBinError(ix, binError/(tof_eff*reco_eff));
    }

    
    //global correction factors (not bin-by-bin)
    double luminosity = 274170;
    double lumi_fraction = .793; //confirmed for my analysis
    double bbc_eff = 0.683; //taken from JDB analysis
    double purity_corrections = 0.975*0.996; //subject to change
    double tpc_eff = 1*1; //100% for each track -- no missing sectors -- this may change after looking at simulation
    double vertex_eff = 0.68; //taken from JDB analysis
    double XnXn_correction = 1/1.1; //still an estimate, need to recalculate

    double total_eff = luminosity*lumi_fraction*bbc_eff*purity_corrections*vertex_eff*XnXn_correction;

    //scale and draw cross sections
    mMass->Scale(1/(total_eff*mMass->GetBinWidth(1)));
    mPt->Scale(1/(total_eff*mPt->GetBinWidth(1)));
    mPt2->Scale(1/(total_eff*mPt2->GetBinWidth(1)));
    mY->Scale(1/(total_eff*mY->GetBinWidth(1)));

    makeCanvas();
    gPad->SetLogy();
    mMass->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
    mMass->GetYaxis()->SetTitle("#frac{d#sigma (#gamma#gamma --> e^{+}e^{-})}{dM} (b)");
    mMass->Draw("PE");
    gPad->Print("note_plots/results_plots/UU_MassXSec.png");

    makeCanvas();
    mPt->GetXaxis()->SetTitle("pT_{ee} (GeV/c)");
    mPt->GetYaxis()->SetTitle("#frac{d#sigma (#gamma#gamma --> e^{+}e^{-})}{dpT} (b)");
    mPt->Draw("PE");
    gPad->Print("note_plots/results_plots/UU_PtXSec.png");

    makeCanvas();
    gPad->SetLogy();
    mPt2->GetXaxis()->SetTitle("pT_{ee}^{2} (GeV/c)^{2}");
    mPt2->GetYaxis()->SetTitle("#frac{d#sigma (#gamma#gamma --> e^{+}e^{-})}{dpT^{2}} (b)");
    mPt2->Draw("PE");
    gPad->Print("note_plots/results_plots/UU_Pt2XSec.png");


    makeCanvas();
    mY->GetXaxis()->SetTitle("y_{ee}");
    mY->GetYaxis()->SetTitle("#frac{d#sigma (#gamma#gamma --> e^{+}e^{-})}{dy} (b)");
    gPad->SetLogy();
    mY->Draw("PE");
    gPad->Print("note_plots/results_plots/UU_YXSec.png");


    //write to root file
    TFile file("output_root_files/crossSectionsUU.root", "RECREATE");
    mMass->Write();
    mPt->Write();
    mPt2->Write();
    mY->Write();


}