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

void crossSectionsAuAu() {

    auto * mMassAuAu = new TH1F("mMassAuAu", "Pair Mass", 20, 0.3, 2.7);
    auto * mPtAuAu = new TH1F("mPtAuAu", "Pair pT", 20, 0, 0.1);
    auto * mPt2AuAu = new TH1F("mPt2AuAu", "Pair pT^{2}", 20, 0, 0.01);
    auto * mYAuAu = new TH1F("mYAuAu", "Pair Rapidity", 20, -1, 1);

    TChain * ch = new TChain("PairDst");
    ch->Add("/Users/Nick/STAR/breit-wheeler/rootFiles/slim_pair_dst_Run10AuAu.root");
    ch->Add("/Users/Nick/STAR/breit-wheeler/rootFiles/slim_pair_dst_Run11AuAu.root");
    TTreeReader myReader(ch);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

    TLorentzVector lv1, lv2, lv, lvn;

    //open efficiency root files  get the 1d histograms for relevant efficiencies
    TFile * trackingEfficiencies = new TFile("/Users/Nick/STAR/breit-wheeler/rootFiles/output_eff_3D_NOMINAL.root");
    TH3D * masspty = (TH3D*)trackingEfficiencies->Get("mMass");  
    TH3D * mc_masspty = (TH3D*)trackingEfficiencies->Get("mc_mMass");
    TH1D * mass = masspty->ProjectionZ();
    TH1D * mc_mass = mc_masspty->ProjectionZ();
    TH1D * pt = masspty->ProjectionY();
    TH1D * mc_pt = mc_masspty->ProjectionY();
    TH1D * y = masspty->ProjectionX();
    TH1D * mc_y = mc_masspty->ProjectionX();

    TH1D * mass_eff = (TH1D*)mass->Clone();
    mass_eff->Divide(mc_mass);
    TH1D * pt_eff = (TH1D*)pt->Clone();
    pt_eff->Divide(mc_pt);
    TH1D * y_eff = (TH1D*)y->Clone();
    y_eff->Divide(mc_y);

    TFile * tofmatchEfficiencies = new TFile("/Users/Nick/STAR/breit-wheeler/rootFiles/TOF_efficiency_y2017_TpcRS_NOMINAL.root");
    TH1D * tofMatch_eff_mass = (TH1D*)tofmatchEfficiencies->Get("tof_eff_mass");
    TH1D * tofMatch_eff_pt = (TH1D*)tofmatchEfficiencies->Get("tof_eff_pt");
    TH1D * tofMatch_eff_pt2 = (TH1D*)tofmatchEfficiencies->Get("tof_eff_pt2");
    TH1D * tofMatch_eff_y = (TH1D*)tofmatchEfficiencies->Get("tof_eff_y");

    

    while (myReader.Next()) {
        //values we will want to use for PID cuts
        double chiee = pow( pair->d1_mNSigmaElectron, 2 ) + pow( pair->d2_mNSigmaElectron, 2 );
        double chipipi = pow( pair -> d1_mNSigmaPion, 2) + pow( pair -> d2_mNSigmaPion, 2);
        double c = 3.0e1; //in cm/ns
    
        //Lorentz vectors for each pair track and get lorentz sum and diff
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
        pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag!=0 && fabs(ddTofVal) < 0.4 && ddTofVal !=0 && chiee < 10 && 3*chiee < chipipi) {
            
            if(lv.M() > 0.4 && lv.M() < 0.76 && fabs(lv.Rapidity()) < 1){
                mPtAuAu->Fill(lv.Pt());
                mPt2AuAu->Fill(pow(lv.Pt(),2)); 
                mYAuAu->Fill(lv.Rapidity());

            }   
            if(lv.Pt() < 0.1 && fabs(lv.Rapidity()) < 1) mMassAuAu->Fill(lv.M());
        }

    }

     //scale the histograms using tofmatch efficiencies and the reco efficiencies

    for(int ix = 0; ix < mMassAuAu->GetNbinsX(); ix++){
        double binContent = mMassAuAu->GetBinContent(ix);
        double binError = mMassAuAu->GetBinError(ix);
        double binCenter = mMassAuAu->GetBinCenter(ix);
        double tofEffBin = tofMatch_eff_mass->FindBin(binCenter);
        double tof_eff = tofMatch_eff_mass->GetBinContent(tofEffBin);
        double recoEffBin = mass_eff->FindBin(binCenter);
        double reco_eff = mass_eff->GetBinContent(recoEffBin);
        if(tof_eff == 0) tof_eff = tofMatch_eff_mass->GetBinContent(5);
        if(reco_eff == 0) reco_eff = mass_eff->GetBinContent(20);
        mMassAuAu->SetBinContent(ix, binContent/(tof_eff*reco_eff));
        mMassAuAu->SetBinError(ix, binError/(tof_eff*reco_eff));
    }

    for(int ix = 0; ix < mPtAuAu->GetNbinsX(); ix++){
        double binContent = mPtAuAu->GetBinContent(ix);
        double binError = mPtAuAu->GetBinError(ix);
        double binCenter = mPtAuAu->GetBinCenter(ix);
        double tofEffBin = tofMatch_eff_pt->FindBin(binCenter);
        double tof_eff = tofMatch_eff_pt->GetBinContent(tofEffBin);
        double recoEffBin = pt_eff->FindBin(binCenter);
        double reco_eff = pt_eff->GetBinContent(recoEffBin);
        if(tof_eff == 0) tof_eff = tofMatch_eff_pt->GetBinContent(1);
        if(reco_eff == 0) reco_eff = pt_eff->GetBinContent(20);
        mPtAuAu->SetBinContent(ix, binContent/(tof_eff*reco_eff));
        mPtAuAu->SetBinError(ix, binError/(tof_eff*reco_eff));

    }

    for(int ix = 0; ix < mPt2AuAu->GetNbinsX(); ix++){
        double binContent = mPt2AuAu->GetBinContent(ix);
        double binError = mPt2AuAu->GetBinError(ix);
        double binCenter = mPt2AuAu->GetBinCenter(ix);
        double tofEffBin = tofMatch_eff_pt2->FindBin(binCenter);
        double tof_eff = tofMatch_eff_pt2->GetBinContent(tofEffBin);
        double recoEffBin = pt_eff->FindBin(binCenter);
        double reco_eff = pt_eff->GetBinContent(recoEffBin);
        if(tof_eff == 0) tof_eff = tofMatch_eff_pt2->GetBinContent(1);
        if(reco_eff == 0) reco_eff = pt_eff->GetBinContent(20);

        mPt2AuAu->SetBinContent(ix, binContent/(tof_eff*reco_eff));
        mPt2AuAu->SetBinError(ix, binError/(tof_eff*reco_eff));

    }

    for(int ix = 0; ix < mYAuAu->GetNbinsX(); ix++){
        double binContent = mYAuAu->GetBinContent(ix);
        double binError = mYAuAu->GetBinError(ix);
        double binCenter = mYAuAu->GetBinCenter(ix);
        double tofEffBin = tofMatch_eff_y->FindBin(binCenter);
        double tof_eff = tofMatch_eff_y->GetBinContent(tofEffBin);
        double recoEffBin = y_eff->FindBin(binCenter);
        double reco_eff = y_eff->GetBinContent(recoEffBin);
        if(tof_eff == 0) tof_eff = tofMatch_eff_y->GetBinContent(5);
        if(reco_eff == 0) reco_eff = y_eff->GetBinContent(20);
        mYAuAu->SetBinContent(ix, binContent/(tof_eff*reco_eff));
        mYAuAu->SetBinError(ix, binError/(tof_eff*reco_eff));

    }

    
    //global correction factors (not bin-by-bin)
    double luminosity = 543000;
    double lumi_fraction = .696; //taken from jdb
    double bbc_eff = 0.683; //taken from jdb
    double vertex_eff = 0.68; //taken from JDB analysis, same as eEvent
    double purity_corrections = 0.975*0.996; //subject to change
    //double tpc_eff = .8*.8; //100% for each track -- no missing sectors -- this may change after looking at simulation
    double XnXn_correction = 1/2.43; //still an estimate, need to recalculate

    double total_eff = luminosity*lumi_fraction*bbc_eff*purity_corrections*vertex_eff*XnXn_correction;


    //scale and draw cross sections
    mMassAuAu->Scale(1/(total_eff*mMassAuAu->GetBinWidth(1)));
    mPtAuAu->Scale(1/(total_eff*mPtAuAu->GetBinWidth(1)));
    mPt2AuAu->Scale(1/(total_eff*mPt2AuAu->GetBinWidth(1)));
    mYAuAu->Scale(1/(total_eff*mYAuAu->GetBinWidth(1)));

    makeCanvas();
    gPad->SetLogy();
    mMassAuAu->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
    mMassAuAu->GetYaxis()->SetTitle("#frac{d#sigma (#gamma#gamma --> e^{+}e^{-})}{dM} (b/(GeV/c^2))");
    mMassAuAu->Draw("PE");
    gPad->Print("note_plots/results_plots/Au_MassXSec.png");

    makeCanvas();
    mPtAuAu->GetXaxis()->SetTitle("pT_{ee} (GeV/c)");
    mPtAuAu->GetYaxis()->SetTitle("#frac{d#sigma (#gamma#gamma --> e^{+}e^{-})}{dpT} (b/(GeV/c))");
    mPtAuAu->Draw("PE");
    gPad->Print("note_plots/results_plots/Au_PtXSec.png");

    makeCanvas();
    gPad->SetLogy();
    mPt2AuAu->GetXaxis()->SetTitle("pT_{ee}^{2} (GeV/c)^{2}");
    mPt2AuAu->GetYaxis()->SetTitle("#frac{d#sigma (#gamma#gamma --> e^{+}e^{-})}{dpT^{2}} (b/(GeV/c)^2)");
    mPt2AuAu->Draw("PE");
    gPad->Print("note_plots/results_plots/Au_Pt2XSec.png");


    makeCanvas();
    mYAuAu->GetXaxis()->SetTitle("y_{ee}");
    mYAuAu->GetYaxis()->SetTitle("#frac{d#sigma (#gamma#gamma --> e^{+}e^{-})}{dy} (b)");
    gPad->SetLogy();
    mYAuAu->Draw("PE");
    gPad->Print("note_plots/results_plots/Au_YXSec.png");

    //write to root file
    TFile file("output_root_files/crossSectionsAuAu.root", "RECREATE");
    mMassAuAu->Write();
    mPtAuAu->Write();
    mPt2AuAu->Write();
    mYAuAu->Write();


}