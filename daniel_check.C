#include <iostream>
#include "FemtoPairFormat.h"

double const pi = 3.14159265358979323;
double c = 3.0e1; //in cm/ns
double me2 = pow(0.00051, 2);

int ican2 = 0;
void makeCanvas() {
    TCanvas * can = new TCanvas(TString::Format("can%d", ican2++), "", 900, 600);
    can->SetTopMargin(0.1);
    can->SetRightMargin(0.1);
    can->SetBottomMargin(0.15);
}

void daniel_check() {

    const Int_t ptbins = 24;
    Double_t edgespt[ptbins + 1] = {0.0, 0.002, 0.004, .006, 0.008, 0.01, 0.013, 0.017, 0.022, 0.025, 0.03, .035, .04, .045, .05, .055, .06, .065, .07, .075, .08, .085, .09, .095, 0.1};
    const Int_t massbins = 19;
    Double_t edgesmass[massbins + 1] = {0.4, .45,  0.5, .55, 0.6, .65, .7, .75, .8, .85, .9, .95, 1, 1.1, 1.2, 1.4, 1.6, 1.8, 2, 2.5};
    const Int_t ybins = 20;
    Double_t edgesy[ybins + 1] = {-1, -.9, -.8, -.7, -.6, -.5, -.4, -.3, -.2, -.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};

    // Corrected histogram
    auto * PtMY = new TH3F("PtMY", "PtMY_Corrected", 24, edgespt, 19, edgesmass, 20, edgesy);
    PtMY->Sumw2(); 

    // Raw (Uncorrected) histogram
    auto * PtMY_raw = new TH3F("PtMY_raw", "PtMY_Uncorrected", 24, edgespt, 19, edgesmass, 20, edgesy);
    PtMY_raw->Sumw2();

    TChain * ch = new TChain("PairDst");
    ch->Add("/Users/Nick/STAR/breit-wheeler/rootFiles/slim_pair_dst_Run10AuAu.root");
    TTreeReader myReader(ch);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

    TLorentzVector lv1, lv2, lv, lvn;

    // Load efficiencies
    TFile * trackingEfficiencies = new TFile("/Users/Nick/STAR/breit-wheeler/rootFiles/output_eff_3D_NOMINAL.root");
    TH3D * masspty = (TH3D*)trackingEfficiencies->Get("mMass");
    TH3D * mc_masspty = (TH3D*)trackingEfficiencies->Get("mc_mMass");

    TEfficiency * ptmy_eff = new TEfficiency(*masspty, *mc_masspty);

    TFile * tofmatchEfficiencies = new TFile("/Users/Nick/STAR/breit-wheeler/rootFiles/TOF_efficiency_y2017_TpcRS_NOMINAL.root");
    TH1D * tofMatch_eff_pt = (TH1D*)tofmatchEfficiencies->Get("tof_eff_pt");

    // Event Loop
    while (myReader.Next()) {
        double chiee = pow(pair->d1_mNSigmaElectron, 2) + pow(pair->d2_mNSigmaElectron, 2);
        double chipipi = pow(pair->d1_mNSigmaPion, 2) + pow(pair->d2_mNSigmaPion, 2);

        lv1.SetPtEtaPhiM(pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.00051);
        lv2.SetPtEtaPhiM(pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.00051);

        lv = lv1 + lv2;
        lvn = lv1 - lv2;

        Float_t p1_2 = pow(lv1.P(), 2);
        Float_t p2_2 = pow(lv2.P(), 2);
        Float_t dTofVal = pair->d1_mTof - pair->d2_mTof;
        Float_t dTofexpVal = pair->d1_mLength/c * sqrt(1 + me2/p1_2) - pair->d2_mLength/c * sqrt(1 + me2/p2_2);
        Float_t ddTofVal = dTofVal - dTofexpVal;

        if(lv1.Pt() < 0.2 || lv2.Pt() < 0.2) continue;
        if(fabs(lv1.Eta()) > 1 || fabs(lv2.Eta()) > 1 || fabs(lv.Rapidity()) > 1) continue;
        if(pair->d1_mNHitsFit < 20 || pair->d2_mNHitsFit < 20) continue;
        if(pair->d1_mNHitsDedx < 15 || pair->d2_mNHitsDedx < 15) continue;

        if(fabs(pair->mVertexZ) < 100 /*&& pair->mGRefMult <= 4*/ && pair->mChargeSum == 0 && pair->d1_mDCA < 1 && pair->d2_mDCA < 1 &&
           pair->d1_mMatchFlag != 0 && pair->d2_mMatchFlag != 0 && fabs(ddTofVal) < 0.4 && ddTofVal != 0 && chiee < 10 && 3*chiee < chipipi) {

            // Calculate efficiencies for corrected histogram
            int eff_bin = ptmy_eff->FindFixBin(lv.Pt(), lv.M(), lv.Rapidity());
            double reco_eff = ptmy_eff->GetEfficiency(eff_bin);
            double tof_eff_val = tofMatch_eff_pt->GetBinContent(tofMatch_eff_pt->FindBin(lv.Pt()));

            if (reco_eff == 0) reco_eff = 1.0;
            if (tof_eff_val == 0) tof_eff_val = 1.0;

            double weight = 1.0 / (reco_eff * tof_eff_val);
            
            // Fill corrected with weight, raw with default 1.0
            PtMY->Fill(lv.Pt(), lv.M(), lv.Rapidity(), weight);
            PtMY_raw->Fill(lv.Pt(), lv.M(), lv.Rapidity());
        }
    }

    // --- PROJECTIONS ---
    
    // Corrected Projections
    TH1D *mass = PtMY->ProjectionY("mass", 
                                   PtMY->GetXaxis()->FindFixBin(0.0001), 
                                   PtMY->GetXaxis()->FindFixBin(0.0999));
                                   
    TH1D *pt = PtMY->ProjectionX("pt", 
                                 PtMY->GetYaxis()->FindFixBin(0.4001), 
                                 PtMY->GetYaxis()->FindFixBin(0.7599));
                                 
    TH1D * mY = PtMY->ProjectionZ("rapidity", 
                                  PtMY->GetXaxis()->FindFixBin(0.0001), 
                                  PtMY->GetXaxis()->FindFixBin(0.0999),
                                  PtMY->GetYaxis()->FindFixBin(0.4001), 
                                  PtMY->GetYaxis()->FindFixBin(0.7599));

    // Raw Projections
    TH1D *pt_raw = PtMY_raw->ProjectionX("pt_raw", 
                                 PtMY_raw->GetYaxis()->FindFixBin(0.4001), 
                                 PtMY_raw->GetYaxis()->FindFixBin(2.6));

    // Project Y (mass) for pT < 0.1
    // The X axis (pT) in this histogram stops exactly at 0.1, so checking bins 0 to 0.0999 perfectly isolates pT < 0.1.
    TH1D *mass_raw = PtMY_raw->ProjectionY("mass_raw", 
                                 PtMY_raw->GetXaxis()->FindFixBin(0.0001), 
                                 PtMY_raw->GetXaxis()->FindFixBin(0.0999));


    // --- SCALING ---
    
    // Global correction factors (applied ONLY to the corrected histograms)
    double luminosity = 679262; // run 10
    double lumi_fraction = 0.69677721;
    double bbc_eff = 0.683;
    double vertex_eff = 0.633;
    double purity_correction = 0.9917 * 0.975;
    double pid_eff = 0.96 * 0.97; // global TOF & chi2
    double XnXn_correction = 1.0 / 2.46;

    double total_eff = luminosity * lumi_fraction * bbc_eff * vertex_eff * XnXn_correction * pid_eff / purity_correction;

    // Scale Corrected by global efficiencies
    mass->Scale(1.0 / total_eff);
    pt->Scale(1.0 / total_eff);
    mY->Scale(1.0 / total_eff);
    
    // Scale all by bin width
    mass->Scale(1.0, "width");
    pt->Scale(1.0, "width");
    mY->Scale(1.0, "width");
    //pt_raw->Scale(1.0, "width"); 
    //mass_raw->Scale(1.0, "width"); // Uncorrected for physics, but normalized for mass bin sizing

    // --- DRAWING ---
    
    makeCanvas();
    mass->SetTitle("Corrected Mass");
    mass->Draw();

    makeCanvas();
    pt->SetTitle("Corrected Pt Cross Section");
    pt->Draw();

    makeCanvas();
    mY->SetTitle("Corrected Rapidity");
    mY->Draw();

    makeCanvas();
    pt_raw->SetTitle("Uncorrected Pt Yield (Bin Width Scaled)");
    pt_raw->SetLineColor(kRed);
    pt_raw->Draw();

    makeCanvas();
    mass_raw->SetTitle("Uncorrected Mass Yield for pT < 0.1 GeV/c");
    mass_raw->SetLineColor(kBlue);
    mass_raw->Draw();
}