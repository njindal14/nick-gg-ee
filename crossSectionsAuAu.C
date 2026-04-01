#include <iostream>
#include "FemtoPairFormat.h"
double const pi = 3.14159265358979323;
double c = 3.0e1; //in cm/ns
double me2 = pow(0.00051,2);

int ican2 = 0;
void makeCanvas()  {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican2++ ), "", 900, 600);
    can->SetTopMargin(0.1);
    can->SetRightMargin(0.1);
    can->SetBottomMargin(0.15);
}

void crossSectionsAuAu() {

    const Int_t ptbins = 24;
    Double_t edgespt[ptbins + 1] = {0.0, 0.002, 0.004, .006, 0.008, 0.01, 0.013, 0.017, 0.022, 0.025, 0.03, .035, .04, .045, .05, .055, .06, .065, .07, .075, .08, .085, .09, .095, 0.1};
    const Int_t massbins = 19;
    Double_t edgesmass[massbins + 1] = {0.4, .45,  0.5, .55, 0.6, .65, .7, .75, .8, .85, .9, .95, 1, 1.1, 1.2, 1.4, 1.6, 1.8, 2, 2.5};
    const Int_t ybins = 20;
    Double_t edgesy[ybins + 1] = {-1, -.9, -.8, -.7, -.6, -.5, -.4, -.3, -.2, -.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
    const Int_t pt2bins = 17;
    Double_t edgespt2[pt2bins + 1] = {0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, .004, .0045, .005, .006, .007, .008, .01, .015 };
    const Int_t cosbins = 14;
    Double_t edgescos[cosbins + 1] = {0, .05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7};
    auto * PtMY = new TH3F("PtMY", "PtMY", 24, edgespt, 19, edgesmass, 20, edgesy);
    auto * mPt2Mcostheta = new TH3F("mPt2Mcostheta", "mPt2Mcostheta",17, edgespt2, 19, edgesmass, 14, edgescos );

    TChain * ch = new TChain("PairDst");
    ch->Add("/Users/Nick/STAR/breit-wheeler/rootFiles/slim_pair_dst_Run10AuAu.root");
    //ch->Add("/Users/Nick/STAR/breit-wheeler/rootFiles/slim_pair_dst_Run11AuAu.root");
    TTreeReader myReader(ch);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

    TLorentzVector lv1, lv2, lv, lvn;

    //open efficiency root files  get the 1d histograms for relevant efficiencies
    TFile * trackingEfficiencies = new TFile("/Users/Nick/STAR/breit-wheeler/rootFiles/output_eff_3D_NOMINAL.root");
    TH3D * masspty = (TH3D*)trackingEfficiencies->Get("mMass");
    TH3D * mc_masspty = (TH3D*)trackingEfficiencies->Get("mc_mMass");

    //trying with my simulation results
    /*TFile * trackingEfficiencies = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/simulation_plots_new_Au.root");
    TH3D * masspty = (TH3D*)trackingEfficiencies->Get("RcPtMY");
    TH3D * mc_masspty = (TH3D*)trackingEfficiencies->Get("McPtMY");*/

    TEfficiency * ptmy_eff = new TEfficiency( * masspty, * mc_masspty);

    TH1D * pt_rc = masspty->ProjectionX();
    TH1D * mc_pt = mc_masspty->ProjectionX();
    TH1D * pt_eff = (TH1D*)pt_rc->Clone();
    pt_eff->Divide(mc_pt);

    TH1D * m_rc = masspty->ProjectionY();
    TH1D * mc_m = mc_masspty->ProjectionY();
    TEfficiency * m_eff = new TEfficiency (* m_rc, * mc_m);




    TFile * tofmatchEfficiencies = new TFile("/Users/Nick/STAR/breit-wheeler/rootFiles/TOF_efficiency_y2017_TpcRS_NOMINAL.root");
    TH1D * tofMatch_eff_mass = (TH1D*)tofmatchEfficiencies->Get("tof_eff_mass");
    TH1D * tofMatch_eff_pt = (TH1D*)tofmatchEfficiencies->Get("tof_eff_pt");
    TH1D * tofMatch_eff_pt2 = (TH1D*)tofmatchEfficiencies->Get("tof_eff_pt2");
    TH1D * tofMatch_eff_y = (TH1D*)tofmatchEfficiencies->Get("tof_eff_y");

    TFile * QED = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/QED_MB_pt_200GeVUPC_WangmeiPro.root");
    TH1D * qed_pt = (TH1D*)QED->Get("ht");
    TH1D * qed_pt2 = (TH1D*)QED->Get("ht2");
    qed_pt->Scale(10);

    TH1F * zVertex11 = new TH1F("zVertex11", "zVertex11", 500, -200, 200);




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

        zVertex11->Fill(pair->mVertexZ);


        if(lv1.Pt() < 0.2 || lv2.Pt() < 0.2) continue;
        if( fabs(lv1.Eta()) > 1 || fabs(lv2.Eta()) > 1 || fabs(lv.Rapidity()) > 1) continue;
        if(pair->d1_mNHitsFit < 20 || pair->d2_mNHitsFit < 20) continue;
        if(pair->d1_mNHitsDedx < 15 || pair->d2_mNHitsDedx < 15) continue;


        TLorentzVector lvbeam;
        lvbeam.SetPxPyPzE(0, 0, sqrt(100*100 - .938*.938), 100 );
        lvbeam.Boost( -(lv.BoostVector()) );
        float costheta = fabs(cos( lv1.Angle(lvbeam.Vect() ) ));



        if( fabs(pair->mVertexZ) < 100 &&  pair->mGRefMult <= 4 && pair->mChargeSum == 0 && pair->d1_mDCA < 1 && pair->d2_mDCA < 1 &&
        pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag!=0 && fabs(lv.Rapidity()) <= 1 && fabs(ddTofVal) < 0.4 && ddTofVal !=0 && chiee < 10 && 3*chiee < chipipi) {

            PtMY->Fill(lv.Pt(), lv.M(), lv.Rapidity());
            mPt2Mcostheta->Fill(pow(lv.Pt(), 2), lv.M(), costheta);


        }

    }


        //3d correction
    for(int ix = 1; ix <= PtMY->GetNbinsX(); ix ++){
        for(int iy = 1; iy <= PtMY->GetNbinsY(); iy ++){
            for(int iz = 1; iz <= PtMY->GetNbinsZ(); iz ++){

                double bincontent = PtMY->GetBinContent(ix, iy, iz);
                double binerror = PtMY->GetBinError(ix, iy, iz);
                std::cout << "bin content: " << bincontent << "\n";
                //std:: cout << "bin error: " << binerror << "\n";

                double pt_val = PtMY->GetXaxis()->GetBinCenter(ix);
                double m_val = PtMY->GetYaxis()->GetBinCenter(iy);
                double y_val = PtMY->GetZaxis()->GetBinCenter(iz);

                std:: cout << "pt val: " << pt_val << "\n";
                std::cout << "m val: " << m_val << "\n";
                std::cout << "y val: " << y_val << "\n";


                int eff_bin = ptmy_eff->FindFixBin(pt_val, m_val, y_val);
                double reco_eff = ptmy_eff->GetEfficiency(eff_bin);

                double tof_eff_pt = tofMatch_eff_pt->GetBinContent(tofMatch_eff_pt->FindBin(pt_val));


                if(reco_eff == 0) reco_eff = 1;
                if(tof_eff_pt == 0) tof_eff_pt = 1;
                PtMY->SetBinContent(ix, iy, iz, bincontent/(reco_eff*tof_eff_pt));
                PtMY->SetBinError(ix, iy, iz, binerror/(reco_eff*tof_eff_pt));
                std::cout << "new bin content: " << PtMY->GetBinContent(ix, iy, iz) << "\n";
            }
        }
    }


        //3d correction
    for(int ix = 1; ix <= mPt2Mcostheta->GetNbinsX(); ix ++){
        for(int iy = 1; iy <= mPt2Mcostheta->GetNbinsY(); iy ++){
            for(int iz = 1; iz <= mPt2Mcostheta->GetNbinsZ(); iz ++){

                double bincontent = mPt2Mcostheta->GetBinContent(ix, iy, iz);
                double binerror = mPt2Mcostheta->GetBinError(ix, iy, iz);
                std::cout << "bin content: " << bincontent << "\n";
                //std:: cout << "bin error: " << binerror << "\n";

                double pt2_val = mPt2Mcostheta->GetXaxis()->GetBinCenter(ix);
                double m_val = mPt2Mcostheta->GetYaxis()->GetBinCenter(iy);
                double costheta_val = mPt2Mcostheta->GetZaxis()->GetBinCenter(iz);

                std:: cout << "pt val: " << pt2_val << "\n";
                std::cout << "m val: " << m_val << "\n";
                std::cout << "y val: " << costheta_val << "\n";

                double reco_eff_m = m_eff->GetEfficiency(m_eff->FindFixBin(m_val));
                double tof_eff_m = tofMatch_eff_mass->GetBinContent(tofMatch_eff_mass->FindBin(m_val));


                if(reco_eff_m == 0) reco_eff_m = 1;
                if(tof_eff_m == 0) tof_eff_m = 1;
                mPt2Mcostheta->SetBinContent(ix, iy, iz, bincontent/(reco_eff_m*tof_eff_m));
                mPt2Mcostheta->SetBinError(ix, iy, iz, binerror/(reco_eff_m*tof_eff_m));
                std::cout << "new bin content: " << mPt2Mcostheta->GetBinContent(ix, iy, iz) << "\n";
            }
        }
    }




    //project onto 1d axes
    TH1D *mass = PtMY->ProjectionY("mass", PtMY->GetXaxis()->FindFixBin(0.), PtMY->GetXaxis()->FindFixBin(0.1) );
    TH1D *pt = PtMY->ProjectionX("pt", PtMY->GetYaxis()->FindFixBin(0.4), PtMY->GetYaxis()->FindFixBin(.76 ));
    TH1D * mY = PtMY->ProjectionZ("rapidity", PtMY->GetXaxis()->FindFixBin(0.), PtMY->GetXaxis()->FindFixBin(0.1),
        PtMY->GetYaxis()->FindFixBin(0.4), PtMY->GetYaxis()->FindFixBin(0.76));


    //get pt2 and cos(theta)
    TH1D * mPt2AuAu = mPt2Mcostheta->ProjectionX("pt2", mPt2Mcostheta->GetYaxis()->FindFixBin(0.4), mPt2Mcostheta->GetYaxis()->FindFixBin(0.76) );
    TH1D * cosTheta = mPt2Mcostheta->ProjectionZ("cosTheta", mPt2Mcostheta->GetYaxis()->FindFixBin(0.4), mPt2Mcostheta->GetYaxis()->FindFixBin(.76) );



    //global correction factors (not bin-by-bin)
    double luminosity = 679262; //for run 10
    //double luminosity = 859435;
    double lumi_fraction = .69677721; //taken from jdb
    double bbc_eff = 0.683; //taken from jdb
    double vertex_eff = 0.633; //taken from JDB analysis
    double purity_correction = 0.9917*.975; //subject to change, from chi2
    double pid_eff = 0.96*0.97; //first is tof eff, second is chi2 eff
    //double tpc_eff = 1.0; //.85; //0.646;//.8*.8; //100% for each track -- no missing sectors -- this may change after looking at simulation
    double XnXn_correction = 1/2.46; //still an estimate, need to recalculate


    double total_eff = luminosity*lumi_fraction*bbc_eff*vertex_eff*XnXn_correction*pid_eff/purity_correction;



      //scale and draw cross sections
    mass->Scale(1/total_eff);
    pt->Scale(1/total_eff);
    mPt2AuAu->Scale(1/total_eff);
    mY->Scale(1/total_eff);
    cosTheta->Scale(1/total_eff);
    mass->Scale(1.0, "width");
    pt->Scale(1.0, "width");
    mPt2AuAu->Scale(1.0, "width");
    mY->Scale(1.0, "width");
    cosTheta->Scale(1.0, "width");




    makeCanvas();
    gPad->SetLogy();
    mass->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
    mass->GetYaxis()->SetTitle("#frac{d#sigma (#gamma#gamma --> e^{+}e^{-})}{dM} (mb/(GeV/c^{2}))");
    mass->SetTitle("M_{ee} Differential Cross Section, Au+Au at 200 GeV");
    mass->SetMarkerStyle(29);
    mass->SetMarkerColor(46);
    mass->SetMarkerSize(1.5);
    mass->SetStats(false);
    gStyle->SetErrorX(0);
    mass->Draw("PE");
    //gPad->Print("note_plots/results_plots/Au_MassXSec.png");

    makeCanvas();
    //pt->Scale(0.000000403);
    pt->GetXaxis()->SetTitle("pT_{ee} (GeV/c)");
    pt->GetYaxis()->SetTitle("#frac{d#sigma (#gamma#gamma --> e^{+}e^{-})}{dpT} (mb/(GeV/c))");
    pt->SetTitle("Differential Cross Section, Au+Au Run 10 at 200 GeV");
    pt->SetMarkerStyle(29);
    pt->SetMarkerColor(46);
    pt->SetMarkerSize(1.5);
    //pt->SetStats(false);
    pt->Draw("PE");
    qed_pt->SetLineColor(kRed);
    //qed_pt->Scale(pt->Integral("width")/qed_pt->Integral("width"));
    qed_pt->Draw("hist;same;c");
    qed_pt->SetLineWidth(3);

    auto * legend = new TLegend(0.8,0.75,1.,.9);
    legend->SetHeader("Legend");
    legend->AddEntry(pt,"AuAu Run10","l");
    legend->AddEntry(qed_pt,"QED Au+Au (XnXn)","l");
    legend->Draw("same");

    //gPad->Print("note_plots/results_plots/Au_PtXSec.png");

    makeCanvas();
    gPad->SetLogy();
    mPt2AuAu->GetXaxis()->SetTitle("pT_{ee}^{2} (GeV/c)^{2}");
    mPt2AuAu->GetYaxis()->SetTitle("#frac{d#sigma (#gamma#gamma --> e^{+}e^{-})}{dpT^{2}} (mb/(GeV/c)^{2})");
    mPt2AuAu->SetTitle("P_{T}^{2} Differential Cross Section, Au+Au at 200 GeV");
    mPt2AuAu->SetMarkerStyle(29);
    mPt2AuAu->SetMarkerColor(46);
    mPt2AuAu->SetMarkerSize(1.5);
    mPt2AuAu->SetStats(false);
    mPt2AuAu->Draw("PE");
    //qed_pt2->SetLineColor(kRed);
    //qed_pt2->Scale(10000000);
    //qed_pt2->Draw("hist;same");
    //gPad->Print("note_plots/results_plots/Au_Pt2XSec.png");


    makeCanvas();
    mY->GetXaxis()->SetTitle("y_{ee}");
    mY->GetYaxis()->SetTitle("#frac{d#sigma (#gamma#gamma --> e^{+}e^{-})}{dy} (mb)");
    mY->SetTitle("y_{ee} Differential Cross Section, Au+Au at 200 GeV");
    mY->SetMarkerStyle(29);
    mY->SetMarkerColor(46);
    mY->SetMarkerSize(1.5);
    mY->SetStats(false);
    gPad->SetLogy();
    mY->Draw("PE");
    //gPad->Print("note_plots/results_plots/Au_YXSec.png");


    makeCanvas();
    cosTheta->GetXaxis()->SetTitle("|cos(#theta')|");
    cosTheta->GetYaxis()->SetTitle("#frac{d#sigma (#gamma#gamma --> e^{+}e^{-})}{d|cos(#theta')|} (mb)");
    cosTheta->SetTitle("|cos(#theta')| Differential Cross Section, Au+Au at 200 GeV");
    cosTheta->SetMarkerStyle(29);
    cosTheta->SetMarkerColor(46);
    cosTheta->SetMarkerSize(1.5);
    cosTheta->SetStats(false);
    cosTheta->Draw("PE");




    Double_t error_m;
    Double_t error_pt;
    Double_t error_pt2;
    Double_t error_y;
    std::cout << "mass total section : " << mass->IntegralAndError(1, 20, error_m, "width") << " +- " << error_m << "\n";
    std::cout << "pT cross section integral: " << pt->IntegralAndError(1, 20, error_pt, "width") << " +- " << error_pt << "\n";
    //std::cout << "pT2 cross section integral: " << mPt2->IntegralAndError(1, 20, error_pt2, "width") << " +- " << error_pt2 << "\n";
    std::cout << "Y cross section integral: " << mY->IntegralAndError(1, 20, error_y, "width") << " +- " << error_y << "\n";
    //Int_t n = 4;
    Double_t x[3]  = {1, 2, 3};
    Double_t y[3]  = {mass->IntegralAndError(1, 20, error_m, "width"), pt->IntegralAndError(1, 20, error_pt, "width"), mY->IntegralAndError(1, 20, error_y, "width")};
    Double_t ex[3] = {0.66, 0.66, 0.66};
    Double_t ey[3] = {error_m, error_pt, error_y};
    makeCanvas();
    TGraph * gr = new TGraphErrors(3,x,y,ex,ey);
    gr->SetTitle("Absolute Cross Sections AuAu");
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->GetXaxis()->SetBinLabel(20,"M_{ee} (GeV/c^{2})");
    gr->GetXaxis()->SetBinLabel(50, "pT (GeV/c)");
    gr->GetXaxis()->SetBinLabel(75, "y_{ee} [1]");
    gr->GetYaxis()->SetTitle("Integrated Cross Section (mb)");
    gr->Draw("AP");

    //write to root file
    TFile file("output_root_files/crossSectionsAuAu.root", "RECREATE");
    mass->Write("mass");
    pt->Write("pt");
    mPt2AuAu->Write("pt2");
    mY->Write("y");
    cosTheta->Write("costheta");

    makeCanvas();
    zVertex11->Draw();


    makeCanvas();
    pt_eff->Draw();

    std::cout << "Integral of vz from -100 to 100 vs total: " << zVertex11->Integral(zVertex11->FindBin(-100), zVertex11->FindBin(100)) << " vs " << zVertex11->Integral() << "\n";

}
