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
    can->SetBottomMargin(0.18);
}

void crossSectionComparisons(){

    TFile * uux = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/crossSectionsUU_NEW.root");
    TH1F * uu_mass = (TH1F*)uux->Get("mass");
    TH1F * uu_pt = (TH1F*)uux->Get("pt");
    TH1F * uu_pt2 = (TH1F*)uux->Get("pt2");
    TH1F * uu_y = (TH1F*)uux->Get("y");
    TH1F * uu_cos = (TH1F*)uux->Get("costheta");

    TFile * aux = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/crossSectionsAuAu.root");
    TH1F * au_mass = (TH1F*)aux->Get("mass");
    TH1F * au_pt = (TH1F*)aux->Get("pt");
    TH1F * au_pt2 = (TH1F*)aux->Get("pt2");
    TH1F * au_y = (TH1F*)aux->Get("y");
    TH1F * au_cos = (TH1F*)aux->Get("costheta");

    TFile * slight_uu = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/SL_plotsBetter.root");
    TH1F * slight_uupt = (TH1F*)slight_uu->Get("mPt");
    TH1F * slight_uuMass = (TH1F*)slight_uu->Get("mMass");

    TFile * slight_au = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/SL_plotsBetter_Au.root");
    TH1F * slight_aupt = (TH1F*)slight_au->Get("mPt");
    TH1F * slight_auMass = (TH1F*)slight_au->Get("mMass");




    //TFile * QED_ratio = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/XnXnRatioAsPt.root");
    //TH1D * QED_pt_ratio = (TH1D*)QED_ratio->Get("hRatio");


   //qed uranium curve
    TFile * QED_pt_UU = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/QED_MB_pt_UU193GeVUPC.root");
    TH1D * QED_pt_UU_hist = (TH1D*)QED_pt_UU->Get("ht");
    QED_pt_UU_hist->Scale(10); //convert from fm^2 to mbarn

    //qed uranium curve
    TFile * UU_deformed_minor = new TFile("/Users/Nick/STAR/breit-wheeler/deformed_UU/QED_pt_UUminorAxis_193GeVUPCXnXn.root");
    TH1D * minor_UU = (TH1D*)UU_deformed_minor->Get("ht");
    minor_UU->Scale(10); //convert from fm^2 to mbarn

    //qed uranium curve
    TFile * UU_deformed_major = new TFile("/Users/Nick/STAR/breit-wheeler/deformed_UU/QED_pt_UUmajorAxis_193GeVUPCXnXn.root");
    TH1D * major_UU = (TH1D*)UU_deformed_major->Get("ht");
    major_UU->Scale(10); //convert from fm^2 to mbarn

    //qed gold curve
    TFile * QED_pt_au = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/QED_MB_pt_200GeVUPC_WangmeiPro.root");
    TH1D * QED_pt_au_hist = (TH1D*)QED_pt_au->Get("ht");
    QED_pt_au_hist->Scale(10); //convert from fm^2 to mbarn

   


    
    makeCanvas();
    uu_mass->SetLineColor(kBlack);
    au_mass->SetLineColor(kRed);
    uu_mass->Draw("PE");
    au_mass->Draw("PE;same");
    uu_mass->SetTitle("Differential Cross Sections, M_{ee}");
    TLegend * mass_leg = new TLegend(0.77,0.3,.97,0.45);
    mass_leg->SetHeader("Legend");
    mass_leg->AddEntry(uu_mass,"U+U Run 12","l");
    mass_leg->AddEntry(au_mass,"Au+Au Run10+11","l");
    mass_leg->Draw("same"); 
    gPad->SetLogy();
    gPad->Print("note_plots/results_plots/massComp.png");

    makeCanvas();
    TH1F * mass_ratio = (TH1F*)uu_mass->Clone();
    mass_ratio->Divide(au_mass);
    mass_ratio->Draw("PE");
    mass_ratio->SetTitle("Mass Differential Cross Section Ratio, UU(XnXn)/AuAu(1-3n)");
    mass_ratio->GetYaxis()->SetTitle("Cross Section Ratio, UU(XnXn)/AuAu(1-3n)");
    TLine *mline = new TLine(0,1,3,1);
    mass_ratio->GetXaxis()->SetRangeUser(0, 2.5);
    mline->SetLineColor(kRed);
    mline->SetLineStyle(4);
    mline->Draw("same");
    gPad->Print("note_plots/results_plots/mass_xsec_ratio.png");

    
    makeCanvas();
    uu_pt->SetLineColor(kBlack);
    au_pt->SetLineColor(kRed);
    uu_pt->Draw("PE");
    uu_pt->GetYaxis()->SetRangeUser(0, 17);
    au_pt->Draw("PE;same");
    uu_pt->SetTitle("Differential Cross Sections, P_{T}");
    QED_pt_UU_hist->SetLineColor(kBlue);
    QED_pt_au_hist->SetLineColor(kGreen);
    QED_pt_UU_hist->Draw("hist;same");
    QED_pt_au_hist->Draw("hist;same");
    minor_UU->SetLineColor(kPink);
    major_UU->SetLineColor(kGreen+2);
    minor_UU->Draw("hist;same");
    major_UU->Draw("hist;same");
    TLegend * pt_leg = new TLegend(0.77,0.3,.97,0.45);
    pt_leg->SetHeader("Legend");
    pt_leg->AddEntry(uu_pt,"U+U Run 12","l");
    pt_leg->AddEntry(au_pt,"Au+Au Run10+11","l");
    pt_leg->AddEntry(QED_pt_UU_hist, "QED code UU Wangmei", "l");
    pt_leg->AddEntry(QED_pt_au_hist, "QED code AuAu Xiaofeng", "l");
    pt_leg->AddEntry(minor_UU, "UU Spherical r = 5.9 fm", "l");
    pt_leg->AddEntry(major_UU, "UU Spherical r = 9.8 fm", "l");
    pt_leg->Draw("same");  
    gPad->Print("note_plots/results_plots/ptComp.png");



    slight_uupt->Scale(1/slight_uupt->Integral("width"));
    slight_aupt->Scale(1/slight_aupt->Integral("width"));

    slight_uuMass->Scale(1/slight_uuMass->Integral("width"));
    slight_auMass->Scale(1/slight_auMass->Integral("width"));

    makeCanvas();
    slight_uupt->Draw("PE");
    slight_aupt->SetLineColor(kRed);
    slight_aupt->Draw("PE;same");

    makeCanvas();
    slight_uuMass->Draw("PE");
    slight_auMass->SetLineColor(kRed);
    slight_auMass->Draw("PE;same");



    makeCanvas();
    TH1F * slight_ratio_pt = (TH1F*)slight_uupt->Clone();
    slight_ratio_pt->Divide(slight_aupt);
    //slight_ratio_pt->Draw();
    //slight_ratio_pt->SetTitle("PT Ratio starlight, UU/AuAu");

    makeCanvas();
    TH1F * qed_ratio_pt = (TH1F*)QED_pt_UU_hist->Clone();
    qed_ratio_pt->Divide(QED_pt_au_hist);


    makeCanvas();
    TH1F * pt_ratio = (TH1F*)uu_pt->Clone();
    //TH1F * qed_ratio = (TH1F*)qed_pt_U->Clone();
    pt_ratio->Divide(au_pt);
    //qed_ratio->Divide(qed_pt_Au);
    pt_ratio->Draw("PE");
    //qed_ratio->SetLineColor(kRed);
    //qed_ratio->Draw("hist;same");
    pt_ratio->GetXaxis()->SetRangeUser(0, 0.1);
    pt_ratio->GetYaxis()->SetTitle("Cross Section Ratio UU(XnXn)/AuAu(1-3n)");
    pt_ratio->GetYaxis()->SetRangeUser(0,3);
    pt_ratio->SetTitle("Differential Cross Section Ratio, UU(XnXn)/AuAu(1-3n)");
    pt_ratio->Fit("pol1", "", "", .015, .05);
    gStyle->SetOptFit();
    TLine *ptline = new TLine(0,1,.1,1);
    ptline->SetLineColor(kRed);
    ptline->SetLineStyle(4);
    ptline->Draw("same");
    slight_ratio_pt->SetLineColor(kGreen);
    slight_ratio_pt->Draw("PE;same");
    qed_ratio_pt->SetLineColor(kRed);
    qed_ratio_pt->Draw("hist;same");

    TLegend * pt_ratio_leg = new TLegend(0.77,0.3,.97,0.45);
    pt_ratio_leg->SetHeader("Legend");
    pt_ratio_leg->AddEntry(pt_ratio,"UU(XnXn)/AuAu(1-3n)","l");
    pt_ratio_leg->AddEntry(qed_ratio_pt,"UU(XnXn)/AuAu(XnXn)","l");
    pt_ratio_leg->Draw("same"); 

    gPad->Print("note_plots/results_plots/pt_xsec_ratio.png");



    makeCanvas();
    uu_pt2->SetLineColor(kBlack);
    au_pt2->SetLineColor(kRed);
    uu_pt2->Draw("PE");
    au_pt2->Draw("PE;same");
    uu_pt2->SetTitle("Differential Cross Sections, P_{T}^{2}");
    TLegend * pt2_leg = new TLegend(0.77,0.3,.97,0.45);
    pt2_leg->SetHeader("Legend");
    pt2_leg->AddEntry(uu_pt2,"U+U Run 12","l");
    pt2_leg->AddEntry(au_pt2,"Au+Au Run10+11","l");
    pt2_leg->Draw("same");  
    gPad->SetLogy();
    gPad->Print("note_plots/results_plots/pt2Comp.png");

    makeCanvas();
    TH1F * pt2_ratio = (TH1F*)uu_pt2->Clone();
    pt2_ratio->Divide(au_pt2);
    pt2_ratio->Draw("PE");
    pt2_ratio->GetYaxis()->SetTitle("Cross Section Ratio, UU(XnXn)/AuAu(1-3n)");
    pt2_ratio->SetTitle("Differential Cross Section Ratio, UU(XnXn)/AuAu(1-3n)");
    pt2_ratio->Fit("pol1", "", "", 0., .005);
    gStyle->SetOptFit();
    pt2_ratio->GetXaxis()->SetRangeUser(0, 0.02);
    TLine *pt2line = new TLine(0,1,.02,1);
    pt2line->SetLineColor(kRed);
    pt2line->SetLineStyle(4);
    pt2line->Draw("same");
    gPad->Print("note_plots/results_plots/pt2_xsec_ratio.png");


    makeCanvas();
    //uu_y->Scale(1/uu_y->Integral("width"));
    //au_y->Scale(1/au_y->Integral("width"));
    uu_y->SetLineColor(kBlack);
    au_y->SetLineColor(kRed);
    uu_y->SetStats(false);
    au_y->SetStats(false);
    uu_y->Draw("PE");
    au_y->Draw("PE;same");
    uu_y->SetTitle("Differential Cross Sections, Rapidity");
    TLegend * y_leg = new TLegend(0.77,0.3,.97,0.45);
    y_leg->SetHeader("Legend");
    y_leg->AddEntry(uu_y,"U+U Run 12","l");
    y_leg->AddEntry(au_y,"Au+Au Run10+11","l");
    y_leg->Draw("same"); 
    gPad->Print("note_plots/results_plots/yComp.png");

    makeCanvas();
    TH1F * y_ratio = (TH1F*)uu_y->Clone();
    y_ratio->Divide(au_y);
    y_ratio->Draw("PE");
    y_ratio->SetTitle("Differential Cross Section Ratio, UU(XnXn)/AuAu(1-3n)");
    y_ratio->GetYaxis()->SetTitle("Cross Section Ratio, UU(XnXn)/AuAu(1-3n)");
    y_ratio->GetXaxis()->SetRangeUser(-0.9, .9);
    TLine *yline = new TLine(-1,1,1,1);
    yline->SetLineColor(kRed);
    yline->SetLineStyle(4);
    yline->Draw("same");
    gPad->Print("note_plots/results_plots/y_xsec_ratio.png");



    makeCanvas();
    uu_cos->SetLineColor(kBlack);
    au_cos->SetLineColor(kRed);
    uu_cos->SetStats(false);
    au_cos->SetStats(false);
    uu_cos->Draw("PE");
    au_cos->Draw("PE;same");
    uu_cos->GetYaxis()->SetRangeUser(0, 0.7);
    uu_cos->SetTitle("Differential Cross Sections, |cos(#theta)|");
    TLegend * cos_leg = new TLegend(0.77,0.3,.97,0.45);
    cos_leg->SetHeader("Legend");
    cos_leg->AddEntry(uu_cos,"U+U Run 12","l");
    cos_leg->AddEntry(au_cos,"Au+Au Run10+11","l");
    cos_leg->Draw("same"); 


    makeCanvas();
    TH1F * cos_ratio = (TH1F*)uu_cos->Clone();
    cos_ratio->Divide(au_cos);
    cos_ratio->Draw("PE");
    cos_ratio->SetTitle("Differential Cross Section Ratio, UU(XnXn)/AuAu(1-3n)");
    cos_ratio->GetYaxis()->SetTitle("Cross Section Ratio, UU(XnXn)/AuAu(1-3n)");
    cos_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
    TLine *cosline = new TLine(-1,1,1,1);
    cosline->SetLineColor(kRed);
    cosline->SetLineStyle(4);
    cosline->Draw("same");

    //cos_ratio->GetXaxis()->SetRangeUser(-0.9, .9);
    //TLine *cosline = new TLine(-1,1,1,1);
    //cosline->SetLineColor(kRed);
    //cosline->SetLineStyle(4);
    //cosline->Draw("same");

      
    Double_t error_m_uu;
    Double_t error_pt_uu;
    Double_t error_y_uu;
    Double_t error_m_au;
    Double_t error_pt_au;
    Double_t error_y_au;
    Double_t error_cos_uu;
    Double_t error_cos_au;

    Double_t x[6]  = {1, 2, 3, 4, 5, 6};
    Double_t y[6]  = {uu_mass->IntegralAndError(1, 20, error_m_uu, "width"), uu_pt->IntegralAndError(1, 20, error_pt_uu, "width"), uu_y->IntegralAndError(1, 20, error_y_uu, "width"), 
    au_mass->IntegralAndError(1, 20, error_m_au, "width"), au_pt->IntegralAndError(1, 20, error_pt_au, "width"), au_y->IntegralAndError(1, 20, error_y_au, "width")};
    Double_t ex[6] = {0.33, 0.33, 0.33, .33, .33, .33};
    Double_t ey[6] = {error_m_uu, error_pt_uu, error_y_uu, error_m_au, error_pt_au, error_y_au};

    makeCanvas();
    TGraph * gr = new TGraphErrors(6,x,y,ex,ey);
    gr->SetTitle("Absolute Cross Sections UU and AuAu");
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);

    gr->GetXaxis()->SetBinLabel(15,"M_{ee} (GeV/c^{2}) UU");
    gr->GetXaxis()->SetBinLabel(28, "pT (GeV/c) UU");
    gr->GetXaxis()->SetBinLabel(42, "y_{ee} [1] UU"); 
    
    gr->GetXaxis()->SetBinLabel(58,"M_{ee} (GeV/c^{2}) AuAu");
    gr->GetXaxis()->SetBinLabel(73, "pT (GeV/c) AuAu");
    gr->GetXaxis()->SetBinLabel(88, "y_{ee} [1] AuAu"); 

    gr->GetYaxis()->SetTitle("Integrated Cross Sections (mb) UU and AuAu");  
    gr->Draw("AP");
    gPad->Print("note_plots/results_plots/UUAbsXSectionsNEW.png");

    //std::cout << "QED UU integral: " << qed_pt_U->Integral("width") << "\n";
   //std::cout << "QED Au integral: " << qed_pt_Au->Integral("width") << "\n";

    std::cout << "UU mass total cross section : " << uu_mass->IntegralAndError(1, 20, error_m_uu, "width") << " +- " << error_m_uu << "\n";
    std::cout << "AuAu mass total cross section : " << au_mass->IntegralAndError(1, 20, error_m_au, "width") << " +- " << error_m_au << "\n";






}