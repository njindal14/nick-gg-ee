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

void crossSectionComparisons(){

    TFile * uux = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/crossSectionsUU.root");
    TH1F * uu_mass = (TH1F*)uux->Get("mMass");
    TH1F * uu_pt = (TH1F*)uux->Get("mPt");
    TH1F * uu_pt2 = (TH1F*)uux->Get("mPt2");
    TH1F * uu_y = (TH1F*)uux->Get("mY");

    TFile * aux = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/crossSectionsAuAu.root");
    TH1F * au_mass = (TH1F*)aux->Get("mMassAuAu");
    TH1F * au_pt = (TH1F*)aux->Get("mPtAuAu");
    TH1F * au_pt2 = (TH1F*)aux->Get("mPt2AuAu");
    TH1F * au_y = (TH1F*)aux->Get("mYAuAu");

    TFile * QEDU = new TFile("/Users/Nick/STAR/QEDCode/QED_MB_pt50100_0.40_0.76Uranium.root");
    TH1F * QEDUpt = (TH1F*)QEDU->Get("ht");

    TFile * QEDAu = new TFile("/Users/Nick/STAR/QEDCode/QED_MB_pt50100_0.40_0.76Gold.root");
    TH1F * QEDAupt = (TH1F*)QEDAu->Get("ht");



    
    makeCanvas();
    uu_mass->Scale(1/uu_mass->Integral());
    au_mass->Scale(1/au_mass->Integral());
    uu_mass->SetLineColor(kBlack);
    au_mass->SetLineColor(kRed);
    au_mass->Draw("PE");
    uu_mass->Draw("PE'same");
    TLegend * mass_leg = new TLegend(0.77,0.3,.97,0.45);
    mass_leg->SetHeader("Legend");
    mass_leg->AddEntry(uu_mass,"U+U Run 12","l");
    mass_leg->AddEntry(au_mass,"Au+Au Run10+11","l");
    mass_leg->Draw("same"); 
    gPad->Print("note_plots/results_plots/massComp.png");

    makeCanvas();
    TH1F * mass_ratio = (TH1F*)uu_mass->Clone();
    mass_ratio->Divide(au_mass);
    mass_ratio->Draw("PE");
    mass_ratio->SetTitle("Mass Differential Cross Section Ratio, UU(12)/AuAu(10+11)");
    TLine *mline = new TLine(0,1,3,1);
    mline->SetLineColor(kRed);
    mline->SetLineStyle(4);
    mline->Draw("same");
    gPad->Print("note_plots/results_plots/mass_xsec_ratio.png");

    
    makeCanvas();
    uu_pt->Scale(1/uu_pt->Integral("width"));
    au_pt->Scale(1/au_pt->Integral("width"));
    uu_pt->SetLineColor(kBlack);
    au_pt->SetLineColor(kRed);
    au_pt->Draw("PE");
    uu_pt->Draw("PE;same");
    TLegend * pt_leg = new TLegend(0.77,0.3,.97,0.45);
    pt_leg->SetHeader("Legend");
    pt_leg->AddEntry(uu_pt,"U+U Run 12","l");
    pt_leg->AddEntry(au_pt,"Au+Au Run10+11","l");
    pt_leg->Draw("same");  
    gPad->Print("note_plots/results_plots/ptComp.png");


    makeCanvas();
    TH1F * pt_ratio = (TH1F*)uu_pt->Clone();
    pt_ratio->Divide(au_pt);
    pt_ratio->Draw("PE");
    pt_ratio->SetTitle("pT Differential Cross Section Ratio, UU(12)/AuAu(10+11)");
    pt_ratio->Fit("pol1", "", "", .01, .07);
    gStyle->SetOptFit(1111);
    TLine *ptline = new TLine(0,1,.1,1);
    ptline->SetLineColor(kRed);
    ptline->SetLineStyle(4);
    ptline->Draw("same");
    gPad->Print("note_plots/results_plots/pt_xsec_ratio.png");



    makeCanvas();
    uu_pt2->Scale(1/uu_pt2->Integral("width"));
    au_pt2->Scale(1/au_pt2->Integral("width"));
    uu_pt2->SetLineColor(kBlack);
    au_pt2->SetLineColor(kRed);
    au_pt2->Draw("PE");
    uu_pt2->Draw("PE;same");
    TLegend * pt2_leg = new TLegend(0.77,0.3,.97,0.45);
    pt2_leg->SetHeader("Legend");
    pt2_leg->AddEntry(uu_pt2,"U+U Run 12","l");
    pt2_leg->AddEntry(au_pt2,"Au+Au Run10+11","l");
    pt2_leg->Draw("same");  
    gPad->Print("note_plots/results_plots/pt2Comp.png");

    makeCanvas();
    TH1F * pt2_ratio = (TH1F*)uu_pt2->Clone();
    pt2_ratio->Divide(au_pt2);
    pt2_ratio->Draw("PE");
    pt2_ratio->SetTitle("pT^{2} Differential Cross Section Ratio, UU(12)/AuAu(10+11)");
    pt2_ratio->Fit("pol1", "", "", .0001, .005);
    gStyle->SetOptFit(1111);
    TLine *pt2line = new TLine(0,1,.01,1);
    pt2line->SetLineColor(kRed);
    pt2line->SetLineStyle(4);
    pt2line->Draw("same");
    gPad->Print("note_plots/results_plots/pt2_xsec_ratio.png");


    makeCanvas();
    uu_y->Scale(1/uu_y->Integral("width"));
    au_y->Scale(1/au_y->Integral("width"));
    uu_y->SetLineColor(kBlack);
    au_y->SetLineColor(kRed);
    au_y->Draw("PE");
    uu_y->Draw("PE;same");
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
    y_ratio->SetTitle("Rap. Differential Cross Section Ratio, UU(12)/AuAu(10+11)");
    TLine *yline = new TLine(-1,1,1,1);
    yline->SetLineColor(kRed);
    yline->SetLineStyle(4);
    yline->Draw("same");
    gPad->Print("note_plots/results_plots/y_xsec_ratio.png");

    makeCanvas();
    //QEDUpt->Scale(1/QEDUpt->Integral());
    //QEDAupt->Scale(1/QEDAupt->Integral());
    QEDUpt->Draw();
    QEDAupt->SetLineColor(kRed);
    QEDAupt->Draw("same");



}