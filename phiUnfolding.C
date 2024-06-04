#include <iostream>
int ican2 = 0;
void makeCanvas()  {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican2++ ), "", 900, 600);
    can->SetTopMargin(0.08);
    can->SetRightMargin(0.15);
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


void phiUnfolding(){
   
    auto * phifit = new TF1("phifit", phiFit, -3.15,3.15,5);


    TFile * sim_out = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/simulation_plots.root");
    TH2F * mcphipt = (TH2F*)sim_out->Get("mcphipt");
    TH2F * recophipt = (TH2F*)sim_out->Get("recophipt");

    TH3F * mcrecophipt = (TH3F*)sim_out->Get("mcrecophipt");
    TH2F * mcrecophi = (TH2F*)sim_out->Get("mcrecophi");


    TH1F * mcphi = (TH1F*)mcphipt->ProjectionX(0,1);
    TH1F * recophi = (TH1F*)recophipt->ProjectionX(0,1);

    makeCanvas();
    mcrecophi->Draw("colz");

    makeCanvas();
    mcphipt->Draw("colz");

    makeCanvas();
    recophipt->Draw("colz");


}