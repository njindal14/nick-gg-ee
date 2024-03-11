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


void phiReweighting(){
   
    auto * phifit = new TF1("phifit", phiFit, -3.15,3.15,5);

    TFile * starlight_out = new TFile("/Users/Nick/STAR/starlight/utils/SL_plotsBetter.root");
    TH2F * slightPhiPt = (TH2F*)starlight_out->Get("mDeltaPhiPt");

    TFile * sim_out = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/simulation_plots.root");
    TH2F * mcphipt = (TH2F*)sim_out->Get("mcphipt");
    TH2F * recophipt = (TH2F*)sim_out->Get("recophipt");


    
    makeCanvas();
    TH1F * slightPhi = (TH1F*)slightPhiPt->ProjectionX(0, 1);
    TH1F * mcphi = (TH1F*)mcphipt->ProjectionX(0,1);
    TH1F * recophi = (TH1F*)recophipt->ProjectionX(0,1);
    slightPhi->Rebin(2);
    slightPhi->Scale(1/slightPhi->Integral());
    mcphi->Scale(1/mcphi->Integral());
    recophi->Scale(1/recophi->Integral());
    int nbins = mcphi->GetNbinsX();
    TH1F * flattened = new TH1F("flattened", "flattened", nbins, -3.1415927, 3.1415927);
    mcphi->GetXaxis()->SetTitle("MC #Delta #phi");
    mcphi->GetYaxis()->SetTitle("Counts");
    //mcphi->Fit("phifit", "", "", -3.15,3.15);
    //gStyle->SetOptFit(1111); 
    double A0 = phifit->GetParameter(0);
    double A1 = phifit->GetParameter(1);
    double A2 = phifit->GetParameter(2);
    double A3 = phifit->GetParameter(3);
    double A4 = phifit->GetParameter(4);

    for(int i = 1; i < nbins+1; ++i){
        double weight = 1/(1 + A1*cos(mcphi->GetBinCenter(i)) + A2*cos(2*mcphi->GetBinCenter(i)) + A3*cos(3*mcphi->GetBinCenter(i)) + A4*cos(4*mcphi->GetBinCenter(i)));
        cout << "weight: " << weight << "\n";
        double newBinContent = weight*mcphi->GetBinContent(i);
        double newBinContentReco = weight*recophi->GetBinContent(i);
        flattened->SetBinContent(i, newBinContent);
        flattened->SetBinError(i, 0.001);
        recophi->SetBinContent(i, newBinContentReco);
        cout << "Set new bin content: " << "Phi = " << mcphi->GetBinCenter(i) << ", content = " << newBinContent << "\n";
    }

    mcphi->SetLineColor(kGreen);
    mcphi->Draw("PE");
    flattened->SetLineColor(kBlue);
    flattened->Draw("PE;same");
    TH1F * embeddedSignal = new TH1F("embeddedSignal", "embeddedSignal", nbins, -3.1415927, 3.1415927);
    double twophistrength = -.2;
    for(int i = 1; i < nbins +1; ++i){
        embeddedSignal->SetBinContent(i, flattened->GetBinContent(i)*(1+twophistrength*cos(2*flattened->GetBinCenter(i))));
        embeddedSignal->SetBinError(i, 0.001);
        recophi->SetBinContent(i, recophi->GetBinContent(i)*(1+twophistrength*cos(2*recophi->GetBinCenter(i))));


    }
    embeddedSignal->SetLineColor(kRed);
    embeddedSignal->Draw("PE;same");

    //recophi->Fit("phifit", "", "", -3.15,3.15);
    gStyle->SetOptFit(1111); 
    recophi->SetLineColor(kBlack);
    recophi->Draw("PE;same");
    mcphi->GetYaxis()->SetRangeUser(0.03,.055);
    TLegend * legend = new TLegend(0.7,0.1,1,0.3);
    legend->AddEntry(mcphi,"MC Phi, pT bin 0");
    legend->AddEntry(flattened, "flattened");
    legend->AddEntry(embeddedSignal, "Embedded 20%");
    legend->AddEntry(recophi, "Weighted reco");
    legend->Draw("same");

    makeCanvas();
    mcphi->SetLineColor(kBlack);
    mcphi->Draw("PE");
    recophi->SetLineColor(kGreen);
    mcphi->SetTitle("Reco vs MC Phi Same weighting");
    recophi->Draw("PE;same");
    TLegend * legend2 = new TLegend(0.7,0.1,1,0.3);
    legend2->AddEntry(mcphi,"MC Phi, pT bin 0");
    legend2->AddEntry(recophi, "Weighted reco");
    legend2->Draw("same");

}