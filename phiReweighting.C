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

    TFile * theory_curves = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/theory_modulations.root");
    TH1F * qed2phi = (TH1F*)theory_curves->Get("QED2phi");
    TH1F * qed4phi = (TH1F*)theory_curves->Get("QED4phi");


    
    //routine to reweight mcphi and see what happens for reco phi
   
    makeCanvas();
    //project mcphi and recophi into first pt bin
    TH1F * mcphi = (TH1F*)mcphipt->ProjectionX(0,1);
    TH1F * recophi = (TH1F*)recophipt->ProjectionX(0,1);
    //normalize
    mcphi->Scale(1/mcphi->Integral());
    recophi->Scale(1/recophi->Integral());
    int nbins = mcphi->GetNbinsX();
    TH1F * flattened = new TH1F("flattened", "flattened", nbins, -3.1415927, 3.1415927);
    mcphi->Fit("phifit", "", "", -3.15,3.15);
    gStyle->SetOptFit(1111); 
    double A0 = phifit->GetParameter(0);
    double A1 = phifit->GetParameter(1);
    double A2 = phifit->GetParameter(2);
    double A3 = phifit->GetParameter(3);
    double A4 = phifit->GetParameter(4);


    //calculate weight to flatten, reweight mcphi and reco according to those weights
    for(int i = 1; i <= nbins; ++i){
        double weight = 1/(1 + A1*cos(mcphi->GetBinCenter(i)) + A2*cos(2*mcphi->GetBinCenter(i)) + A3*cos(3*mcphi->GetBinCenter(i)) + A4*cos(4*mcphi->GetBinCenter(i)));
        cout << "weight: " << weight << "\n";
        double newBinContent = weight*mcphi->GetBinContent(i);
        double newBinError = weight*mcphi->GetBinError(i);
        double newBinContentReco = weight*recophi->GetBinContent(i);
        double newBinErrorReco = weight*recophi->GetBinError(i);
        flattened->SetBinContent(i, newBinContent);
        flattened->SetBinError(i, newBinError);
        recophi->SetBinContent(i, newBinContentReco);
        recophi->SetBinError(i, newBinErrorReco);
        cout << "Set new bin content: " << "Phi = " << mcphi->GetBinCenter(i) << ", content = " << newBinContent << "\n";
    }

    mcphi->SetLineColor(kGreen);
    mcphi->Draw("PE");
    flattened->SetLineColor(kBlue);
    flattened->Draw("PE;same");

    //clone flattened to embed a cos2phi signal into it and into reco
    TH1F * embeddedSignal = (TH1F*)flattened->Clone();
    double twophistrength = -.2;

    for(int i = 1; i <= nbins; ++i){
        embeddedSignal->SetBinContent(i, embeddedSignal->GetBinContent(i)*(1+twophistrength*cos(2*embeddedSignal->GetBinCenter(i))));
        recophi->SetBinContent(i, recophi->GetBinContent(i)*(1+twophistrength*cos(2*recophi->GetBinCenter(i))));
    }
    embeddedSignal->SetLineColor(kRed);
    embeddedSignal->Draw("PE;same");

    //recophi->Fit("phifit", "", "", -3.15,3.15);
    gStyle->SetOptFit(1111); 
    mcphi->GetXaxis()->SetTitle("#Delta #phi");
    mcphi->GetYaxis()->SetTitle("Counts");
    recophi->SetLineColor(kBlack);
    recophi->Draw("PE;same");
    //mcphi->GetYaxis()->SetRangeUser(0.03,.055);
    TLegend * legend = new TLegend(0.7,0.1,1,0.3);
    legend->AddEntry(mcphi,"MC Phi, pT bin 0");
    legend->AddEntry(flattened, "flattened");
    legend->AddEntry(embeddedSignal, "Embedded 20%");
    legend->AddEntry(recophi, "Weighted reco");
    legend->Draw("same");

    
    makeCanvas();
    mcphi->SetLineColor(kBlack);
    mcphi->SetTitle("MC phi with Fit");
    mcphi->Draw("PE");

    makeCanvas();
    recophi->SetLineColor(kGreen);
    recophi->Fit("phifit", "", "", -3.15,3.15);
    gStyle->SetOptFit(1111);
    recophi->SetTitle("Reco phi with Fit");
    recophi->GetXaxis()->SetTitle("#Delta #phi");
    recophi->Draw("PE;same");

    makeCanvas();
    flattened->Fit("phifit", "", "", -3.15,3.15);
    gStyle->SetOptFit(1111);
    flattened->Draw("PE;same");
    
    /*TLegend * legend2 = new TLegend(0.7,0.1,1,0.3);
    legend2->AddEntry(mcphi,"MC Phi, pT bin 0");
    legend2->AddEntry(recophi, "Weighted reco");
    legend2->Draw("same");*/

    makeCanvas();
    qed2phi->Draw("PE");
    qed4phi->Draw("PE;same");

}