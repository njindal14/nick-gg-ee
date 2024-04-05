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


void phiEmbedding(){
   
    TFile * starlight_out = new TFile("/Users/Nick/STAR/starlight/utils/SL_plotsBetter.root");
    TH2F * slightPhiPt = (TH2F*)starlight_out->Get("mDeltaPhiPt");

    TFile * sim_out = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/simulation_plots.root");
    TH2F * mcphipt = (TH2F*)sim_out->Get("mcphipt");
    TH2F * recophipt = (TH2F*)sim_out->Get("recophipt");

    TFile * theory_curves = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/theory_modulations.root");
    TH1F * qed2phi = (TH1F*)theory_curves->Get("QED2phi");
    TH1F * qed4phi = (TH1F*)theory_curves->Get("QED4phi");

    int nbinsmc = mcphipt->GetNbinsY();

    TH1F * reco2phiReweighted = new TH1F("reco2phiReweighted","", nbinsmc, 0, 0.3);
    TH1F * reco4phiReweighted = new TH1F("reco4phiReweighted","", nbinsmc, 0, 0.3);



    for(int i = 0; i <= 15; ++i){
     
        TH1 * mcphi = (TH1*)(mcphipt->ProjectionX( "h", i, i+1 )->Clone( TString::Format( "h_%d", i ) )) ;

        
        //routine to reweight mcphi and see what happens for reco phi
        //project mcphi and recophi into ith pt bin
        TH1 * recophi = (TH1*)(recophipt->ProjectionX( "h", i, i+1 )->Clone( TString::Format( "hr_%d", i ) )) ;

        float pt = mcphipt->GetYaxis()->GetBinCenter(i);

        cout << "Entries MC Phi = " << mcphi->GetEntries() << endl;

        auto * phifit = new TF1("phifit", phiFit, -3.14,3.14,5);
        auto * phifit2 = new TF1("phifit2", phiFit, -3.14,3.14,5);

        cout << "made projections" << "\n";
        //normalize
        mcphi->Scale(1/mcphi->Integral());
        recophi->Scale(1/recophi->Integral());
        int nbins = mcphi->GetNbinsX();
        TH1F * flattened = new TH1F("flattened", "flattened", nbins, -3.1415927, 3.1415927);
        
        mcphi->Fit("phifit", "", "", -3.14,3.14);
        double A0 = phifit->GetParameter(0);
        double A1 = phifit->GetParameter(1);
        double A2 = phifit->GetParameter(2);
        double A3 = phifit->GetParameter(3);
        double A4 = phifit->GetParameter(4);
        cout << "A2: " << A2 << "\n";
        cout << "fit mcphi " << "\n";

        auto c = new TCanvas(TString::Format("c_%d", i), "", 800, 600);
        gStyle->SetOptFit(1111);
        mcphi->Draw();

        //calculate weight to flatten, reweight mcphi and reco according to those weights
        for(int j = 0; j <= nbins; ++j){
            double mcphival = mcphi->GetBinCenter(j);
            double weight = 1/(1 + A1*cos(mcphival) + A2*cos(2*mcphival) + A3*cos(3*mcphival) + A4*cos(4*mcphival));
            cout << "weight: " << weight << "\n";
            cout << "mcphi bin content: " << mcphi->GetBinContent(j) << "\n";
            double newBinContent = weight*mcphi->GetBinContent(j);
            //double newBinError = 1*mcphi->GetBinError(j);
            double newBinContentReco = weight*recophi->GetBinContent(j);
            //double newBinErrorReco = 1*recophi->GetBinError(j);
            flattened->SetBinContent(j, newBinContent);
            //flattened->SetBinError(j, 0.0);
            recophi->SetBinContent(j, newBinContentReco);
            //recophi->SetBinError(j, 0.0);
            cout << "Set new bin content: " << "Phi = " << mcphi->GetBinCenter(j) << ", content = " << newBinContent << "\n";
        }

        //clone flattened to embed a cos2phi signal into it and into reco
        TH1F * embeddedSignal = (TH1F*)flattened->Clone();
        double twophistrength = qed2phi->GetBinContent(qed2phi->FindFixBin(pt));
        double fourphistrength = qed4phi->GetBinContent(qed4phi->FindFixBin(pt));

        for(int k = 0; k <= nbins; ++k){
            embeddedSignal->SetBinContent(k, embeddedSignal->GetBinContent(k)*(1+twophistrength*cos(2*embeddedSignal->GetBinCenter(k)) + fourphistrength*cos(4*embeddedSignal->GetBinCenter(k))));
            recophi->SetBinContent(k, recophi->GetBinContent(k)*(1+twophistrength*cos(2*recophi->GetBinCenter(k))+fourphistrength*(cos(4*recophi->GetBinCenter(k)))));
            //cout << "set reco phi using embedded" << "\n";
        }

        recophi->Fit("phifit2", "", "", -3.15,3.15);
        double recoA2 = phifit2->GetParameter(2);
        double recoA4 = phifit2->GetParameter(4);

        reco2phiReweighted->SetBinContent(i, recoA2);
        reco4phiReweighted->SetBinContent(i, recoA4);


    }

    makeCanvas();
    reco2phiReweighted->Draw("PE");

    makeCanvas();
    reco4phiReweighted->Draw("PE");

    TFile file("output_root_files/signal_embedded_results.root", "RECREATE");
    reco2phiReweighted->Write();
    reco4phiReweighted->Write();

}

