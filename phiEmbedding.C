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

    TFile * sim_out = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/simulation_plots_new.root");
    TH2F * mcphipt = (TH2F*)sim_out->Get("mcphipt");
    TH2F * recophipt = (TH2F*)sim_out->Get("recophipt");

    TFile * theory_curves = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/theory_modulations.root");
    TH1F * qed2phi = (TH1F*)theory_curves->Get("QED2phi");
    TH1F * qed4phi = (TH1F*)theory_curves->Get("QED4phi");

    int nbinsmc = mcphipt->GetNbinsY();

    TH1F * reco2phiReweighted = new TH1F("reco2phiReweighted","", nbinsmc, 0, 0.3);
    TH1F * reco4phiReweighted = new TH1F("reco4phiReweighted","", nbinsmc, 0, 0.3);



    for(int i = 0; i <= 9; ++i){
     
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

    double fourphix[] = {0.00023,0.00211,0.00491,0.00749,0.0101	,0.0131	,0.0164	,0.0201	,0.0239	,0.0304	,0.0356	,0.0416	,0.047	,0.0538	,0.0618	,0.0674	,0.0727	,0.0791	,0.084	,0.0877	,0.0912	,0.0943	,0.098	,0.101	,0.105	,0.108	,0.111	,0.114	,0.117	,0.119	,0.123	,0.126	,0.129	,0.132	,0.137	,0.14	,0.144	,0.149	,0.155	,0.161	,0.168	,0.174	,0.179	,0.184	,0.188	,0.192	,0.196	,0.199	,0.2};
    double fourphiy[] = {0.092-.1,0.068-.1	,0.028-.1	,-0.008-.1,-0.048-.1,-0.06-.1	,-0.06-.1	,-0.052-.1,-0.036-.1,-0.008-.1,0.024-.1	,0.06-.1	,0.104-.1	,0.156-.1	,0.232-.1	,0.288-.1	,0.356	-.1,0.448-.1	,0.528-.1	,0.596-.1	,0.672-.1	,0.744-.1	,0.832-.1	,0.916-.1	,1-.1	,1.08-.1	,1.17-.1	,1.24-.1	,1.3-.1	,1.35-.1	,1.41-.1	,1.47-.1	,1.52-.1	,1.55-.1	,1.6-.1	,1.63-.1	,1.67-.1	,1.7-.1	,1.73-.1	,1.75-.1	,1.76-.1	,1.76-.1	,1.78-.1	,1.79-.1	,1.81-.1	,1.83-.1	,1.85-.1	,1.88-.1	,1.88-.1};
    TGraph * QED4phi = new TGraph(sizeof(fourphix)/sizeof(fourphix[0]), fourphix,fourphiy);

    double twophix[] = {0.00384,0.00565,0.00701,0.00949,0.012,0.0145,0.0172,0.0201,0.0237,0.0271,0.0303,0.0346,0.0384,0.0416,0.045,0.0477,0.0499,0.0524,0.0556,0.059,0.0626,0.0653,0.0676,0.0705,0.0739,0.0764,0.0784,0.0807,0.0832,0.0856,0.0881,0.0908,0.0929,0.0949,0.0976,0.1,0.103,0.105,0.108,0.111,0.114,0.117,0.119,0.121,0.123,0.126,0.129,0.132	,0.136	,0.14,0.145,0.149,0.151,0.155,0.158,0.163,0.168,0.172,0.177,0.182,0.186,0.189,0.193,0.196,0.197};
    double twophiy[] = {0.00765,0.0115,0.0115,0.0153,0.0191,0.0229,0.0344	,0.0421,0.0535,0.065,0.0765,0.0918	,0.115	,0.13	,0.145	,0.164	,0.187	,0.203	,0.241	,0.268	,0.306	,0.333	,0.367	,0.417	,0.467	,0.512	,0.547	,0.593	,0.646	,0.704	,0.78	,0.837	,0.891	,0.956	,1.02	,1.12	,1.19	,1.26	,1.32	,1.41	,1.5	,1.56	,1.6	,1.64	,1.68	,1.72	,1.76	,1.78	,1.81	,1.82	,1.83	,1.83	,1.82	,1.82	,1.82	,1.82	,1.82	,1.83	,1.84	,1.85	,1.86	,1.87	,1.88	,1.9	,1.9};
    TGraph * QED2phi = new TGraph(sizeof(twophix)/sizeof(twophix[0]), twophix,twophiy);



    makeCanvas();
    QED2phi->GetXaxis()->SetTitle("pT (GeV/c)");
    QED2phi->GetYaxis()->SetTitle("A_{2#phi}");
    QED2phi->SetTitle("A_{2#phi} Theory vs Reweighted");
    QED2phi->GetYaxis()->SetRangeUser(-0.5, 2);
    QED2phi->SetLineColor(kBlue);
    QED2phi->SetLineWidth(6);
    QED2phi->Draw();
    reco2phiReweighted->Draw("PE;same");


    makeCanvas();
    QED4phi->GetXaxis()->SetTitle("pT (GeV/c)");
    QED4phi->GetYaxis()->SetTitle("A_{4#phi}");
    QED4phi->SetTitle("A_{4#phi} Theory vs Reweighted");
    QED4phi->GetYaxis()->SetRangeUser(-1, 2);
    QED4phi->SetLineColor(kBlue);
    QED4phi->SetLineWidth(6);
    QED4phi->Draw();
    reco4phiReweighted->Draw("PE;same");

    TFile file("output_root_files/signal_embedded_results.root", "RECREATE");
    reco2phiReweighted->Write();
    reco4phiReweighted->Write();

}

