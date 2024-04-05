



double phiFit(double *x, double *par){
    double fitval;
    fitval = par[0]*(1+par[1]*cos(x[0]) + par[2]*cos(2*x[0]) + par[3]*cos(3*x[0]) + par[4]*cos(4*x[0]));
    return fitval;
}


void daniel_test(){


    TFile * f = new TFile( "/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/simulation_plots.root" );


    TH2 * h2 = (TH2*)f->Get( "mcphipt" );


    for ( int i = 1; i < h2->GetNbinsY(); i++ ){

        auto * phifit = new TF1("phifit", phiFit, -3.14,3.14,5);
        TH1 * h = (TH1*)(h2->ProjectionX( "h", i, i+1 )->Clone( TString::Format( "h_%d", i ) )) ;
        cout << "made projections" << "\n";
        //normalize
        h->Scale(1/h->Integral());
        // int nbins = mcphi->GetNbinsX();
        // TH1F * flattened = new TH1F("flattened", "flattened", nbins, -3.1415927, 3.1415927);
        
        h->Fit("phifit", "", "", -3.14,3.14);

        
        auto c = new TCanvas(TString::Format("c_%d", i), "", 400, 400);
        h->Draw();
    }
    
}