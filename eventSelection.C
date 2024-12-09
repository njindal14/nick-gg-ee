//#include "./FemtoPairFormat.h"


int ican = 0;
void makeCanvas() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican++ ), "", 900, 600 );
    can->SetTopMargin(0.1);
    can->SetRightMargin(0.1);
}

void eventSelection() {

    TFile *Run12FemtoDst = TFile::Open("/Users/Nick/STAR/breit-wheeler/rootFiles/FemtoDst_Run12UU_wZDC.root");

    //TFile *file2 = new TFile("/Users/Nick/STAR/breit-wheeler/rootFiles/pair_dst_Run12UU.root"); // Open the file
   // Navigate to the PairDst folder
    //TDirectory *pairDst = (TDirectory*)file2->Get("PairDst");
    //TH1F *dca = (TH1F*)pairDst->Get("Pairs.d1_mDCA");
    
    TH1F * event_stats = (TH1F*)Run12FemtoDst->Get("event_stats");
    TH2D * n_verts_GRefMult = (TH2D*)Run12FemtoDst->Get("n_verts_GRefMult");
    TH1F * hist_vz = (TH1F*)Run12FemtoDst->Get("vz");
    TH1F * vr = (TH1F*)Run12FemtoDst->Get("vr");
    TH1F * hist_ranking = (TH1F*)Run12FemtoDst->Get("ranking");
    TH1F * hist_delta_vz_wide = (TH1F*)Run12FemtoDst->Get("delta_vz_wide");
    TH1F * hist_RefMult = (TH1F*)Run12FemtoDst->Get("RefMult");
    TH1F * hist_gRefMult = (TH1F*)Run12FemtoDst->Get("gRefMult");

    makeCanvas();
    event_stats->GetYaxis()->SetTitle("N Events");
    event_stats->GetXaxis()->SetTitle("Run 12 U+U");
    event_stats->GetXaxis()->SetBinLabel(1,"All");
    event_stats->GetXaxis()->SetBinLabel(2, "Triggered");
    event_stats->GetXaxis()->SetBinLabel(3, "Vertex");
    event_stats->GetXaxis()->SetBinLabel(4, "V_{z}");
    event_stats->GetXaxis()->SetBinLabel(5, "Analysis");
    event_stats->GetXaxis()->SetRangeUser(0, 6);
    gPad->SetLogy();
    event_stats->Draw();
    gPad->Print("note_plots/event_selection_plots/Run12_event_stats.png");

    makeCanvas();
    gPad->SetLogz();
    n_verts_GRefMult->GetXaxis()->SetRangeUser(0, 700);
    n_verts_GRefMult->Draw("colz");
    gPad->Print("note_plots/event_selection_plots/Run12_n_verts_GRefMult.png");

    makeCanvas();
    gPad->SetLogy();
    hist_vz->Draw();
    gPad->Print("note_plots/event_selection_plots/Run12_vz.png");

    makeCanvas();
    gPad->SetLogy();
    vr->GetXaxis()->SetRangeUser(0,1.);
    vr->Draw();
    gPad->Print("note_plots/event_selection_plots/Run12_vr.png");

    makeCanvas();
    gPad->SetLogy();
    hist_ranking->Draw();
    gPad->Print("note_plots/event_selection_plots/Run12_ranking.png");

    makeCanvas();
    gPad->SetLogy();
    hist_delta_vz_wide->Draw();
    gPad->Print("note_plots/event_selection_plots/Run12_delta_vz_wide.png");


    makeCanvas();
    gPad->SetLogy();
    hist_RefMult->Draw();
    gPad->Print("note_plots/event_selection_plots/Run12_refmult.png");


    makeCanvas();
    gPad->SetLogy();
    hist_gRefMult->Draw();
    gPad->Print("note_plots/event_selection_plots/Run12_grefmult.png");

    TFile file("output_root_files/event_selection.root", "RECREATE");
    event_stats->Write();
    n_verts_GRefMult->Write();
    hist_vz->Write();
    vr->Write();
    hist_ranking->Write();
    hist_delta_vz_wide->Write();
    hist_RefMult->Write();
    hist_gRefMult->Write();

    double lumi_fraction = hist_vz->Integral(hist_vz->FindBin(-100), hist_vz->FindBin(100))/hist_vz->Integral();
    std::cout << "lumi fraction: " << lumi_fraction << "\n";

    double grefmult_eff =  hist_gRefMult->Integral(hist_gRefMult->FindBin(0), hist_gRefMult->FindBin(4))/hist_gRefMult->Integral();
    std::cout << "grefmult_eff: " << grefmult_eff << "\n";

    //get dca cut eff
    //TFile *Run12PairDst = TFile::Open("/Users/Nick/STAR/breit-wheeler/rootFiles/pair_dst_Run12UU.root");
    //TTree *tree = (TTree*)Run12PairDst->GetObject("PairDst");
    //TBranch *branch = (TBranch*)tree->GetBranch("Pairs");

    //TH1F * d1_dca = (TH1F*)Run12PairDst->Get("PairDst.Pairs.d1_mDCA");
    
    //double d1_dca_eff =  d1_dca->Integral(d1_dca->FindBin(0), d1_dca->FindBin(1))/d1_dca->Integral();
    //std::cout << "d1_dca_eff: " << d1_dca_eff << "\n";


}


