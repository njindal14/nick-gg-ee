#include "FemtoPairFormat.h"
double const pi = 3.14159265358979323;
double c = 3.0e1; //in cm/ns
double me2 = pow(0.00051,2);

int ican = 0;
void makeCanvas() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican++ ), "", 900, 600 );
    can->SetTopMargin(0.08);
    can->SetRightMargin(0.3);
}

void TPCCheck() {
    
    TFile *myFile = TFile::Open("/Users/Nick/STAR/breit-wheeler/rootFiles/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");
    TLorentzVector lv1, lv2, lv, lvn;

    TFile * simulation_plots_file = new TFile("/Users/Nick/STAR/breit-wheeler/Code/nick-gg-ee/simulation_plots.root");
    //now import simulation hists here from simulation to compare with data

    TH2D * STARpairPhiVsEta = new TH2D("Pair #phi", "Pair #phi; #eta; #phi", 50, -1, 1, 50, -pi, pi);
    TH2D * STARposPhiVsEta = new TH2D("Positron #phi", "Positron #phi; #eta; #phi", 50, -1, 1, 50, -pi, pi);
    TH2D * STARelPhiVsEta = new TH2D("Electron #phi", "Electron #phi; #eta; #phi", 50, -1, 1, 50, -pi, pi);

   
    while (myReader.Next()) {
        double chiee = pow( pair->d1_mNSigmaElectron, 2 ) + pow( pair->d2_mNSigmaElectron, 2 );
        double chipipi = pow( pair -> d1_mNSigmaPion, 2) + pow( pair -> d2_mNSigmaPion, 2);
        int chargesumval = pair->mChargeSum;

        lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.00051 );
        lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.00051 ); 

        lv = lv1 + lv2;
        lvn = lv1 - lv2;



        double posPhi = lv1.Phi();
        double elPhi = lv2.Phi();
        double pairPhi = lv.Phi();

        Float_t p1_2 = pow(lv1.P(), 2);
        Float_t p2_2 = pow(lv2.P(), 2);
        Float_t mVertexZVal = pair->mVertexZ;
        UShort_t mGRefMultVal = pair->mGRefMult;  
        Float_t dTofVal = pair->d1_mTof - pair->d2_mTof;
        Float_t dTofexpVal = pair->d1_mLength/c * sqrt(1 + me2/p1_2) - pair->d2_mLength/c * sqrt(1 + me2/p2_2);
        Float_t ddTofVal = dTofVal - dTofexpVal;

        if( fabs(mVertexZVal) < 100 &&  mGRefMultVal <= 4 && chargesumval == 0 && pair->d1_mDCA < 1 && pair->d2_mDCA < 1 && 
        pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag!=0 && fabs(ddTofVal < 0.4) && ddTofVal !=0) {
            if(chiee < 10 && 3*chiee < chipipi){
                
                STARpairPhiVsEta->Fill(lv.Eta(), pairPhi);
                STARposPhiVsEta->Fill(lv1.Eta(), posPhi);
                STARelPhiVsEta->Fill(lv2.Eta(), elPhi);
            }
               
        }

    }

    TH1D *pairNegEta = STARpairPhiVsEta->ProjectionY("Pair #phi, #eta < 0", 0, STARpairPhiVsEta->GetXaxis()->FindBin(0.));
    TH1D *pairPosEta = STARpairPhiVsEta->ProjectionY("Pair #phi, #eta > 0", STARpairPhiVsEta->GetXaxis()->FindBin(0.), -1);
    TH1D *pairAllEta = STARpairPhiVsEta->ProjectionY("Pair #phi, all #eta", 0, -1);


    TH1D *posNegEta = STARposPhiVsEta->ProjectionY("Positron #phi, #eta < 0", 0, STARposPhiVsEta->GetXaxis()->FindBin(0.) );
    TH1D *posPosEta = STARposPhiVsEta->ProjectionY("Positron #phi, #eta > 0", STARposPhiVsEta->GetXaxis()->FindBin(0.), -1 );
    TH1D *posAllEta = STARposPhiVsEta->ProjectionY("Positron #phi, all #eta", 0, -1);

    
    TH1D *elNegEta = STARelPhiVsEta->ProjectionY("Electron #phi, #eta < 0", 0, STARelPhiVsEta->GetXaxis()->FindBin(0.) );
    TH1D *elPosEta = STARelPhiVsEta->ProjectionY("Electron #phi, #eta > 0", STARelPhiVsEta->GetXaxis()->FindBin(0.), -1 );
    TH1D *elAllEta = STARelPhiVsEta->ProjectionY("Electron #phi, all #eta", 0, -1 );



    makeCanvas();
    STARpairPhiVsEta->Draw("colz");

    makeCanvas();
    STARposPhiVsEta->Draw("colz");

    makeCanvas();
    STARelPhiVsEta->Draw("colz");

    
    makeCanvas();
    posAllEta->SetLineColor(kBlack);
    posNegEta->SetLineColor(kRed);
    posPosEta->SetLineColor(kBlue);
    posAllEta->GetXaxis()->SetTitle("STAR #phi");
    posAllEta->GetYaxis()->SetTitle("Counts"); 
    posAllEta->GetYaxis()->SetRangeUser(0,450);
    posAllEta->Draw("PE");
    posNegEta->Draw("same;PE");
    posPosEta->Draw("same;PE");
    auto * posLegend = new TLegend(0.77,0.5,.97,0.65);
    posLegend->SetHeader("Legend");
    posLegend->AddEntry(posAllEta,"Positron all rapidity","l");
    posLegend->AddEntry(posPosEta,"Positron #eta > 0","l");
    posLegend->AddEntry(posNegEta,"Positron #eta < 0","l");
    posLegend->Draw("same");

    makeCanvas();
    elAllEta->SetLineColor(kBlack);
    elNegEta->SetLineColor(kRed);
    elPosEta->SetLineColor(kBlue);
    elAllEta->GetXaxis()->SetTitle("STAR #phi");
    elAllEta->GetYaxis()->SetTitle("Counts"); 
    elAllEta->GetYaxis()->SetRangeUser(0,450);
    elAllEta->Draw("PE");
    elNegEta->Draw("same;PE");
    elPosEta->Draw("same;PE");
    auto * elLegend = new TLegend(0.77,0.5,.97,0.65);
    elLegend->SetHeader("Legend");
    elLegend->AddEntry(elAllEta,"Electron all rapidity","l");
    elLegend->AddEntry(elPosEta,"Electron #eta > 0","l");
    elLegend->AddEntry(elNegEta,"Electron #eta < 0","l");
    elLegend->Draw("same");

    makeCanvas();
    pairAllEta->SetLineColor(kBlack);
    pairNegEta->SetLineColor(kRed);
    pairPosEta->SetLineColor(kBlue);
    pairAllEta->GetXaxis()->SetTitle("STAR #phi");
    pairAllEta->GetYaxis()->SetTitle("Counts"); 
    pairAllEta->GetYaxis()->SetRangeUser(0,450);
    pairAllEta->Draw("PE");
    pairNegEta->Draw("same;PE");
    pairPosEta->Draw("same;PE");
    auto * pairLegend = new TLegend(0.77,0.5,.97,0.65);
    pairLegend->SetHeader("Legend");
    pairLegend->AddEntry(pairAllEta,"Pair all rapidity","l");
    pairLegend->AddEntry(pairPosEta,"Pair #eta > 0","l");
    pairLegend->AddEntry(pairNegEta,"Pair #eta < 0","l");
    pairLegend->Draw("same");

    TFile file("tpc_check_plots.root", "RECREATE");
    STARpairPhiVsEta->Write();
    STARposPhiVsEta->Write();
    STARelPhiVsEta->Write();
    pairNegEta->Write();
    pairPosEta->Write();
    posNegEta->Write();
    posPosEta->Write();
    elNegEta->Write();
    elPosEta->Write();
    posAllEta->Write();
    elAllEta->Write();
    pairAllEta->Write();




}