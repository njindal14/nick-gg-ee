#include "FemtoPairFormat.h"

double c = 3.0e1; //in cm/ns
double me2 = pow(0.00051,2);

int ican = 0;
void makeCanvas() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican++ ), "", 900, 600 );
    can->SetTopMargin(0.04);
    can->SetRightMargin(0.3);
}

double chiFit(double *x, double *par){
    double fitval;
    fitval = par[0]*exp(x[0]/par[1]);
    return fitval;
}

double ddToffit(double *x, double *par){
    double fitval;

    fitval = par[0]/(par[2]* sqrt(2*M_PI)) * exp(-0.5*pow((x[0]-par[1])/par[2], 2)) + 
            par[3]/(par[5]* sqrt(2*M_PI)) * exp(-0.5*pow((x[0]-par[4])/par[5], 2)) +
         + par[6];                   
    return fitval;
}

double calc_Phi( TLorentzVector lv1, TLorentzVector lv2) {
    TLorentzVector lvPlus = lv1 + lv2;
    TLorentzVector lvMinus = lv1 - lv2;
    lv1.Boost(-lvPlus.BoostVector());
    lv2.Boost(-lvPlus.BoostVector());
    double Px = lvPlus.Px();
    double Py = lvPlus.Py();
    double Qx = lvMinus.Px();
    double Qy = lvMinus.Py();
    double absPperp = pow((Px*Px)+(Py*Py), 0.5);
    double absQperp = pow((Qx*Qx)+(Qy*Qy), 0.5);
    double PcrossQ = (Px*Qy) - (Py*Qx);
    double PdotQ = (Px*Qx) + (Py*Qy);
    double cosphi = (Px*Qx + Py*Qy) / (absPperp*absQperp);
    double PairPhi = acos(cosphi);
    if ( PcrossQ > 0 ){
        return PairPhi - 3.141592;
    } else {
        return 3.141592 - PairPhi;
    }
}

void eePIDandPurity() { 

    auto * chieeFit = new TF1("chieefit", "[0]", 0, 30);
    chieeFit->SetParameter(0,40);

    auto * mdTof = new TH1F("#DeltaTOF Hist", "#DeltaTOF", 1000, -15, 15);
    auto * mdTofexp = new TH1F("#DeltaTOFExp Hist", "#DeltaTOFexp", 1000, -15, 15);
    auto * mddTof = new TH1F("#Delta#DeltaTOF Hist", "#Delta#DeltaTOF", 1000, -6, 6);
    auto * Xee = new TH1F("#chi_{ee}^{2} Distribution", "#chi_{ee}^{2} Distribution", 100, 0, 30);
    auto * Xee25 = new TH1F("#chi_{ee}^{2} Distribution", "#chi_{ee}^{2} Distribution, #chi_{#pi#pi}^{2} > 25", 200, 0, 15);
    auto * background25 = new TH1F("Background", "#chi_{ee}^{2} Distribution background, #chi_{#pi#pi}^{2} > 25", 200, 0, 15);
    auto * Xee20 = new TH1F("#chi_{ee}^{2} Distribution", "#chi_{ee}^{2} Distribution, #chi_{#pi#pi}^{2} > 20", 200, 0, 15);
    auto * background20 = new TH1F("Background", "#chi_{ee}^{2} Distribution background, #chi_{#pi#pi}^{2} > 20", 200, 0, 15);
    auto * Xee15 = new TH1F("#chi_{ee}^{2} Distribution", "#chi_{ee}^{2} Distribution, #chi_{#pi#pi}^{2} > 15", 200, 0, 15);
    auto * background15 = new TH1F("Background", "#chi_{ee}^{2} Distribution background, #chi_{#pi#pi}^{2} > 15", 200, 0, 15);
    auto * Xee10 = new TH1F("#chi_{ee}^{2} Distribution", "#chi_{ee}^{2} Distribution, #chi_{#pi#pi}^{2} > 10", 200, 0, 15);
    auto * background10 = new TH1F("Background", "#chi_{ee}^{2} Distribution background, #chi_{#pi#pi}^{2} > 10", 200, 0, 15);
    auto * Xee5 = new TH1F("#chi_{ee}^{2} Distribution", "#chi_{ee}^{2} Distribution, #chi_{#pi#pi}^{2} > 5", 200, 0, 15);
    auto * background5 = new TH1F("Background", "#chi_{ee}^{2} Distribution background, #chi_{#pi#pi}^{2} > 5", 200, 0, 15);
    auto * Xee1 = new TH1F("#chi_{ee}^{2} Distribution", "#chi_{ee}^{2} Distribution, #chi_{#pi#pi}^{2} > 1", 200, 0, 15);
    auto * background1 = new TH1F("Background", "#chi_{ee}^{2} Distribution background, #chi_{#pi#pi}^{2} > 1", 200, 0, 15);
    auto * ddTofFit = new TF1("fit", ddToffit, -2, 2, 7);
    ddTofFit->SetParameters(100000.0, 0, 0.2, 50000.0, 0, 0.5, 10);
    ddTofFit->SetParNames("A1", "#lambda1", "#sigma1", "A2","#lambda2", "sigma2", "p0");
    //set some parameter ranges
    ddTofFit->SetParLimits(2,0.05, 0.3);
    ddTofFit->SetParLimits(5,0.1,1);

    //smooth it out
    ddTofFit->SetNpx(1000);
    ddTofFit->SetLineWidth(4);
    
    TFile *myFile = TFile::Open("/Users/Nick/STAR/breit-wheeler/rootFiles/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");
    TLorentzVector lv1, lv2, lv, lvn;

    while (myReader.Next()) {

        double chiee = pow( pair->d1_mNSigmaElectron, 2 ) + pow( pair->d2_mNSigmaElectron, 2 );
        double chipipi = pow( pair -> d1_mNSigmaPion, 2) + pow( pair -> d2_mNSigmaPion, 2);
        int chargesumval = pair->mChargeSum;

        lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.00051 );
        lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.00051 ); 

        lv = lv1 + lv2;
        lvn = lv1 - lv2;

        Float_t p1_2 = pow(lv1.P(), 2);
        Float_t p2_2 = pow(lv2.P(), 2);
        Float_t mVertexZVal = pair->mVertexZ;
        UShort_t mGRefMultVal = pair->mGRefMult;  
        Float_t dTofVal = pair->d1_mTof - pair->d2_mTof;
        Float_t dTofexpVal = pair->d1_mLength/c * sqrt(1 + me2/p1_2) - pair->d2_mLength/c * sqrt(1 + me2/p2_2);
        Float_t ddTofVal = dTofVal - dTofexpVal;

        if( fabs(mVertexZVal) < 100 &&  mGRefMultVal <= 4 && chargesumval == 0 && pair->d1_mDCA < 1 && pair->d2_mDCA < 1 && 
        pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag!=0 && fabs(ddTofVal < 0.4)) {
            
            if(ddTofVal == 0) continue;

            if(chipipi > 30) Xee->Fill(chiee);

            else{
                if (chipipi > 25) {
                    if (3 * chiee < chipipi) {
                        Xee25->Fill(chiee);
                    } else {
                        background25->Fill(chiee);
                    }
                }

                if (chipipi > 20) {
                    if (3 * chiee < chipipi) {
                        Xee20->Fill(chiee);
                    } else {
                        background20->Fill(chiee);
                    }
                }

                if (chipipi > 15) {
                    if (3 * chiee < chipipi) {
                        Xee15->Fill(chiee);
                    } else {
                        background15->Fill(chiee);
                    }
                }

                if (chipipi > 10) {
                    if (3 * chiee < chipipi) {
                        Xee10->Fill(chiee);
                    } else {
                        background10->Fill(chiee);
                    }
                }

                if (chipipi > 5) {
                    if (3 * chiee < chipipi) {
                        Xee5->Fill(chiee);
                    } else {
                        background5->Fill(chiee);
                    }
                }

                if (chipipi > 1) {
                    if (3 * chiee < chipipi) {
                        Xee1->Fill(chiee);
                    } else {
                        background1->Fill(chiee);
                    }
                }
            }
        }
        if(dTofVal != 0 && dTofexpVal !=0 && ddTofVal !=0){
            mdTof->Fill( dTofVal );
            mdTofexp->Fill( dTofexpVal );
            mddTof->Fill( ddTofVal );
        }
    }

    makeCanvas();
    mdTof->SetLineColor(kBlack);
    gPad->SetLogy();
    mdTof->GetXaxis()->SetTitle("#Delta TOF Distrubition (ns)");
    mdTof->GetYaxis()->SetTitle("Counts");
    mdTof->Draw();
    
    mdTofexp->SetLineColor(kRed);
    gPad->SetLogy();
    mdTofexp->Draw("same");

    auto * legend = new TLegend(0.77,0.6,.97,0.75);
    legend->SetHeader("Legend");
    legend->AddEntry(mdTof,"#DeltaTOF","l");
    legend->AddEntry(mdTofexp,"#DeltaTOFexp","l");
    legend->Draw("same");

    
    makeCanvas();
    mddTof->SetLineColor(kBlack);
    gPad->SetLogy();
    mddTof->GetXaxis()->SetTitle("#Delta #Delta TOF Distrubition (ns)");
    mddTof->GetYaxis()->SetTitle("Counts");
    mddTof->Draw();

    mddTof->Fit("fit", "", "", -2,2);
    gStyle->SetOptFit(1111);
    //gPad->Print( "plots/plot_ddTof.png"); 
    double errorddtof;
    double ddtofIntegral = mddTof->IntegralAndError(-34,34,errorddtof,"");
    std::cout << "Integral of ddtof inside cut: " << ddtofIntegral;
    std::cout << "Error of integral: " << errorddtof;   


    
    makeCanvas();
    Xee->SetLineColor(kBlack);
    gPad->SetLogy();
    Xee->GetXaxis()->SetTitle("#chi_{ee}^{2}");
    Xee->GetYaxis()->SetTitle("dN/d(#chi_{ee}^{2})");
    Xee->SetMarkerStyle(20);
    Xee->Draw("PE");
    chieeFit->SetLineWidth(4);
    chieeFit->SetLineColor(kBlue);
    Xee->Fit("chieefit", "", "", 12, 30);
    Xee->Fit("expo", "R+", "", 0., 12.);
    gStyle->SetOptFit(1111);
    chieeFit->Draw("same");
    double background0 = chieeFit->GetParameter(0);
    double binsx0 = Xee->GetNbinsX();
    double totbackground = 0;
    double totsignal = 0;
    for(int ix =1; ix <= Xee->FindFixBin(10); ix++){
        totsignal += Xee->GetBinContent(ix);
        totbackground += background0;
    }
    double purity = totsignal/ (totsignal+totbackground);
    cout << "Purity for chipipi > 30: " << purity*100 << "%\n";
    //gPad->Print("plots/chi2eePlot.png");

    
    makeCanvas();
    Xee25->SetLineColor(kBlack);
    gPad->SetLogy();
    Xee25->GetXaxis()->SetTitle("#chi_{ee}^{2}");
    Xee25->GetYaxis()->SetTitle("dN/d(#chi_{ee}^{2})");
    Xee25->Draw();
    background25->SetLineColor(kRed);
    background25->Draw("same");
    auto * legend25 = new TLegend(0.77,0.5,.97,0.65);
    legend25->SetHeader("Legend");
    legend25->AddEntry(Xee25,"Signal","l");
    legend25->AddEntry(background25,"Background","l");
    legend25->Draw("same");
    TLine * l = new TLine(8.33,0,8.33,60);
    l->SetLineColor(kGreen);
    l->Draw("same");
    int xbins25 = Xee25->GetNbinsX();
    double sig25 = 0;
    double bground25 = 0;
    for(int ix = 1; ix <= xbins25; ix++){
        if(Xee25->GetBinCenter(ix) <= 8.33){
            sig25 += Xee25->GetBinContent(ix);
            bground25 += background25->GetBinContent(ix);
        }
    }
    cout << "Signal 25: " << sig25 << "\n";
    cout << "Background 25: " << bground25 << "\n";
    double purity25 = sig25/(sig25 + bground25);
    cout << "Purity25 : " << purity25 << "\n";
    double sigint25 = Xee25->Integral(Xee25->FindFixBin(0), Xee25->FindFixBin(8.33), "");
    double bkg25 = background25->Integral(background25->FindFixBin(0), background25->FindFixBin(8.33), "");
    double purityIntegral25 = sigint25/(sigint25+bkg25);
    cout << "Purity 25 Integral Method: " << purityIntegral25;
    //gPad->Print("plots/plot_chi25.png");


    
    makeCanvas();
    Xee20->SetLineColor(kBlack);
    gPad->SetLogy();
    Xee20->GetXaxis()->SetTitle("#chi_{ee}^{2}");
    Xee20->GetYaxis()->SetTitle("dN/d(#chi_{ee}^{2})");
    Xee20->Draw();
    background20->SetLineColor(kRed);
    background20->Draw("same");
    auto * legend20 = new TLegend(0.77,0.5,.97,0.65);
    legend20->SetHeader("Legend");
    legend20->AddEntry(Xee20,"Signal","l");
    legend20->AddEntry(background20,"Background","l");
    legend20->Draw("same");
    TLine * l2 = new TLine(6.66,0,6.66,60);
    l2->SetLineColor(kGreen);
    l2->Draw("same");
    
    int xbins20 = Xee20->GetNbinsX();
    double sig20 = 0;
    double bground20 = 0;
    for(int ix = 1; ix <= xbins20; ix++){
        if(Xee20->GetBinCenter(ix) <= 6.66){
            sig20 += Xee20->GetBinContent(ix);
            bground20 += background20->GetBinContent(ix);
        }
    }
    cout << "Signal 20: " << sig20 << "\n";
    cout << "Background 20: " << bground20 << "\n";
    double purity20 = sig20/(sig20 + bground20);
    cout << "Purity20 : " << purity20 << "\n";
    double sigint20 = Xee20->Integral(Xee20->FindFixBin(0), Xee20->FindFixBin(6.66), "");
    double bkg20 = background20->Integral(background20->FindFixBin(0), background20->FindFixBin(6.66), "");
    double purityIntegral20 = sigint20/(sigint20+bkg20);
    cout << "Purity 20 Integral Method: " << purityIntegral20;
    //gPad->Print("plots/plot_chi20.png");


    
    makeCanvas();
    Xee15->SetLineColor(kBlack);
    gPad->SetLogy();
    Xee15->GetXaxis()->SetTitle("#chi_{ee}^{2}");
    Xee15->GetYaxis()->SetTitle("dN/d(#chi_{ee}^{2})");
    Xee15->Draw();
    background15->SetLineColor(kRed);
    background15->Draw("same");
    auto * legend15 = new TLegend(0.77,0.5,.97,0.65);
    legend15->SetHeader("Legend");
    legend15->AddEntry(Xee15,"Signal","l");
    legend15->AddEntry(background15,"Background","l");
    legend15->Draw("same");
    TLine * l3 = new TLine(5,0,5,60);
    l3->SetLineColor(kGreen);
    l3->Draw("same");
    int xbins15 = Xee15->GetNbinsX();
    double sig15 = 0;
    double bground15 = 0;
    for(int ix = 1; ix <= xbins15; ix++){
        if(Xee15->GetBinCenter(ix) <= 5){
            sig15 += Xee15->GetBinContent(ix);
            bground15 += background15->GetBinContent(ix);
        }
    }
    cout << "Signal 15: " << sig15 << "\n";
    cout << "Background 15: " << bground15 << "\n";
    double purity15 = sig15/(sig15 + bground15);
    cout << "Purity15 : " << purity15 << "\n";
    double sigint15 = Xee15->Integral(Xee15->FindFixBin(0), Xee15->FindFixBin(5), "");
    double bkg15 = background15->Integral(background15->FindFixBin(0), background15->FindFixBin(5), "");
    double purityIntegral15 = sigint15/(sigint15+bkg15);
    cout << "Purity 15 Integral Method: " << purityIntegral15;
    //gPad->Print("plots/plot_chi15.png");


    
    makeCanvas();
    Xee10->SetLineColor(kBlack);
    gPad->SetLogy();
    Xee10->GetXaxis()->SetTitle("#chi_{ee}^{2}");
    Xee10->GetYaxis()->SetTitle("dN/d(#chi_{ee}^{2})");
    Xee10->Draw();
    background10->SetLineColor(kRed);
    background10->Draw("same");
    auto * legend10 = new TLegend(0.77,0.5,.97,0.65);
    legend10->SetHeader("Legend");
    legend10->AddEntry(Xee10,"Signal","l");
    legend10->AddEntry(background10,"Background","l");
    legend10->Draw("same");
    TLine * l4 = new TLine(3.33,0,3.33,60);
    l4->SetLineColor(kGreen);
    l4->Draw("same");
    int xbins10 = Xee10->GetNbinsX();
    double sig10 = 0;
    double bground10 = 0;
    for(int ix = 1; ix <= xbins10; ix++){
        if(Xee10->GetBinCenter(ix) <= 3.33){
            sig10 += Xee10->GetBinContent(ix);
            bground10 += background10->GetBinContent(ix);
        }
    }
    cout << "Signal 10: " << sig10 << "\n";
    cout << "Background 10: " << bground10 << "\n";
    double purity10 = sig10/(sig10 + bground10);
    cout << "Purity10 : " << purity10 << "\n";
    double sigint10 = Xee10->Integral(Xee10->FindFixBin(0), Xee10->FindFixBin(3.33), "");
    double bkg10 = background10->Integral(background10->FindFixBin(0), background10->FindFixBin(3.33), "");
    double purityIntegral10 = sigint10/(sigint10+bkg10);
    cout << "Purity 10 Integral Method: " << purityIntegral10;
    //gPad->Print("plots/plot_chi10.png");



    
    makeCanvas();
    Xee5->SetLineColor(kBlack);
    gPad->SetLogy();
    Xee5->GetXaxis()->SetTitle("#chi_{ee}^{2}");
    Xee5->GetYaxis()->SetTitle("dN/d(#chi_{ee}^{2})");
    Xee5->Draw();
    background5->SetLineColor(kRed);
    background5->Draw("same");
    auto * legend5 = new TLegend(0.77,0.5,.97,0.65);
    legend5->SetHeader("Legend");
    legend5->AddEntry(Xee5,"Signal","l");
    legend5->AddEntry(background5,"Background","l");
    legend5->Draw("same");
    TLine * l5 = new TLine(1.66,0,1.66,60);
    l5->SetLineColor(kGreen);
    l5->Draw("same");
    int xbins5 = Xee5->GetNbinsX();
    double sig5 = 0;
    double bground5 = 0;
    for(int ix = 1; ix <= xbins5; ix++){
        if(Xee5->GetBinCenter(ix) <= 1.66){
            sig5 += Xee5->GetBinContent(ix);
            bground5 += background5->GetBinContent(ix);
        }
    }
    cout << "Signal 5: " << sig5 << "\n";
    cout << "Background 5: " << bground5 << "\n";
    double purity5 = sig5/(sig5 + bground5);
    cout << "Purity5 : " << purity5 << "\n";
    double sigint5 = Xee5->Integral(Xee5->FindFixBin(0), Xee5->FindFixBin(1.66), "");
    double bkg5 = background5->Integral(background5->FindFixBin(0), background5->FindFixBin(1.66), "");
    double purityIntegral5 = sigint5/(sigint5+bkg5);
    cout << "Purity 5 Integral Method: " << purityIntegral5;
    //gPad->Print("plots/plot_chi5.png");


    
    makeCanvas();
    Xee1->SetLineColor(kBlack);
    gPad->SetLogy();
    Xee1->GetXaxis()->SetTitle("#chi_{ee}^{2}");
    Xee1->GetYaxis()->SetTitle("dN/d(#chi_{ee}^{2})");
    Xee1->Draw();
    background1->SetLineColor(kRed);
    background1->Draw("same");
    auto * legend1 = new TLegend(0.77,0.5,.97,0.65);
    legend1->SetHeader("Legend");
    legend1->AddEntry(Xee1,"Signal","l");
    legend1->AddEntry(background1,"Background","l");
    legend1->Draw("same");
    TLine * l6 = new TLine(.33,0,.33,60);
    l6->SetLineColor(kGreen);
    l6->Draw("same");
    int xbins1 = Xee1->GetNbinsX();
    double sig1 = 0;
    double bground1 = 0;
    for(int ix = 1; ix <= xbins1; ix++){
        if(Xee1->GetBinCenter(ix) <= .33){
            sig1 += Xee1->GetBinContent(ix);
            bground1 += background1->GetBinContent(ix);
        }
    }
    cout << "Signal 1: " << sig1 << "\n";
    cout << "Background 1: " << bground1 << "\n";
    double purity1 = sig1/(sig1 + bground1);
    cout << "Purity1 : " << purity1 << "\n";
    double sigint1 = Xee1->Integral(Xee1->FindFixBin(0), Xee1->FindFixBin(.33), "");
    double bkg1 = background1->Integral(background1->FindFixBin(0), background1->FindFixBin(.33), "");
    double purityIntegral1 = sigint1/(sigint1+bkg1);
    cout << "Purity 1 Integral Method: " << purityIntegral1;

    //gPad->Print("plots/plot_chi1.png");

    double overallPurity = 0;
    overallPurity += purity1 + purity10 + purity15 + purity5 + purity20 + purity25;
    overallPurity /= 6;
    cout << "AVERAGE SLICE PURITY: " << overallPurity*100 << "%\n";

    double overallPurityIntegral = 0;
    overallPurityIntegral+= purityIntegral1 + purityIntegral5 + purityIntegral10 + purityIntegral15 + purityIntegral20 + purityIntegral25;
    overallPurityIntegral/=6;
    cout << "AVERAGE SLICE PURITY INTEGRAL METHOD: " << overallPurity*100 << "%\n";

    
    TFile file("PurityAndTOFInfo.root", "RECREATE");
    mdTof->Write();
    mdTofexp->Write(); 
    mddTof->Write(); 
    Xee->Write(); 
    Xee25->Write();
    background25->Write();
    Xee20->Write(); 
    background20->Write(); 
    Xee15->Write();
    background15->Write(); 
    Xee10->Write(); 
    background10->Write(); 
    Xee5->Write(); 
    background5->Write();
    Xee1->Write(); 
    background1->Write();
    
    
    


}
      



