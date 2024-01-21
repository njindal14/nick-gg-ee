#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "FemtoPairFormat.h"
#include <vector>
#include "TRandom3.h"
#include "TH3F.h"
#include "vector"
#include "TChain.h"
#include "TTreeReaderArray.h"
#include <iostream>
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TLegend.h"


int ican2 = 0;
void makeCanvas()  {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican2++ ), "", 900, 600);
    can->SetTopMargin(0.08);
    can->SetRightMargin(0.15);
}

double calc_Phi( TLorentzVector lv1, TLorentzVector lv2) {
    TLorentzVector lvPlus = lv1 + lv2;
    lv1.Boost(-lvPlus.BoostVector());
    lv2.Boost(-lvPlus.BoostVector());
    TLorentzVector lvMinus = lv1 - lv2;
    double Px = lvPlus.Px();
    double Py = lvPlus.Py();
    double Qx = lvMinus.Px();
    double Qy = lvMinus.Py();
    double PcrossQ = (Px*Qy) - (Py*Qx);
    double cosphi = (Px*Qx + Py*Qy) / (lvPlus.Pt()*lvMinus.Pt());
    double PairPhi = acos(cosphi);
    if ( PcrossQ > 0 ){
        return PairPhi - 3.141592;
    } else {
        return 3.141592 - PairPhi;
    }
}


double ddTOF ( double mass, double p1, double tof1, double l1,  double p2, double tof2, double l2 ) {

    double c = 3.0e1; //in cm/ns
    double me2 = pow(mass,2);
    double dTOF = tof2 - tof1 ;
    double tof1exp = (l1 / c) * sqrt( 1 + ( (mass*mass)/(p1*p1) ) );
    double tof2exp = (l2 / c) * sqrt( 1 + ( (mass*mass)/(p2*p2) ) );
    double dTOFexp = tof2exp - tof1exp ;
    double ddTOF = dTOF - dTOFexp ;

    return ddTOF;
}

void simAnalysisAndEfficiencies() {

    auto mPtPairReco = new TH1F("mPtPairReco","mPtPairReco", 50, 0, 0.4);
    auto McPtPair = new TH1F("McPtPair", "McPtPair", 50, 0, 0.4);

    auto mMassPairReco = new TH1F("mMassPairReco","mMassPairReco", 50, 0, 3);
    auto McMassPair = new TH1F("McMassPair", "McMassPair", 50, 0, 3);


    auto mEtaPairReco = new TH1F("mEtaPairReco","mEtaPairReco", 100, -1, 1);
    auto McEtaPair = new TH1F("McEtaPair", "McEtaPair", 100, -1, 1);

    auto cos2phivPtMC = new TH2F("Cos2#phivPtMC", "cos2#phi Moments vs P_{T} MC Pairs; 2<cos(2#phi)>; pT (GeV/c)", 100, -2, 2, 15, 0, 0.2);
    auto cos4phivPtMC = new TH2F("Cos4#phivPtMC", "cos4#phi Moments vs P_{T} MC Pairs; 2<cos(4#phi)>; pT (GeV/c)", 100, -2, 2, 15, 0, 0.2);

    auto cos2phivPtReco = new TH2F("Cos2#phivPtReco", "cos2#phi Moments vs P_{T} Reco Pairs; 2<cos(2#phi)>; pT (GeV/c)", 100, -2, 2, 15, 0, 0.2);
    auto cos4phivPtReco = new TH2F("Cos4#phivPtReco", "cos4#phi Moments vs P_{T} Reco Pairs; 2<cos(4#phi)>; pT (GeV/c)", 100, -2, 2, 15, 0, 0.2);




    TChain * ch = new TChain("FemtoDst");
    ch->Add("/Users/Nick/STAR/docker_mount/JobResults/*");
    TTreeReader myReader(ch);

    /*TFile *myFile = TFile::Open("/Users/Nick/STAR/docker_mount/backup/testFemtoDst.root");
    TTreeReader myReader("FemtoDst", myFile);*/



    TTreeReaderArray<float> PtVals(myReader, "Tracks.mPt");
    TTreeReaderArray<float> MCPtVals(myReader, "McTracks.mPt");
    TTreeReaderArray<float> MCEtaVals(myReader, "McTracks.mEta");
    TTreeReaderArray<float> MCPhiVals(myReader, "McTracks.mPhi");


    TTreeReaderArray<float> EtaVals(myReader, "Tracks.mEta");
    TTreeReaderArray<float> PhiVals(myReader, "Tracks.mPhi");
    TTreeReaderArray<unsigned short> PIDs(myReader, "McTracks.mGeantPID");
    TTreeReaderArray<short> MCIdx(myReader, "Tracks.mMcIndex");
    TTreeReaderArray<short> ParentIndex(myReader, "McTracks.mParentIndex");

    int mc_pair_counter = 0;
    int reco_pair_counter = 0;

    while (myReader.Next()) {
        
        TLorentzVector lvp, lvm, lv, Mclvp, Mclvm, Mclv;

        bool mcPair = false;
        bool recoPair = false;

        if(ParentIndex.GetSize() > 1 ){ //check for two tracks
            mcPair = true; //we have an mc pair
            for(int j = 0; j < 2; j++){ //loop through parent indices first 2 tracks for an event
                if(ParentIndex[j] == -1){ //make sure parent index is -1
                    if(PIDs[j] == 2){
                        Mclvp.SetPtEtaPhiM(MCPtVals[j], MCEtaVals[j], MCPhiVals[j], .00051); //set tlorentz vector for positron
                    }
                    if(PIDs[j] == 3){
                        Mclvm.SetPtEtaPhiM(MCPtVals[j], MCEtaVals[j], MCPhiVals[j], .00051); //set tlorentz vector for electron
                    }
                }
            }

            if(mcPair == true) { //check if we have an MCpairm fill histogram
                Mclv = Mclvm + Mclvp;
                McPtPair->Fill(Mclv.Pt());
                McEtaPair->Fill(Mclv.Rapidity());
                McMassPair->Fill(Mclv.M());
                double phival = calc_Phi(Mclvp, Mclvm);
                cos2phivPtMC->Fill(2*cos(2*phival), Mclv.Pt());
                cos4phivPtMC->Fill(2*cos(4*phival), Mclv.Pt());
                mc_pair_counter +=1;
            

                //if we do have an mc pair, now check for reco
                if(MCIdx.GetSize() > 1){ //check if we have more than one MCInx
                    int mcidx_counter = 0;
                    for(int i = 0; i < 2; i++){
                        if(MCIdx[i] != 0 && recoPair != true){
                            if(PIDs[i] == 2){
                                lvp.SetPtEtaPhiM(MCPtVals[i], MCEtaVals[i], MCPhiVals[i], .00051);
                                mcidx_counter +=1;

                            }
                            if(PIDs[i] == 3){
                                lvm.SetPtEtaPhiM(MCPtVals[i], MCEtaVals[i], MCPhiVals[i], .00051);
                                mcidx_counter +=1;

                            }
                        }
                        if(mcidx_counter == 2){
                            recoPair = true;
                            //std::cout << "Reconstructed \n";
                        }
                        //else std::cout << "did not reconstruct \n";
                    }
        
                }
            
                if(recoPair == true){
                    lv = lvp + lvm;
                    mPtPairReco->Fill(lv.Pt());
                    mEtaPairReco->Fill(lv.Rapidity());
                    mMassPairReco->Fill(lv.M());
                    double phival2 = calc_Phi(lvp, lvm);
                    cos2phivPtReco->Fill(2*cos(2*phival2), lv.Pt());
                    cos4phivPtReco->Fill(2*cos(4*phival2), lv.Pt()); 
                    reco_pair_counter  += 1;               
                }
            }
        }
    }

    makeCanvas();
    McPtPair->Draw("PE");
    mPtPairReco->SetLineColor(kRed);
    mPtPairReco->Draw("PE;same");
    //gPad->SetLogy();
    auto ptlegend = new TLegend(0.85,0.1,0.95,0.3);
    ptlegend->SetHeader("Legend","C"); // option "C" allows to center the header
    ptlegend->AddEntry(McPtPair,"MC Pair pT");
    ptlegend->AddEntry(mPtPairReco, "Reco Pair pT (binned in McpT)");
    ptlegend->SetTextSize(.012);
    ptlegend->Draw("same");


    makeCanvas();
    McEtaPair->Draw("PE");
    McEtaPair->GetYaxis()->SetRangeUser(1000, 100000);
    mEtaPairReco->SetLineColor(kRed);
    mEtaPairReco->Draw("PE;same");
    gPad->SetLogy();
    auto etalegend = new TLegend(0.85,0.1,0.95,0.3);
    etalegend->SetHeader("Legend","C"); // option "C" allows to center the header
    etalegend->AddEntry(McEtaPair,"MC Pair #eta");
    etalegend->AddEntry(mEtaPairReco, "Reco Pair #eta (binned in Mc #eta)");
    etalegend->SetTextSize(.012);
    etalegend->Draw("same");


    /*makeCanvas();
    TEfficiency * rapidityPaireff = new TEfficiency(* mEtaPairReco, * McEtaPair);
    rapidityPaireff->SetTitle("Pair Rapidity Efficiency; y_{ee}; Efficiency");
    rapidityPaireff->Draw("PE");


    makeCanvas();
    TEfficiency * ptPaireff = new TEfficiency(* mPtPairReco, * McPtPair);
    ptPaireff->SetTitle("Pair pT Efficiency; pT(GeV/c); Efficiency");
    ptPaireff->Draw("PE");*/


    makeCanvas();
    auto m2Ptcos2phimomentsMC = cos2phivPtMC->ProfileY("m2Ptcos2phimomentsMC",1, -1);
    m2Ptcos2phimomentsMC->SetLineColor(kBlack);
    m2Ptcos2phimomentsMC->GetYaxis()->SetTitle("2<cos(2#phi)>");
    m2Ptcos2phimomentsMC->Draw("PE");
    //gPad->Print("slightEff/plot_Starsimcos2phimomentsMC.png");

    makeCanvas();
    auto m2Ptcos4phimomentsMC = cos4phivPtMC->ProfileY("m2Ptcos4phimomentsMC",1, -1);
    m2Ptcos4phimomentsMC->SetLineColor(kBlack);
    m2Ptcos4phimomentsMC->GetYaxis()->SetTitle("2<cos(4#phi)>");
    m2Ptcos4phimomentsMC->Draw("PE");
    //gPad->Print("slightEff/plot_Starsimcos4phimomentsMC.png");

    makeCanvas();
    auto m2Ptcos2phimomentsReco = cos2phivPtReco->ProfileY("m2Ptcos2phimomentsReco",1, -1);
    m2Ptcos2phimomentsReco->SetLineColor(kBlack);
    m2Ptcos2phimomentsReco->GetYaxis()->SetTitle("2<cos(2#phi)>");
    m2Ptcos2phimomentsReco->Draw("PE");
    //gPad->Print("slightEff/plot_Starsimcos2phimomentsReco.png");

    makeCanvas();
    auto m2Ptcos4phimomentsReco = cos4phivPtReco->ProfileY("m2Ptcos4phimomentsReco",1, -1);
    m2Ptcos4phimomentsReco->SetLineColor(kBlack);
    m2Ptcos4phimomentsReco->GetYaxis()->SetTitle("2<cos(4#phi)>");
    m2Ptcos4phimomentsReco->Draw("PE");
    //gPad->Print("slightEff/plot_Starsimcos4phimomentsReco.png");


    TFile file("sim_plots.root", "RECREATE");
    mMassPairReco->Write();
    McMassPair->Write();

    mPtPairReco->Write();
    McPtPair->Write();

    mEtaPairReco->Write();
    McEtaPair->Write();


    m2Ptcos2phimomentsMC->Write();
    m2Ptcos2phimomentsReco->Write();
    m2Ptcos4phimomentsMC->Write();
    m2Ptcos4phimomentsReco->Write();

    std::cout << "wrote";

    std::cout << "MC Pairs: " << mc_pair_counter << "\n";
    std::cout << "Reco Pairs: " << reco_pair_counter;



}