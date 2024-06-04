#include <iostream>
int ican2 = 0;
void makeCanvas()  {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican2++ ), "", 900, 600);
    can->SetTopMargin(0.08);
    can->SetRightMargin(0.15);
}
/*double calc_Phi( TLorentzVector lv1, TLorentzVector lv2) {
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
}*/
double calc_Phi( TLorentzVector lv1, TLorentzVector lv2) {
    TLorentzVector lvPlus = lv1 + lv2;
    TLorentzVector lvMinus = lv1 - lv2;
    return lvMinus.DeltaPhi(lvPlus);
}


void PairEff(){
    TChain * ch = new TChain("FemtoDst");
    ch->Add("/Users/Nick/STAR/docker_mount/JobResults/*");
    TTreeReader myReader(ch);

    
    TTreeReaderArray<float> MCPtVals(myReader, "McTracks.mPt");
    TTreeReaderArray<float> MCEtaVals(myReader, "McTracks.mEta");
    TTreeReaderArray<float> MCPhiVals(myReader, "McTracks.mPhi");

    TTreeReaderArray<short> MCIdx(myReader, "Tracks.mMcIndex");
    TTreeReaderArray<short> NHitsFit(myReader, "Tracks.mNHitsFit"); 
    TTreeReaderArray<short> NHitsDedx(myReader, "Tracks.mNHitsDedx");
    TTreeReaderArray<float> DCA(myReader, "Tracks.mDCA");

    TTreeReaderArray<float> RCPtVals(myReader, "Tracks.mPt");
    TTreeReaderArray<float> RCEtaVals(myReader, "Tracks.mEta");
    TTreeReaderArray<float> RCPhiVals(myReader, "Tracks.mPhi");





    TTreeReaderArray<unsigned short> PIDs(myReader, "McTracks.mGeantPID");
    TTreeReaderArray<short> ParentIndex(myReader, "McTracks.mParentIndex");

    cout << "Events = " << ch->GetEntries() << endl;
    auto McPtPair = new TH1F("McPtPair", "McPtPair", 50, 0, 0.4);
    auto RcPtPair = new TH1F("RcPtPair", "RcPtPair", 50, 0, 0.4);

    auto McPt2Pair = new TH1F("McPt2Pair", "McPt2Pair", 50, 0, 0.02);
    auto RcPt2Pair = new TH1F("RcPt2Pair", "RcPt2Pair", 50, 0, 0.02);


    auto McMassPair = new TH1F("McMassPair", "McMassPair", 50, 0, 2);
    auto RcMassPair = new TH1F("RcMassPair", "RcMassPair", 50, 0, 2);

    auto McRapidityPair = new TH1F("McRapidityPair", "McRapidityPair", 50, -1, 1);
    auto RcRapidityPair = new TH1F("RcRapidityPair", "RcRapidityPair", 50, -1, 1);

    //for track and pair acceptance checks for tpc, nhitsfit, and nhitsdedx comparisons with data
    //These histograms are used in TPCCheck.C

    //checking tpc acceptance in simulation
    auto McPairEtaVsPhi = new TH2D("McPairEtaVsPhi", "McPairEtaVsPhi", 50, -1, 1, 50, -3.14159, 3.14159);
    auto McPosEtaVsPhi = new TH2D("McPosEtaVsPhi", "McPosEtaVsPhi", 50, -1, 1, 50, -3.14159, 3.14159);
    auto McNegEtaVsPhi = new TH2D("McNegEtaVsPhi", "McNegEtaVsPhi", 50, -1, 1, 50, -3.14159, 3.14159);

    auto RcPairEtaVsPhi = new TH2D("RcPairEtaVsPhi", "RcPairEtaVsPhi", 50, -1, 1, 50, -3.14159, 3.14159);
    auto RcPosEtaVsPhi = new TH2D("RcPosEtaVsPhi", "RcPosEtaVsPhi", 50, -1, 1, 50, -3.14159, 3.14159);
    auto RcNegEtaVsPhi = new TH2D("RcNegEtaVsPhi", "RcNegEtaVsPhi", 50, -1, 1, 50, -3.14159, 3.14159);

    //checking nhitsfit and dedx in simulation
    auto RCPosnhitsfitvpT = new TH2D("RCPosnhitsfitvpT", "RCPosnhitsfitvpT", 50, 0, 50, 80, 0.2, 1);
    auto RCNegnhitsfitvpT = new TH2D("RCNegnhitsfitvpT", "RCNegnhitsfitvpT", 50, 0, 50, 80, 0.2, 1);

    auto RCPosnhitsdedxvpT = new TH2D("RCPosnhitsdedxvpT", "RCPosnhitsdedxvpT", 50, 0, 50, 80, 0.2, 1);
    auto RCNegnhitsdedxvpT = new TH2D("RCNegnhitsdedxvpT", "RCNegnhitsdedxvpT", 50, 0, 50, 80, 0.2, 1);

    auto cos2phivPtMC = new TH2F("Cos2#phivPtMC", "cos2#phi Moments vs P_{T} MC Pairs; 2<cos(2#phi)>; pT (GeV/c)", 100, -2, 2, 7, 0, 0.1);
    auto cos4phivPtMC = new TH2F("Cos4#phivPtMC", "cos4#phi Moments vs P_{T} MC Pairs; 2<cos(4#phi)>; pT (GeV/c)", 100, -2, 2, 7, 0, 0.1);

    auto cos2phivPtReco = new TH2F("Cos2#phivPtReco", "cos2#phi Moments vs P_{T} Reco Pairs; 2<cos(2#phi)>; pT (GeV/c)", 100, -2, 2, 7, 0, 0.1);
    auto cos4phivPtReco = new TH2F("Cos4#phivPtReco", "cos4#phi Moments vs P_{T} Reco Pairs; 2<cos(4#phi)>; pT (GeV/c)", 100, -2, 2, 7, 0, 0.1);

    auto cos2phivPtRecoRc = new TH2F("Cos2#phivPtRecoRc", "cos2#phi Moments vs P_{T} Reco Pairs (binned in RC); 2<cos(2#phi)>; pT (GeV/c)", 100, -2, 2, 20, 0, 0.3);
    auto cos4phivPtRecoRc = new TH2F("Cos4#phivPtRecoRc", "cos4#phi Moments vs P_{T} Reco Pairs (binned in RC); 2<cos(4#phi)>; pT (GeV/c)", 100, -2, 2, 20, 0, 0.3);
    
    auto mcphipt = new TH2F("mcphipt", "mcphipt (binned in MC); #phi; pT (GeV/c)", 25, -3.1415927, 3.1415927, 20, 0, 0.3);
    auto recophipt = new TH2F("recophipt", "recophipt (binned in RC); #phi>; pT (GeV/c)", 25, -3.1415927, 3.1415927, 20, 0, 0.3);

    auto mcrecophipt = new TH3F("mcrecophipt", "mcrecophipt", 30, -3.15, 3.15, 30, -3.15, 3.15, 30, 0, 0.6);

    auto mcrecophi = new TH2F("mcrecophi", "mcrecophi", 30, -3.15, 3.15, 30, -3.15, 3.15);

    auto MCaco = new TH1F("MCaco", "MCaco", 20, 0, .2);
    auto rcaco = new TH1F("rcaco", "rcaco", 20, 0, .2);





    while(myReader.Next()) {
        
        TLorentzVector mcneg, mcpos, mcpair;
        TLorentzVector rcneg, rcpos, rcpair; //(these are still MC values, just from the RC tracks matched MC )

        TLorentzVector rcnegreco, rcposreco, rcpairreco; //rc values

        if ( PIDs.GetSize() == 0){
            cout << "No MC tracks" << endl;
            continue;
        } else if ( PIDs.GetSize() < 2 ){
            cout << "did not find two MC tracks, skipping event" << endl;
            continue;
        }

        if ( ParentIndex[0] != -1 || ParentIndex[1] != -1  ){
            cout << "first two tracks not from STARLight, skipping event" << endl;
            continue;
        }

        if ( MCPtVals[0] < 0.2 || MCPtVals[1] < 0.2 ) continue;

        if ( fabs(MCEtaVals[0]) > 0.9 || fabs(MCEtaVals[1]) > 0.9 ) continue;

        if ( PIDs[0] == 2 && PIDs[1] == 3 ){
            mcpos.SetPtEtaPhiM( MCPtVals[0], MCEtaVals[0], MCPhiVals[0], 0.00051099895000 );
            mcneg.SetPtEtaPhiM( MCPtVals[1], MCEtaVals[1], MCPhiVals[1], 0.00051099895000 );
        } else if ( PIDs[0] == 3 && PIDs[1] == 2 ){
            mcneg.SetPtEtaPhiM( MCPtVals[0], MCEtaVals[0], MCPhiVals[0], 0.00051099895000 );
            mcpos.SetPtEtaPhiM( MCPtVals[1], MCEtaVals[1], MCPhiVals[1], 0.00051099895000 );
        }

        mcpair = mcpos + mcneg;
        McPtPair->Fill( mcpair.Pt() );
        McPt2Pair->Fill(mcpair.Pt()*mcpair.Pt());
        McMassPair->Fill(mcpair.M());
        McRapidityPair->Fill(mcpair.Rapidity());

        McPairEtaVsPhi->Fill(mcpair.Eta(), mcpair.Phi());
        McPosEtaVsPhi->Fill(mcpos.Eta(), mcpos.Phi());
        McNegEtaVsPhi->Fill(mcneg.Eta(), mcneg.Phi());

        double phival = calc_Phi(mcpos, mcneg);
        cos2phivPtMC->Fill(2*cos(2*phival), mcpair.Pt());
        cos4phivPtMC->Fill(2*cos(4*phival), mcpair.Pt());
        mcphipt->Fill(phival, mcpair.Pt());

        if(mcpair.M() > 0.4 && mcpair.M() < 1){
            double aco = 1 - abs(mcpos.Phi() - mcneg.Phi())/M_PI;
            MCaco->Fill(aco);
        }

        //std::cout << "filled mc" << "\n";



        short idxPos = -1, idxNeg = -1;
        for ( size_t i = 0; i < MCIdx.GetSize(); i++  ){
            short idx = MCIdx[i];
            if ( idx == -1 || idx > 1) continue;
            if ( ParentIndex[idx] != -1 ) continue;

            if ( PIDs[idx] == 2 ){
                idxPos = idx;
            } else if ( PIDs[idx] == 3 ){
                idxNeg = idx;
            }
        }


        if ( idxPos <= -1 || idxNeg <= -1 ) continue;

        //code not getting here
        //std::cout << "got past continue";

        rcneg.SetPtEtaPhiM( MCPtVals[idxNeg], MCEtaVals[idxNeg], MCPhiVals[idxNeg], 0.00051099895000 );
        rcpos.SetPtEtaPhiM( MCPtVals[idxPos], MCEtaVals[idxPos], MCPhiVals[idxPos], 0.00051099895000 );
        //std::cout << "Set reco stuff" << "\n";

        rcnegreco.SetPtEtaPhiM( RCPtVals[idxNeg], RCEtaVals[idxNeg], RCPhiVals[idxNeg], 0.00051099895000 );
        rcposreco.SetPtEtaPhiM( RCPtVals[idxPos], RCEtaVals[idxPos], RCPhiVals[idxPos], 0.00051099895000 );
        rcpair = rcneg + rcpos;
        rcpairreco = rcnegreco + rcposreco;
        double phivalreco = calc_Phi(rcnegreco, rcposreco);
        
        if(NHitsDedx[0] >= 15 && NHitsDedx[1] >= 15 && abs(NHitsFit[0]) >= 20 && abs(NHitsFit[1]) >= 20 && DCA[0] < 1 && DCA[1] < 1){
            
            //std::cout << "Filled reco stuff" << "\n";


            RcPtPair->Fill(rcpair.Pt());
            RcPt2Pair->Fill(rcpair.Pt()*rcpair.Pt());
            RcMassPair->Fill(rcpair.M());
            RcRapidityPair->Fill(rcpair.Rapidity());
        

            RcPairEtaVsPhi->Fill(rcpair.Eta(), rcpair.Phi());
            RcPosEtaVsPhi->Fill(rcpos.Eta(), rcpos.Phi());
            RcNegEtaVsPhi->Fill(rcneg.Eta(), rcneg.Phi());

            RCPosnhitsfitvpT->Fill(abs(NHitsFit[0]), rcpos.Pt());
            RCNegnhitsfitvpT->Fill(abs(NHitsFit[1]), rcneg.Pt());
            RCPosnhitsdedxvpT->Fill(NHitsDedx[0], rcpos.Pt());
            RCNegnhitsdedxvpT->Fill(NHitsDedx[1], rcneg.Pt());

            cos2phivPtReco->Fill(2*cos(2*phival), rcpair.Pt());
            cos4phivPtReco->Fill(2*cos(4*phival), rcpair.Pt()); 

            cos2phivPtRecoRc->Fill(2*cos(2*phivalreco), rcpair.Pt());
            cos4phivPtRecoRc->Fill(2*cos(4*phivalreco), rcpair.Pt()); 

            recophipt->Fill(phivalreco, rcpair.Pt());

            mcrecophipt->Fill(phival, phivalreco, mcpair.Pt());

            //in one pt bin
            if(mcpair.Pt() < 0.06 && mcpair.Pt() > 0.04){
                mcrecophi->Fill(phival, phivalreco);
            }
            double rcacoval = 1- abs(rcposreco.Phi() - rcnegreco.Phi())/M_PI;
            rcaco->Fill(rcacoval);



        }
        
    }


   
    makeCanvas();    
    McPtPair->Draw();
    RcPtPair->SetLineColor(kRed);
    RcPtPair->Draw("same");
    McPtPair->GetXaxis()->SetTitle("pT (GeV/c)");
    McPtPair->GetYaxis()->SetTitle("Counts");
    gPad->Print("note_plots/sim_and_eff_plots/simpt.png");


    makeCanvas();    
    TEfficiency * hEff = new TEfficiency(* RcPtPair, * McPtPair);
    hEff->SetTitle("pT Pair Efficiency");
    hEff->Draw("PE");
    gPad->Print("note_plots/sim_and_eff_plots/pteff.png");


    makeCanvas();    
    TEfficiency * mEff = new TEfficiency(* RcMassPair, * McMassPair);
    mEff->SetTitle("Pair Invariant Mass Efficiency");
    mEff->Draw("PE");
    
    gPad->Print("note_plots/sim_and_eff_plots/masseff.png");


    makeCanvas();    
    McMassPair->Draw();
    RcMassPair->SetLineColor(kRed);
    RcMassPair->Draw("same");
    McMassPair->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
    McMassPair->GetYaxis()->SetTitle("Counts");
    gPad->Print("note_plots/sim_and_eff_plots/simmass.png");


    makeCanvas();    
    TEfficiency * YEff = new TEfficiency(* RcRapidityPair, * McRapidityPair);
    YEff->SetTitle("Pair Rapidity Efficiency");
    YEff->Draw("PE");
    gPad->Print("note_plots/sim_and_eff_plots/yeff.png");


    makeCanvas();    
    McRapidityPair->Draw();
    RcRapidityPair->SetLineColor(kRed);
    RcRapidityPair->Draw("same");
    McRapidityPair->GetXaxis()->SetTitle("y_{ee}");
    McRapidityPair->GetYaxis()->SetTitle("Counts");
    gPad->Print("note_plots/sim_and_eff_plots/simy.png");


    makeCanvas();    
    TEfficiency * h2Eff = new TEfficiency(* RcPt2Pair, * McPt2Pair);
    h2Eff->SetTitle("Pair pT^{2} Efficiency");
    h2Eff->Draw("PE");
    gPad->Print("note_plots/sim_and_eff_plots/pt2eff.png");


    makeCanvas();    
    McPt2Pair->Draw();
    RcPt2Pair->SetLineColor(kRed);
    RcPt2Pair->Draw("same");
    McPt2Pair->GetXaxis()->SetTitle("pT^{2} (GeV/c)^{2}");
    McPt2Pair->GetYaxis()->SetTitle("Counts");
    gPad->Print("note_plots/sim_and_eff_plots/simpt2.png");


    makeCanvas();
    auto m2Ptcos2phimomentsMC = cos2phivPtMC->ProfileY("m2Ptcos2phimomentsMC",1, -1);
    m2Ptcos2phimomentsMC->SetLineColor(kBlack);
    m2Ptcos2phimomentsMC->GetYaxis()->SetTitle("2<cos(2#phi)>");
    m2Ptcos2phimomentsMC->Draw("PE");
    gPad->Print("note_plots/sim_and_eff_plots/plot_Starsimcos2phimomentsMC.png");

    makeCanvas();
    auto m2Ptcos4phimomentsMC = cos4phivPtMC->ProfileY("m2Ptcos4phimomentsMC",1, -1);
    m2Ptcos4phimomentsMC->SetLineColor(kBlack);
    m2Ptcos4phimomentsMC->GetYaxis()->SetTitle("2<cos(4#phi)>");
    m2Ptcos4phimomentsMC->Draw("PE");
    gPad->Print("note_plots/sim_and_eff_plots/plot_Starsimcos4phimomentsMC.png");

    makeCanvas();
    auto m2Ptcos2phimomentsReco = cos2phivPtReco->ProfileY("m2Ptcos2phimomentsReco",1, -1);
    m2Ptcos2phimomentsReco->SetLineColor(kBlack);
    m2Ptcos2phimomentsReco->GetYaxis()->SetTitle("2<cos(2#phi)>");
    m2Ptcos2phimomentsReco->Draw("PE");
    gPad->Print("note_plots/sim_and_eff_plots/plot_Starsimcos2phimomentsReco.png");

    makeCanvas();
    auto m2Ptcos4phimomentsReco = cos4phivPtReco->ProfileY("m2Ptcos4phimomentsReco",1, -1);
    m2Ptcos4phimomentsReco->SetLineColor(kBlack);
    m2Ptcos4phimomentsReco->GetYaxis()->SetTitle("2<cos(4#phi)>");
    m2Ptcos4phimomentsReco->Draw("PE");
    gPad->Print("note_plots/sim_and_eff_plots/plot_Starsimcos4phimomentsReco.png");

    makeCanvas();
    auto m2Ptcos2phimomentsRecoRc = cos2phivPtRecoRc->ProfileY("m2Ptcos2phimomentsRecoRc",1, -1);
    m2Ptcos2phimomentsRecoRc->SetLineColor(kBlack);
    m2Ptcos2phimomentsRecoRc->GetYaxis()->SetTitle("2<cos(2#phi)>");
    m2Ptcos2phimomentsRecoRc->Draw("PE");
    gPad->Print("note_plots/sim_and_eff_plots/plot_Starsimcos2phimomentsRecoRc.png");

    makeCanvas();
    auto m2Ptcos4phimomentsRecoRc = cos4phivPtRecoRc->ProfileY("m2Ptcos4phimomentsRecoRc",1, -1);
    m2Ptcos4phimomentsRecoRc->SetLineColor(kBlack);
    m2Ptcos4phimomentsRecoRc->GetYaxis()->SetTitle("2<cos(4#phi)>");
    m2Ptcos4phimomentsRecoRc->Draw("PE");
    gPad->Print("note_plots/sim_and_eff_plots/plot_Starsimcos4phimomentsRecoRc.png");

    makeCanvas();
    MCaco->Scale(1/MCaco->GetEntries());
    rcaco->Scale(1/rcaco->GetEntries());
    for(int i = 0; i < MCaco->GetNbinsX(); ++i){
        double binwidth = MCaco->GetBinWidth(i);
        std::cout << "bin width " << binwidth; 
        std::cout << "bin content: " << MCaco->GetBinContent(i);
        double content = MCaco->GetBinContent(i)/binwidth;
        double rccontent = rcaco->GetBinContent(i)/binwidth;
        MCaco->SetBinContent(i, content);
        rcaco->SetBinContent(i, rccontent);
    }
    MCaco->GetXaxis()->SetTitle("#alpha");
    MCaco->GetYaxis()->SetTitle("1/N * dN/d#alpha");
    MCaco->Draw("PE");
    rcaco->SetLineColor(kRed);
    rcaco->Draw("PE;same");

    TFile file("output_root_files/simulation_plots.root", "RECREATE");
    McPtPair->Write();
    RcPtPair->Write();   
    McMassPair->Write(); 
    RcMassPair->Write(); 
    McRapidityPair->Write(); 
    RcRapidityPair->Write(); 
    hEff->Write("hEff");
    h2Eff->Write("h2Eff");    
    mEff->Write("mEff");
    YEff->Write("YEff");

    McPairEtaVsPhi->Write();
    McPosEtaVsPhi->Write();
    McNegEtaVsPhi->Write();

    RcPairEtaVsPhi->Write();
    RcPosEtaVsPhi->Write();
    RcNegEtaVsPhi->Write();

    RCPosnhitsfitvpT->Write();
    RCNegnhitsfitvpT->Write();
    RCPosnhitsdedxvpT->Write();
    RCNegnhitsdedxvpT->Write();

    m2Ptcos2phimomentsMC->Write();
    m2Ptcos4phimomentsMC->Write();
    m2Ptcos2phimomentsReco->Write();
    m2Ptcos4phimomentsReco->Write();
    m2Ptcos2phimomentsRecoRc->Write();
    m2Ptcos4phimomentsRecoRc->Write();

    recophipt->Write();
    mcphipt->Write();
    mcrecophipt->Write();
    mcrecophi->Write();

    MCaco->Write();
    rcaco->Write();



 
}