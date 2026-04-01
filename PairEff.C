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


void PairEff(){
    TChain * ch = new TChain("FemtoDst");
    //ch->Add("/Users/Nick/STAR/docker_mount/bigJobResults/*");
    ch->Add("/Users/Nick/STAR/docker_mount/JobResultsAu1M/*");
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

    auto McPtM = new TH2F("McMassPt", "McMassPt", 50, 0, 0.2, 50, .3, 3);
    auto RcPtM = new TH2F("RcMassPt", "RcMassPt", 50, 0, 0.2, 50, .3, 3);

    auto McPtMY = new TH3F("McPtMY", "McPtMY", 50, 0, 0.2, 50, .3, 3, 20, -1, 1);
    auto RcPtMY = new TH3F("RcPtMY", "RcPtMY", 50, 0, 0.2, 50, .3, 3, 20, -1, 1);

    auto McPt2Mcostheta = new TH3F("McPt2Mcostheta", "McPt2Mcostheta", 50, 0, 0.02, 50, .3, 3, 20, 0, 1);
    auto RcPt2Mcostheta = new TH3F("RcPt2Mcostheta", "RcPt2Mcostheta", 50, 0, 0.02, 50, .3, 3, 20, 0, 1);




    const Int_t ptbins = 10;
    Double_t edges[ptbins + 1] = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.14, .2};


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

    auto cos2phivPtMC = new TH2F("Cos2#phivPtMC", "cos2#phi Moments vs P_{T} MC Pairs; 2<cos(2#phi)>; pT (GeV/c)", 100, -2, 2, 10, edges);
    auto cos4phivPtMC = new TH2F("Cos4#phivPtMC", "cos4#phi Moments vs P_{T} MC Pairs; 2<cos(4#phi)>; pT (GeV/c)", 100, -2, 2, 10, edges);

    auto cos2phivPtReco = new TH2F("Cos2#phivPtReco", "cos2#phi Moments vs P_{T} Reco Pairs; 2<cos(2#phi)>; pT (GeV/c)", 100, -2, 2, 10, edges);
    auto cos4phivPtReco = new TH2F("Cos4#phivPtReco", "cos4#phi Moments vs P_{T} Reco Pairs; 2<cos(4#phi)>; pT (GeV/c)", 100, -2, 2, 10, edges);

    auto cos2phivPtRecoRc = new TH2F("Cos2#phivPtRecoRc", "cos2#phi Moments vs P_{T} Reco Pairs (binned in RC); 2<cos(2#phi)>; pT (GeV/c)", 100, -2, 2, 10, edges);
    auto cos4phivPtRecoRc = new TH2F("Cos4#phivPtRecoRc", "cos4#phi Moments vs P_{T} Reco Pairs (binned in RC); 2<cos(4#phi)>; pT (GeV/c)", 100, -2, 2, 10, edges);

    auto mcphipt = new TH2F("mcphipt", "mcphipt (binned in MC); #phi; pT (GeV/c)", 25, -3.1415927, 3.1415927, 20, 0, 0.3);
    auto recophipt = new TH2F("recophipt", "recophipt (binned in RC); #phi>; pT (GeV/c)", 25, -3.1415927, 3.1415927, 20, 0, 0.3);

    auto mcrecophipt = new TH3F("mcrecophipt", "mcrecophipt", 30, -3.15, 3.15, 30, -3.15, 3.15, 30, 0, 0.6);

    auto mcrecophi = new TH2F("mcrecophi", "mcrecophi", 30, -3.15, 3.15, 30, -3.15, 3.15);

    auto MCaco = new TH1F("MCaco", "MCaco", 10, 0, .2);
    auto rcaco = new TH1F("rcaco", "rcaco", 10, 0, .2);

    auto ptresolution = new TH2F("ptresolution", "ptresolution", 20, 0, 1, 200, -.1, .1);

    auto ptpairRes = new TH2F("ptpairRes", "ptpairRes", 100, 0, 0.1, 200, -1, 1);

    auto recophiptmass = new TH3F("recophiptmass", "recophiptmass", 20, -3.15, 3.15, 50, 0, 0.2, 50, 0.3, 3);





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

        if ( fabs(MCEtaVals[0]) > 1 || fabs(MCEtaVals[1]) > 1 ) continue;

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

        McPtM->Fill(mcpair.Pt(), mcpair.M());

        McPtMY->Fill(mcpair.Pt(), mcpair.M(), mcpair.Rapidity());

        TLorentzVector lvbeam;
        lvbeam.SetPxPyPzE(0, 0, sqrt(100*100 - .938*.938), 100);
        lvbeam.Boost( -(mcpair.BoostVector()) );
        float costheta = fabs(cos( mcpos.Angle(lvbeam.Vect() ) ));

        McPt2Mcostheta->Fill(pow(mcpair.Pt(), 2), mcpair.M(), costheta);



        double phival = calc_Phi(mcpos, mcneg);
        cos2phivPtMC->Fill(2*cos(2*phival), mcpair.Pt());
        cos4phivPtMC->Fill(2*cos(4*phival), mcpair.Pt());
        mcphipt->Fill(phival, mcpair.Pt());

        if(mcpair.M() > 0.4 && mcpair.M() < .76){
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
        rcpair = rcpos + rcneg;
        rcpairreco = rcposreco + rcnegreco;
        double phivalreco = calc_Phi(rcposreco, rcnegreco);

        if(NHitsDedx[0] >= 15 && NHitsDedx[1] >= 15 && abs(NHitsFit[0]) >= 20 && abs(NHitsFit[1]) >= 20 && DCA[0] < 1 && DCA[1] < 1){

            //std::cout << "Filled reco stuff" << "\n";


            if(rcpair.M() > 0.4 && rcpair.M() < .76 && rcpair.Rapidity() < 1 && rcpair.Rapidity() > -1){
                RcPtPair->Fill(rcpair.Pt());
                RcPtMY->Fill(rcpair.Pt(), rcpair.M(), rcpair.Rapidity());

                

            }
            RcPt2Pair->Fill(rcpair.Pt()*rcpair.Pt());
            if(mcpair.Pt() < 0.1 && fabs(mcpair.Rapidity()) < 1.0) {
                RcMassPair->Fill(rcpair.M());
                RcRapidityPair->Fill(rcpair.Rapidity());
            }

            //RcMassPair->Fill(rcpair.M());

            RcPtM->Fill(rcpair.Pt(), rcpair.M());

            TLorentzVector lvbeamRC;
            lvbeamRC.SetPxPyPzE(0, 0, sqrt(100*100 - .938*.938), 100 );
            lvbeamRC.Boost( -(rcpair.BoostVector()) );
            float costhetaRC = fabs(cos( rcpos.Angle(lvbeam.Vect() ) ));

            RcPt2Mcostheta->Fill(pow(mcpair.Pt(), 2), mcpair.M(), costhetaRC);



            RcPairEtaVsPhi->Fill(rcpair.Eta(), rcpair.Phi());
            RcPosEtaVsPhi->Fill(rcpos.Eta(), rcpos.Phi());
            RcNegEtaVsPhi->Fill(rcneg.Eta(), rcneg.Phi());

            RCPosnhitsfitvpT->Fill(abs(NHitsFit[0]), rcpos.Pt());
            RCNegnhitsfitvpT->Fill(abs(NHitsFit[1]), rcneg.Pt());
            RCPosnhitsdedxvpT->Fill(NHitsDedx[0], rcpos.Pt());
            RCNegnhitsdedxvpT->Fill(NHitsDedx[1], rcneg.Pt());

            if(rcpairreco.M() > 0.4 && rcpairreco.M() < 0.76){
                cos2phivPtReco->Fill(2*cos(2*phival), rcpair.Pt());
                cos4phivPtReco->Fill(2*cos(4*phival), rcpair.Pt());

                cos2phivPtRecoRc->Fill(2*cos(2*phivalreco), rcpair.Pt());
                cos4phivPtRecoRc->Fill(2*cos(4*phivalreco), rcpair.Pt());
            }

            recophipt->Fill(phivalreco, rcpair.Pt());

            mcrecophipt->Fill(phival, phivalreco, mcpair.Pt());

            //in one pt bin
            //if(mcpair.Pt() < 0.06 && mcpair.Pt() > 0.04){
                mcrecophi->Fill(phival, phivalreco);
            //}
            double rcacoval = 1- abs(rcposreco.Phi() - rcnegreco.Phi())/M_PI;
            rcaco->Fill(rcacoval);

            double res_pos = (mcpos.Pt() - rcposreco.Pt())/mcpos.Pt();
            double res_neg = (mcneg.Pt()-rcnegreco.Pt())/mcneg.Pt();
            double res_pair = (mcpair.Pt() - rcpairreco.Pt())/mcpair.Pt();
            ptresolution->Fill(mcneg.Pt(), res_neg);
            ptresolution->Fill(mcpos.Pt(), res_pos);
            ptpairRes->Fill(mcpair.Pt(), res_pair);

            recophiptmass->Fill(phivalreco, rcpairreco.Pt(), rcpairreco.M());

        }

    }



    makeCanvas();
    McPtPair->Draw();
    RcPtPair->SetLineColor(kRed);
    RcPtPair->Draw("same");
    McPtPair->GetXaxis()->SetTitle("pT (GeV/c)");
    McPtPair->GetYaxis()->SetTitle("Counts");
    auto * hlegend = new TLegend(0.55,0.4,.85,0.55);
    hlegend->SetBorderSize(0);
    hlegend->SetTextSize(0.025);
    hlegend->AddEntry(McPtPair,"MC Pairs","l");
    hlegend->AddEntry(RcPtPair,"Reconstructed Pairs","l");
    hlegend->Draw("same");
    //gPad->Print("note_plots/sim_and_eff_plots/simptNEW.png");


    makeCanvas();
    TEfficiency * hEff = new TEfficiency(* RcPtPair, * McPtPair);
    hEff->SetTitle("pT Pair Efficiency; P_{T} (GeV/c); Efficiency");
    hEff->Draw("PE");

    //gPad->Print("note_plots/sim_and_eff_plots/pteffNEW.png");


    makeCanvas();
    TEfficiency * mEff = new TEfficiency(* RcMassPair, * McMassPair);
    mEff->SetTitle("Pair Invariant Mass Efficiency; M_{ee} (GeV/c^{2}); Efficiency");
    mEff->Draw("PE");

    //gPad->Print("note_plots/sim_and_eff_plots/masseffNEW.png");



    makeCanvas();
    McMassPair->Draw();
    RcMassPair->SetLineColor(kRed);
    RcMassPair->Draw("same");
    McMassPair->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
    McMassPair->GetYaxis()->SetTitle("Counts");
    auto * mlegend = new TLegend(0.8,0.75,1.,.9);
    mlegend->SetHeader("Legend");
    mlegend->AddEntry(McMassPair,"MC Pairs","l");
    mlegend->AddEntry(RcMassPair,"Reconstructed Pairs","l");
    mlegend->Draw("same");
    //gPad->Print("note_plots/sim_and_eff_plots/simmassNEW.png");


    makeCanvas();
    TEfficiency * YEff = new TEfficiency(* RcRapidityPair, * McRapidityPair);
    YEff->SetTitle("Pair Rapidity Efficiency; y_{ee}; Efficiency");
    YEff->Draw("PE");
    //gPad->Print("note_plots/sim_and_eff_plots/yeff.png");


    makeCanvas();
    McRapidityPair->Draw();
    RcRapidityPair->SetLineColor(kRed);
    RcRapidityPair->Draw("same");
    McRapidityPair->GetXaxis()->SetTitle("y_{ee}");
    McRapidityPair->GetYaxis()->SetTitle("Counts");

    auto * ylegend = new TLegend(0.8,0.75,1.,.9);
    ylegend->SetHeader("Legend");
    ylegend->AddEntry(McRapidityPair,"MC Pairs","l");
    ylegend->AddEntry(RcRapidityPair,"Reconstructed Pairs","l");
    ylegend->Draw("same");

    //gPad->Print("note_plots/sim_and_eff_plots/simyNEW.png");


    makeCanvas();
    TEfficiency * h2Eff = new TEfficiency(* RcPt2Pair, * McPt2Pair);
    h2Eff->SetTitle("Pair pT^{2} Efficiency; P_{T}^{2} (GeV/c)^{2}; Efficiency");
    h2Eff->Draw("PE");
    //gPad->Print("note_plots/sim_and_eff_plots/pt2effNEW.png");


    makeCanvas();
    McPt2Pair->Draw();
    RcPt2Pair->SetLineColor(kRed);
    RcPt2Pair->Draw("same");
    McPt2Pair->GetXaxis()->SetTitle("pT^{2} (GeV/c)^{2}");
    McPt2Pair->GetYaxis()->SetTitle("Counts");

    auto * pt2legend = new TLegend(0.8,0.75,1.,.9);
    pt2legend->SetHeader("Legend");
    pt2legend->AddEntry(McPt2Pair,"MC Pairs","l");
    pt2legend->AddEntry(RcPt2Pair,"Reconstructed Pairs","l");
    pt2legend->Draw("same");

    //gPad->Print("note_plots/sim_and_eff_plots/simpt2NEW.png");




    makeCanvas();
    auto m2Ptcos2phimomentsMC = cos2phivPtMC->ProfileY("m2Ptcos2phimomentsMC",1, -1);
    m2Ptcos2phimomentsMC->SetLineColor(kBlack);
    m2Ptcos2phimomentsMC->GetYaxis()->SetTitle("2<cos(2#phi)>");
    m2Ptcos2phimomentsMC->Draw("PE");
    //gPad->Print("note_plots/sim_and_eff_plots/plot_Starsimcos2phimomentsMCNEW.png");

    makeCanvas();
    auto m2Ptcos4phimomentsMC = cos4phivPtMC->ProfileY("m2Ptcos4phimomentsMC",1, -1);
    m2Ptcos4phimomentsMC->SetLineColor(kBlack);
    m2Ptcos4phimomentsMC->GetYaxis()->SetTitle("2<cos(4#phi)>");
    m2Ptcos4phimomentsMC->Draw("PE");
    //gPad->Print("note_plots/sim_and_eff_plots/plot_Starsimcos4phimomentsMCNEW.png");

    makeCanvas();
    auto m2Ptcos2phimomentsReco = cos2phivPtReco->ProfileY("m2Ptcos2phimomentsReco",1, -1);
    m2Ptcos2phimomentsReco->SetLineColor(kBlack);
    m2Ptcos2phimomentsReco->GetYaxis()->SetTitle("2<cos(2#phi)>");
    m2Ptcos2phimomentsReco->Draw("PE");
   //gPad->Print("note_plots/sim_and_eff_plots/plot_Starsimcos2phimomentsRecoNEW.png");

    makeCanvas();
    auto m2Ptcos4phimomentsReco = cos4phivPtReco->ProfileY("m2Ptcos4phimomentsReco",1, -1);
    m2Ptcos4phimomentsReco->SetLineColor(kBlack);
    m2Ptcos4phimomentsReco->GetYaxis()->SetTitle("2<cos(4#phi)>");
    m2Ptcos4phimomentsReco->Draw("PE");
    //gPad->Print("note_plots/sim_and_eff_plots/plot_Starsimcos4phimomentsRecoNEW.png");

    makeCanvas();
    auto m2Ptcos2phimomentsRecoRc = cos2phivPtRecoRc->ProfileY("m2Ptcos2phimomentsRecoRc",1, -1);
    m2Ptcos2phimomentsRecoRc->SetLineColor(kBlack);
    m2Ptcos2phimomentsRecoRc->GetYaxis()->SetTitle("2<cos(2#phi)>");
    m2Ptcos2phimomentsRecoRc->Draw("PE");
    //gPad->Print("note_plots/sim_and_eff_plots/plot_Starsimcos2phimomentsRecoRcNEW.png");

    makeCanvas();
    auto m2Ptcos4phimomentsRecoRc = cos4phivPtRecoRc->ProfileY("m2Ptcos4phimomentsRecoRc",1, -1);
    m2Ptcos4phimomentsRecoRc->SetLineColor(kBlack);
    m2Ptcos4phimomentsRecoRc->GetYaxis()->SetTitle("2<cos(4#phi)>");
    m2Ptcos4phimomentsRecoRc->Draw("PE");
    //gPad->Print("note_plots/sim_and_eff_plots/plot_Starsimcos4phimomentsRecoRcNEW.png");

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


    makeCanvas();
    TEfficiency * PtMEff = new TEfficiency(* RcPtM, * McPtM);
    PtMEff->SetTitle("Pair Pt M Efficiency");
    PtMEff->Draw("colz");

    //gPad->Print("note_plots/sim_and_eff_plots/ptmEff.png");



    TEfficiency * PtMYEff = new TEfficiency( * RcPtMY, * McPtMY);
    TEfficiency * Pt2McosthetaEff = new TEfficiency( * RcPt2Mcostheta ,* McPt2Mcostheta);



    makeCanvas();
    ptresolution->GetYaxis()->SetTitle("(MC pt_{e} - reco pt_{e})/MC pt_{e}");
    ptresolution->GetXaxis()->SetTitle("MC pt_{e} (GeV/c)");
    ptresolution->SetTitle("Track Pt Resolution");
    gPad->SetLogz();
    ptresolution->Draw("colz");

    makeCanvas();
    ptpairRes->GetYaxis()->SetTitle("(MC pt_{ee} - reco pt_{ee})/MC pt_{ee}");
    ptpairRes->GetXaxis()->SetTitle("MC pt_{ee} (GeV/c)");
    ptpairRes->SetTitle("Pair Pt Resolution");
    gPad->SetLogz();
    ptpairRes->Draw("colz");



    TFile file("output_root_files/simulation_plots_new_Au.root", "RECREATE");
    McPtPair->Write();
    RcPtPair->Write();
    McMassPair->Write();
    RcMassPair->Write();
    McRapidityPair->Write();
    RcRapidityPair->Write();
    McPtM->Write();
    RcPtM->Write();
    hEff->Write("hEff");
    h2Eff->Write("h2Eff");
    mEff->Write("mEff");
    YEff->Write("YEff");
    PtMEff->Write("PtMEff");
    McPtMY->Write("McPtMY");
    RcPtMY->Write("RcPtMY");
    PtMYEff->Write("PtMYEff");
    Pt2McosthetaEff->Write("Pt2McosthetaEff");


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
    ptresolution->Write();
    ptpairRes->Write();


    recophiptmass->Write();



}
