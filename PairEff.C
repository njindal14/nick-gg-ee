


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



    TTreeReaderArray<unsigned short> PIDs(myReader, "McTracks.mGeantPID");
    TTreeReaderArray<short> ParentIndex(myReader, "McTracks.mParentIndex");

    cout << "Events = " << ch->GetEntries() << endl;
    auto McPtPair = new TH1F("McPtPair", "McPtPair", 50, 0, 0.4);
    auto RcPtPair = new TH1F("RcPtPair", "RcPtPair", 50, 0, 0.4);

    auto McMassPair = new TH1F("McMassPair", "McMassPair", 50, 0, 2);
    auto RcMassPair = new TH1F("RcMassPair", "RcMassPair", 50, 0, 2);

    auto McRapidityPair = new TH1F("McRapidityPair", "McRapidityPair", 50, -1, 1);
    auto RcRapidityPair = new TH1F("RcRapidityPair", "RcRapidityPair", 50, -1, 1);

    //for track and pair acceptance checks for tpc, nhitsfit, and nhitsdedx comparisons with data
    //These histograms are used in TPCCheck.C

    //checking tpc acceptance in simulation
    auto McPairEtaVsPhi = new TH2D("McPairEtaVsPhi", "McPairEtaVsPhi", 50, -1, 1, 50, -3.14159, 3.14159);
    auto McPosEtaVsPhi = new TH2D("McPosEtaVsPhi", "McPosEtaVsPhi", 50, -1, 1, 50, -3.14159, 3.14159);
    auto McNegEtaVsPhi = new TH2D("MCNegEtaVsPhi", "MCNegEtaVsPhi", 50, -1, 1, 50, -3.14159, 3.14159);

    auto RcPairEtaVsPhi = new TH2D("RcPairEtaVsPhi", "RcPairEtaVsPhi", 50, -1, 1, 50, -3.14159, 3.14159);
    auto RcPosEtaVsPhi = new TH2D("RcPosEtaVsPhi", "RcPosEtaVsPhi", 50, -1, 1, 50, -3.14159, 3.14159);
    auto RcNegEtaVsPhi = new TH2D("RcNegEtaVsPhi", "RcNegEtaVsPhi", 50, -1, 1, 50, -3.14159, 3.14159);

    //checking nhitsfit and dedx in simulation
    auto RCPosnhitsfitvpT = new TH2D("RCPosnhitsfitvpT", "RCPosnhitsfitvpT", 50, 0, 50, 80, 0.2, 1);
    auto RCNegnhitsfitvpT = new TH2D("RCNegnhitsfitvpT", "RCNegnhitsfitvpT", 50, 0, 50, 80, 0.2, 1);

    auto RCPosnhitsdedxvpT = new TH2D("RCPosnhitsdedxvpT", "RCPosnhitsdedxvpT", 50, 0, 50, 80, 0.2, 1);
    auto RCNegnhitsdedxvpT = new TH2D("RCNegnhitsdedxvpT", "RCNegnhitsdedxvpT", 50, 0, 50, 80, 0.2, 1);





    while(myReader.Next()) {
        
        TLorentzVector mcneg, mcpos, mcpair;
        TLorentzVector rcneg, rcpos, rcpair; //(these are still MC values, just from the RC tracks matched MC )

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
        if ( fabs(MCEtaVals[0]) > 0.9 || fabs(MCPtVals[1]) > 0.9 ) continue;

        if ( PIDs[0] == 2 && PIDs[1] == 3 ){
            mcpos.SetPtEtaPhiM( MCPtVals[0], MCEtaVals[0], MCPhiVals[0], 0.00051099895000 );
            mcneg.SetPtEtaPhiM( MCPtVals[1], MCEtaVals[1], MCPhiVals[1], 0.00051099895000 );
        } else if ( PIDs[0] == 3 && PIDs[1] == 2 ){
            mcneg.SetPtEtaPhiM( MCPtVals[0], MCEtaVals[0], MCPhiVals[0], 0.00051099895000 );
            mcpos.SetPtEtaPhiM( MCPtVals[1], MCEtaVals[1], MCPhiVals[1], 0.00051099895000 );
        }

        mcpair = mcpos + mcneg;
        McPtPair->Fill( mcpair.Pt() );
        McMassPair->Fill(mcpair.M());
        McRapidityPair->Fill(mcpair.Rapidity());

        McPairEtaVsPhi->Fill(mcpair.Eta(), mcpair.Phi());
        McPosEtaVsPhi->Fill(mcpos.Eta(), mcpos.Phi());
        McNegEtaVsPhi->Fill(mcneg.Eta(), mcneg.Phi());



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

        rcneg.SetPtEtaPhiM( MCPtVals[idxNeg], MCEtaVals[idxNeg], MCPhiVals[idxNeg], 0.00051099895000 );
        rcpos.SetPtEtaPhiM( MCPtVals[idxPos], MCEtaVals[idxPos], MCPhiVals[idxPos], 0.00051099895000 );
        rcpair = rcneg + rcpos;
        RcPtPair->Fill(rcpair.Pt());
        RcMassPair->Fill(rcpair.M());
        RcRapidityPair->Fill(rcpair.Rapidity());

        RcPairEtaVsPhi->Fill(rcpair.Eta(), rcpair.Phi());
        RcPosEtaVsPhi->Fill(rcpos.Eta(), rcpos.Phi());
        RcNegEtaVsPhi->Fill(rcneg.Eta(), rcneg.Phi());

        RCPosnhitsfitvpT->Fill(abs(NHitsFit[0]), rcpos.Pt());
        RCNegnhitsfitvpT->Fill(abs(NHitsFit[1]), rcneg.Pt());
        RCPosnhitsdedxvpT->Fill(NHitsDedx[0], rcpos.Pt());
        RCNegnhitsdedxvpT->Fill(NHitsDedx[1], rcneg.Pt());
        
    }


    McPtPair->Draw();
    RcPtPair->SetLineColor(kRed);
    RcPtPair->Draw("same");

    TCanvas *c2 = new TCanvas("c2");
    TH1* hEff = (TH1*)RcPtPair->Clone( "Eff" );
    hEff->Divide( McPtPair );
    hEff->SetTitle("pT Pair Efficiency");
    hEff->Draw();

    TCanvas * c3 = new TCanvas("c3");
    TH1* mEff = (TH1*)RcMassPair->Clone("mEff");
    mEff->Divide(McMassPair);
    mEff->SetTitle("Pair Invariant Mass Efficiency");
    mEff->Draw();

    TCanvas * c4 = new TCanvas("c4");
    McMassPair->Draw();
    RcMassPair->SetLineColor(kRed);
    RcMassPair->Draw("same");

    TCanvas * c5 = new TCanvas("c5");
    TH1* YEff = (TH1*)RcRapidityPair->Clone("YEff");
    YEff->Divide(McRapidityPair);
    YEff->SetTitle("Pair Rapidity Efficiency");
    YEff->Draw();

    TCanvas * c6 = new TCanvas("c6");
    McRapidityPair->Draw();
    RcRapidityPair->SetLineColor(kRed);
    RcRapidityPair->Draw("same");


    TFile file("simulation_plots.root", "RECREATE");
    McPtPair->Write();
    RcPtPair->Write();   
    McMassPair->Write(); 
    RcMassPair->Write(); 
    McRapidityPair->Write(); 
    RcRapidityPair->Write(); 
    hEff->Write();
    mEff->Write();
    YEff->Write();

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



 
}