


void PairEff(){
    TChain * ch = new TChain("FemtoDst");
    ch->Add("JobResults/*");
    TTreeReader myReader(ch);

    
    TTreeReaderArray<float> MCPtVals(myReader, "McTracks.mPt");
    TTreeReaderArray<float> MCEtaVals(myReader, "McTracks.mEta");
    TTreeReaderArray<float> MCPhiVals(myReader, "McTracks.mPhi");

    TTreeReaderArray<short> MCIdx(myReader, "Tracks.mMcIndex");

    TTreeReaderArray<unsigned short> PIDs(myReader, "McTracks.mGeantPID");
    TTreeReaderArray<short> ParentIndex(myReader, "McTracks.mParentIndex");

    cout << "Events = " << ch->GetEntries() << endl;
    auto McPtPair = new TH1F("McPtPair", "McPtPair", 50, 0, 0.4);
    auto RcPtPair = new TH1F("RcPtPair", "RcPtPair", 50, 0, 0.4);

    while (myReader.Next()) {
        
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
        RcPtPair->Fill( rcpair.Pt() );
        
    }


    McPtPair->Draw();
    RcPtPair->SetLineColor(kRed);
    RcPtPair->Draw("same");

    TCanvas *c2 = new TCanvas("c2");
    TH1* hEff = (TH1*)RcPtPair->Clone( "Eff" );
    hEff->Divide( McPtPair );
    hEff->Draw();
 
}