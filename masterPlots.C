//includes
#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TChain.h"
#include "TH2F.h"
#include "FemtoPairFormat.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TEfficiency.h"



//for plotting
int ican6 = 0;
void makeCanvas2() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican6++ ), "", 900, 600 );
    can->SetTopMargin(0.1);
    can->SetBottomMargin(0.1);
    can->SetRightMargin(0.1);

}

static const double hbarc    = 0.197327053;
static const double Znu       = 79;
static const double RNuc       = 6.38;
static const double Anu       = 197;
double A=Anu;
double c = 3;


double formFactor(double * x, double *par)
{
	//par[0] is R, x is t (pt^2)
    //par[0] = pow(A, 1./3.) * 1.16 * (1. - 1.16 * pow(A, -2. / 3.));
	double q    = sqrt(x[0]/2)*c;
	double arg1 = q * par[0] / hbarc;
	double arg2 = hbarc / (q * 1.16 * (1. - 1.16 * pow(A, -2. / 3.)));
	double sph  = (sin(arg1) - arg1 * cos(arg1)) * 3. * arg2 * arg2 * arg2 / double(A);
	double a0   = 0.70;  // [fm]
	return sph / (1. + (a0 * a0 * x[0]/sqrt(2)*c) / (hbarc * hbarc));
}

double ddToffit(double *x, double *par){
    double fitval;

    fitval = par[0]/(par[2]* sqrt(2*M_PI)) * exp(-0.5*pow((x[0]-par[1])/par[2], 2)) + 
            par[3]/(par[5]* sqrt(2*M_PI)) * exp(-0.5*pow((x[0]-par[4])/par[5], 2)) +
         + par[6];                   
    return fitval;
}

double phiFit(double *x, double *par){
    double fitval;
    fitval = par[0]*(1+par[1]*cos(x[0]) + par[2]*cos(2*x[0]) + par[3]*cos(3*x[0]) + par[4]*cos(4*x[0]));
    return fitval;
}

double chiFit(double *x, double *par){
    double fitval;
    fitval = par[0]*exp(x[0]/par[1]);
    return fitval;
}

double calc_Phi( TLorentzVector lv1, TLorentzVector lv2) {
    TLorentzVector lvPlus = lv1 + lv2;
    TLorentzVector lvMinus = lv1 - lv2;
    return lvMinus.DeltaPhi(lvPlus);
}

double* histMoments( TH2F* hist , int n) {
    int nbinsy = hist->GetNbinsY();
    int nbinsx = hist->GetNbinsX();
    double phivals[nbinsx];
    double weights[nbinsx];
    double *cosnphi_moments = new double[nbinsy];
    for (int i = 1; i < nbinsy+1; i++) {
        auto* onedhist = hist->ProjectionX("1dhist", i, i);
	for (int j = 1; j < nbinsx+1; j++) {
            phivals[j-1] = onedhist->GetXaxis()->GetBinCenter(j);
            weights[j-1] = onedhist->GetBinContent(j);
        }
        cosnphi_moments[i-1] = 0;
        for (int k = 0; k < nbinsx; k++){
            if ( onedhist->GetEntries() > 0) {
                cosnphi_moments[i-1] += weights[k] * n * cos(n*phivals[k]) / (onedhist->GetEntries());
            } else { cosnphi_moments[i-1] += 0; }
        }
    }
    
    return cosnphi_moments;
} 

double* moment_error( TH2F* hist, int n ) {
    int nbinsx = hist->GetNbinsX();
    int nbinsy = hist->GetNbinsY();
    double *moment_errors = new double[nbinsy];
    for (int i = 1; i < nbinsy+1; i++) {
        auto* onedhist = hist->ProjectionX("1dhist", i, i);
        auto* momenthist = new TH1F("momenthist", "ncosnphi moments", nbinsx, -1*n, 1*n);
        for (int j = 1; j < nbinsx+1; j++) {
            for (int k = 0; k < onedhist->GetBinContent(j); k++) {
                momenthist->Fill( n * cos( n * (onedhist->GetXaxis()->GetBinCenter(j))));
            }
        }
        moment_errors[i-1] = momenthist->GetMeanError();
    }
    return moment_errors;
}       

double* histYbins( TH2F* hist) {
    int nbinsy = hist->GetNbinsY();
    double *ybinvals = new double[nbinsy];
    for (int i = 1; i < nbinsy+1; i++) {
        ybinvals[i-1] = hist->GetYaxis()->GetBinCenter(i);
    }
    return ybinvals;
}

void masterPlots() {


    auto * mPt = new TH1F("mPt, PID, M_{ee} + TOF Cuts", "Parent Transverse Momentum (GeV/c)", 50, 0, .4);
    auto * mMass = new TH1F("Parent Mass PID + TOF cuts, P_{T} < 0.2 GeV/c", "Parent Mass (GeV/c^{2})", 50, 0, 3);
    auto * mEta = new TH1F("mEta", "Parent #eta", 100, -1, 1);

    auto * cos4phivPt = new TH2F("Cos4#phivPt", "A_{4#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2}", 400, -2, 2, 20, 0, 0.3);
    auto * cos3phivPt = new TH2F("Cos3#phivPt", "A_{3#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2}", 400, -2, 2, 20, 0, 0.3);
    auto * cos2phivPt = new TH2F("Cos2#phivPt", "A_{2#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2}", 400, -2, 2, 20, 0, 0.3);
    auto * cosphivPt = new TH2F("Cos#phivPt", "A_{1#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2}", 400, -2, 2, 20, 0, 0.3);


    auto * cos4phivPt1n1n = new TH2F("Cos4#phivPt1n1n", "A_{4#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2} 1n1n", 400, -2, 2, 20, 0, 0.3);
    auto * cos3phivPt1n1n = new TH2F("Cos3#phivPt1n1n", "A_{3#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2} 1n1n", 400, -2, 2, 20, 0, 0.3);
    auto * cos2phivPt1n1n = new TH2F("Cos2#phivPt1n1n", "A_{2#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2} 1n1n", 400, -2, 2, 20, 0, 0.3);
    auto * cosphivPt1n1n = new TH2F("Cos#phivPt1n1n", "A_{1#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2} 1n1n", 400, -2, 2, 20, 0, 0.3);

    auto * cos4phivPt2nPlus = new TH2F("Cos4#phivPt2nPlus", "A_{4#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2} 2n+", 400, -2, 2, 20, 0, 0.3);
    auto * cos3phivPt2nPlus = new TH2F("Cos3#phivPt2nPlus", "A_{3#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2} 2n+", 400, -2, 2, 20, 0, 0.3);
    auto * cos2phivPt2nPlus = new TH2F("Cos2#phivPt2nPlus", "A_{2#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2} 2n+", 400, -2, 2, 20, 0, 0.3);
    auto * cosphivPt2nPlus = new TH2F("Cos#phivPt2nPlus", "A_{1#phi}, 0.4 < M_{ee} < 0.76 GeV/c^{2} 2n+", 400, -2, 2, 20, 0, 0.3);




    // Open the file containing the tree.
    /*TFile *myFile = TFile::Open("/Users/Nick/STAR/breit-wheeler/rootFiles/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");*/

    TChain * ch = new TChain("PairDst");
    ch->Add("/Users/Nick/STAR/breit-wheeler/rootFiles/pair_dst_Run12UU.root");
    //ch->Add("/Users/Nick/STAR/breit-wheeler/rootFiles/slim_pair_dst_Run10AuAu.root");
    //ch->Add("/Users/Nick/STAR/breit-wheeler/rootFiles/slim_pair_dst_Run11AuAu.root");
    TTreeReader myReader(ch); 
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");
 

    TFile * starlight_out = new TFile("/Users/Nick/STAR/starlight/utils/SL_plotsBetter.root");
    TH1F * slightPt = (TH1F*)starlight_out->Get("mPt");
    TH1F * slightMass = (TH1F*)starlight_out->Get("mMass");
    TH1F * slightEta = (TH1F*)starlight_out->Get("mEta");
    TH2D* slightcos4phipt = (TH2D*)starlight_out->Get("mCos4phivsPT");
    TH2D* slightcos2phipt = (TH2D*)starlight_out->Get("mCos2phivsPT");

    auto * slight_4phiMoments = slightcos4phipt->ProfileY("slight_4phiMoments", 1, -1);
    auto * slight_2phiMoments = slightcos2phipt->ProfileY("slight_2phiMoments", 1, -1);



    TFile * sim_output = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/simulation_plots.root");

    TH1F * McPt = (TH1F*)sim_output->Get("McPtPair");
    TH1F * McMass = (TH1F*)sim_output->Get("McMassPair");
    TH1F * McRapidity = (TH1F*)sim_output->Get("McRapidityPair");

    TH1F * RcPt = (TH1F*)sim_output->Get("RcPtPair");
    TH1F * RcMass = (TH1F*)sim_output->Get("RcMassPair");
    TH1F * RcRapidity = (TH1F*)sim_output->Get("RcRapidityPair");

    TH1F * mc2phi = (TH1F*)sim_output->Get("m2Ptcos2phimomentsMC");
    TH1F * mc4phi = (TH1F*)sim_output->Get("m2Ptcos4phimomentsMC");
    TH1F * reco2phi = (TH1F*)sim_output->Get("m2Ptcos2phimomentsReco");
    TH1F * reco4phi = (TH1F*)sim_output->Get("m2Ptcos4phimomentsReco");

    TH1F * reco2phiRc = (TH1F*)sim_output->Get("m2Ptcos2phimomentsRecoRc");
    TH1F * reco4phiRc = (TH1F*)sim_output->Get("m2Ptcos4phimomentsRecoRc");

    auto * phifit = new TF1("phifit", phiFit, -3.15,3.15,5);
    auto * mPhi1n1n = new TH1F("mPhi1n1n", "#Delta#phi, 1n1n, pT < 0.15 ", 25, -3.14159, 3.14159);
    auto * mPhi2nPlus = new TH1F("mPhi2nPlus", "#Delta#phi, 2n+, pT < 0.15", 25, -3.14159, 3.14159);


    TLorentzVector lv1, lv2, lv, lvn;

    //loop through events
    while (myReader.Next()) {

        //values we will want to use for PID cuts
        double chiee = pow( pair->d1_mNSigmaElectron, 2 ) + pow( pair->d2_mNSigmaElectron, 2 );
        double chipipi = pow( pair -> d1_mNSigmaPion, 2) + pow( pair -> d2_mNSigmaPion, 2);
        double c = 3.0e1; //in cm/ns
        double me2 = pow(0.00051,2);

        UShort_t mZDCEastVal = pair->mZDCEast;               // ZDC East
        UShort_t mZDCWestVal = pair->mZDCWest;
        Float_t rapidity = pair->mRapidity;
        
        //Lorentz vectors for each pair track and get lorentz sum and diff
        lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.00051 );
        lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.00051 ); 

        lv = lv1 + lv2;
        lvn = lv1 - lv2;
        Float_t p1_2 = pow(lv1.P(),2);
        Float_t p2_2 = pow(lv2.P(),2);
        Float_t mVertexZVal = pair->mVertexZ;
        Float_t Tof1 = pair->d1_mTof;     
        Float_t Tof2 = pair->d2_mTof;  
        Float_t len1 = pair->d1_mLength;
        Float_t len2 = pair->d2_mLength; 
        Float_t mPtVal = pair->mPt;
        Float_t dTofVal = Tof1 - Tof2;
        Float_t texp1 = len1/c * sqrt(1 + me2/p1_2);
        Float_t texp2 = len2/c * sqrt(1 + me2/p2_2);
        Float_t dTofexpVal = texp1 - texp2;
        Float_t ddTofVal = dTofVal - dTofexpVal;
        
        int chargesumval = pair->mChargeSum;

        if( fabs(mVertexZVal) < 100 && pair->mGRefMult <= 4 && chargesumval == 0 && pair->d1_mDCA < 1 && pair->d2_mDCA < 1 && 
        pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag!=0 && fabs(ddTofVal) < 0.4 && chiee < 10 && 3*chiee < chipipi && ddTofVal !=0) {
            
            if(mPtVal < 0.15){
                mMass->Fill(lv.M());

            }
                
            if(lv.M() < 0.76 && lv.M() > .4) {
                double phival;
                TRandom3 rng(123);
                double uniform_double = rng.Uniform(0., 1.);
                if(uniform_double < 0.5){
                    phival = calc_Phi(lv1,lv2);
                }
                else{
                    phival = calc_Phi(lv2,lv1);
                }
                if(lv.Pt() < 0.15){    
                    if(mZDCEastVal < 220 && mZDCWestVal < 220) {mPhi1n1n->Fill(phival);}
                    else if (mZDCEastVal > 220 && mZDCWestVal > 220) {mPhi2nPlus->Fill(phival);}
                }

                mPt->Fill( mPtVal ); 
                mEta->Fill( rapidity );
                cos4phivPt->Fill( 2*cos(4*phival), mPtVal);
                cos3phivPt->Fill( 2*cos(3*phival), mPtVal);
                cos2phivPt->Fill( 2*cos(2*phival), mPtVal);
                cosphivPt->Fill( 2*cos(phival), mPtVal);

                if(mZDCEastVal < 85 && mZDCWestVal < 85){
                    cos4phivPt1n1n->Fill( 2*cos(4*phival), mPtVal);
                    cos3phivPt1n1n->Fill( 2*cos(3*phival), mPtVal);
                    cos2phivPt1n1n->Fill( 2*cos(2*phival), mPtVal);
                    cosphivPt1n1n->Fill( 2*cos(phival), mPtVal);
                }
                else if(mZDCEastVal > 85 && mZDCWestVal > 85){
                    cos4phivPt2nPlus->Fill( 2*cos(4*phival), mPtVal);
                    cos3phivPt2nPlus->Fill( 2*cos(3*phival), mPtVal);
                    cos2phivPt2nPlus->Fill( 2*cos(2*phival), mPtVal);
                    cosphivPt2nPlus->Fill( 2*cos(phival), mPtVal);
                }
            } 
                 
        }

    }


    makeCanvas2();
    mPt->SetTitle("Pair pT from Various Sources");
    mPt->GetXaxis()->SetTitle("Pair pT (GeV/c)");
    mPt->GetYaxis()->SetTitle("Normalized Counts");
    mPt->Scale(1/mPt->Integral(mPt->FindBin(0.04), mPt->FindBin(.1),"Width"));
    mPt->SetMarkerSize(0.7);
    mPt->SetMarkerStyle(kFullDotLarge);
    mPt->Draw("PE");
    mPt->SetStats(false);
    slightPt->SetLineColor(kRed);
    slightPt->Scale(1/slightPt->Integral(slightPt->FindBin(0.04), slightPt->FindBin(.1),"Width"));
    slightPt->Draw("PE;same");
    McPt->SetLineColor(kBlack);
    McPt->Scale(1/McPt->Integral(McPt->FindBin(0.04), McPt->FindBin(.1),"Width"));
    McPt->Draw("PE;same");
    RcPt->SetLineColor(kGreen);
    RcPt->Scale(1/RcPt->Integral(RcPt->FindBin(0.04), RcPt->FindBin(.1),"Width"));
    RcPt->Draw("PE;same");
    auto ptLegend = new TLegend(0.7,0.6,1,0.9);
    ptLegend->SetHeader("Legend","C"); // option "C" allows to center the header
    ptLegend->AddEntry(mPt,"Data");
    ptLegend->AddEntry(slightPt, "slight.out pT");
    ptLegend->AddEntry(McPt, "MC pT");
    ptLegend->AddEntry(RcPt, "Reco pT (binned in MC)");
    ptLegend->SetTextSize(.03);
    ptLegend->Draw("same");
    //gPad->Print("masterPlots/plot_pairpT.png");


    makeCanvas2();
    mMass->SetTitle("Pair Mass from Various Sources");
    mMass->GetXaxis()->SetTitle("Pair Mass (GeV/c^2)");
    mMass->GetYaxis()->SetTitle("Normalized Counts");
    mMass->Scale(1/mMass->Integral(mMass->FindBin(0.5), mMass->FindBin(1), "Width"));
    mMass->SetStats(false);
    mMass->SetMarkerSize(0.7);
    mMass->SetMarkerStyle(kFullDotLarge);
    mMass->Draw("PE");
    slightMass->SetLineColor(kRed);
    slightMass->Scale(1/slightMass->Integral(slightMass->FindBin(0.5), slightMass->FindBin(1),"Width"));
    slightMass->Draw("PE;same");
    McMass->SetLineColor(kBlack);
    McMass->Scale(1/McMass->Integral(McMass->FindBin(0.5), McMass->FindBin(1),"Width"));
    McMass->Draw("PE;same");
    RcMass->SetLineColor(kGreen);
    RcMass->Scale(1/RcMass->Integral(RcMass->FindBin(0.5), RcMass->FindBin(1),"Width"));
    RcMass->Draw("PE;same");
    auto massLegend = new TLegend(0.7,0.6,1,0.9);
    massLegend->SetHeader("Legend","C"); // option "C" allows to center the header
    massLegend->AddEntry(mMass,"Data");
    massLegend->AddEntry(slightMass, "slight.out Mass (1M events)");
    massLegend->AddEntry(McMass, "MC Mass");
    massLegend->AddEntry(RcMass, "Reco Mass (binned in MC)");
    massLegend->SetTextSize(.03);
    massLegend->Draw("same");
    //gPad->Print("masterPlots/plot_pairMass.png");



    makeCanvas2();
    mEta->SetTitle("Pair Rapidity from Various Sources");
    mEta->GetXaxis()->SetTitle("Pair Rapidity");
    mEta->GetYaxis()->SetTitle("Normalized Counts");
    mEta->Scale(1/mEta->Integral(mEta->FindBin(-1), mEta->FindBin(1),"Width"));
    mEta->Draw("PE");
    mEta->SetStats(false);
    mEta->SetMarkerSize(0.7);
    mEta->SetMarkerStyle(kFullDotLarge);
    slightEta->SetLineColor(kRed);
    slightEta->Scale(1/slightEta->Integral(slightEta->FindBin(-1), slightEta->FindBin(1),"Width"));
    slightEta->Draw("PE;same");
    McRapidity->SetLineColor(kBlack);
    McRapidity->Scale(1/McRapidity->Integral(McRapidity->FindBin(-1), McRapidity->FindBin(1),"Width"));
    McRapidity->Draw("PE;same");
    RcRapidity->SetLineColor(kGreen);
    RcRapidity->Scale(1/RcRapidity->Integral(RcRapidity->FindBin(-1), RcRapidity->FindBin(1),"Width"));
    RcRapidity->Draw("PE;same");
    auto etaLegend = new TLegend(0.7,0.6,1,0.9);
    etaLegend->SetHeader("Legend","C"); // option "C" allows to center the header
    etaLegend->AddEntry(mEta,"Data");
    etaLegend->AddEntry(slightEta, "slight.out Rapidity (1M events)");
    etaLegend->AddEntry(McRapidity, "MC Rapidity");
    etaLegend->AddEntry(RcRapidity, "Reco Rapidity (binned in MC)");
    etaLegend->SetTextSize(.03);
    etaLegend->Draw("same");
    //gPad->Print("masterPlots/plot_pairRapidity.png");



    auto *m2Ptcos4phimoments = cos4phivPt->ProfileY("m2Ptcos4phimoments",1, -1);
    auto *m2Ptcos2phimoments = cos2phivPt->ProfileY("m2Ptcos2phimoments",1, -1);

    auto *m2Ptcos4phimoments1n1n = cos4phivPt1n1n->ProfileY("m2Ptcos4phimoments1n1n",1, -1);
    auto *m2Ptcos2phimoments1n1n = cos2phivPt1n1n->ProfileY("m2Ptcos2phimoments1n1n",1, -1);
    
    auto *m2Ptcos4phimoments2nPlus = cos4phivPt2nPlus->ProfileY("m2Ptcos4phimoments2nPlus",1, -1);
    auto *m2Ptcos2phimoments2nPlus = cos2phivPt2nPlus->ProfileY("m2Ptcos2phimoments2nPlus",1, -1);




    double fourphix[] = {0.00023,0.00211,0.00491,0.00749,0.0101	,0.0131	,0.0164	,0.0201	,0.0239	,0.0304	,0.0356	,0.0416	,0.047	,0.0538	,0.0618	,0.0674	,0.0727	,0.0791	,0.084	,0.0877	,0.0912	,0.0943	,0.098	,0.101	,0.105	,0.108	,0.111	,0.114	,0.117	,0.119	,0.123	,0.126	,0.129	,0.132	,0.137	,0.14	,0.144	,0.149	,0.155	,0.161	,0.168	,0.174	,0.179	,0.184	,0.188	,0.192	,0.196	,0.199	,0.2};
    double fourphiy[] = {0.092-.1,0.068-.1	,0.028-.1	,-0.008-.1,-0.048-.1,-0.06-.1	,-0.06-.1	,-0.052-.1,-0.036-.1,-0.008-.1,0.024-.1	,0.06-.1	,0.104-.1	,0.156-.1	,0.232-.1	,0.288-.1	,0.356	-.1,0.448-.1	,0.528-.1	,0.596-.1	,0.672-.1	,0.744-.1	,0.832-.1	,0.916-.1	,1-.1	,1.08-.1	,1.17-.1	,1.24-.1	,1.3-.1	,1.35-.1	,1.41-.1	,1.47-.1	,1.52-.1	,1.55-.1	,1.6-.1	,1.63-.1	,1.67-.1	,1.7-.1	,1.73-.1	,1.75-.1	,1.76-.1	,1.76-.1	,1.78-.1	,1.79-.1	,1.81-.1	,1.83-.1	,1.85-.1	,1.88-.1	,1.88-.1};
    TGraph * QED4phi = new TGraph(sizeof(fourphix)/sizeof(fourphix[0]), fourphix,fourphiy);

    double twophix[] = {0.00384,0.00565,0.00701,0.00949,0.012,0.0145,0.0172,0.0201,0.0237,0.0271,0.0303,0.0346,0.0384,0.0416,0.045,0.0477,0.0499,0.0524,0.0556,0.059,0.0626,0.0653,0.0676,0.0705,0.0739,0.0764,0.0784,0.0807,0.0832,0.0856,0.0881,0.0908,0.0929,0.0949,0.0976,0.1,0.103,0.105,0.108,0.111,0.114,0.117,0.119,0.121,0.123,0.126,0.129,0.132	,0.136	,0.14,0.145,0.149,0.151,0.155,0.158,0.163,0.168,0.172,0.177,0.182,0.186,0.189,0.193,0.196,0.197};
    double twophiy[] = {0.00765,0.0115,0.0115,0.0153,0.0191,0.0229,0.0344	,0.0421,0.0535,0.065,0.0765,0.0918	,0.115	,0.13	,0.145	,0.164	,0.187	,0.203	,0.241	,0.268	,0.306	,0.333	,0.367	,0.417	,0.467	,0.512	,0.547	,0.593	,0.646	,0.704	,0.78	,0.837	,0.891	,0.956	,1.02	,1.12	,1.19	,1.26	,1.32	,1.41	,1.5	,1.56	,1.6	,1.64	,1.68	,1.72	,1.76	,1.78	,1.81	,1.82	,1.83	,1.83	,1.82	,1.82	,1.82	,1.82	,1.82	,1.83	,1.84	,1.85	,1.86	,1.87	,1.88	,1.9	,1.9};
    TGraph * QED2phi = new TGraph(sizeof(twophix)/sizeof(twophix[0]), twophix,twophiy);


    makeCanvas2();
    m2Ptcos2phimoments->SetLineColor(kBlue);
    m2Ptcos2phimoments->Draw("PE");
    m2Ptcos2phimoments->GetYaxis()->SetRangeUser(-1, 2.2);
    m2Ptcos2phimoments->SetStats(false);
    m2Ptcos2phimoments->SetMarkerSize(0.7);
    m2Ptcos2phimoments->SetMarkerStyle(kFullDotLarge);

    m2Ptcos2phimoments1n1n->SetLineColor(kMagenta);
    //m2Ptcos2phimoments1n1n->Draw("PE;same");
    m2Ptcos2phimoments1n1n->SetStats(false);
    m2Ptcos2phimoments1n1n->SetMarkerSize(0.7);
    m2Ptcos2phimoments1n1n->SetMarkerStyle(kFullDotLarge);

    m2Ptcos2phimoments2nPlus->SetLineColor(kOrange+10);
    //m2Ptcos2phimoments2nPlus->Draw("PE;same");
    m2Ptcos2phimoments2nPlus->SetStats(false);
    m2Ptcos2phimoments2nPlus->SetMarkerSize(0.7);
    m2Ptcos2phimoments2nPlus->SetMarkerStyle(kFullDotLarge);
    
    slight_2phiMoments->SetLineColor(kGreen);
    slight_2phiMoments->Draw("PE;same");
    mc2phi->SetLineColor(kBlack);
    //mc2phi->Draw("PE;same");
    reco2phi->SetLineColor(kRed);
    reco2phi->Draw("PE;same");
    reco2phiRc->SetLineColor(kYellow);
    reco2phiRc->Draw("PE;same");
    reco2phi->SetStats(false);
    reco2phiRc->SetStats(false);
    mc2phi->SetStats(false);
    QED2phi->SetLineColor(kBlue);
    QED2phi->SetLineWidth(6);
    QED2phi->Draw("same");
    m2Ptcos2phimoments->GetXaxis()->SetTitle("Pair pT (GeV/c)");
    m2Ptcos2phimoments->GetYaxis()->SetTitle("A_{2#phi}");
    auto twophiLegend = new TLegend(0.75,0.6,1,0.85);
    twophiLegend->SetHeader("Legend","C"); // option "C" allows to center the header
    twophiLegend->AddEntry(m2Ptcos2phimoments,"Raw Data Run 12");
    //twophiLegend->AddEntry(m2Ptcos2phimoments1n1n,"Raw Data Run 12 1n1n");
    //twophiLegend->AddEntry(m2Ptcos2phimoments2nPlus,"Raw Data Run 12 2n+");
    twophiLegend->AddEntry(slight_2phiMoments, "slight.out A_{2#phi}");
    //twophiLegend->AddEntry(mc2phi, "MC A_{2#phi}");
    twophiLegend->AddEntry(reco2phi, "Reco A_{2#phi} (MC bins)");
    twophiLegend->AddEntry(reco4phiRc, "Reco A_{2#phi} (RC bins)");
    twophiLegend->AddEntry(QED2phi, "QED Theory Curve");
    twophiLegend->SetTextSize(.03);
    twophiLegend->Draw("same");
    //gPad->Print("masterPlots/plot_pair2phimomentsRun10to12.png");




    makeCanvas2();
    m2Ptcos4phimoments->SetLineColor(kBlue);
    m2Ptcos4phimoments->Draw("PE");
    m2Ptcos4phimoments->GetYaxis()->SetRangeUser(-1, 2.2);
    m2Ptcos4phimoments->SetStats(false);
    m2Ptcos4phimoments->SetMarkerSize(0.7);
    m2Ptcos4phimoments->SetMarkerStyle(kFullDotLarge);

    m2Ptcos4phimoments1n1n->SetLineColor(kMagenta);
    //m2Ptcos4phimoments1n1n->Draw("PE;same");
    m2Ptcos4phimoments1n1n->SetStats(false);
    m2Ptcos4phimoments1n1n->SetMarkerSize(0.7);
    m2Ptcos4phimoments1n1n->SetMarkerStyle(kFullDotLarge);

    m2Ptcos4phimoments2nPlus->SetLineColor(kOrange+10);
    //m2Ptcos4phimoments2nPlus->Draw("PE;same");
    m2Ptcos4phimoments2nPlus->SetStats(false);
    m2Ptcos4phimoments2nPlus->SetMarkerSize(0.7);
    m2Ptcos4phimoments2nPlus->SetMarkerStyle(kFullDotLarge);
    
    slight_4phiMoments->SetLineColor(kGreen);
    slight_4phiMoments->Draw("PE;same");
    mc4phi->SetLineColor(kBlack);
    //mc4phi->Draw("PE;same");
    reco4phi->SetLineColor(kRed);
    reco4phi->Draw("PE;same");
    reco4phiRc->SetLineColor(kYellow);
    reco4phiRc->Draw("PE;same");
    QED4phi->SetLineColor(kBlue);
    QED4phi->SetLineWidth(6);
    QED4phi->Draw("same");
    m2Ptcos4phimoments->GetXaxis()->SetTitle("Pair pT (GeV/c)");
    m2Ptcos4phimoments->GetYaxis()->SetTitle("A_{4#phi}");
    reco4phi->SetStats(false);
    reco4phiRc->SetStats(false);
    mc4phi->SetStats(false);
    auto fourphiLegend = new TLegend(0.75,0.6,1,0.85);
    fourphiLegend->SetHeader("Legend","C"); // option "C" allows to center the header
    fourphiLegend->AddEntry(m2Ptcos4phimoments,"Raw Data Run 12");
    //fourphiLegend->AddEntry(m2Ptcos4phimoments1n1n,"Raw Data Run 12 1n1n");
    //fourphiLegend->AddEntry(m2Ptcos4phimoments2nPlus,"Raw Data Run 12 2n+");
    fourphiLegend->AddEntry(slight_4phiMoments, "slight.out A_{4#phi}");
    //fourphiLegend->AddEntry(mc4phi, "MC A_{4#phi}");
    fourphiLegend->AddEntry(reco4phi, "Reco A_{4#phi} (MC bins)");
    fourphiLegend->AddEntry(reco4phiRc, "Reco A_{4#phi} (RC bins)");
    fourphiLegend->AddEntry(QED4phi, "QED Theory Curve");
    fourphiLegend->SetTextSize(.03);
    fourphiLegend->Draw("same");
    //gPad->Print("masterPlots/plot_pair4phimomentsRun10to12.png");

    makeCanvas2();
    mPhi1n1n->GetXaxis()->SetTitle("#Delta #phi");
    mPhi1n1n->SetTitle("Modulation for ZDC 1-3n");
    mPhi1n1n->GetYaxis()->SetTitle("Counts");
    mPhi1n1n->Fit("phifit", "", "", -3.15,3.15);
    gStyle->SetOptFit(1111);
    mPhi1n1n->Draw("PE");

    makeCanvas2();
    mPhi2nPlus->GetXaxis()->SetTitle("#Delta #phi");
    mPhi2nPlus->SetTitle("Modulation for ZDC tail (3n+)");
    mPhi2nPlus->Fit("phifit", "", "", -3.15,3.15);
    mPhi2nPlus->GetYaxis()->SetTitle("Counts");
    gStyle->SetOptFit(1111);
    mPhi2nPlus->Draw("PE");

    makeCanvas2();
    slight_4phiMoments->GetXaxis()->SetTitle("pT (GeV/c)");
    slight_4phiMoments->GetYaxis()->SetTitle("A_{4#phi}");
    slight_4phiMoments->Draw("PE");
    

    makeCanvas2();
    slight_2phiMoments->GetXaxis()->SetTitle("pT (GeV/c)");
    slight_2phiMoments->GetYaxis()->SetTitle("A_{2#phi}");
    slight_2phiMoments->Draw("PE");


    TH1F * twophihist = new TH1F("QED2phi", "QED2phi", QED2phi->GetN(), 0, .3); // the histogram (you should set the number of bins, the title etc)
    TH1F * fourphihist = new TH1F("QED4phi", "QED4phi", QED4phi->GetN(), 0, .3); // the histogram (you should set the number of bins, the title etc)

    for(int i=0; i < QED2phi->GetN(); ++i) {
        double x,y;
        twophihist->SetBinContent(i, QED2phi->Eval(twophix[i]));
    }   
    for(int i=0; i < QED4phi->GetN(); ++i) {
        double x,y;
        fourphihist->SetBinContent(i, QED4phi->Eval(fourphix[i]));
    }   

    makeCanvas2();
    twophihist->SetMarkerSize(0.7);
    twophihist->SetMarkerStyle(kFullDotLarge);
    twophihist->GetXaxis()->SetTitle("#Delta #phi");
    twophihist->GetYaxis()->SetTitle("A_{2#phi}");
    twophihist->Draw("P");

    makeCanvas2();
    fourphihist->SetMarkerSize(0.7);
    fourphihist->SetMarkerStyle(kFullDotLarge);
    fourphihist->GetXaxis()->SetTitle("#Delta #phi");
    fourphihist->GetYaxis()->SetTitle("A_{4#phi}");
    fourphihist->Draw("P");
}
