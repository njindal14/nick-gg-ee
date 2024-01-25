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
    can->SetRightMargin(0.05);

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

    auto * cos4phivPt = new TH2F("Cos4#phivPt", "cos4#phi distribution vs P_{T}, 0.5 < M_{ee} < 0.8 GeV/c^{2}", 400, -2, 2, 15, 0, 0.3);
    auto * cos3phivPt = new TH2F("Cos3#phivPt", "cos3#phi distribution vs P_{T}, 0.5 < M_{ee} < 0.8 GeV/c^{2}", 400, -2, 2, 15, 0, 0.3);
    auto * cos2phivPt = new TH2F("Cos2#phivPt", "cos2#phi distribution vs P_{T}, 0.5 < M_{ee} < 0.8 GeV/c^{2}", 400, -2, 2, 15, 0, 0.3);
    auto * cosphivPt = new TH2F("Cos#phivPt", "cos#phi distribution vs P_{T}, 0.5 < M_{ee} < 0.8 GeV/c^{2}", 400, -2, 2, 15, 0, 0.3);



    // Open the file containing the tree.
    TFile *myFile = TFile::Open("/Users/Nick/STAR/breit-wheeler/rootFiles/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

    TFile * starlight_out = new TFile("/Users/Nick/STAR/starlight/utils/SL_plotsBetter.root");
    TH1F * slightPt = (TH1F*)starlight_out->Get("mPt");
    TH1F * slightMass = (TH1F*)starlight_out->Get("mMass");
    TH1F * slightEta = (TH1F*)starlight_out->Get("mEta");
    TH2D* slightcos4phipt = (TH2D*)starlight_out->Get("mCos4phivsPT");
    TH2D* slightcos2phipt = (TH2D*)starlight_out->Get("mCos2phivsPT");
    auto * slight_4phiMoments = slightcos4phipt->ProfileY("slight_4phiMoments", 1, -1);
    auto * slight_2phiMoments = slightcos2phipt->ProfileY("slight_2phiMoments", 1, -1);



    TFile * sim_output = new TFile("/Users/Nick/STAR/breit-wheeler/Code/sim_plots.root");

    TH1F * McPt = (TH1F*)sim_output->Get("McPtPair");
    TH1F * McMass = (TH1F*)sim_output->Get("McMassPair");
    TH1F * McEta = (TH1F*)sim_output->Get("McEtaPair");

    TH1F * RecoPt = (TH1F*)sim_output->Get("mPtPairReco");
    TH1F * RecoMass = (TH1F*)sim_output->Get("mMassPairReco");
    TH1F * RecoEta = (TH1F*)sim_output->Get("mEtaPairReco");

    TH1F * mc2phi = (TH1F*)sim_output->Get("m2Ptcos2phimomentsMC");
    TH1F * mc4phi = (TH1F*)sim_output->Get("m2Ptcos4phimomentsMC");
    TH1F * reco2phi = (TH1F*)sim_output->Get("m2Ptcos2phimomentsReco");
    TH1F * reco4phi = (TH1F*)sim_output->Get("m2Ptcos4phimomentsReco");


    TLorentzVector lv1, lv2, lv, lvn;

    //loop through events
    while (myReader.Next()) {

        //values we will want to use for PID cuts
        double chiee = pow( pair->d1_mNSigmaElectron, 2 ) + pow( pair->d2_mNSigmaElectron, 2 );
        double chipipi = pow( pair -> d1_mNSigmaPion, 2) + pow( pair -> d2_mNSigmaPion, 2);
        double c = 3.0e1; //in cm/ns
        double me2 = pow(0.00051,2);
        double dca1 = pair->d1_mDCA;
        double dca2 = pair->d2_mDCA;
        UShort_t mZDCEastVal = pair->mZDCEast;               // ZDC East
        UShort_t mZDCWestVal = pair->mZDCWest;
        Float_t rapidity = pair->mRapidity;
        
        //Lorentz vectors for each pair track and get lorentz sum and diff
        lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.00051 );
        lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.00051 ); 

        lv = lv1 + lv2;
        lvn = lv1 - lv2;
        double parentMass = lv.M();
        //get parent total momentum for each track
        Float_t p1 = lv1.P();
        Float_t p2 = lv2.P();
        //square it
        Float_t p1_2 = pow(p1,2);
        Float_t p2_2 = pow(p2,2);
        Float_t mVertexZVal = pair->mVertexZ;
        UShort_t mGRefMultVal = pair->mGRefMult;  
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

       if( mVertexZVal < 100 && mVertexZVal > -100 && mGRefMultVal <= 4 && chargesumval == 0 && dca1 < 1 && dca2 < 1 && 
       pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag!=0) {

        if( ddTofVal < 0.4 && ddTofVal > -0.4 && chiee < 10 && 3*chiee < chipipi && ddTofVal !=0) {

            if(mPtVal < 0.2){
                mMass->Fill(lv.M());
            }
                
            if(lv.M() < 0.8 && lv.M() > .5) {
                double phival;
                TRandom3 rng(123);
                double uniform_double = rng.Uniform(0., 1.);
                if(uniform_double < 0.5){
                    phival = calc_Phi(lv1,lv2);
                }
                else{
                    phival = calc_Phi(lv2,lv1);
                }     
                mPt->Fill( mPtVal ); 
                mEta->Fill( rapidity );
                cos4phivPt->Fill( 2*cos(4*phival), mPtVal);
                cos3phivPt->Fill( 2*cos(3*phival), mPtVal);
                cos2phivPt->Fill( 2*cos(2*phival), mPtVal);
                cosphivPt->Fill( 2*cos(phival), mPtVal);
            }
            
            mMass->Fill( lv.M() );
 
                 
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
    RecoPt->SetLineColor(kGreen);
    RecoPt->Scale(1/RecoPt->Integral(RecoPt->FindBin(0.04), RecoPt->FindBin(.1),"Width"));
    RecoPt->Draw("PE;same");
    auto ptLegend = new TLegend(0.7,0.6,1,0.9);
    ptLegend->SetHeader("Legend","C"); // option "C" allows to center the header
    ptLegend->AddEntry(mPt,"Data");
    ptLegend->AddEntry(slightPt, "slight.out pT");
    ptLegend->AddEntry(McPt, "MC pT");
    ptLegend->AddEntry(RecoPt, "Reco pT (binned in MC)");
    ptLegend->SetTextSize(.03);
    ptLegend->Draw("same");
    gPad->Print("masterPlots/plot_pairpT.png");


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
    RecoMass->SetLineColor(kGreen);
    RecoMass->Scale(1/RecoMass->Integral(RecoMass->FindBin(0.5), RecoMass->FindBin(1),"Width"));
    RecoMass->Draw("PE;same");
    auto massLegend = new TLegend(0.7,0.6,1,0.9);
    massLegend->SetHeader("Legend","C"); // option "C" allows to center the header
    massLegend->AddEntry(mMass,"Data");
    massLegend->AddEntry(slightMass, "slight.out Mass (1M events)");
    massLegend->AddEntry(McMass, "MC Mass");
    massLegend->AddEntry(RecoMass, "Reco Mass (binned in MC)");
    massLegend->SetTextSize(.03);
    massLegend->Draw("same");
    gPad->Print("masterPlots/plot_pairMass.png");



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
    McEta->SetLineColor(kBlack);
    McEta->Scale(1/McEta->Integral(McEta->FindBin(-1), McEta->FindBin(1),"Width"));
    McEta->Draw("PE;same");
    RecoEta->SetLineColor(kGreen);
    RecoEta->Scale(1/RecoEta->Integral(RecoEta->FindBin(-1), RecoEta->FindBin(1),"Width"));
    RecoEta->Draw("PE;same");
    auto etaLegend = new TLegend(0.7,0.6,1,0.9);
    etaLegend->SetHeader("Legend","C"); // option "C" allows to center the header
    etaLegend->AddEntry(mEta,"Data");
    etaLegend->AddEntry(slightEta, "slight.out Rapidity (1M events)");
    etaLegend->AddEntry(McEta, "MC Rapidity");
    etaLegend->AddEntry(RecoEta, "Reco Rapidity (binned in MC)");
    etaLegend->SetTextSize(.03);
    etaLegend->Draw("same");
    gPad->Print("masterPlots/plot_pairRapidity.png");



    auto *m2Ptcos4phimoments = cos4phivPt->ProfileY("m2Ptcos4phimoments",1, -1);
    auto *m2Ptcos2phimoments = cos2phivPt->ProfileY("m2Ptcos2phimoments",1, -1);

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
    slight_2phiMoments->SetLineColor(kRed);
    slight_2phiMoments->Draw("PE;same");
    mc2phi->SetLineColor(kBlack);
    mc2phi->Draw("PE;same");
    reco2phi->SetLineColor(kGreen);
    reco2phi->Draw("PE;same");
    reco2phi->SetStats(false);
    mc2phi->SetStats(false);
    QED2phi->SetLineColor(kBlue);
    QED2phi->SetLineWidth(6);
    QED2phi->Draw("same");
    m2Ptcos2phimoments->GetXaxis()->SetTitle("Pair pT (GeV/c)");
    m2Ptcos2phimoments->GetYaxis()->SetTitle("A_{2#phi}");
    auto twophiLegend = new TLegend(0.8,0.7,1,0.9);
    twophiLegend->SetHeader("Legend","C"); // option "C" allows to center the header
    twophiLegend->AddEntry(m2Ptcos2phimoments,"Raw Data");
    twophiLegend->AddEntry(slight_2phiMoments, "slight.out A_{2#phi}");
    twophiLegend->AddEntry(mc2phi, "MC A_{2#phi}");
    twophiLegend->AddEntry(reco2phi, "Reco A_{2#phi}");
    twophiLegend->AddEntry(QED2phi, "QED Theory Curve");
    twophiLegend->SetTextSize(.03);
    twophiLegend->Draw("same");
    gPad->Print("masterPlots/plot_pair2phimoments.png");




    makeCanvas2();
    m2Ptcos4phimoments->SetLineColor(kBlue);
    m2Ptcos4phimoments->Draw("PE");
    m2Ptcos4phimoments->GetYaxis()->SetRangeUser(-1, 2.2);
    m2Ptcos4phimoments->SetStats(false);
    m2Ptcos4phimoments->SetMarkerSize(0.7);
    m2Ptcos4phimoments->SetMarkerStyle(kFullDotLarge);
    slight_4phiMoments->SetLineColor(kRed);
    slight_4phiMoments->Draw("PE;same");
    mc4phi->SetLineColor(kBlack);
    mc4phi->Draw("PE;same");
    reco4phi->SetLineColor(kGreen);
    reco4phi->Draw("PE;same");
    QED4phi->SetLineColor(kBlue);
    QED4phi->SetLineWidth(6);
    QED4phi->Draw("same");
    m2Ptcos4phimoments->GetXaxis()->SetTitle("Pair pT (GeV/c)");
    m2Ptcos4phimoments->GetYaxis()->SetTitle("A_{4#phi}");
    reco4phi->SetStats(false);
    mc4phi->SetStats(false);
    auto fourphiLegend = new TLegend(0.8,0.7,1,0.9);
    fourphiLegend->SetHeader("Legend","C"); // option "C" allows to center the header
    fourphiLegend->AddEntry(m2Ptcos4phimoments,"Raw Data");
    fourphiLegend->AddEntry(slight_4phiMoments, "slight.out A_{4#phi}");
    fourphiLegend->AddEntry(mc4phi, "MC A_{4#phi}");
    fourphiLegend->AddEntry(reco4phi, "Reco A_{4#phi}");
    fourphiLegend->AddEntry(QED4phi, "QED Theory Curve");
    fourphiLegend->SetTextSize(.03);
    fourphiLegend->Draw("same");
    gPad->Print("masterPlots/plot_pair4phimoments.png");

    
    makeCanvas2();
    std::cout << "MADE CANVAS";
    TH1F * pteffBetter = (TH1F*)McPt->Clone();
    pteffBetter->Divide(RecoPt);
    pteffBetter->SetTitle("pT Pair Efficiency; pT (GeV/c); Efficiency");
    pteffBetter->Draw("PE");
    gPad->Print("masterPlots/plot_pTEff.png");









}
