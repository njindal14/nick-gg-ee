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
int ican = 0;
void makeCanvas() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican++ ), "", 900, 600 );
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

void aco() {


    auto * aco = new TH1F("aco", "aco", 10, 0, .2);

    TFile * sim_output = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/simulation_plots_new.root");
    TH1F * MCaco = (TH1F*)sim_output->Get("MCaco");
    TH1F * rcaco = (TH1F*)sim_output->Get("rcaco");

    TFile * mixedeventplots = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/MixedEventplots.root");
    TH1F* mixedAco = (TH1F*)mixedeventplots->Get("acoULSULS");

    TFile * backgroundPlots = new TFile("/Users/Nick/STAR/breit-wheeler/nick-gg-ee/output_root_files/background_plots.root");
    TH1F* acoResample = (TH1F*)backgroundPlots->Get("acoResample");




    TChain * ch = new TChain("PairDst");
    ch->Add("/Users/Nick/STAR/breit-wheeler/rootFiles/pair_dst_Run12UU.root");
    //ch->Add("/Users/Nick/STAR/breit-wheeler/rootFiles/slim_pair_dst_Run10AuAu.root");
    //ch->Add("/Users/Nick/STAR/breit-wheeler/rootFiles/slim_pair_dst_Run11AuAu.root");
    TTreeReader myReader(ch); 
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

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

            if(lv.M() > 0.4 && lv.M() < 1){
                double alpha = 1 - abs(lv1.Phi() - lv2.Phi())/M_PI;
                aco->Fill(alpha);
            }
        }


    }
    
    //normalize mixed event one
    mixedAco->Scale(1/mixedAco->GetEntries());
    for(int i = 0; i < mixedAco->GetNbinsX(); ++i){
        double binwidth = mixedAco->GetBinWidth(i);
        double content = mixedAco->GetBinContent(i)/binwidth;
        mixedAco->SetBinContent(i, content);
    }


    makeCanvas();
    aco->SetLineColor(kBlack);
    aco->SetTitle("#alpha, 0.4 < M < 1");
    aco->Scale(1/aco->GetEntries());
    for(int i = 0; i < aco->GetNbinsX(); ++i){
        double binwidth = aco->GetBinWidth(i);
        double content = aco->GetBinContent(i)/binwidth;
        aco->SetBinContent(i, content);
    }
    aco->GetXaxis()->SetTitle("#alpha");
    aco->GetYaxis()->SetTitle("1/N * dN/d#alpha");
    aco->Draw("PE");
    MCaco->SetLineColor(kBlue);
    rcaco->SetLineColor(kRed);
    MCaco->Draw("PE;same");
    rcaco->Draw("PE;same");
    mixedAco->SetLineColor(kGreen);
    acoResample->SetLineColor(kPink);
    mixedAco->Draw("PE;same");
    acoResample->Draw("PE;same");
    auto acolegend = new TLegend(0.75,0.6,1,0.85);
    acolegend->SetHeader("Legend","C"); // option "C" allows to center the header
    acolegend->AddEntry(aco,"Run 12");
    acolegend->AddEntry(MCaco,"Simulation MC");
    acolegend->AddEntry(rcaco,"Simulation Reco");
    acolegend->AddEntry(mixedAco, "Mixed Event ULSULS");
    acolegend->AddEntry(acoResample, "Starlight resampled pairs");
    acolegend->Draw("same");

    makeCanvas();
    mixedAco->Draw("PE");
    
    makeCanvas();
    TH1F * data_mc_ratio = (TH1F*)aco->Clone();
    data_mc_ratio->Divide(MCaco);
    data_mc_ratio->SetTitle("#alpha ratio, run12/MC");
    data_mc_ratio->Draw("PE");

    makeCanvas();
    MCaco->Draw("PE");




}