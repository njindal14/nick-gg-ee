#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TFile.h>
#include <iostream>
#include <vector>

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

double phiFit(double *x, double *par){
    double fitval;
    fitval = par[0]*(1+par[1]*cos(x[0]) + par[2]*cos(2*x[0]) + par[3]*cos(3*x[0]) + par[4]*cos(4*x[0]));
    return fitval;
}


// Function to unfold phi vs pt histogram using response matrix
TH2D* UnfoldHistogram(TH2D *hMC, const std::vector<TH2D*> &hResponseMatrices) {
    int nBinsPt = hMC->GetNbinsY();
    int nBinsPhi = hMC->GetNbinsX();

    TH2D *hUnfolded = (TH2D*)hMC->Clone("hUnfolded");
    hUnfolded->Reset();

    for (int ptBin = 1; ptBin <= 9; ++ptBin) {
        TH1D *hMC_phi = hMC->ProjectionX(Form("hMC_phi_ptBin%d", ptBin), ptBin, ptBin);
        TH1D *hUnfolded_phi = (TH1D*)hMC_phi->Clone(Form("hUnfolded_phi_ptBin%d", ptBin));
        //hUnfolded_phi->Reset();

        for (int phiBinMC = 1; phiBinMC <= nBinsPhi; ++phiBinMC) {
            double mcValue = hMC_phi->GetBinContent(phiBinMC);
            for (int phiBinReco = 1; phiBinReco <= nBinsPhi; ++phiBinReco) {
                double responseValue = hResponseMatrices[ptBin-1]->GetBinContent(phiBinMC, phiBinReco);
                hUnfolded_phi->AddBinContent(phiBinReco, mcValue * responseValue);
            }
        }

        for (int phiBin = 1; phiBin <= nBinsPhi; ++phiBin) {
            hUnfolded->SetBinContent(phiBin, ptBin, hUnfolded_phi->GetBinContent(phiBin));
        }

        delete hMC_phi;
        delete hUnfolded_phi;
    }

    return hUnfolded;
}


void UnfoldingMatrix(TH2D *hMC, TH2D *hReco, std::vector<TH2D*> &hResponseMatrices) {
    int nBinsPt = hMC->GetNbinsY();
    int nBinsPhi = hMC->GetNbinsX();
    hResponseMatrices.resize(nBinsPt);

    
    // Loop over pt bins
    for (int ptBin = 1; ptBin <= 9; ++ptBin) {
        // Project phi histograms for each pt bin
        TH1D *hMC_phi = hMC->ProjectionX(Form("hMC_phi_ptBin%d", ptBin), ptBin, ptBin);
        TH1D *hReco_phi = hReco->ProjectionX(Form("hReco_phi_ptBin%d", ptBin), ptBin, ptBin);

        // Fit to Fourier series
        TF1 *fitFunc = new TF1("fitFunc", phiFit, -TMath::Pi(), TMath::Pi(), 5);
        fitFunc->SetParNames("a0", "a1", "a2", "a3", "a4");
        hMC_phi->Fit("fitFunc", "", "", -3.14,3.14);


        // Get fitted parameters
        std::vector<double> params(5);
        for (int i = 0; i < 5; ++i) {
            params[i] = fitFunc->GetParameter(i);
        }

        // Flatten MC and reco histograms using fit parameters
        TH1D *hMC_flattened = (TH1D*)hMC_phi->Clone(Form("hMC_flattened_ptBin%d", ptBin));
        TH1D *hReco_flattened = (TH1D*)hReco_phi->Clone(Form("hReco_flattened_ptBin%d", ptBin));

        for (int phiBin = 1; phiBin <= nBinsPhi; ++phiBin) {
            double phi = hMC_flattened->GetXaxis()->GetBinCenter(phiBin);
            double fitValue = phiFit(&phi, params.data());
            if (fitValue != 0) { // Avoid division by zero
                hMC_flattened->SetBinContent(phiBin, hMC_phi->GetBinContent(phiBin) / fitValue);
                hReco_flattened->SetBinContent(phiBin, hReco_phi->GetBinContent(phiBin) / fitValue);
            } else {
                hMC_flattened->SetBinContent(phiBin, 0);
                hReco_flattened->SetBinContent(phiBin, 0);
            }
        }

        // Normalize
        double mcFlattenedIntegral = hMC_flattened->Integral();
        double recoFlattenedIntegral = hReco_flattened->Integral();
        
        if (mcFlattenedIntegral != 0 && recoFlattenedIntegral != 0) {
            hMC_flattened->Scale(1.0 / mcFlattenedIntegral);
            hReco_flattened->Scale(1.0 / recoFlattenedIntegral);
        } else {
            std::cerr << "Error: Integral of flattened histograms is zero in pt bin " << ptBin << std::endl;
            delete hMC_phi;
            delete hReco_phi;
            delete hMC_flattened;
            delete hReco_flattened;
            delete fitFunc;
            continue; // Skip this pt bin if the histogram normalization fails
        }

        // Create a temporary response matrix for the current pt bin
        TMatrixD responseMatrix(nBinsPhi, nBinsPhi);
        responseMatrix.Zero();

        // Reweight flattened MC and reco histograms with embedded signals between -2 and 2 for each harmonic
        for (double a2 = -2; a2 <= 2; a2 += .2) {
            for (double a4 = -2; a4 <= 2; a4 += .2) {
                TH1D *hMC_reweighted = (TH1D*)hMC_flattened->Clone(Form("hMC_reweighted_ptBin%d", ptBin));
                TH1D *hReco_reweighted = (TH1D*)hReco_flattened->Clone(Form("hReco_reweighted_ptBin%d", ptBin));

                for (int phiBin = 1; phiBin <= nBinsPhi; ++phiBin) {
                    double phi = hMC_reweighted->GetXaxis()->GetBinCenter(phiBin);
                    double weight = 1 + a2 * TMath::Cos(2 * phi) + a4 * TMath::Cos(4 * phi);
                    hMC_reweighted->SetBinContent(phiBin, hMC_flattened->GetBinContent(phiBin) * weight);
                    hReco_reweighted->SetBinContent(phiBin, hReco_flattened->GetBinContent(phiBin) * weight);
                }

                /*double mcReweightedIntegral = hMC_reweighted->Integral();
                double recoReweightedIntegral = hReco_reweighted->Integral();

                if (mcReweightedIntegral != 0 && recoReweightedIntegral != 0) {
                    hMC_reweighted->Scale(1.0 / mcReweightedIntegral);
                    hReco_reweighted->Scale(1.0 / recoReweightedIntegral);
                } else {
                    std::cerr << "Error: Integral of reweighted histograms is zero in pt bin " << ptBin << std::endl;
                    delete hMC_reweighted;
                    delete hReco_reweighted;
                    continue; // Skip this combination if the histogram normalization fails
                }*/

                // Compute response and fill response matrix
                for (int phiBinReco = 1; phiBinReco <= nBinsPhi; ++phiBinReco) {
                    for (int phiBinMC = 1; phiBinMC <= nBinsPhi; ++phiBinMC) {
                        responseMatrix(phiBinReco - 1, phiBinMC - 1) += (hReco_reweighted->GetBinContent(phiBinReco) / hMC_reweighted->GetBinContent(phiBinMC));
                    }
                }

                delete hMC_reweighted;
                delete hReco_reweighted;
            }
        }

         // Normalize the response matrix
        for (int i = 0; i < nBinsPhi; ++i) {
            double rowSum = 0;
            for (int j = 0; j < nBinsPhi; ++j) {
                rowSum += responseMatrix(i, j);
            }
            if (rowSum != 0) {
                for (int j = 0; j < nBinsPhi; ++j) {
                    responseMatrix(i, j) /= rowSum;
                }
            }
        }

        // Convert to TH2D for visualization
        hResponseMatrices[ptBin - 1] = new TH2D(Form("hUnfoldingMatrix_ptBin%d", ptBin), 
            Form("Unfolding Matrix for pt Bin %d;#phi_{MC};#phi_{Reco}", ptBin), 
            nBinsPhi, -TMath::Pi(), TMath::Pi(), 
            nBinsPhi, -TMath::Pi(), TMath::Pi());
        for (int i = 1; i <= nBinsPhi; ++i) {
            for (int j = 1; j <= nBinsPhi; ++j) {
                hResponseMatrices[ptBin - 1]->SetBinContent(i, j, responseMatrix(i-1, j-1));
            }
        }
        
        delete hMC_phi;
        delete hReco_phi;
        delete hMC_flattened;
        delete hReco_flattened;
        delete fitFunc;
    }
}

// Function to run the unfolding matrix calculation on input histograms
void RunUnfoldingMatrix(const char* inputFileName, const char* mcHistName, const char* recoHistName) {
    // Open the input file
    TFile *inputFile = TFile::Open(inputFileName);
    if (!inputFile) {
        std::cerr << "Error: Could not open input file " << inputFileName << std::endl;
        return;
    }

    // Get the histograms
    TH2D *hMC = (TH2D*)inputFile->Get(mcHistName);
    TH2D *hReco = (TH2D*)inputFile->Get(recoHistName);
    if (!hMC || !hReco) {
        std::cerr << "Error: Could not find histograms in the input file" << std::endl;
        inputFile->Close();
        return;
    }

    std::vector<TH2D*> hResponseMatrices;
    UnfoldingMatrix(hMC, hReco, hResponseMatrices);

    TFile *outputFile = new TFile("UnfoldingMatrices.root", "RECREATE");
    for (size_t i = 0; i < hResponseMatrices.size(); ++i) {
        if (hResponseMatrices[i]) {
            hResponseMatrices[i]->Write();
        }
    }
    outputFile->Close();

    // Apply the response matrices to unfold the input histogram
    TH2D *hUnfolded = UnfoldHistogram(hMC, hResponseMatrices);

    outputFile = new TFile("UnfoldedHistogram.root", "RECREATE");
    hUnfolded->Write();
    outputFile->Close();

    inputFile->Close();

    for (TH2D* hist : hResponseMatrices) {
        delete hist;
    }
    delete hUnfolded;

    std::cout << "Unfolding matrices and unfolded histogram have been saved to UnfoldingMatrices.root and UnfoldedHistogram.root" << std::endl;
}

// Example usage
int phiUnfolding2() {
    const char* inputFileName = "output_root_files/simulation_plots.root";
    const char* mcHistName = "mcphipt";
    const char* recoHistName = "recophipt";

    RunUnfoldingMatrix(inputFileName, mcHistName, recoHistName);

    return 0;
}