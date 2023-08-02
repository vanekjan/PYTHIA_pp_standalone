// main91.cc is a part of the PYTHIA event generator.
// Copyright (C) 2022 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: analysis; root;

// This is a simple test program.
// It studies the charged multiplicity distribution at the LHC.
// Modified by Rene Brun, Axel Naumann and Bernhard Meirose
// to use ROOT for histogramming.

// Stdlib header file for input and output.
#include <iostream>
#include <fstream>
#include <sstream>

// ROOT, for histogramming.
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TF1.h"
//#include "TRandom3.h"
//#include "TString.h"
//#include "TMath.h"


// ROOT, for saving file.
//#include "TFile.h"

using namespace std;


int findBinPt( TVector3 Mom, float const *bins, const int nBins )
{
  int bin = -1;

  //find pT bin of Lambda
  for(int j = 0; j < nBins; j++) //loop over pT bins
  {
    if(Mom.Pt() > bins[j] && Mom.Pt() <= bins[j+1])
    {
      bin = j;
      break; //stop after bin is found
    }
  }
  
  return bin;
}

int findBinEta( TVector3 Mom, float const *bins, const int nBins )
{
  int bin = -1;

  //find pT bin of Lambda
  for(int j = 0; j < nBins; j++) //loop over pT bins
  {
    if(Mom.Eta() > bins[j] && Mom.Eta() <= bins[j+1])
    {
      bin = j;
      break; //stop after bin is found
    }
  }
  
  return bin;
}


struct L_MC {
  // general information
  // lab frame
  enum {
    nL_MAX_MC=10 //check that thi is enough
  };

  //number of L
  Int_t nL_MC;

  //L momentum for eta and pT bins and for Delta eta and Delta phi
  Float_t L_px_MC[nL_MAX_MC];
  Float_t L_py_MC[nL_MAX_MC];
  Float_t L_pz_MC[nL_MAX_MC];

  //proton momntum in L rest frame
  Float_t p_pxStar_MC[nL_MAX_MC];
  Float_t p_pyStar_MC[nL_MAX_MC];
  Float_t p_pzStar_MC[nL_MAX_MC];

  //for cuts
  Float_t L_decayL_MC[nL_MAX_MC];
  
  Float_t p_pT_MC[nL_MAX_MC];
  Float_t p_eta_MC[nL_MAX_MC];
  Float_t p_phi_MC[nL_MAX_MC];
  
  Float_t pi_pT_MC[nL_MAX_MC];
  Float_t pi_eta_MC[nL_MAX_MC];
  Float_t pi_phi_MC[nL_MAX_MC];
  
  //L charge
  Int_t L_charge_MC[nL_MAX_MC]; //1 - L, -1 - Lbar
};


struct L_from_pairs {
  // general information
  // lab frame
  enum {
    nL_MAX=100 //check that thi is enough
  };

  //number of L
  Int_t nL;

  //L momentum for eta and pT bins and for Delta eta and Delta phi
  Float_t L_px[nL_MAX];
  Float_t L_py[nL_MAX];
  Float_t L_pz[nL_MAX];
  Float_t L_Minv[nL_MAX];

  //proton momntum in L rest frame
  Float_t p_pxStar[nL_MAX];
  Float_t p_pyStar[nL_MAX];
  Float_t p_pzStar[nL_MAX];

  //for cuts
  Float_t L_decayL[nL_MAX];
  Float_t L_theta[nL_MAX]; //pointing angle
  
  Float_t p_pT[nL_MAX];
  Float_t p_eta[nL_MAX];
  Float_t p_phi[nL_MAX];
  
  Float_t pi_pT[nL_MAX];
  Float_t pi_eta[nL_MAX];
  Float_t pi_phi[nL_MAX];
  
  //L charge
  Int_t L_charge[nL_MAX]; //1 - L, -1 - Lbar
  Int_t L_US_LS_flag[nL_MAX]; // 1 - US, 0 - LS
};


void Read_PYTHIA_tree( int mEnergy = 200 )
{
  cout<<"Start"<<endl;

  ifstream fileList;

  if(mEnergy == 510)
  {
    //MB without strict TOF matching
    fileList.open("/gpfs/mnt/gpfs02/eic/janvanek/PYTHIA_pp/Lambda_MB/output/pp_510/tree_1B_events/fileList.list");

  }
  else if(mEnergy == 200)
  {
    //MB without strict TOF matching
    fileList.open("/gpfs/mnt/gpfs02/eic/janvanek/PYTHIA_pp/Lambda_MB/output/pp_200/tree_1B_events/fileList.list");
    
  }
  else if (mEnergy == 0) //testing file
  {
    fileList.open("/gpfs/mnt/gpfs02/eic/janvanek/PYTHIA_pp/Lambda_MB/output/pp_200/tree_1B_events/testList.list"); //one root file for quick tests
  }
  else
  {
    cout<<"Not a valid colliison energy! Aborting!"<<endl;
    return;
  }

  
  TChain *L_MC_chain = new TChain("L_MC_tree");
  
  TChain *L_from_pairs_chain = new TChain("L_tree");

  string fileFromList;


  while(getline(fileList, fileFromList))
  {
    L_MC_chain->Add(fileFromList.c_str());

    L_from_pairs_chain->Add(fileFromList.c_str());
  }


  L_MC Lambda_MC;
  
  //number of L
  L_MC_chain->SetBranchAddress( "nL_MC", &Lambda_MC.nL_MC);

  //L momentum for eta and pT bins and for Delta eta and Delta phi
  L_MC_chain->SetBranchAddress( "L_px_MC", &Lambda_MC.L_px_MC );
  L_MC_chain->SetBranchAddress( "L_py_MC", &Lambda_MC.L_py_MC );
  L_MC_chain->SetBranchAddress( "L_pz_MC", &Lambda_MC.L_pz_MC );

  //proton momntum in L rest frame
  L_MC_chain->SetBranchAddress( "p_pxStar_MC", &Lambda_MC.p_pxStar_MC );
  L_MC_chain->SetBranchAddress( "p_pyStar_MC", &Lambda_MC.p_pyStar_MC );
  L_MC_chain->SetBranchAddress( "p_pzStar_MC", &Lambda_MC.p_pzStar_MC );

  //for cuts
  L_MC_chain->SetBranchAddress( "L_decayL_MC", &Lambda_MC.L_decayL_MC );
  
  L_MC_chain->SetBranchAddress( "p_pT_MC", &Lambda_MC.p_pT_MC );
  L_MC_chain->SetBranchAddress( "p_eta_MC", &Lambda_MC.p_eta_MC );
  L_MC_chain->SetBranchAddress( "p_phi_MC", &Lambda_MC.p_phi_MC );
  
  L_MC_chain->SetBranchAddress( "pi_pT_MC", &Lambda_MC.pi_pT_MC );
  L_MC_chain->SetBranchAddress( "pi_eta_MC", &Lambda_MC.pi_eta_MC );
  L_MC_chain->SetBranchAddress( "pi_phi_MC", &Lambda_MC.pi_phi_MC );
  
  //L charge
  L_MC_chain->SetBranchAddress( "L_charge_MC", &Lambda_MC.L_charge_MC ); //1 - L, -1 - Lbar
  
  //-------------------------------------------------------------------------------
  
  L_from_pairs Lambda_from_pair;  
  
  //number of L
  L_from_pairs_chain->SetBranchAddress( "nL", &Lambda_from_pair.nL );

  //L momentum for eta and pT bins and for Delta eta and Delta phi
  L_from_pairs_chain->SetBranchAddress( "L_px", &Lambda_from_pair.L_px );
  L_from_pairs_chain->SetBranchAddress( "L_py", &Lambda_from_pair.L_py );
  L_from_pairs_chain->SetBranchAddress( "L_pz", &Lambda_from_pair.L_pz );
  L_from_pairs_chain->SetBranchAddress( "L_Minv", &Lambda_from_pair.L_Minv );

  //proton momntum in L rest frame
  L_from_pairs_chain->SetBranchAddress( "p_pxStar", &Lambda_from_pair.p_pxStar );
  L_from_pairs_chain->SetBranchAddress( "p_pyStar", &Lambda_from_pair.p_pyStar );
  L_from_pairs_chain->SetBranchAddress( "p_pzStar", &Lambda_from_pair.p_pzStar );

  //for cuts
  L_from_pairs_chain->SetBranchAddress( "L_decayL", &Lambda_from_pair.L_decayL );
  L_from_pairs_chain->SetBranchAddress( "L_theta", &Lambda_from_pair.L_theta );
  
  L_from_pairs_chain->SetBranchAddress( "p_pT", &Lambda_from_pair.p_pT );
  L_from_pairs_chain->SetBranchAddress( "p_eta", &Lambda_from_pair.p_eta );
  L_from_pairs_chain->SetBranchAddress( "p_phi", &Lambda_from_pair.p_phi );
  
  L_from_pairs_chain->SetBranchAddress( "pi_pT", &Lambda_from_pair.pi_pT );
  L_from_pairs_chain->SetBranchAddress( "pi_eta", &Lambda_from_pair.pi_eta );
  L_from_pairs_chain->SetBranchAddress( "pi_phi", &Lambda_from_pair.pi_phi );
  
  //L charge
  L_from_pairs_chain->SetBranchAddress( "L_charge", &Lambda_from_pair.L_charge ); //1 - L, -1 - Lbar
  L_from_pairs_chain->SetBranchAddress( "L_US_LS_flag", &Lambda_from_pair.L_US_LS_flag ); // 1 - US, 0 - LS
  
  //____________________________________________________________________
  
  

  // Create file on which histogram(s) can be saved.
  TFile* outFile = new TFile("/gpfs/mnt/gpfs02/eic/janvanek/PYTHIA_pp/Lambda_MB/output/pp_200/new_out_hist/output_Lambda_pp_200_MB_1B_events_hists_work.root", "RECREATE");

  // Book histogram.
  //TH1F *mult = new TH1F("mult","charged multiplicity", 100, -0.5, 799.5);
  
  const int nPtBins = 8;
  float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5.};
  
  const int nPtBins_corr = 2;
  float const pT_bins_corr[nPtBins_corr+1] = { 0.5, 1.5, 5.};

  const int nEtaBins = 3;
  float const eta_bins[nEtaBins+1] = { -1, -0.4, 0.4, 1 };


  //histograms  
  
  //true L and Lbar
  TH1F *L0_L0bar_cosThetaProdPlane = new TH1F("L0_L0bar_cosThetaProdPlane", "L0_L0bar_cosThetaProdPlane", 10, -1, 1);
  
  TH1F *L0_L0bar_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0bar_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane = new TH1F("L0_L0_cosThetaProdPlane", "L0_L0_cosThetaProdPlane", 10, -1, 1);
  
  TH1F *L0_L0_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane = new TH1F("L0bar_L0bar_cosThetaProdPlane", "L0bar_L0bar_cosThetaProdPlane", 10, -1, 1);
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];
 
  TH1F *L0bar_L0bar_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];
  
  //delta eta vs. delta phi histograms for re-weighing of ME
  //TH2F *L0_L0bar_delta_eta_vs_delta_phi_hist = new TH2F("L0_L0bar_delta_eta_vs_delta_phi_hist", "L0_L0bar_delta_eta_vs_delta_phi_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0_L0_delta_eta_vs_delta_phi_hist = new TH2F("L0_L0_delta_eta_vs_delta_phi_hist", "L0_L0_delta_eta_vs_delta_phi_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0bar_L0bar_delta_eta_vs_delta_phi_hist = new TH2F("L0bar_L0bar_delta_eta_vs_delta_phi_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  
  //TH3F *L0_L0bar_delta_eta_vs_delta_phi_vs_delta_pT_hist = new TH3F("L0_L0bar_delta_eta_vs_delta_phi_vs_delta_pT_hist", "L0_L0bar_delta_eta_vs_delta_phi_vs_delta_pT_hist", 100, 0, 2, 100, 0, TMath::TwoPi(), 100, 0, 5);
  //TH3F *L0_L0_delta_eta_vs_delta_phi_vs_delta_pT_hist = new TH3F("L0_L0_delta_eta_vs_delta_phi_vs_delta_pT_hist", "L0_L0_delta_eta_vs_delta_phi_vs_delta_pT_hist", 100, 0, 2, 100, 0, TMath::TwoPi(), 100, 0, 5);
  //TH3F *L0bar_L0bar_delta_eta_vs_delta_phi_vs_delta_pT_hist = new TH3F("L0bar_L0bar_delta_eta_vs_delta_phi_vs_delta_pT_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_vs_delta_pT_hist", 100, 0, 2, 100, 0, TMath::TwoPi(), 100, 0, 5);
  
  //TH1F *L0_L0bar_delta_pT = new TH1F("L0_L0bar_delta_pT", "L0_L0bar_delta_pT", 100, 0, 5);
  
  TH2F *L0_L0bar_eta1_vs_eta2_hist = new TH2F("L0_L0bar_eta1_vs_eta2_hist", "L0_L0bar_eta1_vs_eta2_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_phi1_vs_phi2_hist = new TH2F("L0_L0bar_phi1_vs_phi2_hist", "L0_L0bar_phi1_vs_phi2_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pT1_vs_pT2_hist = new TH2F("L0_L0bar_pT1_vs_pT2_hist", "L0_L0bar_pT1_vs_pT2_hist", 20, 0, 5, 20, 0, 5);
  
  
  //-------------------------------------------------------------------------------------------------------------------------------
  //mixed event histograms
/*  
  TH1F *L0_L0bar_cosThetaProdPlane_ME = new TH1F("L0_L0bar_cosThetaProdPlane_ME", "L0_L0bar_cosThetaProdPlane_ME", 10, -1, 1);
  
  TH1F *L0_L0bar_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0bar_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_ME = new TH1F("L0_L0_cosThetaProdPlane_ME", "L0_L0_cosThetaProdPlane_ME", 10, -1, 1);
  
  TH1F *L0_L0_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME = new TH1F("L0bar_L0bar_cosThetaProdPlane_ME", "L0bar_L0bar_cosThetaProdPlane_ME", 10, -1, 1);
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr];
 
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];
*/  
  //-----------------------------------------------------------------------------------------------
  
  //mixed event histograms with delta eta vs. delta phi re-weighing
  TH1F *L0_L0bar_cosThetaProdPlane_ME_weight = new TH1F("L0_L0bar_cosThetaProdPlane_ME_weight", "L0_L0bar_cosThetaProdPlane_ME_weight", 10, -1, 1);
  
  TH1F *L0_L0bar_cosThetaProdPlane_ME_weight_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0bar_cosThetaProdPlane_ME_weight_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_ME_weight = new TH1F("L0_L0_cosThetaProdPlane_ME_weight", "L0_L0_cosThetaProdPlane_ME_weight", 10, -1, 1);
  
  TH1F *L0_L0_cosThetaProdPlane_ME_weight_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_ME_weight_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_weight = new TH1F("L0bar_L0bar_cosThetaProdPlane_ME_weight", "L0bar_L0bar_cosThetaProdPlane_ME_weight", 10, -1, 1);
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_weight_pT_hist[nPtBins_corr][nPtBins_corr];
 
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_weight_eta_hist[nEtaBins][nEtaBins];
  
  //delta eta vs. delta phi histograms from ME for re-weighing of ME
  //TH2F *L0_L0bar_delta_eta_vs_delta_phi_ME_hist = new TH2F("L0_L0bar_delta_eta_vs_delta_phi_ME_hist", "L0_L0bar_delta_eta_vs_delta_phi_ME_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0_L0_delta_eta_vs_delta_phi_ME_hist = new TH2F("L0_L0_delta_eta_vs_delta_phi_ME_hist", "L0_L0_delta_eta_vs_delta_phi_ME_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0bar_L0bar_delta_eta_vs_delta_phi_ME_hist = new TH2F("L0bar_L0bar_delta_eta_vs_delta_phi_ME_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_ME_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  
  //TH3F *L0_L0bar_delta_eta_vs_delta_phi_vs_delta_pT_ME_hist = new TH3F("L0_L0bar_delta_eta_vs_delta_phi_vs_delta_pT_ME_hist", "L0_L0bar_delta_eta_vs_delta_phi_vs_delta_pT_ME_hist", 100, 0, 2, 100, 0, TMath::TwoPi(), 100, 0, 5);
  //TH3F *L0_L0_delta_eta_vs_delta_phi_vs_delta_pT_ME_hist = new TH3F("L0_L0_delta_eta_vs_delta_phi_vs_delta_pT_ME_hist", "L0_L0_delta_eta_vs_delta_phi_vs_delta_pT_ME_hist", 100, 0, 2, 100, 0, TMath::TwoPi(), 100, 0, 5);
  //TH3F *L0bar_L0bar_delta_eta_vs_delta_phi_vs_delta_pT_ME_hist = new TH3F("L0bar_L0bar_delta_eta_vs_delta_phi_vs_delta_pT_ME_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_vs_delta_pT_ME_hist", 100, 0, 2, 100, 0, TMath::TwoPi(), 100, 0, 5);
  
  //TH1F *L0_L0bar_delta_pT_ME = new TH1F("L0_L0bar_delta_pT_ME", "L0_L0bar_delta_pT_ME", 100, 0, 5);
  
  TH2F *L0_L0bar_eta1_vs_eta2_ME_hist = new TH2F("L0_L0bar_eta1_vs_eta2_ME_hist", "L0_L0bar_eta1_vs_eta2_ME_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_phi1_vs_phi2_ME_hist = new TH2F("L0_L0bar_phi1_vs_phi2_ME_hist", "L0_L0bar_phi1_vs_phi2_ME_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pT1_vs_pT2_ME_hist = new TH2F("L0_L0bar_pT1_vs_pT2_ME_hist", "L0_L0bar_pT1_vs_pT2_ME_hist", 20, 0, 5, 20, 0, 5);
    
  //_____________________________________________________________________________________________________________________________________________________________________
  
  //p pi pairs - to simulate combinatorial background
  //invariant mass histograms  
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US_LS[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs

  TH2F *L0_inv_mass_vs_L0_inv_mass_US[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_US_LS[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs

  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs
  
  
  //unlike-sign (signal + background)
  TH1F *L0_L0bar_cosThetaProdPlane_US = new TH1F("L0_L0bar_cosThetaProdPlane_US", "L0_L0bar_cosThetaProdPlane_US", 10, -1, 1);
  
  TH1F *L0_L0bar_cosThetaProdPlane_US_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0bar_cosThetaProdPlane_US_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_US = new TH1F("L0_L0_cosThetaProdPlane_US", "L0_L0_cosThetaProdPlane_US", 10, -1, 1);
  
  TH1F *L0_L0_cosThetaProdPlane_US_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_US_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US = new TH1F("L0bar_L0bar_cosThetaProdPlane_US", "L0bar_L0bar_cosThetaProdPlane_US", 10, -1, 1);
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_pT_hist[nPtBins_corr][nPtBins_corr];
 
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_eta_hist[nEtaBins][nEtaBins];
  
  //delta eta vs. delta phi histograms for re-weighing of ME
  TH2F *L0_L0bar_delta_eta_vs_delta_phi_US_hist = new TH2F("L0_L0bar_delta_eta_vs_delta_phi_US_hist", "L0_L0bar_delta_eta_vs_delta_phi_US_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0_L0_delta_eta_vs_delta_phi_US_hist = new TH2F("L0_L0_delta_eta_vs_delta_phi_US_hist", "L0_L0_delta_eta_vs_delta_phi_US_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0bar_L0bar_delta_eta_vs_delta_phi_US_hist = new TH2F("L0bar_L0bar_delta_eta_vs_delta_phi_US_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_US_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  
  //------------------------------------------------------------------------------------------
  
  //US matched to LS (background)
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS = new TH1F("L0_L0bar_cosThetaProdPlane_US_LS", "L0_L0bar_cosThetaProdPlane_US_LS", 10, -1, 1);
  
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_US_LS = new TH1F("L0_L0_cosThetaProdPlane_US_LS", "L0_L0_cosThetaProdPlane_US_LS", 10, -1, 1);
  
  TH1F *L0_L0_cosThetaProdPlane_US_LS_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_US_LS_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS = new TH1F("L0bar_L0bar_cosThetaProdPlane_US_LS", "L0bar_L0bar_cosThetaProdPlane_US_LS", 10, -1, 1);
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_pT_hist[nPtBins_corr][nPtBins_corr];
 
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_eta_hist[nEtaBins][nEtaBins];
  
  //delta eta vs. delta phi histograms for re-weighing of ME
  TH2F *L0_L0bar_delta_eta_vs_delta_phi_US_LS_hist = new TH2F("L0_L0bar_delta_eta_vs_delta_phi_US_LS_hist", "L0_L0bar_delta_eta_vs_delta_phi_US_LS_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0_L0_delta_eta_vs_delta_phi_US_LS_hist = new TH2F("L0_L0_delta_eta_vs_delta_phi_US_LS_hist", "L0_L0_delta_eta_vs_delta_phi_US_LS_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_hist = new TH2F("L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  
  //-------------------------------------------------------------------------------------------------------------------------------
  
/*
  //mixed event  
  //US
  TH1F *L0_L0bar_cosThetaProdPlane_US_ME = new TH1F("L0_L0bar_cosThetaProdPlane_US_ME", "L0_L0bar_cosThetaProdPlane_US_ME", 10, -1, 1);
  
  TH1F *L0_L0bar_cosThetaProdPlane_US_ME_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0bar_cosThetaProdPlane_US_ME_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_US_ME = new TH1F("L0_L0_cosThetaProdPlane_US_ME", "L0_L0_cosThetaProdPlane_US_ME", 10, -1, 1);
  
  TH1F *L0_L0_cosThetaProdPlane_US_ME_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_US_ME_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_ME = new TH1F("L0bar_L0bar_cosThetaProdPlane_US_ME", "L0bar_L0bar_cosThetaProdPlane_US_ME", 10, -1, 1);
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_ME_pT_hist[nPtBins_corr][nPtBins_corr];
 
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_ME_eta_hist[nEtaBins][nEtaBins];
  
  //-------------------------------------------------------------------------------------------
  
  //US matched to LS (background)
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_ME = new TH1F("L0_L0bar_cosThetaProdPlane_US_LS_ME", "L0_L0bar_cosThetaProdPlane_US_LS_ME", 10, -1, 1);
  
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_ME_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_ME_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_US_LS_ME = new TH1F("L0_L0_cosThetaProdPlane_US_LS_ME", "L0_L0_cosThetaProdPlane_US_LS_ME", 10, -1, 1);
  
  TH1F *L0_L0_cosThetaProdPlane_US_LS_ME_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_US_LS_ME_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_ME = new TH1F("L0bar_L0bar_cosThetaProdPlane_US_LS_ME", "L0bar_L0bar_cosThetaProdPlane_US_LS_ME", 10, -1, 1);
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_ME_pT_hist[nPtBins_corr][nPtBins_corr];
 
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_ME_eta_hist[nEtaBins][nEtaBins];
  
  //--------------------------------------------------------------------------------
*/  
  //mixed event after re-weighing
  //US
  TH1F *L0_L0bar_cosThetaProdPlane_US_ME_weight = new TH1F("L0_L0bar_cosThetaProdPlane_US_ME_weight", "L0_L0bar_cosThetaProdPlane_US_ME_weight", 10, -1, 1);
  
  TH1F *L0_L0bar_cosThetaProdPlane_US_ME_weight_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_US_ME_weight = new TH1F("L0_L0_cosThetaProdPlane_US_ME_weight", "L0_L0_cosThetaProdPlane_US_ME_weight", 10, -1, 1);
  
  TH1F *L0_L0_cosThetaProdPlane_US_ME_weight_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_US_ME_weight_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_ME_weight = new TH1F("L0bar_L0bar_cosThetaProdPlane_US_ME_weight", "L0bar_L0bar_cosThetaProdPlane_US_ME_weight", 10, -1, 1);
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT_hist[nPtBins_corr][nPtBins_corr];
 
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[nEtaBins][nEtaBins];
  
  //delta eta vs. delta phi histograms for re-weighing of ME
  TH2F *L0_L0bar_delta_eta_vs_delta_phi_US_ME_hist = new TH2F("L0_L0bar_delta_eta_vs_delta_phi_US_ME_hist", "L0_L0bar_delta_eta_vs_delta_phi_US_ME_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0_L0_delta_eta_vs_delta_phi_US_ME_hist = new TH2F("L0_L0_delta_eta_vs_delta_phi_US_ME_hist", "L0_L0_delta_eta_vs_delta_phi_US_ME_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0bar_L0bar_delta_eta_vs_delta_phi_US_ME_hist = new TH2F("L0bar_L0bar_delta_eta_vs_delta_phi_US_ME_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_US_ME_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  
  //-----------------------------------------------------------------------------------
  //US matched to LS (background)
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_ME_weight = new TH1F("L0_L0bar_cosThetaProdPlane_US_LS_ME_weight", "L0_L0bar_cosThetaProdPlane_US_LS_ME_weight", 10, -1, 1);
  
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_US_LS_ME_weight = new TH1F("L0_L0_cosThetaProdPlane_US_LS_ME_weight", "L0_L0_cosThetaProdPlane_US_LS_ME_weight", 10, -1, 1);
  
  TH1F *L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight = new TH1F("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight", "L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight", 10, -1, 1);
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_hist[nPtBins_corr][nPtBins_corr];
 
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_hist[nEtaBins][nEtaBins];
    
  //delta eta vs. delta phi histograms for re-weighing of ME
  TH2F *L0_L0bar_delta_eta_vs_delta_phi_US_LS_ME_hist = new TH2F("L0_L0bar_delta_eta_vs_delta_phi_US_LS_ME_hist", "L0_L0bar_delta_eta_vs_delta_phi_US_LS_ME_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0_L0_delta_eta_vs_delta_phi_US_LS_ME_hist = new TH2F("L0_L0_delta_eta_vs_delta_phi_US_LS_ME_hist", "L0_L0_delta_eta_vs_delta_phi_US_LS_ME_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_ME_hist = new TH2F("L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_ME_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_ME_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
    
 
  //_____________________________________________________________________________________________________________________________________________________________________

  //true MC histograms after cuts	  
  TH1F *L0_L0bar_cosThetaProdPlane_cuts = new TH1F("L0_L0bar_cosThetaProdPlane_cuts", "L0_L0bar_cosThetaProdPlane_cuts", 10, -1, 1);
  
  TH1F *L0_L0bar_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0bar_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_cuts = new TH1F("L0_L0_cosThetaProdPlane_cuts", "L0_L0_cosThetaProdPlane_cuts", 10, -1, 1);
  
  TH1F *L0_L0_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_cuts = new TH1F("L0bar_L0bar_cosThetaProdPlane_cuts", "L0bar_L0bar_cosThetaProdPlane_cuts", 10, -1, 1);  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  //delta eta vs. delta phi histograms from MC for re-weighing of ME
  //TH2F *L0_L0bar_delta_eta_vs_delta_phi_cuts_hist = new TH2F("L0_L0bar_delta_eta_vs_delta_phi_cuts_hist", "L0_L0bar_delta_eta_vs_delta_phi_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0_L0_delta_eta_vs_delta_phi_cuts_hist = new TH2F("L0_L0_delta_eta_vs_delta_phi_cuts_hist", "L0_L0_delta_eta_vs_delta_phi_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0bar_L0bar_delta_eta_vs_delta_phi_cuts_hist = new TH2F("L0bar_L0bar_delta_eta_vs_delta_phi_cuts_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  
  //TH3F *L0_L0bar_delta_eta_vs_delta_phi_delta_pT_cuts_hist = new TH3F("L0_L0bar_delta_eta_vs_delta_phi_delta_pT_cuts_hist", "L0_L0bar_delta_eta_vs_delta_phi_delta_pT_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi(), 100, 0, 5);
  //TH3F *L0_L0_delta_eta_vs_delta_phi_delta_pT_cuts_hist = new TH3F("L0_L0_delta_eta_vs_delta_phi_delta_pT_cuts_hist", "L0_L0_delta_eta_vs_delta_phi_delta_pT_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi(), 100, 0, 5);
  //TH3F *L0bar_L0bar_delta_eta_vs_delta_delta_pT_phi_cuts_hist = new TH3F("L0bar_L0bar_delta_eta_vs_delta_delta_pT_phi_cuts_hist", "L0bar_L0bar_delta_eta_vs_delta_delta_pT_phi_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi(), 100, 0, 5);
  
  //TH1F *L0_L0bar_delta_pT_cuts = new TH1F("L0_L0bar_delta_pT_cuts", "L0_L0bar_delta_pT_cuts", 100, 0, 5);
  
  TH2F *L0_L0bar_eta1_vs_eta2_cuts_hist = new TH2F("L0_L0bar_eta1_vs_eta2_cuts_hist", "L0_L0bar_eta1_vs_eta2_cuts_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_phi1_vs_phi2_cuts_hist = new TH2F("L0_L0bar_phi1_vs_phi2_cuts_hist", "L0_L0bar_phi1_vs_phi2_cuts_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pT1_vs_pT2_cuts_hist = new TH2F("L0_L0bar_pT1_vs_pT2_cuts_hist", "L0_L0bar_pT1_vs_pT2_cuts_hist", 20, 0, 5, 20, 0, 5);
  

  TH2F *L0_L0bar_p1_pT1_vs_p2_pT2_cuts_hist = new TH2F("L0_L0bar_p1_pT1_vs_p2_pT2_cuts_hist", "L0_L0bar_p1_pT1_vs_p2_pT2_cuts_hist", 100, 0, 5, 100, 0, 5);
  TH2F *L0_L0bar_pi1_pT1_vs_pi2_pT2_cuts_hist = new TH2F("L0_L0bar_pi1_pT1_vs_pi2_pT2_cuts_hist", "L0_L0bar_pi1_pT1_vs_pi2_pT2_cuts_hist", 100, 0, 1.5, 100, 0, 1.5);
   
  TH2F *L0_L0bar_p1_eta1_vs_p2_eta2_cuts_hist = new TH2F("L0_L0bar_p1_eta1_vs_p2_eta2_cuts_hist", "L0_L0bar_p1_eta1_vs_p2_eta2_cuts_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_pi1_eta1_vs_pi2_eta2_cuts_hist = new TH2F("L0_L0bar_pi1_eta1_vs_pi2_eta2_cuts_hist", "L0_L0bar_pi1_eta1_vs_pi2_eta2_cuts_hist", 20, -1, 1, 20, -1, 1);
  
  TH2F *L0_L0bar_p1_phi1_vs_p2_phi2_cuts_hist = new TH2F("L0_L0bar_p1_phi1_vs_p2_phi2_cuts_hist", "L0_L0bar_p1_phi1_vs_p2_phi2_cuts_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pi1_phi1_vs_pi2_phi2_cuts_hist = new TH2F("L0_L0bar_pi1_phi1_vs_pi2_phi2_cuts_hist", "L0_L0bar_pi1_phi1_vs_pi2_phi2_cuts_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());

  
  //--------------------------------------------------------------------------------------------------------------------------------------
/*  
  //mixed event histograms
  TH1F *L0_L0bar_cosThetaProdPlane_ME_cuts = new TH1F("L0_L0bar_cosThetaProdPlane_ME_cuts", "L0_L0bar_cosThetaProdPlane_ME_cuts", 10, -1, 1);
  
  TH1F *L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_ME_cuts = new TH1F("L0_L0_cosThetaProdPlane_ME_cuts", "L0_L0_cosThetaProdPlane_ME_cuts", 10, -1, 1);
  
  TH1F *L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_cuts = new TH1F("L0bar_L0bar_cosThetaProdPlane_ME_cuts", "L0bar_L0bar_cosThetaProdPlane_ME_cuts", 10, -1, 1);
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
 
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[nEtaBins][nEtaBins];
*/  
  //--------------------------------------------------------------------------------------------
  
  //mixed event histograms after re-weighing
  TH1F *L0_L0bar_cosThetaProdPlane_ME_weight_cuts = new TH1F("L0_L0bar_cosThetaProdPlane_ME_weight_cuts", "L0_L0bar_cosThetaProdPlane_ME_weight_cuts", 10, -1, 1);
  
  TH1F *L0_L0bar_cosThetaProdPlane_ME_weight_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0bar_cosThetaProdPlane_ME_weight_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_ME_weight_cuts = new TH1F("L0_L0_cosThetaProdPlane_ME_weight_cuts", "L0_L0_cosThetaProdPlane_ME_weight_cuts", 10, -1, 1);
  
  TH1F *L0_L0_cosThetaProdPlane_ME_weight_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_ME_weight_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_weight_cuts = new TH1F("L0bar_L0bar_cosThetaProdPlane_ME_weight_cuts", "L0bar_L0bar_cosThetaProdPlane_ME_weight_cuts", 10, -1, 1);
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_weight_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
 
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_weight_eta_cuts_hist[nEtaBins][nEtaBins];
  
  //delta eta vs. delta phi histograms from ME for re-weighing of ME
  //TH2F *L0_L0bar_delta_eta_vs_delta_phi_ME_cuts_hist = new TH2F("L0_L0bar_delta_eta_vs_delta_phi_ME_cuts_hist", "L0_L0bar_delta_eta_vs_delta_phi_ME_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0_L0_delta_eta_vs_delta_phi_ME_cuts_hist = new TH2F("L0_L0_delta_eta_vs_delta_phi_ME_cuts_hist", "L0_L0_delta_eta_vs_delta_phi_ME_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0bar_L0bar_delta_eta_vs_delta_phi_ME_cuts_hist = new TH2F("L0bar_L0bar_delta_eta_vs_delta_phi_ME_cuts_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_ME_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  
  //TH3F *L0_L0bar_delta_eta_vs_delta_phi_delta_pT_ME_cuts_hist = new TH3F("L0_L0bar_delta_eta_vs_delta_phi_delta_pT_ME_cuts_hist", "L0_L0bar_delta_eta_vs_delta_phi_delta_pT_ME_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi(), 100, 0, 5);
  //TH3F *L0_L0_delta_eta_vs_delta_phi_delta_pT_ME_cuts_hist = new TH3F("L0_L0_delta_eta_vs_delta_phi_delta_pT_ME_cuts_hist", "L0_L0_delta_eta_vs_delta_phi_delta_pT_ME_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi(), 100, 0, 5);
  //TH3F *L0bar_L0bar_delta_eta_vs_delta_phi_delta_pT_ME_cuts_hist = new TH3F("L0bar_L0bar_delta_eta_vs_delta_phi_delta_pT_ME_cuts_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_delta_pT_ME_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi(), 100, 0, 5);
  
  //TH1F *L0_L0bar_delta_pT_ME_cuts = new TH1F("L0_L0bar_delta_pT_ME_cuts", "L0_L0bar_delta_pT_ME_cuts", 100, 0, 5);
  
  TH2F *L0_L0bar_eta1_vs_eta2_ME_cuts_hist = new TH2F("L0_L0bar_eta1_vs_eta2_ME_cuts_hist", "L0_L0bar_eta1_vs_eta2_ME_cuts_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_phi1_vs_phi2_ME_cuts_hist = new TH2F("L0_L0bar_phi1_vs_phi2_ME_cuts_hist", "L0_L0bar_phi1_vs_phi2_ME_cuts_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pT1_vs_pT2_ME_cuts_hist = new TH2F("L0_L0bar_pT1_vs_pT2_ME_cuts_hist", "L0_L0bar_pT1_vs_pT2_ME_cuts_hist", 20, 0, 5, 20, 0, 5);
  
  TH2F *L0_L0bar_eta1_vs_eta2_ME_cuts_weight_hist = new TH2F("L0_L0bar_eta1_vs_eta2_ME_cuts_weight_hist", "L0_L0bar_eta1_vs_eta2_ME_cuts_weight_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_phi1_vs_phi2_ME_cuts_weight_hist = new TH2F("L0_L0bar_phi1_vs_phi2_ME_cuts_weight_hist", "L0_L0bar_phi1_vs_phi2_ME_cuts_weight_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pT1_vs_pT2_ME_cuts_weight_hist = new TH2F("L0_L0bar_pT1_vs_pT2_ME_cuts_weight_hist", "L0_L0bar_pT1_vs_pT2_ME_cuts_weight_hist", 20, 0, 5, 20, 0, 5);
  
  TH2F *L0_L0bar_p1_pT1_vs_p2_pT2_ME_cuts_hist = new TH2F("L0_L0bar_p1_pT1_vs_p2_pT2_ME_cuts_hist", "L0_L0bar_p1_pT1_vs_p2_pT2_ME_cuts_hist", 100, 0, 5, 100, 0, 5);
  TH2F *L0_L0bar_pi1_pT1_vs_pi2_pT2_ME_cuts_hist = new TH2F("L0_L0bar_pi1_pT1_vs_pi2_pT2_ME_cuts_hist", "L0_L0bar_pi1_pT1_vs_pi2_pT2_ME_cuts_hist", 100, 0, 1.5, 100, 0, 1.5);
  
  TH2F *L0_L0bar_p1_eta1_vs_p2_eta2_ME_cuts_hist = new TH2F("L0_L0bar_p1_eta1_vs_p2_eta2_ME_cuts_hist", "L0_L0bar_p1_eta1_vs_p2_eta2_ME_cuts_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_pi1_eta1_vs_pi2_eta2_ME_cuts_hist = new TH2F("L0_L0bar_pi1_eta1_vs_pi2_eta2_ME_cuts_hist", "L0_L0bar_pi1_eta1_vs_pi2_eta2_ME_cuts_hist", 20, -1, 1, 20, -1, 1);
  
  TH2F *L0_L0bar_p1_phi1_vs_p2_phi2_ME_cuts_hist = new TH2F("L0_L0bar_p1_phi1_vs_p2_phi2_ME_cuts_hist", "L0_L0bar_p1_phi1_vs_p2_phi2_ME_cuts_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pi1_phi1_vs_pi2_phi2_ME_cuts_hist = new TH2F("L0_L0bar_pi1_phi1_vs_pi2_phi2_ME_cuts_hist", "L0_L0bar_pi1_phi1_vs_pi2_phi2_ME_cuts_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  
  
  TH2F *L0_L0bar_p1_pT1_vs_p2_pT2_ME_cuts_hist_weight = new TH2F("L0_L0bar_p1_pT1_vs_p2_pT2_ME_cuts_hist_weight", "L0_L0bar_p1_pT1_vs_p2_pT2_ME_cuts_hist_weight", 100, 0, 5, 100, 0, 5);
  TH2F *L0_L0bar_pi1_pT1_vs_pi2_pT2_ME_cuts_hist_weight = new TH2F("L0_L0bar_pi1_pT1_vs_pi2_pT2_ME_cuts_hist_weight", "L0_L0bar_pi1_pT1_vs_pi2_pT2_ME_cuts_hist_weight", 100, 0, 1.5, 100, 0, 1.5);
  
  TH2F *L0_L0bar_p1_eta1_vs_p2_eta2_ME_cuts_hist_weight = new TH2F("L0_L0bar_p1_eta1_vs_p2_eta2_ME_cuts_hist_weight", "L0_L0bar_p1_eta1_vs_p2_eta2_ME_cuts_hist_weight", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_pi1_eta1_vs_pi2_eta2_ME_cuts_hist_weight = new TH2F("L0_L0bar_pi1_eta1_vs_pi2_eta2_ME_cuts_hist_weight", "L0_L0bar_pi1_eta1_vs_pi2_eta2_ME_cuts_hist_weight", 20, -1, 1, 20, -1, 1);
  
  TH2F *L0_L0bar_p1_phi1_vs_p2_phi2_ME_cuts_hist_weight = new TH2F("L0_L0bar_p1_phi1_vs_p2_phi2_ME_cuts_hist_weight", "L0_L0bar_p1_phi1_vs_p2_phi2_ME_cuts_hist_weight", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pi1_phi1_vs_pi2_phi2_ME_cuts_hist_weight = new TH2F("L0_L0bar_pi1_phi1_vs_pi2_phi2_ME_cuts_hist_weight", "L0_L0bar_pi1_phi1_vs_pi2_phi2_ME_cuts_hist_weight", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  
  
  //_____________________________________________________________________________________________________________________________________________________________________
  
  //p pi pairs after cuts - to simulate combinatorial background
  //invariant mass histograms  
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US_cuts[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs

  TH2F *L0_inv_mass_vs_L0_inv_mass_US_cuts[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_US_LS_cuts[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs

  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs
  
  //unlike-sign (signal + background)
  TH1F *L0_L0bar_cosThetaProdPlane_US_cuts = new TH1F("L0_L0bar_cosThetaProdPlane_US_cuts", "L0_L0bar_cosThetaProdPlane_US_cuts", 10, -1, 1);
  
  TH1F *L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_US_cuts = new TH1F("L0_L0_cosThetaProdPlane_US_cuts", "L0_L0_cosThetaProdPlane_US_cuts", 10, -1, 1);
  
  TH1F *L0_L0_cosThetaProdPlane_US_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_US_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_cuts = new TH1F("L0bar_L0bar_cosThetaProdPlane_US_cuts", "L0bar_L0bar_cosThetaProdPlane_US_cuts", 10, -1, 1);
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
 
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  //delta eta vs. delta phi histograms for re-weighing of ME
  TH2F *L0_L0bar_delta_eta_vs_delta_phi_US_cuts_hist = new TH2F("L0_L0bar_delta_eta_vs_delta_phi_US_cuts_hist", "L0_L0bar_delta_eta_vs_delta_phi_US_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0_L0_delta_eta_vs_delta_phi_US_cuts_hist = new TH2F("L0_L0_delta_eta_vs_delta_phi_US_cuts_hist", "L0_L0_delta_eta_vs_delta_phi_US_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0bar_L0bar_delta_eta_vs_delta_phi_US_cuts_hist = new TH2F("L0bar_L0bar_delta_eta_vs_delta_phi_US_cuts_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_US_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  
  //-----------------------------------------------------------------------------------------------
  
  //US matched to LS (background)
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_cuts = new TH1F("L0_L0bar_cosThetaProdPlane_US_LS_cuts", "L0_L0bar_cosThetaProdPlane_US_LS_cuts", 10, -1, 1);
  
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_US_LS_cuts = new TH1F("L0_L0_cosThetaProdPlane_US_LS_cuts", "L0_L0_cosThetaProdPlane_US_LS_cuts", 10, -1, 1);
  
  TH1F *L0_L0_cosThetaProdPlane_US_LS_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_US_LS_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_cuts = new TH1F("L0bar_L0bar_cosThetaProdPlane_US_LS_cuts", "L0bar_L0bar_cosThetaProdPlane_US_LS_cuts", 10, -1, 1);
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
 
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[nEtaBins][nEtaBins];
  
  //delta eta vs. delta phi histograms for re-weighing of ME
  TH2F *L0_L0bar_delta_eta_vs_delta_phi_US_LS_cuts_hist = new TH2F("L0_L0bar_delta_eta_vs_delta_phi_US_LS_cuts_hist", "L0_L0bar_delta_eta_vs_delta_phi_US_LS_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0_L0_delta_eta_vs_delta_phi_US_LS_cuts_hist = new TH2F("L0_L0_delta_eta_vs_delta_phi_US_LS_cuts_hist", "L0_L0_delta_eta_vs_delta_phi_US_LS_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_cuts_hist = new TH2F("L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_cuts_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  
  //----------------------------------------------------------------------------------------
/*  
  //mixed event
  //US
  TH1F *L0_L0bar_cosThetaProdPlane_US_ME_cuts = new TH1F("L0_L0bar_cosThetaProdPlane_US_ME_cuts", "L0_L0bar_cosThetaProdPlane_US_ME_cuts", 10, -1, 1);
  
  TH1F *L0_L0bar_cosThetaProdPlane_US_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0bar_cosThetaProdPlane_US_ME_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_US_ME_cuts = new TH1F("L0_L0_cosThetaProdPlane_US_ME_cuts", "L0_L0_cosThetaProdPlane_US_ME_cuts", 10, -1, 1);
  
  TH1F *L0_L0_cosThetaProdPlane_US_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_US_ME_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_ME_cuts = new TH1F("L0bar_L0bar_cosThetaProdPlane_US_ME_cuts", "L0bar_L0bar_cosThetaProdPlane_US_ME_cuts", 10, -1, 1);
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
 
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_ME_eta_cuts_hist[nEtaBins][nEtaBins];
  
  //-----------------------------------------------------------------------------------------
  
  //US matched to LS (background)
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_ME_cuts = new TH1F("L0_L0bar_cosThetaProdPlane_US_LS_ME_cuts", "L0_L0bar_cosThetaProdPlane_US_LS_ME_cuts", 10, -1, 1);
  
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_US_LS_ME_cuts = new TH1F("L0_L0_cosThetaProdPlane_US_LS_ME_cuts", "L0_L0_cosThetaProdPlane_US_LS_ME_cuts", 10, -1, 1);
  
  TH1F *L0_L0_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_ME_cuts = new TH1F("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_cuts", "L0bar_L0bar_cosThetaProdPlane_US_LS_ME_cuts", 10, -1, 1);
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
 
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[nEtaBins][nEtaBins];
*/  
  //------------------------------------------------------------------------------------------------
  
  //mixed event after re-weighing
  //US
  TH1F *L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts = new TH1F("L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts", "L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts", 10, -1, 1);
  
  TH1F *L0_L0bar_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_US_ME_weight_cuts = new TH1F("L0_L0_cosThetaProdPlane_US_ME_weight_cuts", "L0_L0_cosThetaProdPlane_US_ME_weight_cuts", 10, -1, 1);
  
  TH1F *L0_L0_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_ME_weight_cuts = new TH1F("L0bar_L0bar_cosThetaProdPlane_US_ME_weight_cuts", "L0bar_L0bar_cosThetaProdPlane_US_ME_weight_cuts", 10, -1, 1);
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
 
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[nEtaBins][nEtaBins];
  
  //delta eta vs. delta phi histograms for re-weighing of ME
  TH2F *L0_L0bar_delta_eta_vs_delta_phi_US_ME_cuts_hist = new TH2F("L0_L0bar_delta_eta_vs_delta_phi_US_ME_cuts_hist", "L0_L0bar_delta_eta_vs_delta_phi_US_ME_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0_L0_delta_eta_vs_delta_phi_US_ME_cuts_hist = new TH2F("L0_L0_delta_eta_vs_delta_phi_US_ME_cuts_hist", "L0_L0_delta_eta_vs_delta_phi_US_ME_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0bar_L0bar_delta_eta_vs_delta_phi_US_ME_cuts_hist = new TH2F("L0bar_L0bar_delta_eta_vs_delta_phi_US_ME_cuts_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_US_ME_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  
  //----------------------------------------------------------------------------------------- 
  
  //US matched to LS (background)
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts = new TH1F("L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts", "L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts", 10, -1, 1);
  
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_US_LS_ME_weight_cuts = new TH1F("L0_L0_cosThetaProdPlane_US_LS_ME_weight_cuts", "L0_L0_cosThetaProdPlane_US_LS_ME_weight_cuts", 10, -1, 1);
  
  TH1F *L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts = new TH1F("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts", "L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts", 10, -1, 1);
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
 
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[nEtaBins][nEtaBins];
  
  //delta eta vs. delta phi histograms for re-weighing of ME
  TH2F *L0_L0bar_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist = new TH2F("L0_L0bar_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist", "L0_L0bar_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0_L0_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist = new TH2F("L0_L0_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist", "L0_L0_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
  TH2F *L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist = new TH2F("L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());  

  //_________________________________________________________________________________________________________
    
  
  //L-L correlation histograms in bins
  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      //true L and Lbar
      L0_L0bar_cosThetaProdPlane_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
/*      
      //mixed event
      L0_L0bar_cosThetaProdPlane_ME_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_ME_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_ME_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
*/      
      //mixed event after re-weighing
      L0_L0bar_cosThetaProdPlane_ME_weight_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_ME_weight_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_ME_weight_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_ME_weight_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_ME_weight_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_ME_weight_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
      //---------------------------------------------------------
      
      //pi p pairs
      //US
      L0_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_US_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_US_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_US_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
      //US matched to LS
      L0_L0bar_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_US_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_LS_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_US_LS_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_US_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_LS_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_US_LS_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_US_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_LS_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_US_LS_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
/*      
      //mixed event
      //US
      L0_L0bar_cosThetaProdPlane_US_ME_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_US_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_US_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_US_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_US_ME_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_US_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_US_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_US_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_US_ME_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_US_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_US_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_US_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
      //US matched to LS
      L0_L0bar_cosThetaProdPlane_US_LS_ME_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_LS_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_US_LS_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_LS_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_US_LS_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_US_LS_ME_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_LS_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_US_LS_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_LS_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_US_LS_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_US_LS_ME_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
*/     
      
      //mixed event after re-weighing
      //US
      L0_L0bar_cosThetaProdPlane_US_ME_weight_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_US_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_US_ME_weight_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_US_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_US_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_US_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
      //US matched to LS
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_LS_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_US_LS_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
      
      //invariant mass
      
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2), 80, 1, 1.2, 180, 1, 1.2); //for US-US Lambda pairs
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), 80, 1, 1.2, 180, 1, 1.2); //for US-LS Lambda pairs

      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2), 80, 1, 1.2, 180, 1, 1.2); //for US-US Lambda pairs
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), 80, 1, 1.2, 180, 1, 1.2); //for US-LS Lambda pairs

      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2] = new TH2F(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2), 80, 1, 1.2, 180, 1, 1.2); //for US-US Lambda pairs
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2] = new TH2F(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), 80, 1, 1.2, 180, 1, 1.2); //for US-LS Lambda pairs
      
      
      L0_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0bar_inv_mass_US_cuts_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_cuts_pT1_%i_pT2_%i", pTbin1, pTbin2), 80, 1, 1.2, 180, 1, 1.2); //for US-US Lambda pairs
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts_pT1_%i_pT2_%i", pTbin1, pTbin2), 80, 1, 1.2, 180, 1, 1.2); //for US-LS Lambda pairs

      L0_inv_mass_vs_L0_inv_mass_US_cuts[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0_inv_mass_US_cuts_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_cuts_pT1_%i_pT2_%i", pTbin1, pTbin2), 80, 1, 1.2, 180, 1, 1.2); //for US-US Lambda pairs
      L0_inv_mass_vs_L0_inv_mass_US_LS_cuts[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0_inv_mass_US_LS_cuts_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_LS_cuts_pT1_%i_pT2_%i", pTbin1, pTbin2), 80, 1, 1.2, 180, 1, 1.2); //for US-LS Lambda pairs

      L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2] = new TH2F(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts_pT1_%i_pT2_%i", pTbin1, pTbin2), 80, 1, 1.2, 180, 1, 1.2); //for US-US Lambda pairs
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2] = new TH2F(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts_pT1_%i_pT2_%i", pTbin1, pTbin2), 80, 1, 1.2, 180, 1, 1.2); //for US-LS Lambda pairs
      
      
    }
  }

  for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
  {
    for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
    {
      //true L and Lbar
      L0_L0bar_cosThetaProdPlane_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
/*      
      //mixed event
      L0_L0bar_cosThetaProdPlane_ME_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_ME_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_ME_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
*/      
      //mixed event after re-weighing
      L0_L0bar_cosThetaProdPlane_ME_weight_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_ME_weight_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_ME_weight_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_ME_weight_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_ME_weight_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_ME_weight_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      
      //---------------------------------------------------------
      
      //pi p pairs
      
      //US
      L0_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_US_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_US_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_US_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      
      //US matched to LS
      L0_L0bar_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_US_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_LS_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_US_LS_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_US_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_LS_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_US_LS_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_US_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_LS_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_US_LS_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      
/*      
      //mixed event
      //US
      L0_L0bar_cosThetaProdPlane_US_ME_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_US_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_US_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_US_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_US_ME_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_US_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_US_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_US_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_US_ME_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_US_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_US_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_US_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      
      //US matched to LS
      L0_L0bar_cosThetaProdPlane_US_LS_ME_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_LS_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_US_LS_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_LS_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_US_LS_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_US_LS_ME_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_LS_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_US_LS_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_LS_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_US_LS_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_US_LS_ME_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
*/      
      //mixed event after re-weighing
      //US
      L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_US_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_US_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_US_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_US_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      
      //US matched to LS
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_LS_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_US_LS_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      
    }
  }
  //_____________________________________________________________________________________________________________________________________________________________________
  
  
  //mixed event vectors, true L and Lbar
  vector<TVector3> L_vector_ME;
  vector<TVector3> L_cuts_vector_ME;
  
  vector<float> L_decayL_ME_vector;
  vector<float> L_decayL_ME_cuts_vector;
 
  vector<TVector3> p_star_vector_ME;
  vector<TVector3> p_star_cuts_vector_ME;
  
  vector<int> L_pT_bin_vector_ME;
  vector<int> L_pT_bin_cuts_vector_ME;
  
  vector<int> L_eta_bin_vector_ME;
  vector<int> L_eta_bin_cuts_vector_ME;
  
  vector<float> L_p_pT_cuts_vector_ME;
  vector<float> L_pi_pT_cuts_vector_ME;
  
  vector<float> L_p_eta_cuts_vector_ME;
  vector<float> L_pi_eta_cuts_vector_ME;
  
  vector<float> L_p_phi_cuts_vector_ME;
  vector<float> L_pi_phi_cuts_vector_ME;

  
  vector<TVector3> Lbar_vector_ME;
  vector<TVector3> Lbar_cuts_vector_ME;
  
  vector<float> Lbar_decayL_ME_vector;
  vector<float> Lbar_decayL_ME_cuts_vector;
  
  vector<TVector3> pBar_star_vector_ME;
  vector<TVector3> pBar_star_cuts_vector_ME;  
 
  vector<int> Lbar_pT_bin_vector_ME;
  vector<int> Lbar_pT_bin_cuts_vector_ME;
  
  vector<int> Lbar_eta_bin_vector_ME;
  vector<int> Lbar_eta_bin_cuts_vector_ME;
  
  vector<float> Lbar_p_pT_cuts_vector_ME;
  vector<float> Lbar_pi_pT_cuts_vector_ME;
  
  vector<float> Lbar_p_eta_cuts_vector_ME;
  vector<float> Lbar_pi_eta_cuts_vector_ME;
  
  vector<float> Lbar_p_phi_cuts_vector_ME;
  vector<float> Lbar_pi_phi_cuts_vector_ME;
  
  //---------------------------------------
  
  //mixed event of pi p
  //US
  vector<TLorentzVector> L_vector_US_ME;  
  vector<TVector3> p_star_vector_US_ME;
  
  vector<int> L_pT_bin_vector_US_ME;  
  vector<int> L_eta_bin_vector_US_ME;

  
  vector<TLorentzVector> Lbar_vector_US_ME;
  vector<TVector3> pBar_star_vector_US_ME;
 
  vector<int> Lbar_pT_bin_vector_US_ME;  
  vector<int> Lbar_eta_bin_vector_US_ME;
  
  //after cuts
  vector<TLorentzVector> L_vector_US_ME_cuts; 
  vector<TVector3> p_star_vector_US_ME_cuts;
  
  vector<int> L_pT_bin_vector_US_ME_cuts;  
  vector<int> L_eta_bin_vector_US_ME_cuts;

  
  vector<TLorentzVector> Lbar_vector_US_ME_cuts;
  vector<TVector3> pBar_star_vector_US_ME_cuts;
 
  vector<int> Lbar_pT_bin_vector_US_ME_cuts;  
  vector<int> Lbar_eta_bin_vector_US_ME_cuts;

  
  
  //US matched to LS
  vector<TLorentzVector> L_vector_LS_ME;
  vector<TVector3> p_star_vector_LS_ME;
  
  vector<int> L_pT_bin_vector_LS_ME;  
  vector<int> L_eta_bin_vector_LS_ME;

  
  vector<TLorentzVector> Lbar_vector_LS_ME;
  vector<TVector3> pBar_star_vector_LS_ME;
 
  vector<int> Lbar_pT_bin_vector_LS_ME;  
  vector<int> Lbar_eta_bin_vector_LS_ME;
  
  //after cuts
  vector<TLorentzVector> L_vector_LS_ME_cuts;
  vector<TVector3> p_star_vector_LS_ME_cuts;
  
  vector<int> L_pT_bin_vector_LS_ME_cuts;  
  vector<int> L_eta_bin_vector_LS_ME_cuts;

  
  vector<TLorentzVector> Lbar_vector_LS_ME_cuts;
  vector<TVector3> pBar_star_vector_LS_ME_cuts;
 
  vector<int> Lbar_pT_bin_vector_LS_ME_cuts;  
  vector<int> Lbar_eta_bin_vector_LS_ME_cuts;
  
  //_____________________________________________________________________________________________________________________________________________________________________
  

  
  Long64_t nEntries_MC = L_MC_chain->GetEntries();
  
  cout<<"nEntries MC: "<<nEntries_MC<<endl;

  //loop over MC Lambda chain
  //one entry is one PYTHIA event with arrays of Lambdas and Lambda-bars
  for (Long64_t iEntry_MC = 0; iEntry_MC < nEntries_MC; iEntry_MC++)
  { 
    if(iEntry_MC % 10000000 == 0)
    {
      cout<<"Working on MC event:"<<iEntry_MC<<endl;
      
      //cout<<Lambda_MC.nL_MC<<endl;
      //cout<<Lambda_MC.L_charge_MC[0]<<endl;
    }
  
    L_MC_chain->GetEntry(iEntry_MC);
   
    //------------------------------------------------------------------
   
    //vectors of (anti-)proton 4-momenta boosted into mother rest frame
    vector<TVector3> L_vector;
    vector<TVector3> L_cuts_vector;
    
    vector<float> L_decayL_MC_vector;
    vector<float> L_decayL_MC_cuts_vector;
    
    vector<int> L_pT_bin_vector;
    vector<int> L_pT_bin_cuts_vector;
    
    vector<int> L_eta_bin_vector;
    vector<int> L_eta_bin_cuts_vector;
    
    
    vector<float> L_p_pT_cuts_vector;
    vector<float> L_pi_pT_cuts_vector;

    vector<float> L_p_eta_cuts_vector;
    vector<float> L_pi_eta_cuts_vector;
    
    vector<float> L_p_phi_cuts_vector;
    vector<float> L_pi_phi_cuts_vector;
    
    
    vector<TVector3> p_star_vector;
    vector<TVector3> p_star_cuts_vector;
    
    
    vector<TVector3> Lbar_vector;
    vector<TVector3> Lbar_cuts_vector;
    
    vector<float> Lbar_decayL_MC_vector;
    vector<float> Lbar_decayL_MC_cuts_vector;
    
    vector<int> Lbar_pT_bin_vector;
    vector<int> Lbar_pT_bin_cuts_vector;
    
    vector<int> Lbar_eta_bin_vector;
    vector<int> Lbar_eta_bin_cuts_vector;
    

    vector<float> Lbar_p_pT_cuts_vector;
    vector<float> Lbar_pi_pT_cuts_vector;

    vector<float> Lbar_p_eta_cuts_vector;
    vector<float> Lbar_pi_eta_cuts_vector;
    
    vector<float> Lbar_p_phi_cuts_vector;
    vector<float> Lbar_pi_phi_cuts_vector;

    vector<TVector3> pBar_star_vector;
    vector<TVector3> pBar_star_cuts_vector;    
    
    //cout<<Lambda_MC.nL_MC<<endl;
    //cout<<Lambda_MC.L_charge_MC[0]<<endl;

    //loop over MC Lambda chain
    //one entry is one PYTHIA event with arrays of Lambdas and Lambda-bars
    //sort L and Lbar and L (Lbar) before and after cut
    //L pairs created after this for loop
    for (int iLambda_MC = 0; iLambda_MC < Lambda_MC.nL_MC; iLambda_MC++)
    {  
      //Lambda momentum
      TVector3 L_mom(Lambda_MC.L_px_MC[iLambda_MC], Lambda_MC.L_py_MC[iLambda_MC], Lambda_MC.L_pz_MC[iLambda_MC]);    
      
      //if(fabs(L_mom.Eta()) > 0.3 ) continue; //for testing of acceptance effect

      //p momentum boosted to mother rest frame
      TVector3 p_mom_star(Lambda_MC.p_pxStar_MC[iLambda_MC], Lambda_MC.p_pyStar_MC[iLambda_MC], Lambda_MC.p_pzStar_MC[iLambda_MC]);
           
        
      //find bins               
      
      int pT_bin_corr = findBinPt( L_mom, pT_bins_corr, nPtBins_corr );
      
      if( pT_bin_corr == -1 ) continue;

      
      int eta_bin = findBinEta( L_mom, eta_bins, nEtaBins );

      if( eta_bin == -1 ) continue;


      //Lambda
      if( Lambda_MC.L_charge_MC[iLambda_MC] > 0 )
      {
        L_vector.push_back(L_mom);
        
        L_decayL_MC_vector.push_back(Lambda_MC.L_decayL_MC[iLambda_MC]);
        
        L_pT_bin_vector.push_back(pT_bin_corr);
        L_eta_bin_vector.push_back(eta_bin);
      
        p_star_vector.push_back(p_mom_star);
      }

      //Lambda-bar
      if( Lambda_MC.L_charge_MC[iLambda_MC] < 0 )
      {                  
        Lbar_vector.push_back(L_mom);
        
        Lbar_decayL_MC_vector.push_back(Lambda_MC.L_decayL_MC[iLambda_MC]);
        
        Lbar_pT_bin_vector.push_back(pT_bin_corr);
        Lbar_eta_bin_vector.push_back(eta_bin);

      
        pBar_star_vector.push_back(p_mom_star);      
      }
      
      
      //Lambda cuts
      //decay length cuts in mm (default PYTHIA units)
      if(Lambda_MC.L_decayL_MC[iLambda_MC] < 20 || Lambda_MC.L_decayL_MC[iLambda_MC] > 250) continue;

      
      //daughter cuts
      if( fabs(Lambda_MC.p_eta_MC[iLambda_MC]) >= 1. || fabs(Lambda_MC.pi_eta_MC[iLambda_MC]) >= 1. ) continue;       
      if( Lambda_MC.p_pT_MC[iLambda_MC] < 0.15 || Lambda_MC.p_pT_MC[iLambda_MC] > 20. ) continue;
      if( Lambda_MC.pi_pT_MC[iLambda_MC] < 0.15 || Lambda_MC.pi_pT_MC[iLambda_MC] > 20. ) continue;  
      

      //Lambda
      if( Lambda_MC.L_charge_MC[iLambda_MC] > 0 )
      {
        L_cuts_vector.push_back(L_mom);
        
        L_decayL_MC_cuts_vector.push_back(Lambda_MC.L_decayL_MC[iLambda_MC]);
        
        L_pT_bin_cuts_vector.push_back(pT_bin_corr);
        L_eta_bin_cuts_vector.push_back(eta_bin);
      
        p_star_cuts_vector.push_back(p_mom_star);
        
        L_p_pT_cuts_vector.push_back(Lambda_MC.p_pT_MC[iLambda_MC]); 
        L_pi_pT_cuts_vector.push_back(Lambda_MC.pi_pT_MC[iLambda_MC]);
        
        L_p_eta_cuts_vector.push_back(Lambda_MC.p_eta_MC[iLambda_MC]); 
        L_pi_eta_cuts_vector.push_back(Lambda_MC.pi_eta_MC[iLambda_MC]);
        
        L_p_phi_cuts_vector.push_back(Lambda_MC.p_phi_MC[iLambda_MC]); 
        L_pi_phi_cuts_vector.push_back(Lambda_MC.pi_phi_MC[iLambda_MC]);
     
      }

      //Lambda-bar
      if( Lambda_MC.L_charge_MC[iLambda_MC] < 0 )
      {
        Lbar_cuts_vector.push_back(L_mom);
        
        Lbar_decayL_MC_cuts_vector.push_back(Lambda_MC.L_decayL_MC[iLambda_MC]);
        
        Lbar_pT_bin_cuts_vector.push_back(pT_bin_corr);
        Lbar_eta_bin_cuts_vector.push_back(eta_bin);
      
        pBar_star_cuts_vector.push_back(p_mom_star);
        
        Lbar_p_pT_cuts_vector.push_back(Lambda_MC.p_pT_MC[iLambda_MC]); 
        Lbar_pi_pT_cuts_vector.push_back(Lambda_MC.pi_pT_MC[iLambda_MC]);
        
        Lbar_p_eta_cuts_vector.push_back(Lambda_MC.p_eta_MC[iLambda_MC]);
        Lbar_pi_eta_cuts_vector.push_back(Lambda_MC.pi_eta_MC[iLambda_MC]);
        
        Lbar_p_phi_cuts_vector.push_back(Lambda_MC.p_phi_MC[iLambda_MC]); 
        Lbar_pi_phi_cuts_vector.push_back(Lambda_MC.pi_phi_MC[iLambda_MC]);
      
      }         
        
    } //end loop over MC Lambdas in event
    //_________________________________________________________________________________________________________
    
     
    //true L and Lbar
    
    //fill L-L correlation histograms before cuts
    
    //L0-L0bar
    if(p_star_vector.size() > 0 && pBar_star_vector.size() > 0)
    {
      for(unsigned int iLambda = 0; iLambda < p_star_vector.size(); iLambda++)
      {
        for(unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector.size(); iLambdaBar++)
        {
          double theta_star = p_star_vector.at(iLambda).Angle(pBar_star_vector.at(iLambdaBar));
         
          L0_L0bar_cosThetaProdPlane->Fill(TMath::Cos(theta_star));          
          L0_L0bar_cosThetaProdPlane_pT_hist[L_pT_bin_vector.at(iLambda)][Lbar_pT_bin_vector.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
          L0_L0bar_cosThetaProdPlane_eta_hist[L_eta_bin_vector.at(iLambda)][Lbar_eta_bin_vector.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
          
          L0_L0bar_eta1_vs_eta2_hist->Fill(L_vector.at(iLambda).Eta() , Lbar_vector.at(iLambdaBar).Eta());
          L0_L0bar_phi1_vs_phi2_hist->Fill(L_vector.at(iLambda).Phi() , Lbar_vector.at(iLambdaBar).Phi());
          L0_L0bar_pT1_vs_pT2_hist->Fill(L_vector.at(iLambda).Pt() , Lbar_vector.at(iLambdaBar).Pt());
          
          
          //delta eta vs. delta phi for ME re-weighing
          /*
          float delta_eta = fabs( L_vector.at(iLambda).Eta() - Lbar_vector.at(iLambdaBar).Eta() );
          float delta_phi = fabs( L_vector.at(iLambda).Phi() - Lbar_vector.at(iLambdaBar).Phi() );
          float delta_pT = fabs( L_vector.at(iLambda).Pt() - Lbar_vector.at(iLambdaBar).Pt() );
          
          L0_L0bar_delta_eta_vs_delta_phi_hist->Fill(delta_eta, delta_phi);
          //L0_L0bar_delta_eta_vs_delta_phi_vs_delta_pT_hist->Fill(delta_eta, delta_phi, delta_pT);
          L0_L0bar_delta_pT->Fill(delta_pT);
          */
        }
      }      
    }
    
    //L0-L0
    if(p_star_vector.size() > 1)
    {
      for(unsigned int iLambda1 = 0; iLambda1 < p_star_vector.size(); iLambda1++)
      {
        for(unsigned int iLambda2 = iLambda1+1; iLambda2 < p_star_vector.size(); iLambda2++)
        {
          double theta_star = p_star_vector.at(iLambda1).Angle(p_star_vector.at(iLambda2));
          
          L0_L0_cosThetaProdPlane->Fill(TMath::Cos(theta_star));
          L0_L0_cosThetaProdPlane_pT_hist[L_pT_bin_vector.at(iLambda1)][L_pT_bin_vector.at(iLambda2)]->Fill(TMath::Cos(theta_star));
          L0_L0_cosThetaProdPlane_eta_hist[L_eta_bin_vector.at(iLambda1)][L_eta_bin_vector.at(iLambda2)]->Fill(TMath::Cos(theta_star));
          
          
          //delta eta vs. delta phi for ME re-weighing
          
          float delta_eta = fabs( L_vector.at(iLambda1).Eta() - L_vector.at(iLambda2).Eta() );
          float delta_phi = fabs( L_vector.at(iLambda1).Phi() - L_vector.at(iLambda2).Phi() );
          
          L0_L0_delta_eta_vs_delta_phi_hist->Fill(delta_eta, delta_phi);
                      
        }
      }    
    }
    
    //L0bar-L0bar
    if(pBar_star_vector.size() > 1)
    {
      for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pBar_star_vector.size(); iLambdaBar1++)
      {
        for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < pBar_star_vector.size(); iLambdaBar2++)
        {
          double theta_star = pBar_star_vector.at(iLambdaBar1).Angle(pBar_star_vector.at(iLambdaBar2));
          
          L0bar_L0bar_cosThetaProdPlane->Fill(TMath::Cos(theta_star));
          L0bar_L0bar_cosThetaProdPlane_pT_hist[Lbar_pT_bin_vector.at(iLambdaBar1)][Lbar_pT_bin_vector.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
          L0bar_L0bar_cosThetaProdPlane_eta_hist[Lbar_eta_bin_vector.at(iLambdaBar1)][Lbar_eta_bin_vector.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
          
          
          //delta eta vs. delta phi for ME re-weighing
          
          float delta_eta = fabs( Lbar_vector.at(iLambdaBar1).Eta() - Lbar_vector.at(iLambdaBar2).Eta() );
          float delta_phi = fabs( Lbar_vector.at(iLambdaBar1).Phi() - Lbar_vector.at(iLambdaBar2).Phi() );
          
          L0bar_L0bar_delta_eta_vs_delta_phi_hist->Fill(delta_eta, delta_phi);    
        }     
      }     
    }
    //_____________________________________________________________________________________________________________________________________________________________________

    //mixed event before cuts
    if( p_star_vector.size() >= 1 && pBar_star_vector.size() == 0 && p_star_vector_ME.size() < 1e3)
    {
      L_vector_ME.push_back(L_vector.at(0));
      
      L_decayL_ME_vector.push_back(L_decayL_MC_vector.at(0));
    
      p_star_vector_ME.push_back(p_star_vector.at(0));

      L_pT_bin_vector_ME.push_back(L_pT_bin_vector.at(0));
      L_eta_bin_vector_ME.push_back(L_eta_bin_vector.at(0));
    }

    if( p_star_vector.size() == 0 && pBar_star_vector.size() >= 1 && pBar_star_vector_ME.size() < 1e3)
    {    
      Lbar_vector_ME.push_back(Lbar_vector.at(0));
      
      Lbar_decayL_ME_vector.push_back(Lbar_decayL_MC_vector.at(0));
      
      pBar_star_vector_ME.push_back(pBar_star_vector.at(0));

      Lbar_pT_bin_vector_ME.push_back(Lbar_pT_bin_vector.at(0));
      Lbar_eta_bin_vector_ME.push_back(Lbar_eta_bin_vector.at(0));
    }
    //_________________________________________________________________________________________________________
    
     
    
    //fill L-L correlation histograms afrer cuts
    
    //L0-L0bar
    if(p_star_cuts_vector.size() > 0 && pBar_star_cuts_vector.size() > 0)
    {
      for(unsigned int iLambda = 0; iLambda < p_star_cuts_vector.size(); iLambda++)
      {
        for(unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_cuts_vector.size(); iLambdaBar++)
        {
          //if( fabs( L_cuts_vector.at(iLambda).Eta() - Lbar_cuts_vector.at(iLambdaBar).Eta()) < 0.3 ) continue;
        
          double theta_star = p_star_cuts_vector.at(iLambda).Angle(pBar_star_cuts_vector.at(iLambdaBar));
          
          L0_L0bar_cosThetaProdPlane_cuts->Fill(TMath::Cos(theta_star));
          L0_L0bar_cosThetaProdPlane_pT_cuts_hist[L_pT_bin_cuts_vector.at(iLambda)][Lbar_pT_bin_cuts_vector.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
          L0_L0bar_cosThetaProdPlane_eta_cuts_hist[L_eta_bin_cuts_vector.at(iLambda)][Lbar_eta_bin_cuts_vector.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));  
          
          
          L0_L0bar_eta1_vs_eta2_cuts_hist->Fill(L_cuts_vector.at(iLambda).Eta() , Lbar_cuts_vector.at(iLambdaBar).Eta());
          L0_L0bar_phi1_vs_phi2_cuts_hist->Fill(L_cuts_vector.at(iLambda).Phi() , Lbar_cuts_vector.at(iLambdaBar).Phi());
          L0_L0bar_pT1_vs_pT2_cuts_hist->Fill(L_cuts_vector.at(iLambda).Pt() , Lbar_cuts_vector.at(iLambdaBar).Pt());
          
          L0_L0bar_p1_pT1_vs_p2_pT2_cuts_hist->Fill(L_p_pT_cuts_vector.at(iLambda) , Lbar_p_pT_cuts_vector.at(iLambdaBar));
          L0_L0bar_pi1_pT1_vs_pi2_pT2_cuts_hist->Fill(L_pi_pT_cuts_vector.at(iLambda) , Lbar_pi_pT_cuts_vector.at(iLambdaBar));       
          
          L0_L0bar_p1_eta1_vs_p2_eta2_cuts_hist->Fill(L_p_eta_cuts_vector.at(iLambda) , Lbar_p_eta_cuts_vector.at(iLambdaBar));
          L0_L0bar_pi1_eta1_vs_pi2_eta2_cuts_hist->Fill(L_pi_eta_cuts_vector.at(iLambda) , Lbar_pi_eta_cuts_vector.at(iLambdaBar));
          
          L0_L0bar_p1_phi1_vs_p2_phi2_cuts_hist->Fill(L_p_phi_cuts_vector.at(iLambda) , Lbar_p_phi_cuts_vector.at(iLambdaBar));
          L0_L0bar_pi1_phi1_vs_pi2_phi2_cuts_hist->Fill(L_pi_phi_cuts_vector.at(iLambda) , Lbar_pi_phi_cuts_vector.at(iLambdaBar));
          
          //delta eta vs. delta phi for ME re-weighing
          /*
          float delta_eta = fabs( L_cuts_vector.at(iLambda).Eta() - Lbar_cuts_vector.at(iLambdaBar).Eta() );
          float delta_phi = fabs( L_cuts_vector.at(iLambda).Phi() - Lbar_cuts_vector.at(iLambdaBar).Phi() );
          float delta_pT = fabs( L_cuts_vector.at(iLambda).Pt() - Lbar_cuts_vector.at(iLambdaBar).Pt() );
          
          L0_L0bar_delta_eta_vs_delta_phi_cuts_hist->Fill(delta_eta, delta_phi);
          //L0_L0bar_delta_eta_vs_delta_phi_delta_pT_cuts_hist->Fill(delta_eta, delta_phi, delta_pT); 
          L0_L0bar_delta_pT_cuts->Fill(delta_pT);
          */
        }   
      
      }
    }
    

   
    //L0-L0
    if(p_star_cuts_vector.size() > 0)
    {
      for(unsigned int iLambda1 = 0; iLambda1 < p_star_cuts_vector.size(); iLambda1++)
      {
        for(unsigned int iLambda2 = iLambda1+1; iLambda2 < p_star_cuts_vector.size(); iLambda2++)
        {
          double theta_star = p_star_cuts_vector.at(iLambda1).Angle(p_star_cuts_vector.at(iLambda2));
          
          L0_L0_cosThetaProdPlane_cuts->Fill(TMath::Cos(theta_star));
          L0_L0_cosThetaProdPlane_pT_cuts_hist[L_pT_bin_cuts_vector.at(iLambda1)][L_pT_bin_cuts_vector.at(iLambda2)]->Fill(TMath::Cos(theta_star));
          L0_L0_cosThetaProdPlane_eta_cuts_hist[L_eta_bin_cuts_vector.at(iLambda1)][L_eta_bin_cuts_vector.at(iLambda2)]->Fill(TMath::Cos(theta_star));           
          
          
          //delta eta vs. delta phi for ME re-weighing          
          float delta_eta = fabs( L_cuts_vector.at(iLambda1).Eta() - L_cuts_vector.at(iLambda2).Eta() );
          float delta_phi = fabs( L_cuts_vector.at(iLambda1).Phi() - L_cuts_vector.at(iLambda2).Phi() );
          
          L0_L0_delta_eta_vs_delta_phi_cuts_hist->Fill(delta_eta, delta_phi);
        }   
      
      }
    
    }
    
    //L0bar-L0bar
    if(pBar_star_cuts_vector.size() > 0)
    {
      for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pBar_star_cuts_vector.size(); iLambdaBar1++)
      {
        for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < pBar_star_cuts_vector.size(); iLambdaBar2++)
        {
          double theta_star = pBar_star_cuts_vector.at(iLambdaBar1).Angle(pBar_star_cuts_vector.at(iLambdaBar2));
          
          L0bar_L0bar_cosThetaProdPlane_cuts->Fill(TMath::Cos(theta_star));
          L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[Lbar_pT_bin_cuts_vector.at(iLambdaBar1)][Lbar_pT_bin_cuts_vector.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
          L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[Lbar_eta_bin_cuts_vector.at(iLambdaBar1)][Lbar_eta_bin_cuts_vector.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
          
          //delta eta vs. delta phi for ME re-weighing          
          float delta_eta = fabs( Lbar_cuts_vector.at(iLambdaBar1).Eta() - Lbar_cuts_vector.at(iLambdaBar2).Eta() );
          float delta_phi = fabs( Lbar_cuts_vector.at(iLambdaBar1).Phi() - Lbar_cuts_vector.at(iLambdaBar2).Phi() );
          
          L0bar_L0bar_delta_eta_vs_delta_phi_cuts_hist->Fill(delta_eta, delta_phi);
        }   
      
      } 
    
    }
    
    //mixed event after cuts
    if( p_star_cuts_vector.size() >= 1 && pBar_star_cuts_vector.size() == 0 && p_star_cuts_vector_ME.size() < 3e3)
    {
      L_cuts_vector_ME.push_back(L_cuts_vector.at(0));
      
      L_decayL_ME_cuts_vector.push_back(L_decayL_MC_cuts_vector.at(0));
   
      p_star_cuts_vector_ME.push_back(p_star_cuts_vector.at(0));

      L_pT_bin_cuts_vector_ME.push_back(L_pT_bin_cuts_vector.at(0));
      L_eta_bin_cuts_vector_ME.push_back(L_eta_bin_cuts_vector.at(0));
      
      L_p_pT_cuts_vector_ME.push_back(L_p_pT_cuts_vector.at(0));
      L_pi_pT_cuts_vector_ME.push_back(L_pi_pT_cuts_vector.at(0));
      
      L_p_eta_cuts_vector_ME.push_back(L_p_eta_cuts_vector.at(0));
      L_pi_eta_cuts_vector_ME.push_back(L_pi_eta_cuts_vector.at(0));
      
      L_p_phi_cuts_vector_ME.push_back(L_p_phi_cuts_vector.at(0));
      L_pi_phi_cuts_vector_ME.push_back(L_pi_phi_cuts_vector.at(0));

    }

    if( p_star_cuts_vector.size() == 0 && pBar_star_cuts_vector.size() >= 1 && pBar_star_cuts_vector_ME.size() < 3e3)
    {    
      Lbar_cuts_vector_ME.push_back(Lbar_cuts_vector.at(0));
      
      Lbar_decayL_ME_cuts_vector.push_back(Lbar_decayL_MC_cuts_vector.at(0));
      
      pBar_star_cuts_vector_ME.push_back(pBar_star_cuts_vector.at(0));

      Lbar_pT_bin_cuts_vector_ME.push_back(Lbar_pT_bin_cuts_vector.at(0));
      Lbar_eta_bin_cuts_vector_ME.push_back(Lbar_eta_bin_cuts_vector.at(0));
      
      Lbar_p_pT_cuts_vector_ME.push_back(Lbar_p_pT_cuts_vector.at(0));
      Lbar_pi_pT_cuts_vector_ME.push_back(Lbar_pi_pT_cuts_vector.at(0));
      
      Lbar_p_eta_cuts_vector_ME.push_back(Lbar_p_eta_cuts_vector.at(0));
      Lbar_pi_eta_cuts_vector_ME.push_back(Lbar_pi_eta_cuts_vector.at(0));
      
      Lbar_p_phi_cuts_vector_ME.push_back(Lbar_p_phi_cuts_vector.at(0));
      Lbar_pi_phi_cuts_vector_ME.push_back(Lbar_pi_phi_cuts_vector.at(0));

    }
    
    //---------------------------------------------------------------------------------------------------------------------
    
    L_vector.clear();
    L_cuts_vector.clear();
    
    L_decayL_MC_vector.clear();
    L_decayL_MC_cuts_vector.clear();
    
    L_pT_bin_vector.clear();
    L_pT_bin_cuts_vector.clear();
    
    L_eta_bin_vector.clear();
    L_eta_bin_cuts_vector.clear();
       
    p_star_vector.clear();
    p_star_cuts_vector.clear();
    
    L_p_pT_cuts_vector.clear();
    L_pi_pT_cuts_vector.clear();
    
    L_p_eta_cuts_vector.clear();
    L_pi_eta_cuts_vector.clear();
    
    L_p_phi_cuts_vector.clear();
    L_pi_phi_cuts_vector.clear();
    
    
    Lbar_vector.clear();
    Lbar_cuts_vector.clear();
    
    Lbar_decayL_MC_vector.clear();
    Lbar_decayL_MC_cuts_vector.clear();
    
    Lbar_pT_bin_vector.clear();
    Lbar_pT_bin_cuts_vector.clear();
    
    Lbar_eta_bin_vector.clear();
    Lbar_eta_bin_cuts_vector.clear();
    
    pBar_star_vector.clear();
    pBar_star_cuts_vector.clear();   
    
    Lbar_p_pT_cuts_vector.clear();
    Lbar_pi_pT_cuts_vector.clear();
    
    Lbar_p_eta_cuts_vector.clear();
    Lbar_pi_eta_cuts_vector.clear();
    
    Lbar_p_phi_cuts_vector.clear();
    Lbar_pi_phi_cuts_vector.clear();
    
    
  }//end for loop over MC Lambdas

  
  //______________________________________________________________________________________________________________________________________________________________________________________
  
  Long64_t nEntries_RC = L_from_pairs_chain->GetEntries();

  cout<<"nEntries RC: "<<nEntries_RC<<endl;

  //loop over MC Lambda chain
  //one entry is one PYTHIA event with arrays of Lambdas and Lambda-bars
  for (Long64_t iEntry_RC = 0; iEntry_RC < nEntries_RC; iEntry_RC++)
  { 
    if(iEntry_RC % 10000000 == 0)
    {
      cout<<"Working on RC event:"<<iEntry_RC<<endl;
    }
  
    L_from_pairs_chain->GetEntry(iEntry_RC);
    
    // analyze pi p pairs
    
    //L and Lbar "candidates" from pi p pairs
    vector<TLorentzVector> L_vector_US;
    vector<int> L_pT_bin_vector_US;
    vector<int> L_eta_bin_vector_US;
    vector<int> L_cuts_flag_US; //flag if pair passed analysis cuts
    
    vector<TVector3> p_star_vector_US; //to store p fourmomentum in L (pi p pair) rest frame
    
    vector<float> p_pT_vector_US; //for auto-correlation check
    vector<float> pi_pT_vector_US;
   
    
    vector<TLorentzVector> Lbar_vector_US;
    vector<int> Lbar_pT_bin_vector_US;
    vector<int> Lbar_eta_bin_vector_US;
    vector<int> Lbar_cuts_flag_US; //flag if pair passed analysis cuts
    
    vector<TVector3> pBar_star_vector_US; //to store p fourmomentum in Lbar (pi pBar pair) rest frame    
    
    vector<float> pBar_pT_vector_US; //for auto-correlation check
    vector<float> piBar_pT_vector_US;
   
    
    vector<TLorentzVector> L_vector_LS;
    vector<int> L_pT_bin_vector_LS;
    vector<int> L_eta_bin_vector_LS;
    vector<int> L_cuts_flag_LS; //flag if pair passed analysis cuts

    vector<TVector3> p_star_vector_LS;
    
    vector<float> p_pT_vector_LS; //for auto-correlation check
    vector<float> pi_pT_vector_LS;
  
    
    vector<TLorentzVector> Lbar_vector_LS;
    vector<int> Lbar_pT_bin_vector_LS;
    vector<int> Lbar_eta_bin_vector_LS;
    vector<int> Lbar_cuts_flag_LS; //flag if pair passed analysis cuts
    
    vector<TVector3> pBar_star_vector_LS;
    
    vector<float> pBar_pT_vector_LS; //for auto-correlation check
    vector<float> piBar_pT_vector_LS;
    
    //cout<<Lambda_from_pair.nL<<endl;
  
    //analyze L from p-pi pairs
    //separate L and Lbar, before and after cuts, and US vs. LS p-pi pairs

    for(unsigned int iLambda_RC = 0; iLambda_RC < Lambda_from_pair.nL; iLambda_RC++)
    {
      TVector3 L_mom_RC(Lambda_from_pair.L_px[iLambda_RC], Lambda_from_pair.L_py[iLambda_RC], Lambda_from_pair.L_pz[iLambda_RC]);
      TLorentzVector L_fourmom_RC(L_mom_RC, Lambda_from_pair.L_Minv[iLambda_RC]);
     
      
      int pT_bin_L_RC = findBinPt( L_mom_RC, pT_bins_corr, nPtBins_corr );
      int eta_bin_L_RC = findBinEta( L_mom_RC, eta_bins, nEtaBins );
      
      if( pT_bin_L_RC < 0 || eta_bin_L_RC < 0 ) continue;        
      
      
      TVector3 p_star_RC(Lambda_from_pair.p_pxStar[iLambda_RC], Lambda_from_pair.p_pyStar[iLambda_RC], Lambda_from_pair.p_pzStar[iLambda_RC]);
     
      
      //cuts        
      int cuts_flag = 1;      
      
      if(Lambda_from_pair.L_decayL[iLambda_RC] < 20 || Lambda_from_pair.L_decayL[iLambda_RC] > 250) cuts_flag = 0;
      
      //cos of pointing angle
      if(cos(Lambda_from_pair.L_theta[iLambda_RC]) < 0.996) cuts_flag = 0;
              
      //if(L_fourmom.Pt() < 0.5) continue; //added pT cut - pT integrated in data starts at 0.5 GeV/c
      //if(L_fourmom.Rapidity() >= 1) continue;
      
      //daughter cuts
      if( fabs(Lambda_from_pair.p_eta[iLambda_RC]) >= 1. || fabs(Lambda_from_pair.pi_eta[iLambda_RC]) >= 1. )  cuts_flag = 0;       
      if( Lambda_from_pair.p_pT[iLambda_RC] < 0.15 || Lambda_from_pair.p_pT[iLambda_RC] > 20. )  cuts_flag = 0;
      if( Lambda_from_pair.pi_pT[iLambda_RC] < 0.15 || Lambda_from_pair.pi_pT[iLambda_RC] > 20. )  cuts_flag = 0;
      
      //L, US
      if( Lambda_from_pair.L_charge[iLambda_RC] > 0 && Lambda_from_pair.L_US_LS_flag[iLambda_RC] == 1 )
      {
        L_vector_US.push_back(L_fourmom_RC);
        L_pT_bin_vector_US.push_back(pT_bin_L_RC);
        L_eta_bin_vector_US.push_back(eta_bin_L_RC);
        p_star_vector_US.push_back(p_star_RC);
        
        L_cuts_flag_US.push_back(cuts_flag);
        
        p_pT_vector_US.push_back(Lambda_from_pair.p_pT[iLambda_RC]);
        piBar_pT_vector_US.push_back(Lambda_from_pair.pi_pT[iLambda_RC]);

      }
      
      //L, LS
      if( Lambda_from_pair.L_charge[iLambda_RC] > 0 && Lambda_from_pair.L_US_LS_flag[iLambda_RC] == 0 )
      {
        L_vector_LS.push_back(L_fourmom_RC);
        L_pT_bin_vector_LS.push_back(pT_bin_L_RC);
        L_eta_bin_vector_LS.push_back(eta_bin_L_RC);
        p_star_vector_LS.push_back(p_star_RC);
        
        L_cuts_flag_LS.push_back(cuts_flag);
        
        p_pT_vector_LS.push_back(Lambda_from_pair.p_pT[iLambda_RC]);
        pi_pT_vector_LS.push_back(Lambda_from_pair.pi_pT[iLambda_RC]);
      }
      
      //Lbar, US
      if( Lambda_from_pair.L_charge[iLambda_RC] < 0 && Lambda_from_pair.L_US_LS_flag[iLambda_RC] == 1 )
      {
        Lbar_vector_US.push_back(L_fourmom_RC);
        Lbar_pT_bin_vector_US.push_back(pT_bin_L_RC);
        Lbar_eta_bin_vector_US.push_back(eta_bin_L_RC);
        pBar_star_vector_US.push_back(p_star_RC);
        
        Lbar_cuts_flag_US.push_back(cuts_flag);
        
        pBar_pT_vector_US.push_back(Lambda_from_pair.p_pT[iLambda_RC]);
        pi_pT_vector_US.push_back(Lambda_from_pair.pi_pT[iLambda_RC]);

      }
      
      //Lbar, LS
      if( Lambda_from_pair.L_charge[iLambda_RC] < 0 && Lambda_from_pair.L_US_LS_flag[iLambda_RC] == 0 )
      {
        Lbar_vector_LS.push_back(L_fourmom_RC);
        Lbar_pT_bin_vector_LS.push_back(pT_bin_L_RC);
        Lbar_eta_bin_vector_LS.push_back(eta_bin_L_RC);
        pBar_star_vector_LS.push_back(p_star_RC);
        
        Lbar_cuts_flag_LS.push_back(cuts_flag);
        
        pBar_pT_vector_LS.push_back(Lambda_from_pair.p_pT[iLambda_RC]);
        piBar_pT_vector_LS.push_back(Lambda_from_pair.pi_pT[iLambda_RC]);
      }
      
            
    }//end for loop over L from p pi pairs
      
    
    //---------------------------------------------------------------------------------------------    
    
    //analyze L and Lbar candidates created from pi p pairs

    //US paired with US
    //L-Lbar
    if( L_vector_US.size() > 0 && Lbar_vector_US.size() > 0 )
    {
      for( unsigned int iLambda = 0; iLambda < L_vector_US.size(); iLambda++)
      {
        for( unsigned int iLambdaBar = 0; iLambdaBar < Lbar_vector_US.size(); iLambdaBar++)
        {
          L0_inv_mass_vs_L0bar_inv_mass_US[L_pT_bin_vector_US.at(iLambda)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill( L_vector_US.at(iLambda).M(), Lbar_vector_US.at(iLambdaBar).M() );
          
          float theta_star = p_star_vector_US.at(iLambda).Angle(pBar_star_vector_US.at(iLambdaBar));
          
          float delta_eta = fabs( L_vector_US.at(iLambda).Eta() - Lbar_vector_US.at(iLambdaBar).Eta() );
          float delta_phi = fabs( L_vector_US.at(iLambda).Phi() - Lbar_vector_US.at(iLambdaBar).Phi() );
          
          
          L0_L0bar_cosThetaProdPlane_US->Fill(TMath::Cos(theta_star));
          L0_L0bar_cosThetaProdPlane_US_pT_hist[L_pT_bin_vector_US.at(iLambda)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
          L0_L0bar_cosThetaProdPlane_US_eta_hist[L_eta_bin_vector_US.at(iLambda)][Lbar_eta_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
          
          //delta eta vs. delta phi for ME re-weighing   
          L0_L0bar_delta_eta_vs_delta_phi_US_hist->Fill(delta_eta, delta_phi);

          if( L_cuts_flag_US.at(iLambda) == 1 && Lbar_cuts_flag_US.at(iLambdaBar) == 1) //both L in pair passed cuts
          {
            L0_inv_mass_vs_L0bar_inv_mass_US_cuts[L_pT_bin_vector_US.at(iLambda)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill( L_vector_US.at(iLambda).M(), Lbar_vector_US.at(iLambdaBar).M() );
            
                       
            L0_L0bar_cosThetaProdPlane_US_cuts->Fill(TMath::Cos(theta_star));
            L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[L_pT_bin_vector_US.at(iLambda)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
            L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[L_eta_bin_vector_US.at(iLambda)][Lbar_eta_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
            
            //delta eta vs. delta phi for ME re-weighing              
            L0_L0bar_delta_eta_vs_delta_phi_US_cuts_hist->Fill(delta_eta, delta_phi);
            
            
          }
        
        }
      
      }
      
    }
    
    //L-L
    if( L_vector_US.size() > 1 )
    {
      for( unsigned int iLambda1 = 0; iLambda1 < L_vector_US.size(); iLambda1++ )
      {
        for( unsigned int iLambda2 = iLambda1+1; iLambda2 < L_vector_US.size(); iLambda2++ )
        {
          //check auto-correlation
          if( p_pT_vector_US.at(iLambda1) == p_pT_vector_US.at(iLambda2) ) continue;
          if( piBar_pT_vector_US.at(iLambda1) == piBar_pT_vector_US.at(iLambda2) ) continue;
          
          L0_inv_mass_vs_L0_inv_mass_US[L_pT_bin_vector_US.at(iLambda1)][L_pT_bin_vector_US.at(iLambda2)]->Fill( L_vector_US.at(iLambda1).M(), L_vector_US.at(iLambda2).M() );         
          
          float theta_star = p_star_vector_US.at(iLambda1).Angle(p_star_vector_US.at(iLambda2));
          
          float delta_eta = fabs( L_vector_US.at(iLambda1).Eta() - L_vector_US.at(iLambda2).Eta() );
          float delta_phi = fabs( L_vector_US.at(iLambda1).Phi() - L_vector_US.at(iLambda2).Phi() );
          
          L0_L0_cosThetaProdPlane_US->Fill(TMath::Cos(theta_star));
          L0_L0_cosThetaProdPlane_US_pT_hist[L_pT_bin_vector_US.at(iLambda1)][L_pT_bin_vector_US.at(iLambda2)]->Fill(TMath::Cos(theta_star));
          L0_L0_cosThetaProdPlane_US_eta_hist[L_eta_bin_vector_US.at(iLambda1)][L_eta_bin_vector_US.at(iLambda2)]->Fill(TMath::Cos(theta_star));
          
          //delta eta vs. delta phi for ME re-weighing   
          L0_L0_delta_eta_vs_delta_phi_US_hist->Fill(delta_eta, delta_phi);
          
          
          if( L_cuts_flag_US.at(iLambda1) == 1 && L_cuts_flag_US.at(iLambda2) == 1) //both L in pair passed cuts
          {
            L0_inv_mass_vs_L0_inv_mass_US_cuts[L_pT_bin_vector_US.at(iLambda1)][L_pT_bin_vector_US.at(iLambda2)]->Fill( L_vector_US.at(iLambda1).M(), L_vector_US.at(iLambda2).M() );
            
             
            L0_L0_cosThetaProdPlane_US_cuts->Fill(TMath::Cos(theta_star));
            L0_L0_cosThetaProdPlane_US_pT_cuts_hist[L_pT_bin_vector_US.at(iLambda1)][L_pT_bin_vector_US.at(iLambda2)]->Fill(TMath::Cos(theta_star));
            L0_L0_cosThetaProdPlane_US_eta_cuts_hist[L_eta_bin_vector_US.at(iLambda1)][L_eta_bin_vector_US.at(iLambda2)]->Fill(TMath::Cos(theta_star));
            
            //delta eta vs. delta phi for ME re-weighing   
            L0_L0_delta_eta_vs_delta_phi_US_cuts_hist->Fill(delta_eta, delta_phi);
            
            
          }
        
        }
      
      }
        
    }
    
    
    //Lbar-Lbar
    if( Lbar_vector_US.size() > 1 )
    {
      for( unsigned int iLambdaBar1 = 0; iLambdaBar1 < Lbar_vector_US.size(); iLambdaBar1++ )
      {
        for( unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < Lbar_vector_US.size(); iLambdaBar2++ )
        {
          //check auto-correlation
          if( pBar_pT_vector_US.at(iLambdaBar1) == pBar_pT_vector_US.at(iLambdaBar2) ) continue;
          if( pi_pT_vector_US.at(iLambdaBar1) == pi_pT_vector_US.at(iLambdaBar2) ) continue;
          
          L0bar_inv_mass_vs_L0bar_inv_mass_US[Lbar_pT_bin_vector_US.at(iLambdaBar1)][Lbar_pT_bin_vector_US.at(iLambdaBar2)]->Fill( Lbar_vector_US.at(iLambdaBar1).M(), Lbar_vector_US.at(iLambdaBar2).M() );
          
          float theta_star = pBar_star_vector_US.at(iLambdaBar1).Angle(pBar_star_vector_US.at(iLambdaBar2));
          
          float delta_eta = fabs( Lbar_vector_US.at(iLambdaBar1).Eta() - Lbar_vector_US.at(iLambdaBar2).Eta() );
          float delta_phi = fabs( Lbar_vector_US.at(iLambdaBar1).Phi() - Lbar_vector_US.at(iLambdaBar2).Phi() );
          
          
          L0bar_L0bar_cosThetaProdPlane_US->Fill(TMath::Cos(theta_star));
          L0bar_L0bar_cosThetaProdPlane_US_pT_hist[Lbar_pT_bin_vector_US.at(iLambdaBar1)][Lbar_pT_bin_vector_US.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
          L0bar_L0bar_cosThetaProdPlane_US_eta_hist[Lbar_eta_bin_vector_US.at(iLambdaBar1)][Lbar_eta_bin_vector_US.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
          
          //delta eta vs. delta phi for ME re-weighing   
          L0bar_L0bar_delta_eta_vs_delta_phi_US_hist->Fill(delta_eta, delta_phi);
          
          
          if( Lbar_cuts_flag_US.at(iLambdaBar1) == 1 && Lbar_cuts_flag_US.at(iLambdaBar2) == 1) //both L in pair passed cuts
          {
            L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts[Lbar_pT_bin_vector_US.at(iLambdaBar1)][Lbar_pT_bin_vector_US.at(iLambdaBar2)]->Fill( Lbar_vector_US.at(iLambdaBar1).M(), Lbar_vector_US.at(iLambdaBar2).M() );
            
            
            L0bar_L0bar_cosThetaProdPlane_US_cuts->Fill(TMath::Cos(theta_star));
            L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[Lbar_pT_bin_vector_US.at(iLambdaBar1)][Lbar_pT_bin_vector_US.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
            L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[Lbar_eta_bin_vector_US.at(iLambdaBar1)][Lbar_eta_bin_vector_US.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
            
            //delta eta vs. delta phi for ME re-weighing   
            L0bar_L0bar_delta_eta_vs_delta_phi_US_cuts_hist->Fill(delta_eta, delta_phi);
            
            
          }
        
        }
      
      }    
    
    }
    //--------------------------------------------------------------------

    //US paired with LS
    //L-Lbar
    //L - US, Lbar - LS
    if( L_vector_US.size() > 0 && Lbar_vector_LS.size() > 0 )
    {
      for( unsigned int iLambda = 0; iLambda < L_vector_US.size(); iLambda++)
      {
        for( unsigned int iLambdaBar = 0; iLambdaBar < Lbar_vector_LS.size(); iLambdaBar++)
        {
          //check auto-correlation (prevent daughter sharing in the L-Lbar pair)
          if( piBar_pT_vector_US.at(iLambda) == piBar_pT_vector_LS.at(iLambdaBar) ) continue;
          
          L0_inv_mass_vs_L0bar_inv_mass_US_LS[L_pT_bin_vector_US.at(iLambda)][Lbar_pT_bin_vector_LS.at(iLambdaBar)]->Fill( L_vector_US.at(iLambda).M(), Lbar_vector_LS.at(iLambdaBar).M() );          
        
          float theta_star = p_star_vector_US.at(iLambda).Angle(pBar_star_vector_LS.at(iLambdaBar));
          
          float delta_eta = fabs( L_vector_US.at(iLambda).Eta() - Lbar_vector_LS.at(iLambdaBar).Eta() );
          float delta_phi = fabs( L_vector_US.at(iLambda).Phi() - Lbar_vector_LS.at(iLambdaBar).Phi() );
          
          
          L0_L0bar_cosThetaProdPlane_US_LS->Fill(TMath::Cos(theta_star));
          L0_L0bar_cosThetaProdPlane_US_LS_pT_hist[L_pT_bin_vector_US.at(iLambda)][Lbar_pT_bin_vector_LS.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
          L0_L0bar_cosThetaProdPlane_US_LS_eta_hist[L_eta_bin_vector_US.at(iLambda)][Lbar_eta_bin_vector_LS.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
          
          //delta eta vs. delta phi for ME re-weighing   
          L0_L0bar_delta_eta_vs_delta_phi_US_LS_hist->Fill(delta_eta, delta_phi);
          
          
          if( L_cuts_flag_US.at(iLambda) == 1 && Lbar_cuts_flag_LS.at(iLambdaBar) == 1) //both L in pair passed cuts
          {
            L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[L_pT_bin_vector_US.at(iLambda)][Lbar_pT_bin_vector_LS.at(iLambdaBar)]->Fill( L_vector_US.at(iLambda).M(), Lbar_vector_LS.at(iLambdaBar).M() );
          
            
            L0_L0bar_cosThetaProdPlane_US_LS_cuts->Fill(TMath::Cos(theta_star));
            L0_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[L_pT_bin_vector_US.at(iLambda)][Lbar_pT_bin_vector_LS.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
            L0_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[L_eta_bin_vector_US.at(iLambda)][Lbar_eta_bin_vector_LS.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));    
            
            //delta eta vs. delta phi for ME re-weighing   
            L0_L0bar_delta_eta_vs_delta_phi_US_LS_cuts_hist->Fill(delta_eta, delta_phi);     
            
            
          }
        
        }
      
      }
    
    }
    
    //L - LS, Lbar - US
    if( L_vector_LS.size() > 0 && Lbar_vector_US.size() > 0 )
    {
      for( unsigned int iLambda = 0; iLambda < L_vector_LS.size(); iLambda++)
      {
        for( unsigned int iLambdaBar = 0; iLambdaBar < Lbar_vector_US.size(); iLambdaBar++)
        {
          //check auto-correlation (prevent daughter sharing in the L-Lbar pair)
          if( pi_pT_vector_LS.at(iLambda) == pi_pT_vector_US.at(iLambdaBar) ) continue;
          
          L0_inv_mass_vs_L0bar_inv_mass_US_LS[L_pT_bin_vector_LS.at(iLambda)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill( L_vector_LS.at(iLambda).M(), Lbar_vector_US.at(iLambdaBar).M() );
        
          float theta_star = p_star_vector_LS.at(iLambda).Angle(pBar_star_vector_US.at(iLambdaBar));
          
          float delta_eta = fabs( L_vector_LS.at(iLambda).Eta() - Lbar_vector_US.at(iLambdaBar).Eta() );
          float delta_phi = fabs( L_vector_LS.at(iLambda).Phi() - Lbar_vector_US.at(iLambdaBar).Phi() );
          
          
          L0_L0bar_cosThetaProdPlane_US_LS->Fill(TMath::Cos(theta_star));
          L0_L0bar_cosThetaProdPlane_US_LS_pT_hist[L_pT_bin_vector_LS.at(iLambda)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
          L0_L0bar_cosThetaProdPlane_US_LS_eta_hist[L_eta_bin_vector_LS.at(iLambda)][Lbar_eta_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
          
          //delta eta vs. delta phi for ME re-weighing   
          L0_L0bar_delta_eta_vs_delta_phi_US_LS_hist->Fill(delta_eta, delta_phi);
          
          
          if( L_cuts_flag_LS.at(iLambda) == 1 && Lbar_cuts_flag_US.at(iLambdaBar) == 1) //both L in pair passed cuts
          {
            L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[L_pT_bin_vector_LS.at(iLambda)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill( L_vector_LS.at(iLambda).M(), Lbar_vector_US.at(iLambdaBar).M() );
            
            
            L0_L0bar_cosThetaProdPlane_US_LS_cuts->Fill(TMath::Cos(theta_star));
            L0_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[L_pT_bin_vector_LS.at(iLambda)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
            L0_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[L_eta_bin_vector_LS.at(iLambda)][Lbar_eta_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));       
            
            //delta eta vs. delta phi for ME re-weighing   
            L0_L0bar_delta_eta_vs_delta_phi_US_LS_cuts_hist->Fill(delta_eta, delta_phi);   
            
            
          }
        
        }
      
      }
    
    }
    
    //L-L
    int nFills_L_L_US_LS = 0;
    
    if( L_vector_US.size() > 0 && L_vector_LS.size() > 0 )
    {
      for( unsigned int iLambda = 0; iLambda < L_vector_US.size(); iLambda++ )
      {
        for( unsigned int iLambdaBck = 0; iLambdaBck < L_vector_LS.size(); iLambdaBck++ )
        {
          //check auto-correlation (prevent daughter sharing in the L-L pair)
          if( p_pT_vector_US.at(iLambda) == p_pT_vector_LS.at(iLambdaBck) ) continue;
          
          //interate order of US and LS for L
          if( nFills_L_L_US_LS % 2 == 0 )
          {
            L0_inv_mass_vs_L0_inv_mass_US_LS[L_pT_bin_vector_US.at(iLambda)][L_pT_bin_vector_LS.at(iLambdaBck)]->Fill( L_vector_US.at(iLambda).M(), L_vector_LS.at(iLambdaBck).M() );
          
            float theta_star = p_star_vector_US.at(iLambda).Angle(p_star_vector_LS.at(iLambdaBck));
            
            float delta_eta = fabs( L_vector_US.at(iLambda).Eta() - L_vector_LS.at(iLambdaBck).Eta() );
            float delta_phi = fabs( L_vector_US.at(iLambda).Phi() - L_vector_LS.at(iLambdaBck).Phi() );
            
            
            L0_L0_cosThetaProdPlane_US_LS->Fill(TMath::Cos(theta_star));
            L0_L0_cosThetaProdPlane_US_LS_pT_hist[L_pT_bin_vector_US.at(iLambda)][L_pT_bin_vector_LS.at(iLambdaBck)]->Fill(TMath::Cos(theta_star));
            L0_L0_cosThetaProdPlane_US_LS_eta_hist[L_eta_bin_vector_US.at(iLambda)][L_eta_bin_vector_LS.at(iLambdaBck)]->Fill(TMath::Cos(theta_star));
            
            //delta eta vs. delta phi for ME re-weighing   
            L0_L0_delta_eta_vs_delta_phi_US_LS_hist->Fill(delta_eta, delta_phi);
            
            
            if( L_cuts_flag_US.at(iLambda) == 1 && L_cuts_flag_LS.at(iLambdaBck) == 1) //both L in pair passed cuts
            {
              L0_inv_mass_vs_L0_inv_mass_US_LS_cuts[L_pT_bin_vector_US.at(iLambda)][L_pT_bin_vector_LS.at(iLambdaBck)]->Fill( L_vector_US.at(iLambda).M(), L_vector_LS.at(iLambdaBck).M() );
              
              
              L0_L0_cosThetaProdPlane_US_LS_cuts->Fill(TMath::Cos(theta_star));
              L0_L0_cosThetaProdPlane_US_LS_pT_cuts_hist[L_pT_bin_vector_US.at(iLambda)][L_pT_bin_vector_LS.at(iLambdaBck)]->Fill(TMath::Cos(theta_star));
              L0_L0_cosThetaProdPlane_US_LS_eta_cuts_hist[L_eta_bin_vector_US.at(iLambda)][L_eta_bin_vector_LS.at(iLambdaBck)]->Fill(TMath::Cos(theta_star));  
              
              //delta eta vs. delta phi for ME re-weighing   
              L0_L0_delta_eta_vs_delta_phi_US_LS_cuts_hist->Fill(delta_eta, delta_phi);        
              
              
            }          
          
          }
          else
          {
            L0_inv_mass_vs_L0_inv_mass_US_LS[L_pT_bin_vector_LS.at(iLambdaBck)][L_pT_bin_vector_US.at(iLambda)]->Fill( L_vector_LS.at(iLambdaBck).M(), L_vector_US.at(iLambda).M() );
          
            float theta_star = p_star_vector_LS.at(iLambdaBck).Angle(p_star_vector_US.at(iLambda));
            
            float delta_eta = fabs( L_vector_US.at(iLambda).Eta() - L_vector_LS.at(iLambdaBck).Eta() );
            float delta_phi = fabs( L_vector_US.at(iLambda).Phi() - L_vector_LS.at(iLambdaBck).Phi() );
            
                
            L0_L0_cosThetaProdPlane_US_LS->Fill(TMath::Cos(theta_star));
            L0_L0_cosThetaProdPlane_US_LS_pT_hist[L_pT_bin_vector_LS.at(iLambdaBck)][L_pT_bin_vector_US.at(iLambda)]->Fill(TMath::Cos(theta_star));
            L0_L0_cosThetaProdPlane_US_LS_eta_hist[L_eta_bin_vector_LS.at(iLambdaBck)][L_eta_bin_vector_US.at(iLambda)]->Fill(TMath::Cos(theta_star));
            
            //delta eta vs. delta phi for ME re-weighing   
            L0_L0_delta_eta_vs_delta_phi_US_LS_hist->Fill(delta_eta, delta_phi);
            
            
            if( L_cuts_flag_LS.at(iLambdaBck) == 1 && L_cuts_flag_US.at(iLambda) == 1) //both L in pair passed cuts
            {
              L0_inv_mass_vs_L0_inv_mass_US_LS_cuts[L_pT_bin_vector_LS.at(iLambdaBck)][L_pT_bin_vector_US.at(iLambda)]->Fill( L_vector_LS.at(iLambdaBck).M(), L_vector_US.at(iLambda).M() );
              
             
              L0_L0_cosThetaProdPlane_US_LS_cuts->Fill(TMath::Cos(theta_star));
              L0_L0_cosThetaProdPlane_US_LS_pT_cuts_hist[L_pT_bin_vector_LS.at(iLambdaBck)][L_pT_bin_vector_US.at(iLambda)]->Fill(TMath::Cos(theta_star));
              L0_L0_cosThetaProdPlane_US_LS_eta_cuts_hist[L_eta_bin_vector_LS.at(iLambdaBck)][L_eta_bin_vector_US.at(iLambda)]->Fill(TMath::Cos(theta_star));   
              
              //delta eta vs. delta phi for ME re-weighing   
              L0_L0_delta_eta_vs_delta_phi_US_LS_cuts_hist->Fill(delta_eta, delta_phi);       
              
              
            }
            
          }//end else
          
          nFills_L_L_US_LS++;
        
        }
      
      }
    
    }
    
    
    //Lbar-Lbar
    int nFills_Lbar_Lbar_US_LS = 0;
    
    if( Lbar_vector_US.size() > 0 && Lbar_vector_LS.size() > 0 )
    {
      for( unsigned int iLambdaBar = 0; iLambdaBar < Lbar_vector_US.size(); iLambdaBar++ )
      {
        for( unsigned int iLambdaBarBck = 0; iLambdaBarBck < Lbar_vector_LS.size(); iLambdaBarBck++ )
        {
          //check auto-correlation (prevent daughter sharing in the L-L pair)
          if( pBar_pT_vector_US.at(iLambdaBar) == pBar_pT_vector_LS.at(iLambdaBarBck) ) continue;
          
          //interate order of US and LS for L
          if( nFills_Lbar_Lbar_US_LS % 2 == 0 )
          {
            L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[Lbar_pT_bin_vector_US.at(iLambdaBar)][Lbar_pT_bin_vector_LS.at(iLambdaBarBck)]->Fill( Lbar_vector_US.at(iLambdaBar).M(), Lbar_vector_LS.at(iLambdaBarBck).M() );
          
            float theta_star = pBar_star_vector_US.at(iLambdaBar).Angle(pBar_star_vector_LS.at(iLambdaBarBck));
            
            float delta_eta = fabs( Lbar_vector_US.at(iLambdaBar).Eta() - Lbar_vector_LS.at(iLambdaBarBck).Eta() );
            float delta_phi = fabs( Lbar_vector_US.at(iLambdaBar).Phi() - Lbar_vector_LS.at(iLambdaBarBck).Phi() );
            
                 
            L0bar_L0bar_cosThetaProdPlane_US_LS->Fill(TMath::Cos(theta_star));
            L0bar_L0bar_cosThetaProdPlane_US_LS_pT_hist[Lbar_pT_bin_vector_US.at(iLambdaBar)][Lbar_pT_bin_vector_LS.at(iLambdaBarBck)]->Fill(TMath::Cos(theta_star));
            L0bar_L0bar_cosThetaProdPlane_US_LS_eta_hist[Lbar_eta_bin_vector_US.at(iLambdaBar)][Lbar_eta_bin_vector_LS.at(iLambdaBarBck)]->Fill(TMath::Cos(theta_star));
            
            //delta eta vs. delta phi for ME re-weighing   
            L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_hist->Fill(delta_eta, delta_phi);
            
            
            if( Lbar_cuts_flag_US.at(iLambdaBar) == 1 && Lbar_cuts_flag_LS.at(iLambdaBarBck) == 1) //both L in pair passed cuts
            {
              L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[Lbar_pT_bin_vector_US.at(iLambdaBar)][Lbar_pT_bin_vector_LS.at(iLambdaBarBck)]->Fill( Lbar_vector_US.at(iLambdaBar).M(), Lbar_vector_LS.at(iLambdaBarBck).M() );
              
              
              L0bar_L0bar_cosThetaProdPlane_US_LS_cuts->Fill(TMath::Cos(theta_star));
              L0bar_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[Lbar_pT_bin_vector_US.at(iLambdaBar)][Lbar_pT_bin_vector_LS.at(iLambdaBarBck)]->Fill(TMath::Cos(theta_star));
              L0bar_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[Lbar_eta_bin_vector_US.at(iLambdaBar)][Lbar_eta_bin_vector_LS.at(iLambdaBarBck)]->Fill(TMath::Cos(theta_star));     
              
              //delta eta vs. delta phi for ME re-weighing   
              L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_cuts_hist->Fill(delta_eta, delta_phi);     
              
              
            }          
          
          }
          else
          {
            L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[Lbar_pT_bin_vector_LS.at(iLambdaBarBck)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill( Lbar_vector_LS.at(iLambdaBarBck).M(), Lbar_vector_US.at(iLambdaBar).M() );
          
            float theta_star = pBar_star_vector_LS.at(iLambdaBarBck).Angle(pBar_star_vector_US.at(iLambdaBar));
            
            float delta_eta = fabs( Lbar_vector_US.at(iLambdaBar).Eta() - Lbar_vector_LS.at(iLambdaBarBck).Eta() );
            float delta_phi = fabs( Lbar_vector_US.at(iLambdaBar).Phi() - Lbar_vector_LS.at(iLambdaBarBck).Phi() );
                   
               
            L0bar_L0bar_cosThetaProdPlane_US_LS->Fill(TMath::Cos(theta_star));
            L0bar_L0bar_cosThetaProdPlane_US_LS_pT_hist[Lbar_pT_bin_vector_LS.at(iLambdaBarBck)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
            L0bar_L0bar_cosThetaProdPlane_US_LS_eta_hist[Lbar_eta_bin_vector_LS.at(iLambdaBarBck)][Lbar_eta_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
            
            //delta eta vs. delta phi for ME re-weighing   
            L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_hist->Fill(delta_eta, delta_phi);
            
            
            if( Lbar_cuts_flag_LS.at(iLambdaBarBck) == 1 && Lbar_cuts_flag_US.at(iLambdaBar) == 1) //both L in pair passed cuts
            {
              L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[Lbar_pT_bin_vector_LS.at(iLambdaBarBck)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill( Lbar_vector_LS.at(iLambdaBarBck).M(), Lbar_vector_US.at(iLambdaBar).M() );
              
              
              L0bar_L0bar_cosThetaProdPlane_US_LS_cuts->Fill(TMath::Cos(theta_star));
              L0bar_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[Lbar_pT_bin_vector_LS.at(iLambdaBarBck)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
              L0bar_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[Lbar_eta_bin_vector_LS.at(iLambdaBarBck)][Lbar_eta_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));          
              
              //delta eta vs. delta phi for ME re-weighing   
              L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_cuts_hist->Fill(delta_eta, delta_phi);
              
              
            }
            
          }//end else
          
          nFills_Lbar_Lbar_US_LS++;
        
        }
      
      }
    
    }
    
    //----------------------------------------------------------------
    
    //select L and Lbar for mixed event from pi p pairs
    //before cuts
    if( L_vector_US.size() == 1 && Lbar_vector_US.size() == 0 && p_star_vector_US_ME.size() < 1e3)
    {
      //if( L_Minv_flag_US.at(0) == 0 ) continue;
      
      L_vector_US_ME.push_back(L_vector_US.at(0));
    
      p_star_vector_US_ME.push_back(p_star_vector_US.at(0));

      L_pT_bin_vector_US_ME.push_back(L_pT_bin_vector_US.at(0));
      L_eta_bin_vector_US_ME.push_back(L_eta_bin_vector_US.at(0));     
    }

    if( L_vector_US.size() == 0 && Lbar_vector_US.size() == 1 && pBar_star_vector_US_ME.size() < 1e3)
    {
      //if( Lbar_Minv_flag_US.at(0) == 0 ) continue;
      
      Lbar_vector_US_ME.push_back(Lbar_vector_US.at(0));
      
      pBar_star_vector_US_ME.push_back(pBar_star_vector_US.at(0));

      Lbar_pT_bin_vector_US_ME.push_back(Lbar_pT_bin_vector_US.at(0));
      Lbar_eta_bin_vector_US_ME.push_back(Lbar_eta_bin_vector_US.at(0));

    }
    
    //after cuts
    if( L_vector_US.size() == 1 && Lbar_vector_US.size() == 0 && p_star_vector_US_ME_cuts.size() < 1e3)
    {
      if( L_cuts_flag_US.at(0) == 0 ) continue;
      //if( L_Minv_flag_US.at(0) == 0 ) continue;
      
      L_vector_US_ME_cuts.push_back(L_vector_US.at(0));
      
      p_star_vector_US_ME_cuts.push_back(p_star_vector_US.at(0));

      L_pT_bin_vector_US_ME_cuts.push_back(L_pT_bin_vector_US.at(0));
      L_eta_bin_vector_US_ME_cuts.push_back(L_eta_bin_vector_US.at(0));     
    }

    if( L_vector_US.size() == 0 && Lbar_vector_US.size() == 1 && pBar_star_vector_US_ME_cuts.size() < 1e3)
    {
      //cout<<"fill Lbar ME"<<endl;
      if( Lbar_cuts_flag_US.at(0) == 0 ) continue;
      //if( Lbar_Minv_flag_US.at(0) == 0 ) continue;
      
      Lbar_vector_US_ME_cuts.push_back(Lbar_vector_US.at(0));
      
      pBar_star_vector_US_ME_cuts.push_back(pBar_star_vector_US.at(0));

      Lbar_pT_bin_vector_US_ME_cuts.push_back(Lbar_pT_bin_vector_US.at(0));
      Lbar_eta_bin_vector_US_ME_cuts.push_back(Lbar_eta_bin_vector_US.at(0));

    }

    //_________________________________________________________________________________________________________
    
    //before cuts
    if( L_vector_LS.size() == 1 && Lbar_vector_LS.size() == 0 && p_star_vector_LS_ME.size() < 1e3)
    {
      //if( L_Minv_flag_LS.at(0) == 0 ) continue;
      
      L_vector_LS_ME.push_back(L_vector_LS.at(0));
      
      p_star_vector_LS_ME.push_back(p_star_vector_LS.at(0));

      L_pT_bin_vector_LS_ME.push_back(L_pT_bin_vector_LS.at(0));
      L_eta_bin_vector_LS_ME.push_back(L_eta_bin_vector_LS.at(0));
    }

    if( L_vector_LS.size() == 0 && Lbar_vector_LS.size() == 1 && pBar_star_vector_LS_ME.size() < 1e3)
    {
      //if( Lbar_Minv_flag_LS.at(0) == 0 ) continue;
    
      Lbar_vector_LS_ME.push_back(Lbar_vector_LS.at(0));
      
      pBar_star_vector_LS_ME.push_back(pBar_star_vector_LS.at(0));

      Lbar_pT_bin_vector_LS_ME.push_back(Lbar_pT_bin_vector_LS.at(0));
      Lbar_eta_bin_vector_LS_ME.push_back(Lbar_eta_bin_vector_LS.at(0));
    }
    
    //after cuts
    if( L_vector_LS.size() == 1 && Lbar_vector_LS.size() == 0 && p_star_vector_LS_ME_cuts.size() < 1e3)
    {
      if( L_cuts_flag_LS.at(0) == 0 ) continue;
      //if( L_Minv_flag_LS.at(0) == 0 ) continue;
      
      L_vector_LS_ME_cuts.push_back(L_vector_LS.at(0));
      
      p_star_vector_LS_ME_cuts.push_back(p_star_vector_LS.at(0));

      L_pT_bin_vector_LS_ME_cuts.push_back(L_pT_bin_vector_LS.at(0));
      L_eta_bin_vector_LS_ME_cuts.push_back(L_eta_bin_vector_LS.at(0));     
    }

    if( L_vector_LS.size() == 0 && Lbar_vector_LS.size() == 1 && pBar_star_vector_LS_ME_cuts.size() < 1e3)
    {
      //cout<<"fill Lbar ME"<<endl;
      if( Lbar_cuts_flag_LS.at(0) == 0 ) continue;
      //if( Lbar_Minv_flag_LS.at(0) == 0 ) continue;
      
      Lbar_vector_LS_ME_cuts.push_back(Lbar_vector_LS.at(0));
      
      pBar_star_vector_LS_ME_cuts.push_back(pBar_star_vector_LS.at(0));

      Lbar_pT_bin_vector_LS_ME_cuts.push_back(Lbar_pT_bin_vector_LS.at(0));
      Lbar_eta_bin_vector_LS_ME_cuts.push_back(Lbar_eta_bin_vector_LS.at(0));
    }

  }//end RC Lambda loop
  //_____________________________________________________________________________________________________________________________________________________________________
  
  //cout<<p_star_vector_ME.size()<<endl;
  //cout<<pBar_star_vector_ME.size()<<endl;
  
  cout<<"Start MC ME"<<endl;  
  
  //cout<<p_star_vector_ME.size()<<endl;
  //cout<<pBar_star_vector_ME.size()<<endl;
  
  //mixed event before cuts
  //fill delta eta vs. delta phi for ME first
  for(unsigned int iLambda = 0; iLambda < p_star_vector_ME.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_ME.size(); iLambdaBar++)
    {
      /* 
      if(iLambda >= pBar_star_vector_ME.size()) break; //use each L and Lbar only once
      
      L0_L0bar_eta1_vs_eta2_ME_hist->Fill(L_vector_ME.at(iLambda).Eta() , Lbar_vector_ME.at(iLambda).Eta());
      L0_L0bar_phi1_vs_phi2_ME_hist->Fill(L_vector_ME.at(iLambda).Phi() , Lbar_vector_ME.at(iLambda).Phi());
      L0_L0bar_pT1_vs_pT2_ME_hist->Fill(L_vector_ME.at(iLambda).Pt() , Lbar_vector_ME.at(iLambda).Pt());
      

      */
    
      L0_L0bar_eta1_vs_eta2_ME_hist->Fill(L_vector_ME.at(iLambda).Eta() , Lbar_vector_ME.at(iLambdaBar).Eta());
      L0_L0bar_phi1_vs_phi2_ME_hist->Fill(L_vector_ME.at(iLambda).Phi() , Lbar_vector_ME.at(iLambdaBar).Phi());
      L0_L0bar_pT1_vs_pT2_ME_hist->Fill(L_vector_ME.at(iLambda).Pt() , Lbar_vector_ME.at(iLambdaBar).Pt());
          
      
      /*
      //re-weight ME      
      float delta_eta = fabs( L_vector_ME.at(iLambda).Eta() - Lbar_vector_ME.at(iLambdaBar).Eta() );
      float delta_phi = fabs( L_vector_ME.at(iLambda).Phi() - Lbar_vector_ME.at(iLambdaBar).Phi() );
      float delta_pT = fabs( L_vector_ME.at(iLambda).Pt() - Lbar_vector_ME.at(iLambdaBar).Pt() );
      
      L0_L0bar_delta_eta_vs_delta_phi_ME_hist->Fill(delta_eta, delta_phi);
      //L0_L0bar_delta_eta_vs_delta_phi_vs_delta_pT_ME_hist->Fill(delta_eta, delta_phi, delta_pT);
      L0_L0bar_delta_pT_ME->Fill(delta_pT);
      */
      
    }   
  
  }
  
  
  for(unsigned int iLambda = 0; iLambda < p_star_vector_ME.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_ME.size(); iLambdaBar++)
    {
      /*
      if(iLambda >= pBar_star_vector_ME.size()) break; 
    
      double theta_star = p_star_vector_ME.at(iLambda).Angle(pBar_star_vector_ME.at(iLambda));
      
      //using eta1 vs. eta2, phi2 vs. phi2 and pT1 VS. pT2 distributions
      int eta1_bin = L0_L0bar_eta1_vs_eta2_hist->GetXaxis()->FindBin(L_vector_ME.at(iLambda).Eta());
      int eta2_bin = L0_L0bar_eta1_vs_eta2_hist->GetYaxis()->FindBin(Lbar_vector_ME.at(iLambda).Eta());
      
      int phi1_bin = L0_L0bar_phi1_vs_phi2_hist->GetXaxis()->FindBin(L_vector_ME.at(iLambda).Phi());
      int phi2_bin = L0_L0bar_phi1_vs_phi2_hist->GetYaxis()->FindBin(Lbar_vector_ME.at(iLambda).Phi());
      
      int pT1_bin = L0_L0bar_pT1_vs_pT2_hist->GetXaxis()->FindBin(L_vector_ME.at(iLambda).Pt());
      int pT2_bin = L0_L0bar_pT1_vs_pT2_hist->GetYaxis()->FindBin(Lbar_vector_ME.at(iLambda).Pt());
      
      float weight = 0;
      
      if( L0_L0bar_eta1_vs_eta2_ME_hist->GetBinContent(eta1_bin, eta2_bin) != 0 && L0_L0bar_phi1_vs_phi2_ME_hist->GetBinContent(phi1_bin, phi2_bin) != 0 && L0_L0bar_pT1_vs_pT2_ME_hist->GetBinContent(pT1_bin, pT2_bin) != 0)
      {
        float weight_eta = L0_L0bar_eta1_vs_eta2_hist->GetBinContent(eta1_bin, eta2_bin)/L0_L0bar_eta1_vs_eta2_ME_hist->GetBinContent(eta1_bin, eta2_bin);
        float weight_phi = L0_L0bar_phi1_vs_phi2_hist->GetBinContent(phi1_bin, phi2_bin)/L0_L0bar_phi1_vs_phi2_ME_hist->GetBinContent(phi1_bin, phi2_bin);
        float weight_pT = L0_L0bar_pT1_vs_pT2_hist->GetBinContent(pT1_bin, pT2_bin)/L0_L0bar_pT1_vs_pT2_ME_hist->GetBinContent(pT1_bin, pT2_bin);
        
        weight = weight_eta*weight_phi*weight_pT; 
        //weight = weight_eta;      
      }
      
    
      L0_L0bar_cosThetaProdPlane_ME_weight->Fill(TMath::Cos(theta_star), weight);          
      L0_L0bar_cosThetaProdPlane_ME_weight_pT_hist[L_pT_bin_vector_ME.at(iLambda)][Lbar_pT_bin_vector_ME.at(iLambda)]->Fill(TMath::Cos(theta_star), weight);
      L0_L0bar_cosThetaProdPlane_ME_weight_eta_hist[L_eta_bin_vector_ME.at(iLambda)][Lbar_eta_bin_vector_ME.at(iLambda)]->Fill(TMath::Cos(theta_star), weight);
      */
      
      
      double theta_star = p_star_vector_ME.at(iLambda).Angle(pBar_star_vector_ME.at(iLambdaBar));
      
      //using eta1 vs. eta2, phi2 vs. phi2 and pT1 VS. pT2 distributions
      int eta1_bin = L0_L0bar_eta1_vs_eta2_hist->GetXaxis()->FindBin(L_vector_ME.at(iLambda).Eta());
      int eta2_bin = L0_L0bar_eta1_vs_eta2_hist->GetYaxis()->FindBin(Lbar_vector_ME.at(iLambdaBar).Eta());
      
      int phi1_bin = L0_L0bar_phi1_vs_phi2_hist->GetXaxis()->FindBin(L_vector_ME.at(iLambda).Phi());
      int phi2_bin = L0_L0bar_phi1_vs_phi2_hist->GetYaxis()->FindBin(Lbar_vector_ME.at(iLambdaBar).Phi());
      
      int pT1_bin = L0_L0bar_pT1_vs_pT2_hist->GetXaxis()->FindBin(L_vector_ME.at(iLambda).Pt());
      int pT2_bin = L0_L0bar_pT1_vs_pT2_hist->GetYaxis()->FindBin(Lbar_vector_ME.at(iLambdaBar).Pt());
      
      //float weight = 0;
      float weight = 1;
      
      if( L0_L0bar_eta1_vs_eta2_ME_hist->GetBinContent(eta1_bin, eta2_bin) != 0 && L0_L0bar_phi1_vs_phi2_ME_hist->GetBinContent(phi1_bin, phi2_bin) != 0 && L0_L0bar_pT1_vs_pT2_ME_hist->GetBinContent(pT1_bin, pT2_bin) != 0)
      {
        float weight_eta = L0_L0bar_eta1_vs_eta2_hist->GetBinContent(eta1_bin, eta2_bin)/L0_L0bar_eta1_vs_eta2_ME_hist->GetBinContent(eta1_bin, eta2_bin);
        float weight_phi = L0_L0bar_phi1_vs_phi2_hist->GetBinContent(phi1_bin, phi2_bin)/L0_L0bar_phi1_vs_phi2_ME_hist->GetBinContent(phi1_bin, phi2_bin);
        float weight_pT = L0_L0bar_pT1_vs_pT2_hist->GetBinContent(pT1_bin, pT2_bin)/L0_L0bar_pT1_vs_pT2_ME_hist->GetBinContent(pT1_bin, pT2_bin);
        
        //weight = weight_eta*weight_phi*weight_pT; 
        //weight = weight_eta;      
      }
      
    
      L0_L0bar_cosThetaProdPlane_ME_weight->Fill(TMath::Cos(theta_star), weight);          
      L0_L0bar_cosThetaProdPlane_ME_weight_pT_hist[L_pT_bin_vector_ME.at(iLambda)][Lbar_pT_bin_vector_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star), weight);
      L0_L0bar_cosThetaProdPlane_ME_weight_eta_hist[L_eta_bin_vector_ME.at(iLambda)][Lbar_eta_bin_vector_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star), weight);
      
    }   
  
  }
  
  //---------------------------------------------------------------------------------------------------------------
  
  //L0-L0
  for(unsigned int iLambda1 = 0; iLambda1 < p_star_vector_ME.size(); iLambda1++)
  {
    //if(iLambda1 > 1e3 ) break;
    
    for(unsigned int iLambda2 = iLambda1+1; iLambda2 < p_star_vector_ME.size(); iLambda2++)
    {
      //if(iLambda2 > 1e3 ) break;
      //re-weight ME      
      float delta_eta = fabs( L_vector_ME.at(iLambda1).Eta() - L_vector_ME.at(iLambda2).Eta() );
      float delta_phi = fabs( L_vector_ME.at(iLambda1).Phi() - L_vector_ME.at(iLambda2).Phi() );
      
      L0_L0_delta_eta_vs_delta_phi_ME_hist->Fill(delta_eta, delta_phi);
    }   
  
  }
  
  for(unsigned int iLambda1 = 0; iLambda1 < p_star_vector_ME.size(); iLambda1++)
  {
    //if(iLambda1 > 1e3 ) break;
    
    for(unsigned int iLambda2 = iLambda1+1; iLambda2 < p_star_vector_ME.size(); iLambda2++)
    {
      //if(iLambda2 > 1e3 ) break;
      
      double theta_star = p_star_vector_ME.at(iLambda1).Angle(p_star_vector_ME.at(iLambda2));
/*      
      L0_L0_cosThetaProdPlane_ME->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_ME_pT_hist[L_pT_bin_vector_ME.at(iLambda1)][L_pT_bin_vector_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_ME_eta_hist[L_eta_bin_vector_ME.at(iLambda1)][L_eta_bin_vector_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star));
*/      
      //re-weight ME      
      float delta_eta = fabs( L_vector_ME.at(iLambda1).Eta() - L_vector_ME.at(iLambda2).Eta() );
      float delta_phi = fabs( L_vector_ME.at(iLambda1).Phi() - L_vector_ME.at(iLambda2).Phi() );
      
      int delta_eta_bin = L0_L0_delta_eta_vs_delta_phi_hist->GetXaxis()->FindBin(delta_eta);
      int delta_phi_bin = L0_L0_delta_eta_vs_delta_phi_hist->GetYaxis()->FindBin(delta_phi);
      
      //weight = (delta eta vs. delta phi)_true/(delta eta vs. delta phi)_ME
      //done manually, bin by bin, to preserve original delta eta vs. delta phi histograms
      //float weight = 0; 
      float weight = 1;
      
      //check that denominator is not 0    
      if(L0_L0_delta_eta_vs_delta_phi_ME_hist->GetBinContent(delta_eta_bin, delta_phi_bin) != 0)
      {
        //weight = L0_L0_delta_eta_vs_delta_phi_hist->GetBinContent(delta_eta_bin, delta_phi_bin)/L0_L0_delta_eta_vs_delta_phi_ME_hist->GetBinContent(delta_eta_bin, delta_phi_bin);
      }       
      
      L0_L0_cosThetaProdPlane_ME_weight->Fill(TMath::Cos(theta_star), weight);          
      L0_L0_cosThetaProdPlane_ME_weight_pT_hist[L_pT_bin_vector_ME.at(iLambda1)][L_pT_bin_vector_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star), weight);
      L0_L0_cosThetaProdPlane_ME_weight_eta_hist[L_eta_bin_vector_ME.at(iLambda1)][L_eta_bin_vector_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star), weight);
    }   
  
  }
  
  //---------------------------------------------------------------------------------------------------------------  

  //L0bar-L0bar
  for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pBar_star_vector_ME.size(); iLambdaBar1++)
  {
    //if(iLambdaBar1 > 1e3 ) break;
    
    for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < pBar_star_vector_ME.size(); iLambdaBar2++)
    {
      //if(iLambdaBar2 > 1e3 ) break;
      
      //re-weight ME      
      float delta_eta = fabs( Lbar_vector_ME.at(iLambdaBar1).Eta() - Lbar_vector_ME.at(iLambdaBar2).Eta() );
      float delta_phi = fabs( Lbar_vector_ME.at(iLambdaBar1).Phi() - Lbar_vector_ME.at(iLambdaBar2).Phi() );
      
      L0bar_L0bar_delta_eta_vs_delta_phi_ME_hist->Fill(delta_eta, delta_phi);     
    }   
  
  }
  
  for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pBar_star_vector_ME.size(); iLambdaBar1++)
  {
    //if(iLambdaBar1 > 1e3 ) break;
    
    for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < pBar_star_vector_ME.size(); iLambdaBar2++)
    {
      //if(iLambdaBar2 > 1e3 ) break;
      
      double theta_star = pBar_star_vector_ME.at(iLambdaBar1).Angle(pBar_star_vector_ME.at(iLambdaBar2));
/*      
      L0bar_L0bar_cosThetaProdPlane_ME->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_ME_pT_hist[Lbar_pT_bin_vector_ME.at(iLambdaBar1)][Lbar_pT_bin_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_ME_eta_hist[Lbar_eta_bin_vector_ME.at(iLambdaBar1)][Lbar_eta_bin_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));     
*/      
      //re-weight ME      
      float delta_eta = fabs( Lbar_vector_ME.at(iLambdaBar1).Eta() - Lbar_vector_ME.at(iLambdaBar2).Eta() );
      float delta_phi = fabs( Lbar_vector_ME.at(iLambdaBar1).Phi() - Lbar_vector_ME.at(iLambdaBar2).Phi() );
      
      int delta_eta_bin = L0bar_L0bar_delta_eta_vs_delta_phi_hist->GetXaxis()->FindBin(delta_eta);
      int delta_phi_bin = L0bar_L0bar_delta_eta_vs_delta_phi_hist->GetYaxis()->FindBin(delta_phi);
          
      //float weight = 0;
      float weight = 1;
      
      if(L0bar_L0bar_delta_eta_vs_delta_phi_ME_hist->GetBinContent(delta_eta_bin, delta_phi_bin) != 0)
      {
        //weight = L0bar_L0bar_delta_eta_vs_delta_phi_hist->GetBinContent(delta_eta_bin, delta_phi_bin)/L0bar_L0bar_delta_eta_vs_delta_phi_ME_hist->GetBinContent(delta_eta_bin, delta_phi_bin);
      }
       
      
      L0bar_L0bar_cosThetaProdPlane_ME_weight->Fill(TMath::Cos(theta_star), weight);          
      L0bar_L0bar_cosThetaProdPlane_ME_weight_pT_hist[Lbar_pT_bin_vector_ME.at(iLambdaBar1)][Lbar_pT_bin_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star), weight);
      L0bar_L0bar_cosThetaProdPlane_ME_weight_eta_hist[Lbar_eta_bin_vector_ME.at(iLambdaBar1)][Lbar_eta_bin_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star), weight);       
    }   
  
  }
  //_____________________________________________________________________________________________________________________________________________________________________
  

  //mixed event after cuts
  for(unsigned int iLambda = 0; iLambda < p_star_cuts_vector_ME.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_cuts_vector_ME.size(); iLambdaBar++)
    {
      //if( fabs( L_cuts_vector_ME.at(iLambda).Eta() - Lbar_cuts_vector_ME.at(iLambdaBar).Eta()) < 0.3 ) continue;
            
      L0_L0bar_eta1_vs_eta2_ME_cuts_hist->Fill(L_cuts_vector_ME.at(iLambda).Eta() , Lbar_cuts_vector_ME.at(iLambdaBar).Eta());
      L0_L0bar_phi1_vs_phi2_ME_cuts_hist->Fill(L_cuts_vector_ME.at(iLambda).Phi() , Lbar_cuts_vector_ME.at(iLambdaBar).Phi());
      L0_L0bar_pT1_vs_pT2_ME_cuts_hist->Fill(L_cuts_vector_ME.at(iLambda).Pt() , Lbar_cuts_vector_ME.at(iLambdaBar).Pt());
      
      L0_L0bar_p1_pT1_vs_p2_pT2_ME_cuts_hist->Fill(L_p_pT_cuts_vector_ME.at(iLambda) , Lbar_p_pT_cuts_vector_ME.at(iLambdaBar));
      L0_L0bar_pi1_pT1_vs_pi2_pT2_ME_cuts_hist->Fill(L_pi_pT_cuts_vector_ME.at(iLambda) , Lbar_pi_pT_cuts_vector_ME.at(iLambdaBar));
      
      L0_L0bar_p1_eta1_vs_p2_eta2_ME_cuts_hist->Fill(L_p_eta_cuts_vector_ME.at(iLambda) , Lbar_p_eta_cuts_vector_ME.at(iLambdaBar));
      L0_L0bar_pi1_eta1_vs_pi2_eta2_ME_cuts_hist->Fill(L_pi_eta_cuts_vector_ME.at(iLambda) , Lbar_pi_eta_cuts_vector_ME.at(iLambdaBar));
      
      L0_L0bar_p1_phi1_vs_p2_phi2_ME_cuts_hist->Fill(L_p_phi_cuts_vector_ME.at(iLambda) , Lbar_p_phi_cuts_vector_ME.at(iLambdaBar));
      L0_L0bar_pi1_phi1_vs_pi2_phi2_ME_cuts_hist->Fill(L_pi_phi_cuts_vector_ME.at(iLambda) , Lbar_pi_phi_cuts_vector_ME.at(iLambdaBar));
      
    }   
  
  }
  
/*
  for(unsigned int iLambda = 0; iLambda < p_star_cuts_vector_ME.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_cuts_vector_ME.size(); iLambdaBar++)
    {

      int phi1_bin = L0_L0bar_phi1_vs_phi2_hist->GetXaxis()->FindBin(L_cuts_vector_ME.at(iLambda).Phi());
      int phi2_bin = L0_L0bar_phi1_vs_phi2_hist->GetYaxis()->FindBin(Lbar_cuts_vector_ME.at(iLambdaBar).Phi());
      
      
            
      float weight = 0;
      
      if(  L0_L0bar_phi1_vs_phi2_ME_cuts_hist->GetBinContent(phi1_bin, phi2_bin) != 0 )
      {    
        
        float weight_phi = L0_L0bar_phi1_vs_phi2_cuts_hist->GetBinContent(phi1_bin, phi2_bin)/L0_L0bar_phi1_vs_phi2_ME_cuts_hist->GetBinContent(phi1_bin, phi2_bin);
         
        weight = weight_phi;
        
      }
      
      
      L0_L0bar_eta1_vs_eta2_ME_cuts_weight_hist->Fill(L_cuts_vector_ME.at(iLambda).Eta() , Lbar_cuts_vector_ME.at(iLambdaBar).Eta(), weight );
            
    }   
  
  }

  for(unsigned int iLambda = 0; iLambda < p_star_cuts_vector_ME.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_cuts_vector_ME.size(); iLambdaBar++)
    {

      int phi1_bin = L0_L0bar_phi1_vs_phi2_hist->GetXaxis()->FindBin(L_cuts_vector_ME.at(iLambda).Phi());
      int phi2_bin = L0_L0bar_phi1_vs_phi2_hist->GetYaxis()->FindBin(Lbar_cuts_vector_ME.at(iLambdaBar).Phi());
      
      int eta1_bin = L0_L0bar_eta1_vs_eta2_hist->GetXaxis()->FindBin(L_cuts_vector_ME.at(iLambda).Eta());
      int eta2_bin = L0_L0bar_eta1_vs_eta2_hist->GetYaxis()->FindBin(Lbar_cuts_vector_ME.at(iLambdaBar).Eta());      
      
            
      float weight = 0;
      
      if(  L0_L0bar_phi1_vs_phi2_ME_cuts_hist->GetBinContent(phi1_bin, phi2_bin) != 0 )
      {    
        
        float weight_phi = L0_L0bar_phi1_vs_phi2_cuts_hist->GetBinContent(phi1_bin, phi2_bin)/L0_L0bar_phi1_vs_phi2_ME_cuts_hist->GetBinContent(phi1_bin, phi2_bin);
        float weight_eta = L0_L0bar_eta1_vs_eta2_cuts_hist->GetBinContent(eta1_bin, eta2_bin)/L0_L0bar_eta1_vs_eta2_ME_cuts_weight_hist->GetBinContent(eta1_bin, eta2_bin);
        
         
        weight = weight_phi*weight_eta;
        
      }
      
      
      L0_L0bar_pT1_vs_pT2_ME_cuts_weight_hist->Fill(L_cuts_vector_ME.at(iLambda).Pt() , Lbar_cuts_vector_ME.at(iLambdaBar).Pt(), weight );
      
    }   
  
  }
*/
  for(unsigned int iLambda = 0; iLambda < p_star_cuts_vector_ME.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_cuts_vector_ME.size(); iLambdaBar++)
    {
      //if( fabs( L_cuts_vector_ME.at(iLambda).Eta() - Lbar_cuts_vector_ME.at(iLambdaBar).Eta()) < 0.3 ) continue;
                       
      //using eta1 vs. eta2, phi2 vs. phi2 and pT1 VS. pT2 distributions
      int eta1_bin = L0_L0bar_eta1_vs_eta2_hist->GetXaxis()->FindBin(L_cuts_vector_ME.at(iLambda).Eta());
      int eta2_bin = L0_L0bar_eta1_vs_eta2_hist->GetYaxis()->FindBin(Lbar_cuts_vector_ME.at(iLambdaBar).Eta());
      
      int phi1_bin = L0_L0bar_phi1_vs_phi2_hist->GetXaxis()->FindBin(L_cuts_vector_ME.at(iLambda).Phi());
      int phi2_bin = L0_L0bar_phi1_vs_phi2_hist->GetYaxis()->FindBin(Lbar_cuts_vector_ME.at(iLambdaBar).Phi());
      
      int pT1_bin = L0_L0bar_pT1_vs_pT2_hist->GetXaxis()->FindBin(L_cuts_vector_ME.at(iLambda).Pt());
      int pT2_bin = L0_L0bar_pT1_vs_pT2_hist->GetYaxis()->FindBin(Lbar_cuts_vector_ME.at(iLambdaBar).Pt());
      
      
            
      float weight = 0;
      
      if( L0_L0bar_eta1_vs_eta2_ME_cuts_hist->GetBinContent(eta1_bin, eta2_bin) != 0 && L0_L0bar_phi1_vs_phi2_ME_cuts_hist->GetBinContent(phi1_bin, phi2_bin) != 0 && L0_L0bar_pT1_vs_pT2_ME_cuts_hist->GetBinContent(pT1_bin, pT2_bin) != 0 )
      //if( L0_L0bar_eta1_vs_eta2_ME_cuts_weight_hist->GetBinContent(eta1_bin, eta2_bin) != 0 && L0_L0bar_phi1_vs_phi2_ME_cuts_hist->GetBinContent(phi1_bin, phi2_bin) != 0 && L0_L0bar_pT1_vs_pT2_ME_cuts_weight_hist->GetBinContent(pT1_bin, pT2_bin) != 0 )
      {
      
        float weight_eta = L0_L0bar_eta1_vs_eta2_cuts_hist->GetBinContent(eta1_bin, eta2_bin)/L0_L0bar_eta1_vs_eta2_ME_cuts_hist->GetBinContent(eta1_bin, eta2_bin);
        float weight_phi = L0_L0bar_phi1_vs_phi2_cuts_hist->GetBinContent(phi1_bin, phi2_bin)/L0_L0bar_phi1_vs_phi2_ME_cuts_hist->GetBinContent(phi1_bin, phi2_bin); //phi is first - use basic ME before weight
        float weight_pT = L0_L0bar_pT1_vs_pT2_cuts_hist->GetBinContent(pT1_bin, pT2_bin)/L0_L0bar_pT1_vs_pT2_ME_cuts_hist->GetBinContent(pT1_bin, pT2_bin);
                
        
        weight = weight_eta*weight_phi*weight_pT;
        
      }
      
      
      L0_L0bar_p1_pT1_vs_p2_pT2_ME_cuts_hist_weight->Fill(L_p_pT_cuts_vector_ME.at(iLambda) , Lbar_p_pT_cuts_vector_ME.at(iLambdaBar), weight);
      L0_L0bar_pi1_pT1_vs_pi2_pT2_ME_cuts_hist_weight->Fill(L_pi_pT_cuts_vector_ME.at(iLambda) , Lbar_pi_pT_cuts_vector_ME.at(iLambdaBar), weight);
      
      L0_L0bar_p1_eta1_vs_p2_eta2_ME_cuts_hist_weight->Fill(L_p_eta_cuts_vector_ME.at(iLambda) , Lbar_p_eta_cuts_vector_ME.at(iLambdaBar), weight);
      L0_L0bar_pi1_eta1_vs_pi2_eta2_ME_cuts_hist_weight->Fill(L_pi_eta_cuts_vector_ME.at(iLambda) , Lbar_pi_eta_cuts_vector_ME.at(iLambdaBar), weight);
      
      L0_L0bar_p1_phi1_vs_p2_phi2_ME_cuts_hist_weight->Fill(L_p_phi_cuts_vector_ME.at(iLambda) , Lbar_p_phi_cuts_vector_ME.at(iLambdaBar), weight);
      L0_L0bar_pi1_phi1_vs_pi2_phi2_ME_cuts_hist_weight->Fill(L_pi_phi_cuts_vector_ME.at(iLambda) , Lbar_pi_phi_cuts_vector_ME.at(iLambdaBar), weight);
      
    }   
  
  }
  
  L0_L0bar_p1_pT1_vs_p2_pT2_ME_cuts_hist_weight->Scale(L0_L0bar_p1_pT1_vs_p2_pT2_cuts_hist->Integral()/L0_L0bar_p1_pT1_vs_p2_pT2_ME_cuts_hist_weight->Integral());
  L0_L0bar_pi1_pT1_vs_pi2_pT2_ME_cuts_hist_weight->Scale(L0_L0bar_pi1_pT1_vs_pi2_pT2_cuts_hist->Integral()/L0_L0bar_pi1_pT1_vs_pi2_pT2_ME_cuts_hist_weight->Integral());
  
  L0_L0bar_p1_eta1_vs_p2_eta2_ME_cuts_hist_weight->Scale(L0_L0bar_p1_eta1_vs_p2_eta2_cuts_hist->Integral()/L0_L0bar_p1_eta1_vs_p2_eta2_ME_cuts_hist_weight->Integral());
  L0_L0bar_pi1_eta1_vs_pi2_eta2_ME_cuts_hist_weight->Scale(L0_L0bar_pi1_eta1_vs_pi2_eta2_cuts_hist->Integral()/L0_L0bar_pi1_eta1_vs_pi2_eta2_ME_cuts_hist_weight->Integral());
  
  L0_L0bar_p1_phi1_vs_p2_phi2_ME_cuts_hist_weight->Scale(L0_L0bar_p1_phi1_vs_p2_phi2_cuts_hist->Integral()/L0_L0bar_p1_phi1_vs_p2_phi2_ME_cuts_hist_weight->Integral());
  L0_L0bar_pi1_phi1_vs_pi2_phi2_ME_cuts_hist_weight->Scale(L0_L0bar_pi1_phi1_vs_pi2_phi2_cuts_hist->Integral()/L0_L0bar_pi1_phi1_vs_pi2_phi2_ME_cuts_hist_weight->Integral());
  
  
  for(unsigned int iLambda = 0; iLambda < p_star_cuts_vector_ME.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_cuts_vector_ME.size(); iLambdaBar++)
    {
      //if( fabs( L_cuts_vector_ME.at(iLambda).Eta() - Lbar_cuts_vector_ME.at(iLambdaBar).Eta()) < 0.3 ) continue;
            
      double theta_star = p_star_cuts_vector_ME.at(iLambda).Angle(pBar_star_cuts_vector_ME.at(iLambdaBar));
           
      //using eta1 vs. eta2, phi2 vs. phi2 and pT1 VS. pT2 distributions
      int eta1_bin = L0_L0bar_eta1_vs_eta2_hist->GetXaxis()->FindBin(L_cuts_vector_ME.at(iLambda).Eta());
      int eta2_bin = L0_L0bar_eta1_vs_eta2_hist->GetYaxis()->FindBin(Lbar_cuts_vector_ME.at(iLambdaBar).Eta());
      
      int phi1_bin = L0_L0bar_phi1_vs_phi2_hist->GetXaxis()->FindBin(L_cuts_vector_ME.at(iLambda).Phi());
      int phi2_bin = L0_L0bar_phi1_vs_phi2_hist->GetYaxis()->FindBin(Lbar_cuts_vector_ME.at(iLambdaBar).Phi());
      
      int pT1_bin = L0_L0bar_pT1_vs_pT2_hist->GetXaxis()->FindBin(L_cuts_vector_ME.at(iLambda).Pt());
      int pT2_bin = L0_L0bar_pT1_vs_pT2_hist->GetYaxis()->FindBin(Lbar_cuts_vector_ME.at(iLambdaBar).Pt());
      
      //bins of daughters
      int p_phi1_bin = L0_L0bar_p1_phi1_vs_p2_phi2_cuts_hist->GetXaxis()->FindBin(L_p_phi_cuts_vector_ME.at(iLambda));
      int p_phi2_bin = L0_L0bar_p1_phi1_vs_p2_phi2_cuts_hist->GetYaxis()->FindBin(Lbar_p_phi_cuts_vector_ME.at(iLambdaBar));
      
      int pi_phi1_bin = L0_L0bar_pi1_phi1_vs_pi2_phi2_cuts_hist->GetXaxis()->FindBin(L_pi_phi_cuts_vector_ME.at(iLambda));
      int pi_phi2_bin = L0_L0bar_pi1_phi1_vs_pi2_phi2_cuts_hist->GetYaxis()->FindBin(Lbar_pi_phi_cuts_vector_ME.at(iLambdaBar));
      
      int p_pT1_bin = L0_L0bar_p1_pT1_vs_p2_pT2_cuts_hist->GetXaxis()->FindBin(L_p_pT_cuts_vector_ME.at(iLambda));
      int p_pT2_bin = L0_L0bar_p1_pT1_vs_p2_pT2_cuts_hist->GetYaxis()->FindBin(Lbar_p_pT_cuts_vector_ME.at(iLambdaBar));
      
      int pi_pT1_bin = L0_L0bar_pi1_pT1_vs_pi2_pT2_cuts_hist->GetXaxis()->FindBin(L_pi_pT_cuts_vector_ME.at(iLambda));
      int pi_pT2_bin = L0_L0bar_pi1_pT1_vs_pi2_pT2_cuts_hist->GetYaxis()->FindBin(Lbar_pi_pT_cuts_vector_ME.at(iLambdaBar));
      
      int p_eta1_bin = L0_L0bar_p1_eta1_vs_p2_eta2_cuts_hist->GetXaxis()->FindBin(L_p_eta_cuts_vector_ME.at(iLambda));
      int p_eta2_bin = L0_L0bar_p1_eta1_vs_p2_eta2_cuts_hist->GetYaxis()->FindBin(Lbar_p_eta_cuts_vector_ME.at(iLambdaBar));
      
      int pi_eta1_bin = L0_L0bar_pi1_eta1_vs_pi2_eta2_cuts_hist->GetXaxis()->FindBin(L_pi_eta_cuts_vector_ME.at(iLambda));
      int pi_eta2_bin = L0_L0bar_pi1_eta1_vs_pi2_eta2_cuts_hist->GetYaxis()->FindBin(Lbar_pi_eta_cuts_vector_ME.at(iLambdaBar));

            
      //float weight = 0;
      float weight = 1;
      
      //if( L0_L0bar_eta1_vs_eta2_ME_cuts_hist->GetBinContent(eta1_bin, eta2_bin) != 0 && L0_L0bar_phi1_vs_phi2_ME_cuts_hist->GetBinContent(phi1_bin, phi2_bin) != 0 && L0_L0bar_pT1_vs_pT2_ME_cuts_hist->GetBinContent(pT1_bin, pT2_bin) != 0 )
      if( L0_L0bar_eta1_vs_eta2_ME_cuts_hist->GetBinContent(eta1_bin, eta2_bin) != 0 && L0_L0bar_phi1_vs_phi2_ME_cuts_hist->GetBinContent(phi1_bin, phi2_bin) != 0 && L0_L0bar_pT1_vs_pT2_ME_cuts_hist->GetBinContent(pT1_bin, pT2_bin) != 0 && L0_L0bar_pi1_pT1_vs_pi2_pT2_ME_cuts_hist_weight->GetBinContent(pi_pT1_bin, pi_pT2_bin) != 0 && L0_L0bar_pi1_eta1_vs_pi2_eta2_ME_cuts_hist_weight->GetBinContent(pi_eta1_bin, pi_eta2_bin) != 0 && L0_L0bar_pi1_phi1_vs_pi2_phi2_ME_cuts_hist_weight->GetBinContent(pi_phi1_bin, pi_phi2_bin) != 0 && L0_L0bar_p1_pT1_vs_p2_pT2_ME_cuts_hist_weight->GetBinContent(p_pT1_bin, p_pT2_bin) != 0 && L0_L0bar_p1_eta1_vs_p2_eta2_ME_cuts_hist_weight->GetBinContent(p_eta1_bin, p_eta2_bin) != 0 && L0_L0bar_p1_phi1_vs_p2_phi2_ME_cuts_hist_weight->GetBinContent(p_phi1_bin, p_phi2_bin) != 0 )
      //if( L0_L0bar_eta1_vs_eta2_ME_hist->GetBinContent(eta1_bin, eta2_bin) != 0 && L0_L0bar_phi1_vs_phi2_ME_hist->GetBinContent(phi1_bin, phi2_bin) != 0 && L0_L0bar_pT1_vs_pT2_ME_hist->GetBinContent(pT1_bin, pT2_bin) != 0 && L0_L0bar_pi1_pT1_vs_pi2_pT2_ME_cuts_hist_weight->GetBinContent(pi_pT1_bin, pi_pT2_bin) != 0 && L0_L0bar_pi1_eta1_vs_pi2_eta2_ME_cuts_hist_weight->GetBinContent(pi_eta1_bin, pi_eta2_bin) != 0 && L0_L0bar_pi1_phi1_vs_pi2_phi2_ME_cuts_hist_weight->GetBinContent(pi_phi1_bin, pi_phi2_bin) != 0 && L0_L0bar_p1_pT1_vs_p2_pT2_ME_cuts_hist_weight->GetBinContent(p_pT1_bin, p_pT2_bin) != 0 && L0_L0bar_p1_eta1_vs_p2_eta2_ME_cuts_hist_weight->GetBinContent(p_eta1_bin, p_eta2_bin) != 0 && L0_L0bar_p1_phi1_vs_p2_phi2_ME_cuts_hist_weight->GetBinContent(p_phi1_bin, p_phi2_bin) != 0 )
      //if( L0_L0bar_eta1_vs_eta2_ME_hist->GetBinContent(eta1_bin, eta2_bin) != 0 && L0_L0bar_phi1_vs_phi2_ME_hist->GetBinContent(phi1_bin, phi2_bin) != 0 && L0_L0bar_pT1_vs_pT2_ME_hist->GetBinContent(pT1_bin, pT2_bin) != 0)
      {
      
        float weight_eta = L0_L0bar_eta1_vs_eta2_cuts_hist->GetBinContent(eta1_bin, eta2_bin)/L0_L0bar_eta1_vs_eta2_ME_cuts_hist->GetBinContent(eta1_bin, eta2_bin);
        float weight_phi = L0_L0bar_phi1_vs_phi2_cuts_hist->GetBinContent(phi1_bin, phi2_bin)/L0_L0bar_phi1_vs_phi2_ME_cuts_hist->GetBinContent(phi1_bin, phi2_bin);
        float weight_pT = L0_L0bar_pT1_vs_pT2_cuts_hist->GetBinContent(pT1_bin, pT2_bin)/L0_L0bar_pT1_vs_pT2_ME_cuts_hist->GetBinContent(pT1_bin, pT2_bin);
        
        float weight_p_pT = L0_L0bar_p1_pT1_vs_p2_pT2_cuts_hist->GetBinContent(p_pT1_bin, p_pT2_bin)/L0_L0bar_p1_pT1_vs_p2_pT2_ME_cuts_hist_weight->GetBinContent(p_pT1_bin, p_pT2_bin);
        float weight_pi_pT = L0_L0bar_pi1_pT1_vs_pi2_pT2_cuts_hist->GetBinContent(pi_pT1_bin, pi_pT2_bin)/L0_L0bar_pi1_pT1_vs_pi2_pT2_ME_cuts_hist_weight->GetBinContent(pi_pT1_bin, pi_pT2_bin);
        
        float weight_p_eta = L0_L0bar_p1_eta1_vs_p2_eta2_cuts_hist->GetBinContent(p_eta1_bin, p_eta2_bin)/L0_L0bar_p1_eta1_vs_p2_eta2_ME_cuts_hist_weight->GetBinContent(p_eta1_bin, p_eta2_bin);
        float weight_pi_eta = L0_L0bar_pi1_eta1_vs_pi2_eta2_cuts_hist->GetBinContent(pi_eta1_bin, pi_eta2_bin)/L0_L0bar_pi1_eta1_vs_pi2_eta2_ME_cuts_hist_weight->GetBinContent(pi_eta1_bin, pi_eta2_bin);
        
        float weight_p_phi = L0_L0bar_p1_phi1_vs_p2_phi2_cuts_hist->GetBinContent(p_phi1_bin, p_phi2_bin)/L0_L0bar_p1_phi1_vs_p2_phi2_ME_cuts_hist_weight->GetBinContent(p_phi1_bin, p_phi2_bin);
        float weight_pi_phi = L0_L0bar_pi1_phi1_vs_pi2_phi2_cuts_hist->GetBinContent(pi_phi1_bin, pi_phi2_bin)/L0_L0bar_pi1_phi1_vs_pi2_phi2_ME_cuts_hist_weight->GetBinContent(pi_phi1_bin, pi_phi2_bin);
               
        //float weight_eta = L0_L0bar_eta1_vs_eta2_hist->GetBinContent(eta1_bin, eta2_bin)/L0_L0bar_eta1_vs_eta2_ME_hist->GetBinContent(eta1_bin, eta2_bin);
        //float weight_phi = L0_L0bar_phi1_vs_phi2_hist->GetBinContent(phi1_bin, phi2_bin)/L0_L0bar_phi1_vs_phi2_ME_hist->GetBinContent(phi1_bin, phi2_bin);
        //float weight_pT = L0_L0bar_pT1_vs_pT2_hist->GetBinContent(pT1_bin, pT2_bin)/L0_L0bar_pT1_vs_pT2_ME_hist->GetBinContent(pT1_bin, pT2_bin);

        
        //weight = weight_eta*weight_phi*weight_pT;
        //weight = weight_eta*weight_phi*weight_pT*weight_pi_pT;
        
        //weight = weight_eta*weight_phi*weight_pT*weight_pi_pT*weight_pi_eta*weight_pi_phi;
        //weight = weight_eta*weight_phi*weight_pT*weight_pi_pT*weight_pi_eta*weight_pi_phi*weight_p_pT*weight_p_eta*weight_p_phi;
        
        
        //weight = weight_pi_eta*weight_pi_phi*weight_pi_pT;
        //weight = weight_pi_pT;
        //weight = weight_phi;
      }
      
      if(isnan(weight)) continue;
      
      L0_L0bar_cosThetaProdPlane_ME_weight_cuts->Fill(TMath::Cos(theta_star), weight);          
      L0_L0bar_cosThetaProdPlane_ME_weight_pT_cuts_hist[L_pT_bin_cuts_vector_ME.at(iLambda)][Lbar_pT_bin_cuts_vector_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star), weight);
      L0_L0bar_cosThetaProdPlane_ME_weight_eta_cuts_hist[L_eta_bin_cuts_vector_ME.at(iLambda)][Lbar_eta_bin_cuts_vector_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star), weight);
      
      L0_L0bar_pT1_vs_pT2_ME_cuts_weight_hist->Fill(L_cuts_vector_ME.at(iLambda).Pt() , Lbar_cuts_vector_ME.at(iLambdaBar).Pt(), weight );
      L0_L0bar_eta1_vs_eta2_ME_cuts_weight_hist->Fill(L_cuts_vector_ME.at(iLambda).Eta() , Lbar_cuts_vector_ME.at(iLambdaBar).Eta(), weight );
      L0_L0bar_phi1_vs_phi2_ME_cuts_weight_hist->Fill(L_cuts_vector_ME.at(iLambda).Phi(), Lbar_cuts_vector_ME.at(iLambdaBar).Phi(), weight);
      /*
      L0_L0bar_p1_pT1_vs_p2_pT2_ME_cuts_hist_weight->Fill(L_p_pT_cuts_vector_ME.at(iLambda) , Lbar_p_pT_cuts_vector_ME.at(iLambdaBar), weight);
      L0_L0bar_pi1_pT1_vs_pi2_pT2_ME_cuts_hist_weight->Fill(L_pi_pT_cuts_vector_ME.at(iLambda) , Lbar_pi_pT_cuts_vector_ME.at(iLambdaBar), weight);
      
      L0_L0bar_p1_eta1_vs_p2_eta2_ME_cuts_hist_weight->Fill(L_p_eta_cuts_vector_ME.at(iLambda) , Lbar_p_eta_cuts_vector_ME.at(iLambdaBar), weight);
      L0_L0bar_pi1_eta1_vs_pi2_eta2_ME_cuts_hist_weight->Fill(L_pi_eta_cuts_vector_ME.at(iLambda) , Lbar_pi_eta_cuts_vector_ME.at(iLambdaBar), weight);
      
      L0_L0bar_p1_phi1_vs_p2_phi2_ME_cuts_hist_weight->Fill(L_p_phi_cuts_vector_ME.at(iLambda) , Lbar_p_phi_cuts_vector_ME.at(iLambdaBar), weight);
      L0_L0bar_pi1_phi1_vs_pi2_phi2_ME_cuts_hist_weight->Fill(L_pi_phi_cuts_vector_ME.at(iLambda) , Lbar_pi_phi_cuts_vector_ME.at(iLambdaBar), weight);
      */
    }   
  
  }
  
  //---------------------------------------------------------------------------------------------------------------------------

  //L0-L0
  for(unsigned int iLambda1 = 0; iLambda1 < p_star_cuts_vector_ME.size(); iLambda1++)
  {
    //if(iLambda1 > 1e3 ) break;
    for(unsigned int iLambda2 = iLambda1+1; iLambda2 < p_star_cuts_vector_ME.size(); iLambda2++)
    {
      //if(iLambda2 > 1e3 ) break;
      
      //re-weight ME      
      float delta_eta = fabs( L_cuts_vector_ME.at(iLambda1).Eta() - L_cuts_vector_ME.at(iLambda2).Eta() );
      float delta_phi = fabs( L_cuts_vector_ME.at(iLambda1).Phi() - L_cuts_vector_ME.at(iLambda2).Phi() );
      
      L0_L0_delta_eta_vs_delta_phi_ME_cuts_hist->Fill(delta_eta, delta_phi);

    }   
  
  }
  
  for(unsigned int iLambda1 = 0; iLambda1 < p_star_cuts_vector_ME.size(); iLambda1++)
  {
    //if(iLambda1 > 1e3 ) break;
    for(unsigned int iLambda2 = iLambda1+1; iLambda2 < p_star_cuts_vector_ME.size(); iLambda2++)
    {
      //if(iLambda2 > 1e3 ) break;
      double theta_star = p_star_cuts_vector_ME.at(iLambda1).Angle(p_star_cuts_vector_ME.at(iLambda2));
/*      
      L0_L0_cosThetaProdPlane_ME_cuts->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[L_pT_bin_cuts_vector_ME.at(iLambda1)][L_pT_bin_cuts_vector_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[L_eta_bin_cuts_vector_ME.at(iLambda1)][L_eta_bin_cuts_vector_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star));
*/      
      
      //re-weight ME      
      float delta_eta = fabs( L_cuts_vector_ME.at(iLambda1).Eta() - L_cuts_vector_ME.at(iLambda2).Eta() );
      float delta_phi = fabs( L_cuts_vector_ME.at(iLambda1).Phi() - L_cuts_vector_ME.at(iLambda2).Phi() );
      
      int delta_eta_bin = L0_L0_delta_eta_vs_delta_phi_cuts_hist->GetXaxis()->FindBin(delta_eta);
      int delta_phi_bin = L0_L0_delta_eta_vs_delta_phi_cuts_hist->GetYaxis()->FindBin(delta_phi);
          
      float weight = 0;
      
      if( L0_L0_delta_eta_vs_delta_phi_ME_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin) != 0 )
      {
        weight = L0_L0_delta_eta_vs_delta_phi_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin)/L0_L0_delta_eta_vs_delta_phi_ME_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin);
      }      
      
      L0_L0_cosThetaProdPlane_ME_weight_cuts->Fill(TMath::Cos(theta_star), weight);          
      L0_L0_cosThetaProdPlane_ME_weight_pT_cuts_hist[L_pT_bin_cuts_vector_ME.at(iLambda1)][L_pT_bin_cuts_vector_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star), weight);
      L0_L0_cosThetaProdPlane_ME_weight_eta_cuts_hist[L_eta_bin_cuts_vector_ME.at(iLambda1)][L_eta_bin_cuts_vector_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star), weight);
    }   
  
  }    

  //---------------------------------------------------------------------------------------------------------------------------

  //L0bar-L0bar
  for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pBar_star_cuts_vector_ME.size(); iLambdaBar1++)
  {
    //if(iLambdaBar1 > 1e3 ) break;
    for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < pBar_star_cuts_vector_ME.size(); iLambdaBar2++)
    {     
      //if(iLambdaBar2 > 1e3 ) break;
      //re-weight ME      
      float delta_eta = fabs( Lbar_cuts_vector_ME.at(iLambdaBar1).Eta() - Lbar_cuts_vector_ME.at(iLambdaBar2).Eta() );
      float delta_phi = fabs( Lbar_cuts_vector_ME.at(iLambdaBar1).Phi() - Lbar_cuts_vector_ME.at(iLambdaBar2).Phi() );
      
      L0bar_L0bar_delta_eta_vs_delta_phi_ME_cuts_hist->Fill(delta_eta, delta_phi);

    }   
  
  }  
  
  for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pBar_star_cuts_vector_ME.size(); iLambdaBar1++)
  {
    //if(iLambdaBar1 > 1e3 ) break;
    
    for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < pBar_star_cuts_vector_ME.size(); iLambdaBar2++)
    {
      //if(iLambdaBar2 > 1e3 ) break;
      
      double theta_star = pBar_star_cuts_vector_ME.at(iLambdaBar1).Angle(pBar_star_cuts_vector_ME.at(iLambdaBar2));
/*      
      L0bar_L0bar_cosThetaProdPlane_ME_cuts->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[Lbar_pT_bin_cuts_vector_ME.at(iLambdaBar1)][Lbar_pT_bin_cuts_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[Lbar_eta_bin_cuts_vector_ME.at(iLambdaBar1)][Lbar_eta_bin_cuts_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
*/      
      //re-weight ME      
      float delta_eta = fabs( Lbar_cuts_vector_ME.at(iLambdaBar1).Eta() - Lbar_cuts_vector_ME.at(iLambdaBar2).Eta() );
      float delta_phi = fabs( Lbar_cuts_vector_ME.at(iLambdaBar1).Phi() - Lbar_cuts_vector_ME.at(iLambdaBar2).Phi() );
      
      int delta_eta_bin = L0bar_L0bar_delta_eta_vs_delta_phi_cuts_hist->GetXaxis()->FindBin(delta_eta);
      int delta_phi_bin = L0bar_L0bar_delta_eta_vs_delta_phi_cuts_hist->GetYaxis()->FindBin(delta_phi);
          
      float weight = 0;
      
      if( L0bar_L0bar_delta_eta_vs_delta_phi_ME_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin) != 0 )
      {
        weight = L0bar_L0bar_delta_eta_vs_delta_phi_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin)/L0bar_L0bar_delta_eta_vs_delta_phi_ME_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin);
      }       
      
      L0bar_L0bar_cosThetaProdPlane_ME_weight_cuts->Fill(TMath::Cos(theta_star), weight);          
      L0bar_L0bar_cosThetaProdPlane_ME_weight_pT_cuts_hist[Lbar_pT_bin_cuts_vector_ME.at(iLambdaBar1)][Lbar_pT_bin_cuts_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star), weight);
      L0bar_L0bar_cosThetaProdPlane_ME_weight_eta_cuts_hist[Lbar_eta_bin_cuts_vector_ME.at(iLambdaBar1)][Lbar_eta_bin_cuts_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star), weight);  
    }   
  
  }  
  //_____________________________________________________________________________________________________________________________________________________________________
  
  cout<<"Start pair ME"<<endl;
  
  //mixed event with pi p pairs
  //US paired with US

  //L-Lbar before cuts
  for(unsigned int iLambda = 0; iLambda < p_star_vector_US_ME.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_US_ME.size(); iLambdaBar++)
    {      
      //re-weight ME      
      float delta_eta = fabs( L_vector_US_ME.at(iLambda).Eta() - Lbar_vector_US_ME.at(iLambdaBar).Eta() );
      float delta_phi = fabs( L_vector_US_ME.at(iLambda).Phi() - Lbar_vector_US_ME.at(iLambdaBar).Phi() );
      
      L0_L0bar_delta_eta_vs_delta_phi_US_ME_hist->Fill(delta_eta, delta_phi);
      
    }   
  
  }
  
  for(unsigned int iLambda = 0; iLambda < p_star_vector_US_ME.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_US_ME.size(); iLambdaBar++)
    {
      double theta_star = p_star_vector_US_ME.at(iLambda).Angle(pBar_star_vector_US_ME.at(iLambdaBar));
/*           
      L0_L0bar_cosThetaProdPlane_US_ME->Fill(TMath::Cos(theta_star));          
      L0_L0bar_cosThetaProdPlane_US_ME_pT_hist[L_pT_bin_vector_US_ME.at(iLambda)][Lbar_pT_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_US_ME_eta_hist[L_eta_bin_vector_US_ME.at(iLambda)][Lbar_eta_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
*/      
      //re-weight ME      
      float delta_eta = fabs( L_vector_US_ME.at(iLambda).Eta() - Lbar_vector_US_ME.at(iLambdaBar).Eta() );
      float delta_phi = fabs( L_vector_US_ME.at(iLambda).Phi() - Lbar_vector_US_ME.at(iLambdaBar).Phi() );
      
      int delta_eta_bin = L0_L0bar_delta_eta_vs_delta_phi_US_hist->GetXaxis()->FindBin(delta_eta);
      int delta_phi_bin = L0_L0bar_delta_eta_vs_delta_phi_US_hist->GetYaxis()->FindBin(delta_phi);
          
      float weight = 0;
      
      if( L0_L0bar_delta_eta_vs_delta_phi_US_ME_hist->GetBinContent(delta_eta_bin, delta_phi_bin) != 0 )
      {
        weight = L0_L0bar_delta_eta_vs_delta_phi_US_hist->GetBinContent(delta_eta_bin, delta_phi_bin)/L0_L0bar_delta_eta_vs_delta_phi_US_ME_hist->GetBinContent(delta_eta_bin, delta_phi_bin);
      }     
      
      L0_L0bar_cosThetaProdPlane_US_ME_weight->Fill(TMath::Cos(theta_star), weight);          
      L0_L0bar_cosThetaProdPlane_US_ME_weight_pT_hist[L_pT_bin_vector_US_ME.at(iLambda)][Lbar_pT_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star), weight);
      L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[L_eta_bin_vector_US_ME.at(iLambda)][Lbar_eta_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star), weight);
      
    }   
  
  }

  //-----------------------------------------------------------------------------------------------------------

  //L-Lbar after cuts
  for(unsigned int iLambda = 0; iLambda < p_star_vector_US_ME_cuts.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_US_ME_cuts.size(); iLambdaBar++)
    {      
      //re-weight ME      
      float delta_eta = fabs( L_vector_US_ME_cuts.at(iLambda).Eta() - Lbar_vector_US_ME_cuts.at(iLambdaBar).Eta() );
      float delta_phi = fabs( L_vector_US_ME_cuts.at(iLambda).Phi() - Lbar_vector_US_ME_cuts.at(iLambdaBar).Phi() );
   
      L0_L0bar_delta_eta_vs_delta_phi_US_ME_cuts_hist->Fill(delta_eta, delta_phi);
    }   
  
  }
  
  for(unsigned int iLambda = 0; iLambda < p_star_vector_US_ME_cuts.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_US_ME_cuts.size(); iLambdaBar++)
    {
      double theta_star = p_star_vector_US_ME_cuts.at(iLambda).Angle(pBar_star_vector_US_ME_cuts.at(iLambdaBar));
/*
      L0_L0bar_cosThetaProdPlane_US_ME_cuts->Fill(TMath::Cos(theta_star));          
      L0_L0bar_cosThetaProdPlane_US_ME_pT_cuts_hist[L_pT_bin_vector_US_ME_cuts.at(iLambda)][Lbar_pT_bin_vector_US_ME_cuts.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_US_ME_eta_cuts_hist[L_eta_bin_vector_US_ME_cuts.at(iLambda)][Lbar_eta_bin_vector_US_ME_cuts.at(iLambdaBar)]->Fill(TMath::Cos(theta_star)); 
*/      
      //re-weight ME      
      float delta_eta = fabs( L_vector_US_ME_cuts.at(iLambda).Eta() - Lbar_vector_US_ME_cuts.at(iLambdaBar).Eta() );
      float delta_phi = fabs( L_vector_US_ME_cuts.at(iLambda).Phi() - Lbar_vector_US_ME_cuts.at(iLambdaBar).Phi() );
      
      int delta_eta_bin = L0_L0bar_delta_eta_vs_delta_phi_US_cuts_hist->GetXaxis()->FindBin(delta_eta);
      int delta_phi_bin = L0_L0bar_delta_eta_vs_delta_phi_US_cuts_hist->GetYaxis()->FindBin(delta_phi);
          
      float weight = 0;
      
      if( L0_L0bar_delta_eta_vs_delta_phi_US_ME_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin) != 0 )
      {
        weight = L0_L0bar_delta_eta_vs_delta_phi_US_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin)/L0_L0bar_delta_eta_vs_delta_phi_US_ME_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin);
      }       
      
      L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts->Fill(TMath::Cos(theta_star), weight);          
      L0_L0bar_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[L_pT_bin_vector_US_ME_cuts.at(iLambda)][Lbar_pT_bin_vector_US_ME_cuts.at(iLambdaBar)]->Fill(TMath::Cos(theta_star), weight);
      L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[L_eta_bin_vector_US_ME_cuts.at(iLambda)][Lbar_eta_bin_vector_US_ME_cuts.at(iLambdaBar)]->Fill(TMath::Cos(theta_star), weight);     
    
    }   
  
  }
  
  //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  //L0-L0 before cuts
  for(unsigned int iLambda1 = 0; iLambda1 < p_star_vector_US_ME.size(); iLambda1++)
  {
    for(unsigned int iLambda2 = iLambda1+1; iLambda2 < p_star_vector_US_ME.size(); iLambda2++)
    {      
      //re-weight ME      
      float delta_eta = fabs( L_vector_US_ME.at(iLambda1).Eta() - L_vector_US_ME.at(iLambda2).Eta() );
      float delta_phi = fabs( L_vector_US_ME.at(iLambda1).Phi() - L_vector_US_ME.at(iLambda2).Phi() );

      L0_L0_delta_eta_vs_delta_phi_US_ME_hist->Fill(delta_eta, delta_phi);
    }   
  
  }
  
  for(unsigned int iLambda1 = 0; iLambda1 < p_star_vector_US_ME.size(); iLambda1++)
  {
    for(unsigned int iLambda2 = iLambda1+1; iLambda2 < p_star_vector_US_ME.size(); iLambda2++)
    {
      double theta_star = p_star_vector_US_ME.at(iLambda1).Angle(p_star_vector_US_ME.at(iLambda2));
 /*     
      L0_L0_cosThetaProdPlane_US_ME->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_US_ME_pT_hist[L_pT_bin_vector_US_ME.at(iLambda1)][L_pT_bin_vector_US_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_US_ME_eta_hist[L_eta_bin_vector_US_ME.at(iLambda1)][L_eta_bin_vector_US_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star));
 */     
      //re-weight ME      
      float delta_eta = fabs( L_vector_US_ME.at(iLambda1).Eta() - L_vector_US_ME.at(iLambda2).Eta() );
      float delta_phi = fabs( L_vector_US_ME.at(iLambda1).Phi() - L_vector_US_ME.at(iLambda2).Phi() );
      
      int delta_eta_bin = L0_L0_delta_eta_vs_delta_phi_US_hist->GetXaxis()->FindBin(delta_eta);
      int delta_phi_bin = L0_L0_delta_eta_vs_delta_phi_US_hist->GetYaxis()->FindBin(delta_phi);
          
      float weight = 0;
      
      if( L0_L0_delta_eta_vs_delta_phi_US_ME_hist->GetBinContent(delta_eta_bin, delta_phi_bin) != 0 )
      {
        weight = L0_L0_delta_eta_vs_delta_phi_US_hist->GetBinContent(delta_eta_bin, delta_phi_bin)/L0_L0_delta_eta_vs_delta_phi_US_ME_hist->GetBinContent(delta_eta_bin, delta_phi_bin);
      }     
      
      L0_L0_cosThetaProdPlane_US_ME_weight->Fill(TMath::Cos(theta_star), weight);          
      L0_L0_cosThetaProdPlane_US_ME_weight_pT_hist[L_pT_bin_vector_US_ME.at(iLambda1)][L_pT_bin_vector_US_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star), weight);
      L0_L0_cosThetaProdPlane_US_ME_weight_eta_hist[L_eta_bin_vector_US_ME.at(iLambda1)][L_eta_bin_vector_US_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star), weight);

    }   
  
  }
  
  //-----------------------------------------------------------------------------------------------------------
    
  //L0-L0 after cuts
  for(unsigned int iLambda1 = 0; iLambda1 < p_star_vector_US_ME_cuts.size(); iLambda1++)
  {
    for(unsigned int iLambda2 = iLambda1+1; iLambda2 < p_star_vector_US_ME_cuts.size(); iLambda2++)
    {     
      //re-weight ME      
      float delta_eta = fabs( L_vector_US_ME_cuts.at(iLambda1).Eta() - L_vector_US_ME_cuts.at(iLambda2).Eta() );
      float delta_phi = fabs( L_vector_US_ME_cuts.at(iLambda1).Phi() - L_vector_US_ME_cuts.at(iLambda2).Phi() );
      
      L0_L0_delta_eta_vs_delta_phi_US_ME_cuts_hist->Fill(delta_eta, delta_phi);
     
    }   
  
  }
  
  for(unsigned int iLambda1 = 0; iLambda1 < p_star_vector_US_ME_cuts.size(); iLambda1++)
  {
    for(unsigned int iLambda2 = iLambda1+1; iLambda2 < p_star_vector_US_ME_cuts.size(); iLambda2++)
    {
      double theta_star = p_star_vector_US_ME_cuts.at(iLambda1).Angle(p_star_vector_US_ME_cuts.at(iLambda2));
/*           
      L0_L0_cosThetaProdPlane_US_ME_cuts->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_US_ME_pT_cuts_hist[L_pT_bin_vector_US_ME_cuts.at(iLambda1)][L_pT_bin_vector_US_ME_cuts.at(iLambda2)]->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_US_ME_eta_cuts_hist[L_eta_bin_vector_US_ME_cuts.at(iLambda1)][L_eta_bin_vector_US_ME_cuts.at(iLambda2)]->Fill(TMath::Cos(theta_star));
*/      
      //re-weight ME      
      float delta_eta = fabs( L_vector_US_ME_cuts.at(iLambda1).Eta() - L_vector_US_ME_cuts.at(iLambda2).Eta() );
      float delta_phi = fabs( L_vector_US_ME_cuts.at(iLambda1).Phi() - L_vector_US_ME_cuts.at(iLambda2).Phi() );
      
      int delta_eta_bin = L0_L0_delta_eta_vs_delta_phi_US_cuts_hist->GetXaxis()->FindBin(delta_eta);
      int delta_phi_bin = L0_L0_delta_eta_vs_delta_phi_US_cuts_hist->GetYaxis()->FindBin(delta_phi);
          
      float weight = 0;
      
      if( L0_L0_delta_eta_vs_delta_phi_US_ME_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin) != 0 )
      {
        weight = L0_L0_delta_eta_vs_delta_phi_US_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin)/L0_L0_delta_eta_vs_delta_phi_US_ME_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin);
      }     
      
      L0_L0_cosThetaProdPlane_US_ME_weight_cuts->Fill(TMath::Cos(theta_star), weight);          
      L0_L0_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[L_pT_bin_vector_US_ME_cuts.at(iLambda1)][L_pT_bin_vector_US_ME_cuts.at(iLambda2)]->Fill(TMath::Cos(theta_star), weight);
      L0_L0_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[L_eta_bin_vector_US_ME_cuts.at(iLambda1)][L_eta_bin_vector_US_ME_cuts.at(iLambda2)]->Fill(TMath::Cos(theta_star), weight);
      
    }   
  
  }
  
  //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  //L0bar-L0bar before cuts
  for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pBar_star_vector_US_ME.size(); iLambdaBar1++)
  {
    for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < pBar_star_vector_US_ME.size(); iLambdaBar2++)
    {     
      //re-weight ME      
      float delta_eta = fabs( Lbar_vector_US_ME.at(iLambdaBar1).Eta() - Lbar_vector_US_ME.at(iLambdaBar2).Eta() );
      float delta_phi = fabs( Lbar_vector_US_ME.at(iLambdaBar1).Phi() - Lbar_vector_US_ME.at(iLambdaBar2).Phi() );

      L0bar_L0bar_delta_eta_vs_delta_phi_US_ME_hist->Fill(delta_eta, delta_phi); 
    }   
  
  }
  
  for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pBar_star_vector_US_ME.size(); iLambdaBar1++)
  {
    for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < pBar_star_vector_US_ME.size(); iLambdaBar2++)
    {
      double theta_star = pBar_star_vector_US_ME.at(iLambdaBar1).Angle(pBar_star_vector_US_ME.at(iLambdaBar2));
/*      
      L0bar_L0bar_cosThetaProdPlane_US_ME->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_US_ME_pT_hist[Lbar_pT_bin_vector_US_ME.at(iLambdaBar1)][Lbar_pT_bin_vector_US_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_US_ME_eta_hist[Lbar_eta_bin_vector_US_ME.at(iLambdaBar1)][Lbar_eta_bin_vector_US_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
*/      
      //re-weight ME      
      float delta_eta = fabs( Lbar_vector_US_ME.at(iLambdaBar1).Eta() - Lbar_vector_US_ME.at(iLambdaBar2).Eta() );
      float delta_phi = fabs( Lbar_vector_US_ME.at(iLambdaBar1).Phi() - Lbar_vector_US_ME.at(iLambdaBar2).Phi() );
      
      int delta_eta_bin = L0bar_L0bar_delta_eta_vs_delta_phi_US_hist->GetXaxis()->FindBin(delta_eta);
      int delta_phi_bin = L0bar_L0bar_delta_eta_vs_delta_phi_US_hist->GetYaxis()->FindBin(delta_phi);
          
      float weight = 0;
      
      if( L0bar_L0bar_delta_eta_vs_delta_phi_US_ME_hist->GetBinContent(delta_eta_bin, delta_phi_bin) != 0 )
      {
        weight = L0bar_L0bar_delta_eta_vs_delta_phi_US_hist->GetBinContent(delta_eta_bin, delta_phi_bin)/L0bar_L0bar_delta_eta_vs_delta_phi_US_ME_hist->GetBinContent(delta_eta_bin, delta_phi_bin);
      } 
      
      L0bar_L0bar_cosThetaProdPlane_US_ME_weight->Fill(TMath::Cos(theta_star), weight);          
      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT_hist[Lbar_pT_bin_vector_US_ME.at(iLambdaBar1)][Lbar_pT_bin_vector_US_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star), weight);
      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[Lbar_eta_bin_vector_US_ME.at(iLambdaBar1)][Lbar_eta_bin_vector_US_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star), weight); 
      
    }   
  
  }
  
  //-----------------------------------------------------------------------------------------------------------
  
  //L0bar-L0bar after cuts
  for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pBar_star_vector_US_ME_cuts.size(); iLambdaBar1++)
  {
    for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < pBar_star_vector_US_ME_cuts.size(); iLambdaBar2++)
    {       
      //re-weight ME      
      float delta_eta = fabs( Lbar_vector_US_ME_cuts.at(iLambdaBar1).Eta() - Lbar_vector_US_ME_cuts.at(iLambdaBar2).Eta() );
      float delta_phi = fabs( Lbar_vector_US_ME_cuts.at(iLambdaBar1).Phi() - Lbar_vector_US_ME_cuts.at(iLambdaBar2).Phi() );
      
      L0bar_L0bar_delta_eta_vs_delta_phi_US_ME_cuts_hist->Fill(delta_eta, delta_phi);       
    }   
  
  }
  
  for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pBar_star_vector_US_ME_cuts.size(); iLambdaBar1++)
  {
    for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < pBar_star_vector_US_ME_cuts.size(); iLambdaBar2++)
    {
      double theta_star = pBar_star_vector_US_ME_cuts.at(iLambdaBar1).Angle(pBar_star_vector_US_ME_cuts.at(iLambdaBar2));
/*
      L0bar_L0bar_cosThetaProdPlane_US_ME_cuts->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_US_ME_pT_cuts_hist[Lbar_pT_bin_vector_US_ME_cuts.at(iLambdaBar1)][Lbar_pT_bin_vector_US_ME_cuts.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_US_ME_eta_cuts_hist[Lbar_eta_bin_vector_US_ME_cuts.at(iLambdaBar1)][Lbar_eta_bin_vector_US_ME_cuts.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));  
*/      
      //re-weight ME      
      float delta_eta = fabs( Lbar_vector_US_ME_cuts.at(iLambdaBar1).Eta() - Lbar_vector_US_ME_cuts.at(iLambdaBar2).Eta() );
      float delta_phi = fabs( Lbar_vector_US_ME_cuts.at(iLambdaBar1).Phi() - Lbar_vector_US_ME_cuts.at(iLambdaBar2).Phi() );
      
      int delta_eta_bin = L0bar_L0bar_delta_eta_vs_delta_phi_US_cuts_hist->GetXaxis()->FindBin(delta_eta);
      int delta_phi_bin = L0bar_L0bar_delta_eta_vs_delta_phi_US_cuts_hist->GetYaxis()->FindBin(delta_phi);
          
      float weight = 0;
      
      if( L0bar_L0bar_delta_eta_vs_delta_phi_US_ME_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin) != 0 )
      {
        weight = L0bar_L0bar_delta_eta_vs_delta_phi_US_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin)/L0bar_L0bar_delta_eta_vs_delta_phi_US_ME_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin);
      }      
      
      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_cuts->Fill(TMath::Cos(theta_star), weight);          
      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[Lbar_pT_bin_vector_US_ME_cuts.at(iLambdaBar1)][Lbar_pT_bin_vector_US_ME_cuts.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star), weight);
      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[Lbar_eta_bin_vector_US_ME_cuts.at(iLambdaBar1)][Lbar_eta_bin_vector_US_ME_cuts.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star), weight); 
      
    }   
  
  }
  //_____________________________________________________________________________________________________________________________________________________________________
  
  //US paired with LS
  //L-Lbar before cuts
  //L - US, Lbar - LS   
  for( unsigned int iLambda = 0; iLambda < p_star_vector_US_ME.size(); iLambda++)
  {
    for( unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_LS_ME.size(); iLambdaBar++)
    {     
      //re-weight ME      
      float delta_eta = fabs( L_vector_US_ME.at(iLambda).Eta() - Lbar_vector_LS_ME.at(iLambdaBar).Eta() );
      float delta_phi = fabs( L_vector_US_ME.at(iLambda).Phi() - Lbar_vector_LS_ME.at(iLambdaBar).Phi() );
      
      L0_L0bar_delta_eta_vs_delta_phi_US_LS_ME_hist->Fill(delta_eta, delta_phi);
    }
  
  } 
  
  //L - LS, Lbar - US  
  for( unsigned int iLambda = 0; iLambda < p_star_vector_LS_ME.size(); iLambda++)
  {
    for( unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_US_ME.size(); iLambdaBar++)
    {      
      //re-weight ME      
      float delta_eta = fabs( L_vector_LS_ME.at(iLambda).Eta() - Lbar_vector_US_ME.at(iLambdaBar).Eta() );
      float delta_phi = fabs( L_vector_LS_ME.at(iLambda).Phi() - Lbar_vector_US_ME.at(iLambdaBar).Phi() );

      L0_L0bar_delta_eta_vs_delta_phi_US_LS_ME_hist->Fill(delta_eta, delta_phi);
    }
  
  }
  
  
  for( unsigned int iLambda = 0; iLambda < p_star_vector_US_ME.size(); iLambda++)
  {
    for( unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_LS_ME.size(); iLambdaBar++)
    {     
      float theta_star = p_star_vector_US_ME.at(iLambda).Angle(pBar_star_vector_LS_ME.at(iLambdaBar));
 /*     
      L0_L0bar_cosThetaProdPlane_US_LS_ME->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_US_LS_ME_pT_hist[L_pT_bin_vector_US_ME.at(iLambda)][Lbar_pT_bin_vector_LS_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_US_LS_ME_eta_hist[L_eta_bin_vector_US_ME.at(iLambda)][Lbar_eta_bin_vector_LS_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
 */     
      //re-weight ME      
      float delta_eta = fabs( L_vector_US_ME.at(iLambda).Eta() - Lbar_vector_LS_ME.at(iLambdaBar).Eta() );
      float delta_phi = fabs( L_vector_US_ME.at(iLambda).Phi() - Lbar_vector_LS_ME.at(iLambdaBar).Phi() );
      
      int delta_eta_bin = L0_L0bar_delta_eta_vs_delta_phi_US_LS_hist->GetXaxis()->FindBin(delta_eta);
      int delta_phi_bin = L0_L0bar_delta_eta_vs_delta_phi_US_LS_hist->GetYaxis()->FindBin(delta_phi);
          
      float weight = 0;
      
      if( L0_L0bar_delta_eta_vs_delta_phi_US_LS_ME_hist->GetBinContent(delta_eta_bin, delta_phi_bin) != 0 )
      {
        weight = L0_L0bar_delta_eta_vs_delta_phi_US_LS_hist->GetBinContent(delta_eta_bin, delta_phi_bin)/L0_L0bar_delta_eta_vs_delta_phi_US_LS_ME_hist->GetBinContent(delta_eta_bin, delta_phi_bin);
      }
      
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight->Fill(TMath::Cos(theta_star), weight);          
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_hist[L_pT_bin_vector_US_ME.at(iLambda)][Lbar_pT_bin_vector_LS_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star), weight);
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_hist[L_eta_bin_vector_US_ME.at(iLambda)][Lbar_eta_bin_vector_LS_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star), weight);
    
    }
  
  } 
  
  //L - LS, Lbar - US  
  for( unsigned int iLambda = 0; iLambda < p_star_vector_LS_ME.size(); iLambda++)
  {
    for( unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_US_ME.size(); iLambdaBar++)
    {    
      float theta_star = p_star_vector_LS_ME.at(iLambda).Angle(pBar_star_vector_US_ME.at(iLambdaBar));
/*      
      L0_L0bar_cosThetaProdPlane_US_LS_ME->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_US_LS_ME_pT_hist[L_pT_bin_vector_LS_ME.at(iLambda)][Lbar_pT_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_US_LS_ME_eta_hist[L_eta_bin_vector_LS_ME.at(iLambda)][Lbar_eta_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
*/      
      //re-weight ME      
      float delta_eta = fabs( L_vector_LS_ME.at(iLambda).Eta() - Lbar_vector_US_ME.at(iLambdaBar).Eta() );
      float delta_phi = fabs( L_vector_LS_ME.at(iLambda).Phi() - Lbar_vector_US_ME.at(iLambdaBar).Phi() );
      
      int delta_eta_bin = L0_L0bar_delta_eta_vs_delta_phi_US_LS_hist->GetXaxis()->FindBin(delta_eta);
      int delta_phi_bin = L0_L0bar_delta_eta_vs_delta_phi_US_LS_hist->GetYaxis()->FindBin(delta_phi);
          
      float weight = 0;
      
      if( L0_L0bar_delta_eta_vs_delta_phi_US_LS_ME_hist->GetBinContent(delta_eta_bin, delta_phi_bin) != 0 )
      {
        weight = L0_L0bar_delta_eta_vs_delta_phi_US_LS_hist->GetBinContent(delta_eta_bin, delta_phi_bin)/L0_L0bar_delta_eta_vs_delta_phi_US_LS_ME_hist->GetBinContent(delta_eta_bin, delta_phi_bin);
      }
      
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight->Fill(TMath::Cos(theta_star), weight);          
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_hist[L_pT_bin_vector_LS_ME.at(iLambda)][Lbar_pT_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star), weight);
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_hist[L_eta_bin_vector_LS_ME.at(iLambda)][Lbar_eta_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star), weight);
    
    }
  
  }
  
  //-----------------------------------------------------------------------------------------------------------
    
  //L-Lbar after cuts
  //L - US, Lbar - LS
  for( unsigned int iLambda = 0; iLambda < p_star_vector_US_ME_cuts.size(); iLambda++)
  {
    for( unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_LS_ME_cuts.size(); iLambdaBar++)
    {      
      //re-weight ME      
      float delta_eta = fabs( L_vector_US_ME_cuts.at(iLambda).Eta() - Lbar_vector_LS_ME_cuts.at(iLambdaBar).Eta() );
      float delta_phi = fabs( L_vector_US_ME_cuts.at(iLambda).Phi() - Lbar_vector_LS_ME_cuts.at(iLambdaBar).Phi() );
      
      L0_L0bar_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist->Fill(delta_eta, delta_phi); 
    }
  
  } 
  
  //L - LS, Lbar - US  
  for( unsigned int iLambda = 0; iLambda < p_star_vector_LS_ME_cuts.size(); iLambda++)
  {
    for( unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_US_ME_cuts.size(); iLambdaBar++)
    {      
      //re-weight ME      
      float delta_eta = fabs( L_vector_LS_ME_cuts.at(iLambda).Eta() - Lbar_vector_US_ME_cuts.at(iLambdaBar).Eta() );
      float delta_phi = fabs( L_vector_LS_ME_cuts.at(iLambda).Phi() - Lbar_vector_US_ME_cuts.at(iLambdaBar).Phi() );
      
      L0_L0bar_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist->Fill(delta_eta, delta_phi);    
    }
  
  }
  
  for( unsigned int iLambda = 0; iLambda < p_star_vector_US_ME_cuts.size(); iLambda++)
  {
    for( unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_LS_ME_cuts.size(); iLambdaBar++)
    {     
      float theta_star = p_star_vector_US_ME_cuts.at(iLambda).Angle(pBar_star_vector_LS_ME_cuts.at(iLambdaBar));
/*
      L0_L0bar_cosThetaProdPlane_US_LS_ME_cuts->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[L_pT_bin_vector_US_ME_cuts.at(iLambda)][Lbar_pT_bin_vector_LS_ME_cuts.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[L_eta_bin_vector_US_ME_cuts.at(iLambda)][Lbar_eta_bin_vector_LS_ME_cuts.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
*/      
      //re-weight ME      
      float delta_eta = fabs( L_vector_US_ME_cuts.at(iLambda).Eta() - Lbar_vector_LS_ME_cuts.at(iLambdaBar).Eta() );
      float delta_phi = fabs( L_vector_US_ME_cuts.at(iLambda).Phi() - Lbar_vector_LS_ME_cuts.at(iLambdaBar).Phi() );
      
      int delta_eta_bin = L0_L0bar_delta_eta_vs_delta_phi_US_LS_cuts_hist->GetXaxis()->FindBin(delta_eta);
      int delta_phi_bin = L0_L0bar_delta_eta_vs_delta_phi_US_LS_cuts_hist->GetYaxis()->FindBin(delta_phi);
          
      float weight = 0;
      
      if( L0_L0bar_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin) )
      {
        weight = L0_L0bar_delta_eta_vs_delta_phi_US_LS_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin)/L0_L0bar_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin);
      }
      
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts->Fill(TMath::Cos(theta_star), weight);          
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[L_pT_bin_vector_US_ME_cuts.at(iLambda)][Lbar_pT_bin_vector_LS_ME_cuts.at(iLambdaBar)]->Fill(TMath::Cos(theta_star), weight);
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[L_eta_bin_vector_US_ME_cuts.at(iLambda)][Lbar_eta_bin_vector_LS_ME_cuts.at(iLambdaBar)]->Fill(TMath::Cos(theta_star), weight);         
 
    }
  
  } 
  
  //L - LS, Lbar - US  
  for( unsigned int iLambda = 0; iLambda < p_star_vector_LS_ME_cuts.size(); iLambda++)
  {
    for( unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_US_ME_cuts.size(); iLambdaBar++)
    {    
      float theta_star = p_star_vector_LS_ME_cuts.at(iLambda).Angle(pBar_star_vector_US_ME_cuts.at(iLambdaBar));
/*
      L0_L0bar_cosThetaProdPlane_US_LS_ME_cuts->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[L_pT_bin_vector_LS_ME_cuts.at(iLambda)][Lbar_pT_bin_vector_US_ME_cuts.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[L_eta_bin_vector_LS_ME_cuts.at(iLambda)][Lbar_eta_bin_vector_US_ME_cuts.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));  
*/      
      //re-weight ME      
      float delta_eta = fabs( L_vector_LS_ME_cuts.at(iLambda).Eta() - Lbar_vector_US_ME_cuts.at(iLambdaBar).Eta() );
      float delta_phi = fabs( L_vector_LS_ME_cuts.at(iLambda).Phi() - Lbar_vector_US_ME_cuts.at(iLambdaBar).Phi() );
      
      int delta_eta_bin = L0_L0bar_delta_eta_vs_delta_phi_US_LS_cuts_hist->GetXaxis()->FindBin(delta_eta);
      int delta_phi_bin = L0_L0bar_delta_eta_vs_delta_phi_US_LS_cuts_hist->GetYaxis()->FindBin(delta_phi);
          
      float weight = 0;
      
      if( L0_L0bar_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin) )
      {
        weight = L0_L0bar_delta_eta_vs_delta_phi_US_LS_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin)/L0_L0bar_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin);
      }
      
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts->Fill(TMath::Cos(theta_star), weight);          
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[L_pT_bin_vector_LS_ME_cuts.at(iLambda)][Lbar_pT_bin_vector_US_ME_cuts.at(iLambdaBar)]->Fill(TMath::Cos(theta_star), weight);
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[L_eta_bin_vector_LS_ME_cuts.at(iLambda)][Lbar_eta_bin_vector_US_ME_cuts.at(iLambdaBar)]->Fill(TMath::Cos(theta_star), weight);      
    
    }
  
  } 
  
  //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  //L-L before cuts
  for( unsigned int iLambda = 0; iLambda < p_star_vector_US_ME.size(); iLambda++ )
  {
    for( unsigned int iLambdaBck = 0; iLambdaBck < p_star_vector_LS_ME.size(); iLambdaBck++ )
    { 
      //re-weight ME      
      float delta_eta = fabs( L_vector_US_ME.at(iLambda).Eta() - L_vector_LS_ME.at(iLambdaBck).Eta() );
      float delta_phi = fabs( L_vector_US_ME.at(iLambda).Phi() - L_vector_LS_ME.at(iLambdaBck).Phi() );
          
      L0_L0_delta_eta_vs_delta_phi_US_LS_ME_hist->Fill(delta_eta, delta_phi);   
    }
  
  }
  
  
  int nFills_L_L_US_LS_ME = 0;  
  
  for( unsigned int iLambda = 0; iLambda < p_star_vector_US_ME.size(); iLambda++ )
  {
    for( unsigned int iLambdaBck = 0; iLambdaBck < p_star_vector_LS_ME.size(); iLambdaBck++ )
    { 
      //re-weight ME      
      float delta_eta = fabs( L_vector_US_ME.at(iLambda).Eta() - L_vector_LS_ME.at(iLambdaBck).Eta() );
      float delta_phi = fabs( L_vector_US_ME.at(iLambda).Phi() - L_vector_LS_ME.at(iLambdaBck).Phi() );
      
      int delta_eta_bin = L0_L0_delta_eta_vs_delta_phi_US_LS_hist->GetXaxis()->FindBin(delta_eta);
      int delta_phi_bin = L0_L0_delta_eta_vs_delta_phi_US_LS_hist->GetYaxis()->FindBin(delta_phi);
          
      float weight = 0;
      
      if( L0_L0_delta_eta_vs_delta_phi_US_LS_ME_hist->GetBinContent(delta_eta_bin, delta_phi_bin) != 0 )
      {
        weight = L0_L0_delta_eta_vs_delta_phi_US_LS_hist->GetBinContent(delta_eta_bin, delta_phi_bin)/L0_L0_delta_eta_vs_delta_phi_US_LS_ME_hist->GetBinContent(delta_eta_bin, delta_phi_bin);
      }
      

      //interate order of US and LS for L
      if( nFills_L_L_US_LS_ME % 2 == 0 )
      {
        float theta_star = p_star_vector_US_ME.at(iLambda).Angle(p_star_vector_LS_ME.at(iLambdaBck));
/*                  
        L0_L0_cosThetaProdPlane_US_LS_ME->Fill(TMath::Cos(theta_star));
        L0_L0_cosThetaProdPlane_US_LS_ME_pT_hist[L_pT_bin_vector_US_ME.at(iLambda)][L_pT_bin_vector_LS_ME.at(iLambdaBck)]->Fill(TMath::Cos(theta_star));
        L0_L0_cosThetaProdPlane_US_LS_ME_eta_hist[L_eta_bin_vector_US_ME.at(iLambda)][L_eta_bin_vector_LS_ME.at(iLambdaBck)]->Fill(TMath::Cos(theta_star));
*/        
        
        L0_L0_cosThetaProdPlane_US_LS_ME_weight->Fill(TMath::Cos(theta_star), weight);          
        L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT_hist[L_pT_bin_vector_US_ME.at(iLambda)][L_pT_bin_vector_LS_ME.at(iLambdaBck)]->Fill(TMath::Cos(theta_star), weight);
        L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta_hist[L_eta_bin_vector_US_ME.at(iLambda)][L_eta_bin_vector_LS_ME.at(iLambdaBck)]->Fill(TMath::Cos(theta_star), weight);
     
      
      }
      else
      {
        float theta_star = p_star_vector_LS_ME.at(iLambdaBck).Angle(p_star_vector_US_ME.at(iLambda));
/*                  
        L0_L0_cosThetaProdPlane_US_LS_ME->Fill(TMath::Cos(theta_star));
        L0_L0_cosThetaProdPlane_US_LS_ME_pT_hist[L_pT_bin_vector_LS_ME.at(iLambdaBck)][L_pT_bin_vector_US_ME.at(iLambda)]->Fill(TMath::Cos(theta_star));
        L0_L0_cosThetaProdPlane_US_LS_ME_eta_hist[L_eta_bin_vector_LS_ME.at(iLambdaBck)][L_eta_bin_vector_US_ME.at(iLambda)]->Fill(TMath::Cos(theta_star));
*/        
        
        L0_L0_cosThetaProdPlane_US_LS_ME_weight->Fill(TMath::Cos(theta_star), weight);          
        L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT_hist[L_pT_bin_vector_LS_ME.at(iLambdaBck)][L_pT_bin_vector_US_ME.at(iLambda)]->Fill(TMath::Cos(theta_star), weight);
        L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta_hist[L_eta_bin_vector_LS_ME.at(iLambdaBck)][L_eta_bin_vector_US_ME.at(iLambda)]->Fill(TMath::Cos(theta_star), weight);
        
      }//end else
      
      nFills_L_L_US_LS_ME++;
    
    }
  
  }
  
  //-----------------------------------------------------------------------------------------------------------
  
  //L-L after cuts
  for( unsigned int iLambda = 0; iLambda < p_star_vector_US_ME_cuts.size(); iLambda++ )
  {
    for( unsigned int iLambdaBck = 0; iLambdaBck < p_star_vector_LS_ME_cuts.size(); iLambdaBck++ )
    {
      //re-weight ME      
      float delta_eta = fabs( L_vector_US_ME_cuts.at(iLambda).Eta() - L_vector_LS_ME_cuts.at(iLambdaBck).Eta() );
      float delta_phi = fabs( L_vector_US_ME_cuts.at(iLambda).Phi() - L_vector_LS_ME_cuts.at(iLambdaBck).Phi() );
      
      L0_L0_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist->Fill(delta_eta, delta_phi); 

    }  
  }
  
  
  int nFills_L_L_US_LS_ME_cuts = 0;  
  
  for( unsigned int iLambda = 0; iLambda < p_star_vector_US_ME_cuts.size(); iLambda++ )
  {
    for( unsigned int iLambdaBck = 0; iLambdaBck < p_star_vector_LS_ME_cuts.size(); iLambdaBck++ )
    {
      //re-weight ME      
      float delta_eta = fabs( L_vector_US_ME_cuts.at(iLambda).Eta() - L_vector_LS_ME_cuts.at(iLambdaBck).Eta() );
      float delta_phi = fabs( L_vector_US_ME_cuts.at(iLambda).Phi() - L_vector_LS_ME_cuts.at(iLambdaBck).Phi() );
      
      int delta_eta_bin = L0_L0_delta_eta_vs_delta_phi_US_LS_cuts_hist->GetXaxis()->FindBin(delta_eta);
      int delta_phi_bin = L0_L0_delta_eta_vs_delta_phi_US_LS_cuts_hist->GetYaxis()->FindBin(delta_phi);
          
      float weight = 0;
      
      if( L0_L0_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin) != 0 )
      {
        weight = L0_L0_delta_eta_vs_delta_phi_US_LS_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin)/L0_L0_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin);
      }
       
            
      //interate order of US and LS for L
      if( nFills_L_L_US_LS_ME_cuts % 2 == 0 )
      {
        float theta_star = p_star_vector_US_ME_cuts.at(iLambda).Angle(p_star_vector_LS_ME_cuts.at(iLambdaBck));
/*
        L0_L0_cosThetaProdPlane_US_LS_ME_cuts->Fill(TMath::Cos(theta_star));
        L0_L0_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[L_pT_bin_vector_US_ME_cuts.at(iLambda)][L_pT_bin_vector_LS_ME_cuts.at(iLambdaBck)]->Fill(TMath::Cos(theta_star));
        L0_L0_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[L_eta_bin_vector_US_ME_cuts.at(iLambda)][L_eta_bin_vector_LS_ME_cuts.at(iLambdaBck)]->Fill(TMath::Cos(theta_star));
        
*/        
        L0_L0_cosThetaProdPlane_US_LS_ME_weight_cuts->Fill(TMath::Cos(theta_star), weight);          
        L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[L_pT_bin_vector_US_ME_cuts.at(iLambda)][L_pT_bin_vector_LS_ME_cuts.at(iLambdaBck)]->Fill(TMath::Cos(theta_star), weight);
        L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[L_eta_bin_vector_US_ME_cuts.at(iLambda)][L_eta_bin_vector_LS_ME_cuts.at(iLambdaBck)]->Fill(TMath::Cos(theta_star), weight);          
               
      
      }
      else
      {
        float theta_star = p_star_vector_LS_ME_cuts.at(iLambdaBck).Angle(p_star_vector_US_ME_cuts.at(iLambda));
/*
        L0_L0_cosThetaProdPlane_US_LS_ME_cuts->Fill(TMath::Cos(theta_star));
        L0_L0_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[L_pT_bin_vector_LS_ME_cuts.at(iLambdaBck)][L_pT_bin_vector_US_ME_cuts.at(iLambda)]->Fill(TMath::Cos(theta_star));
        L0_L0_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[L_eta_bin_vector_LS_ME_cuts.at(iLambdaBck)][L_eta_bin_vector_US_ME_cuts.at(iLambda)]->Fill(TMath::Cos(theta_star)); 
*/        
        L0_L0_cosThetaProdPlane_US_LS_ME_weight_cuts->Fill(TMath::Cos(theta_star), weight);          
        L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[L_pT_bin_vector_LS_ME_cuts.at(iLambdaBck)][L_pT_bin_vector_US_ME_cuts.at(iLambda)]->Fill(TMath::Cos(theta_star), weight);
        L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[L_eta_bin_vector_LS_ME_cuts.at(iLambdaBck)][L_eta_bin_vector_US_ME_cuts.at(iLambda)]->Fill(TMath::Cos(theta_star), weight);        
        
      }//end else
      
      nFills_L_L_US_LS_ME_cuts++;
    
    }
  
  }
  
  //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  //Lbar-Lbar before cuts
  for( unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_US_ME.size(); iLambdaBar++ )
  {
    for( unsigned int iLambdaBarBck = 0; iLambdaBarBck < pBar_star_vector_LS_ME.size(); iLambdaBarBck++ )
    {  
      //re-weight ME      
      float delta_eta = fabs( Lbar_vector_US_ME.at(iLambdaBar).Eta() - Lbar_vector_LS_ME.at(iLambdaBarBck).Eta() );
      float delta_phi = fabs( Lbar_vector_US_ME.at(iLambdaBar).Phi() - Lbar_vector_LS_ME.at(iLambdaBarBck).Phi() );
      
      L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_ME_hist->Fill(delta_eta, delta_phi);  
    
    }
  
  }
  
  
  int nFills_Lbar_Lbar_US_LS_ME = 0;
  
  for( unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_US_ME.size(); iLambdaBar++ )
  {
    for( unsigned int iLambdaBarBck = 0; iLambdaBarBck < pBar_star_vector_LS_ME.size(); iLambdaBarBck++ )
    {  
      //re-weight ME      
      float delta_eta = fabs( Lbar_vector_US_ME.at(iLambdaBar).Eta() - Lbar_vector_LS_ME.at(iLambdaBarBck).Eta() );
      float delta_phi = fabs( Lbar_vector_US_ME.at(iLambdaBar).Phi() - Lbar_vector_LS_ME.at(iLambdaBarBck).Phi() );
      
      int delta_eta_bin = L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_hist->GetXaxis()->FindBin(delta_eta);
      int delta_phi_bin = L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_hist->GetYaxis()->FindBin(delta_phi);
          
      float weight = 0;
      
      if( L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_ME_hist->GetBinContent(delta_eta_bin, delta_phi_bin) != 0 )
      {
        weight = L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_hist->GetBinContent(delta_eta_bin, delta_phi_bin)/L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_ME_hist->GetBinContent(delta_eta_bin, delta_phi_bin);
      }
      
       
      //interate order of US and LS for L
      if( nFills_Lbar_Lbar_US_LS_ME % 2 == 0 )
      {
        float theta_star = pBar_star_vector_US_ME.at(iLambdaBar).Angle(pBar_star_vector_LS_ME.at(iLambdaBarBck));
/*                  
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME->Fill(TMath::Cos(theta_star));
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_pT_hist[Lbar_pT_bin_vector_US_ME.at(iLambdaBar)][Lbar_pT_bin_vector_LS_ME.at(iLambdaBarBck)]->Fill(TMath::Cos(theta_star));
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_eta_hist[Lbar_eta_bin_vector_US_ME.at(iLambdaBar)][Lbar_eta_bin_vector_LS_ME.at(iLambdaBarBck)]->Fill(TMath::Cos(theta_star));
*/        
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight->Fill(TMath::Cos(theta_star), weight);          
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_hist[Lbar_pT_bin_vector_US_ME.at(iLambdaBar)][Lbar_pT_bin_vector_LS_ME.at(iLambdaBarBck)]->Fill(TMath::Cos(theta_star), weight);
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_hist[Lbar_eta_bin_vector_US_ME.at(iLambdaBar)][Lbar_eta_bin_vector_LS_ME.at(iLambdaBarBck)]->Fill(TMath::Cos(theta_star), weight);
        
      
      }
      else
      {
        float theta_star = pBar_star_vector_LS_ME.at(iLambdaBarBck).Angle(pBar_star_vector_US_ME.at(iLambdaBar));
/*                  
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME->Fill(TMath::Cos(theta_star));
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_pT_hist[Lbar_pT_bin_vector_LS_ME.at(iLambdaBarBck)][Lbar_pT_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_eta_hist[Lbar_eta_bin_vector_LS_ME.at(iLambdaBarBck)][Lbar_eta_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
*/        
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight->Fill(TMath::Cos(theta_star), weight);          
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_hist[Lbar_pT_bin_vector_LS_ME.at(iLambdaBarBck)][Lbar_pT_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star), weight);
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_hist[Lbar_eta_bin_vector_LS_ME.at(iLambdaBarBck)][Lbar_eta_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star), weight);
        
      }//end else
      
      nFills_Lbar_Lbar_US_LS_ME++;
    
    }
  
  }
  
  //-----------------------------------------------------------------------------------------------------------
  
  //Lbar-Lbar after cuts
  for( unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_US_ME_cuts.size(); iLambdaBar++ )
  {
    for( unsigned int iLambdaBarBck = 0; iLambdaBarBck < pBar_star_vector_LS_ME_cuts.size(); iLambdaBarBck++ )
    {     
      //re-weight ME      
      float delta_eta = fabs( Lbar_vector_US_ME_cuts.at(iLambdaBar).Eta() - Lbar_vector_LS_ME_cuts.at(iLambdaBarBck).Eta() );
      float delta_phi = fabs( Lbar_vector_US_ME_cuts.at(iLambdaBar).Phi() - Lbar_vector_LS_ME_cuts.at(iLambdaBarBck).Phi() );

      L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist->Fill(delta_eta, delta_phi);  
    }
  
  }
  
  
  int nFills_Lbar_Lbar_US_LS_ME_cuts = 0;
  
  for( unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_US_ME_cuts.size(); iLambdaBar++ )
  {
    for( unsigned int iLambdaBarBck = 0; iLambdaBarBck < pBar_star_vector_LS_ME_cuts.size(); iLambdaBarBck++ )
    {     
      //re-weight ME      
      float delta_eta = fabs( Lbar_vector_US_ME_cuts.at(iLambdaBar).Eta() - Lbar_vector_LS_ME_cuts.at(iLambdaBarBck).Eta() );
      float delta_phi = fabs( Lbar_vector_US_ME_cuts.at(iLambdaBar).Phi() - Lbar_vector_LS_ME_cuts.at(iLambdaBarBck).Phi() );
      
      int delta_eta_bin = L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_cuts_hist->GetXaxis()->FindBin(delta_eta);
      int delta_phi_bin = L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_cuts_hist->GetYaxis()->FindBin(delta_phi);
          
      float weight = 0;
      
      if( L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin) != 0 )
      {
        weight = L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin)/L0bar_L0bar_delta_eta_vs_delta_phi_US_LS_ME_cuts_hist->GetBinContent(delta_eta_bin, delta_phi_bin);
      }      
      
      //interate order of US and LS for L
      if( nFills_Lbar_Lbar_US_LS_ME_cuts % 2 == 0 )
      {
        float theta_star = pBar_star_vector_US_ME_cuts.at(iLambdaBar).Angle(pBar_star_vector_LS_ME_cuts.at(iLambdaBarBck));        
/*
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_cuts->Fill(TMath::Cos(theta_star));
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[Lbar_pT_bin_vector_US_ME_cuts.at(iLambdaBar)][Lbar_pT_bin_vector_LS_ME_cuts.at(iLambdaBarBck)]->Fill(TMath::Cos(theta_star));
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[Lbar_eta_bin_vector_US_ME_cuts.at(iLambdaBar)][Lbar_eta_bin_vector_LS_ME_cuts.at(iLambdaBarBck)]->Fill(TMath::Cos(theta_star));  
*/        
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts->Fill(TMath::Cos(theta_star), weight);          
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[Lbar_pT_bin_vector_US_ME_cuts.at(iLambdaBar)][Lbar_pT_bin_vector_LS_ME_cuts.at(iLambdaBarBck)]->Fill(TMath::Cos(theta_star), weight);
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[Lbar_eta_bin_vector_US_ME_cuts.at(iLambdaBar)][Lbar_eta_bin_vector_LS_ME_cuts.at(iLambdaBarBck)]->Fill(TMath::Cos(theta_star), weight);            
      
      }
      else
      {
        float theta_star = pBar_star_vector_LS_ME_cuts.at(iLambdaBarBck).Angle(pBar_star_vector_US_ME_cuts.at(iLambdaBar));
/*                  
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_cuts->Fill(TMath::Cos(theta_star));
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[Lbar_pT_bin_vector_LS_ME_cuts.at(iLambdaBarBck)][Lbar_pT_bin_vector_US_ME_cuts.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[Lbar_eta_bin_vector_LS_ME_cuts.at(iLambdaBarBck)][Lbar_eta_bin_vector_US_ME_cuts.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));       
*/        
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts->Fill(TMath::Cos(theta_star), weight);          
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[Lbar_pT_bin_vector_LS_ME_cuts.at(iLambdaBarBck)][Lbar_pT_bin_vector_US_ME_cuts.at(iLambdaBar)]->Fill(TMath::Cos(theta_star), weight);
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[Lbar_eta_bin_vector_LS_ME_cuts.at(iLambdaBarBck)][Lbar_eta_bin_vector_US_ME_cuts.at(iLambdaBar)]->Fill(TMath::Cos(theta_star), weight);   
        
      }//end else
      
      nFills_Lbar_Lbar_US_LS_ME_cuts++;
    
    }
  
  }
  
  //_________________________________________________________________________________________________________

  
  
  outFile->cd();
  outFile->Write();
  outFile->Close();

  
  //outFile->Close();

  // Done.
  return;
}
