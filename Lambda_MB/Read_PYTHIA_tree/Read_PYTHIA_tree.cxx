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
#include "TRandom3.h"
//#include "TString.h"
#include "TMath.h"


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

TLorentzVector smearMom(TLorentzVector const& b, TF1 const * const fMomResolution)
{ 
  float const pt = b.Perp();
  float const sPt = gRandom->Gaus(pt, pt * fMomResolution->Eval(pt));
  

  TLorentzVector sMom;
  sMom.SetXYZM(sPt * cos(b.Phi()), sPt * sin(b.Phi()), sPt * sinh(b.Eta()), b.M());
  return sMom;
}

TLorentzVector smearMomEtaPhi(TLorentzVector const& b, TF1 const * const fMomResolution, TF1 const * const fEtaResolution, TF1 const * const fPhiResolution)
{
  float const pt = b.Perp();
  float const sPt = gRandom->Gaus(pt, pt * fMomResolution->Eval(pt) ); //additional factor of 10 for testing of smearing effect
  
  float const eta = b.Eta();
  float const sEta = gRandom->Gaus(eta, fEtaResolution->Eval(eta) );
  
  float const phi = b.Phi();
  float const sPhi = gRandom->Gaus(phi, fPhiResolution->Eval(phi) );

  TLorentzVector sMom;
  sMom.SetXYZM(sPt * cos(sPhi), sPt * sin(sPhi), sPt * sinh(sEta), b.M());
  return sMom;
}

TLorentzVector smearPhi(TLorentzVector const& b, TF1 const * const fPhiResolution)
{
  
  float const phi = b.Phi();
  float const sPhi = gRandom->Gaus(phi, fPhiResolution->Eval(phi) );

  TLorentzVector sMom;
  //sMom.SetXYZM(sPt * cos(sPhi), sPt * sin(sPhi), sPt * sinh(sEta), b.M());
  sMom.SetPtEtaPhiM(b.Pt(), b.Eta(), sPhi, b.M());
  return sMom;
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



const double pi_mass_PDG = 0.13957039; //p mass on GeV/c^2 from latest PDG
const double p_mass_PDG = 0.93827208816; //p mass on GeV/c^2 from latest PDG
const double L_mass_PDG = 1.115683; //mass in GeV/c^2 from latest PDG

//mEnergy is collision energy - 200 or 510, 0 - for testing on small sample
//daughterSmearFlag = 0 - no smearing, 1 - daughter pT smearing, 2 - daughter pT and eta smearing
void Read_PYTHIA_tree( int mEnergy = 200 , int daughterSmearFlag = 0)
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
  
  //-------------------------------------------------------------------------------------
  
  //input file with momentum resolution

  TFile *momResFile = new TFile("./input/Momentum_resolution_Run16_SL16j_new_02.root", "read");
   
  TH1F *PiPlusMomRes = (TH1F*)momResFile->Get("PiPlusMomRes");
  TH1F *PiMinusMomRes = (TH1F*)momResFile->Get("PiMinusMomRes");
  
  TH1F *KPlusMomRes = (TH1F*)momResFile->Get("KPlusMomRes");
  TH1F *KMinusMomRes = (TH1F*)momResFile->Get("KMinusMomRes"); 
    
  TF1 *PiPlusMomResFit = new TF1("PiPlusMomResFit", "[0]+[1]/x+[2]*x+[3]*x*x+[4]/x/x", 0.01, 12); 
  TF1 *PiMinusMomResFit = new TF1("PiPlusMomResFit", "[0]+[1]/x+[2]*x+[3]*x*x+[4]/x/x", 0.01, 12); 
  
  TF1 *KPlusMomResFit = new TF1("PiPlusMomResFit", "[0]+[1]/x+[2]*x+[3]*x*x+[4]/x/x", 0.01, 12); 
  TF1 *KMinusMomResFit = new TF1("PiPlusMomResFit", "[0]+[1]/x+[2]*x+[3]*x*x+[4]/x/x", 0.01, 12);
  
  PiPlusMomRes->Fit(PiPlusMomResFit, "0 r");
  PiMinusMomRes->Fit(PiMinusMomResFit, "0 r");

  KPlusMomRes->Fit(KPlusMomResFit, "0 r", "", 0.3, 12.);
  KMinusMomRes->Fit(KMinusMomResFit, "0 r", "", 0.3, 12.);
  
  //-----------------------------------------------------------------------------
  
  TFile *etaResFile = new TFile("./input/PiPlus_eta_phi_resolution_Run16_SL16j_physics_work.root");
  
  TH1F *PiPlusEtaRes = (TH1F*)etaResFile->Get("PiPlusEtaRes");
  
  TF1 *FitEtaRes = new TF1("FitEtaRes", "[0]+[1]*x", -1, 1);
  FitEtaRes->SetParameters(0.045, 0.01);
  
  PiPlusEtaRes->Fit(FitEtaRes, "0 r");
  
  
  TH1F *PiPlusPhiRes = (TH1F*)etaResFile->Get("PiPlusPhiRes");
  
  TF1 *FitPhiRes = new TF1("FitPhiRes", "[0]+[1]*x", -TMath::Pi(), TMath::Pi());
  FitPhiRes->SetParameters(0.045, 0.01);
  
  PiPlusPhiRes->Fit(FitPhiRes, "0 r");
  
  //-----------------------------------------------------------------------------

  
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
  TH1F *L_p_pT = new TH1F("L_p_pT", "L_p_pT", 100, 0,5);
  TH1F *L_pi_pT = new TH1F("L_pi_pT", "L_pi_pT", 100, 0,5);

  TH1F *Lbar_p_pT = new TH1F("Lbar_p_pT", "Lbar_p_pT", 100, 0,5);
  TH1F *Lbar_pi_pT = new TH1F("Lbar_pi_pT", "Lbar_pi_pT", 100, 0,5);
  
  TH1F *L_Minv_hist = new TH1F("L_Minv_hist", "L_Minv_hist", 180, 1., 1.2);
  TH1F *Lbar_Minv_hist = new TH1F("Lbar_Minv_hist", "Lbar_Minv_hist", 180, 1., 1.2);
  
  TH1F *L_cos_theta = new TH1F("L_cos_theta", "L_cos_theta", 100, 0.98, 1.);
  
  
  TH1F *L0_L0bar_cosThetaProdPlane = new TH1F("L0_L0bar_cosThetaProdPlane", "L0_L0bar_cosThetaProdPlane", 10, -1, 1);  
  TH1F *L0_L0bar_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane = new TH1F("L0_L0_cosThetaProdPlane", "L0_L0_cosThetaProdPlane", 10, -1, 1);  
  TH1F *L0_L0_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane = new TH1F("L0bar_L0bar_cosThetaProdPlane", "L0bar_L0bar_cosThetaProdPlane", 10, -1, 1);  
  TH1F *L0bar_L0bar_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr]; 
  TH1F *L0bar_L0bar_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];
    
  //QA histograms with pair kinematics
  TH2F *L0_L0bar_eta1_vs_eta2_hist = new TH2F("L0_L0bar_eta1_vs_eta2_hist", "L0_L0bar_eta1_vs_eta2_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_phi1_vs_phi2_hist = new TH2F("L0_L0bar_phi1_vs_phi2_hist", "L0_L0bar_phi1_vs_phi2_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pT1_vs_pT2_hist = new TH2F("L0_L0bar_pT1_vs_pT2_hist", "L0_L0bar_pT1_vs_pT2_hist", 20, 0, 5, 20, 0, 5);
  
  TH2F *L0_L0_eta1_vs_eta2_hist = new TH2F("L0_L0_eta1_vs_eta2_hist", "L0_L0_eta1_vs_eta2_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0_phi1_vs_phi2_hist = new TH2F("L0_L0_phi1_vs_phi2_hist", "L0_L0_phi1_vs_phi2_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0_pT1_vs_pT2_hist = new TH2F("L0_L0_pT1_vs_pT2_hist", "L0_L0_pT1_vs_pT2_hist", 20, 0, 5, 20, 0, 5);
  
  TH2F *L0bar_L0bar_eta1_vs_eta2_hist = new TH2F("L0bar_L0bar_eta1_vs_eta2_hist", "L0bar_L0bar_eta1_vs_eta2_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0bar_L0bar_phi1_vs_phi2_hist = new TH2F("L0bar_L0bar_phi1_vs_phi2_hist", "L0bar_L0bar_phi1_vs_phi2_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0bar_L0bar_pT1_vs_pT2_hist = new TH2F("L0bar_L0bar_pT1_vs_pT2_hist", "L0bar_L0bar_pT1_vs_pT2_hist", 20, 0, 5, 20, 0, 5);
  
  //--------------------------------------------
  
  //pi kinematics
  TH2F *L0_L0bar_pi_pT1_vs_pi_pT2_hist = new TH2F("L0_L0bar_pi_pT1_vs_pi_pT2_hist", "L0_L0bar_pi_pT1_vs_pi_pT2_hist", 100, 0, 2, 100, 0, 2);

  
  
  //mother vs. daughter
  //from L-Lbar pairs
  
  //mother y vs. daughter eta
  TH2F *L0_y_vs_p_eta[nPtBins+1];
  TH2F *L0_y_vs_pi_eta[nPtBins+1];
  
  TH2F *L0bar_y_vs_p_eta[nPtBins+1];
  TH2F *L0bar_y_vs_pi_eta[nPtBins+1];
  
  //mother y vs. daughter pT
  TH2F *L0_y_vs_p_pT[nPtBins+1];
  TH2F *L0_y_vs_pi_pT[nPtBins+1];
  
  TH2F *L0bar_y_vs_p_pT[nPtBins+1];
  TH2F *L0bar_y_vs_pi_pT[nPtBins+1];
  
  //-------------------------------------------------------------------------------------------------------------------------------
  //mixed event histograms
  
  //mixed event histograms with re-weight after cuts
  TH1F *L0_L0bar_cosThetaProdPlane_ME = new TH1F("L0_L0bar_cosThetaProdPlane_ME", "L0_L0bar_cosThetaProdPlane_ME", 10, -1, 1);  
  TH1F *L0_L0bar_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_ME = new TH1F("L0_L0_cosThetaProdPlane_ME", "L0_L0_cosThetaProdPlane_ME", 10, -1, 1);  
  TH1F *L0_L0_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME = new TH1F("L0bar_L0bar_cosThetaProdPlane_ME", "L0bar_L0bar_cosThetaProdPlane_ME", 10, -1, 1);  
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr]; 
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];
  
 
  //---------------------------------------
  
  //for ME reweight
  TH2F *L0_L0bar_eta1_vs_eta2_ME_hist = new TH2F("L0_L0bar_eta1_vs_eta2_ME_hist", "L0_L0bar_eta1_vs_eta2_ME_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_phi1_vs_phi2_ME_hist = new TH2F("L0_L0bar_phi1_vs_phi2_ME_hist", "L0_L0bar_phi1_vs_phi2_ME_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pT1_vs_pT2_ME_hist = new TH2F("L0_L0bar_pT1_vs_pT2_ME_hist", "L0_L0bar_pT1_vs_pT2_ME_hist", 20, 0, 5, 20, 0, 5);
  
  //----------------------------------------------
  
  TH2F *L0_L0_eta1_vs_eta2_ME_hist = new TH2F("L0_L0_eta1_vs_eta2_ME_hist", "L0_L0_eta1_vs_eta2_ME_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0_phi1_vs_phi2_ME_hist = new TH2F("L0_L0_phi1_vs_phi2_ME_hist", "L0_L0_phi1_vs_phi2_ME_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0_pT1_vs_pT2_ME_hist = new TH2F("L0_L0_pT1_vs_pT2_ME_hist", "L0_L0_pT1_vs_pT2_ME_hist", 20, 0, 5, 20, 0, 5);
  
  //----------------------------------------------
  
  TH2F *L0bar_L0bar_eta1_vs_eta2_ME_hist = new TH2F("L0bar_L0bar_eta1_vs_eta2_ME_hist", "L0bar_L0bar_eta1_vs_eta2_ME_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0bar_L0bar_phi1_vs_phi2_ME_hist = new TH2F("L0bar_L0bar_phi1_vs_phi2_ME_hist", "L0bar_L0bar_phi1_vs_phi2_ME_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0bar_L0bar_pT1_vs_pT2_ME_hist = new TH2F("L0bar_L0bar_pT1_vs_pT2_ME_hist", "L0bar_L0bar_pT1_vs_pT2_ME_hist", 20, 0, 5, 20, 0, 5);
  
  //--------------------------------------------
  
  //pi kinematics - with ME reweight
  TH2F *L0_L0bar_pi_pT1_vs_pi_pT2_ME_hist = new TH2F("L0_L0bar_pi_pT1_vs_pi_pT2_ME_hist", "L0_L0bar_pi_pT1_vs_pi_pT2_ME_hist", 100, 0, 2, 100, 0, 2);

    
  //_____________________________________________________________________________________________________________________________________________________________________
  
  
    
 
  //____________________________________________________________________________________________________________________________________________________________________________________________________________________________

  //true MC histograms after cuts
  TH1F *L_p_pT_cuts = new TH1F("L_p_pT_cuts", "L_p_pT_cuts", 100, 0,5);
  TH1F *L_pi_pT_cuts = new TH1F("L_pi_pT_cuts", "L_pi_pT_cuts", 100, 0,5);

  TH1F *Lbar_p_pT_cuts = new TH1F("Lbar_p_pT_cuts", "Lbar_p_pT_cuts", 100, 0,5);
  TH1F *Lbar_pi_pT_cuts = new TH1F("Lbar_pi_pT_cuts", "Lbar_pi_pT_cuts", 100, 0,5);
  
  TH1F *L_Minv_hist_cuts = new TH1F("L_Minv_hist_cuts", "L_Minv_hist_cuts", 180, 1., 1.2);
  TH1F *Lbar_Minv_hist_cuts = new TH1F("Lbar_Minv_hist_cuts", "Lbar_Minv_hist_cuts", 180, 1., 1.2);
  
  //--------------------------------------------------------------------------------------
  
  TH1F *L0_L0bar_cosThetaProdPlane_cuts = new TH1F("L0_L0bar_cosThetaProdPlane_cuts", "L0_L0bar_cosThetaProdPlane_cuts", 10, -1, 1);  
  TH1F *L0_L0bar_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_cuts = new TH1F("L0_L0_cosThetaProdPlane_cuts", "L0_L0_cosThetaProdPlane_cuts", 10, -1, 1);  
  TH1F *L0_L0_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_cuts = new TH1F("L0bar_L0bar_cosThetaProdPlane_cuts", "L0bar_L0bar_cosThetaProdPlane_cuts", 10, -1, 1);  
  TH1F *L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];
  
  //for ME reweight
  TH2F *L0_L0bar_eta1_vs_eta2_cuts_hist = new TH2F("L0_L0bar_eta1_vs_eta2_cuts_hist", "L0_L0bar_eta1_vs_eta2_cuts_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_phi1_vs_phi2_cuts_hist = new TH2F("L0_L0bar_phi1_vs_phi2_cuts_hist", "L0_L0bar_phi1_vs_phi2_cuts_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pT1_vs_pT2_cuts_hist = new TH2F("L0_L0bar_pT1_vs_pT2_cuts_hist", "L0_L0bar_pT1_vs_pT2_cuts_hist", 20, 0, 5, 20, 0, 5);
    
  //--------------------------------------------
  
  TH2F *L0_L0_eta1_vs_eta2_cuts_hist = new TH2F("L0_L0_eta1_vs_eta2_cuts_hist", "L0_L0_eta1_vs_eta2_cuts_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0_phi1_vs_phi2_cuts_hist = new TH2F("L0_L0_phi1_vs_phi2_cuts_hist", "L0_L0_phi1_vs_phi2_cuts_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0_pT1_vs_pT2_cuts_hist = new TH2F("L0_L0_pT1_vs_pT2_cuts_hist", "L0_L0_pT1_vs_pT2_cuts_hist", 20, 0, 5, 20, 0, 5);
  
  //--------------------------------------------
  
  TH2F *L0bar_L0bar_eta1_vs_eta2_cuts_hist = new TH2F("L0bar_L0bar_eta1_vs_eta2_cuts_hist", "L0bar_L0bar_eta1_vs_eta2_cuts_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0bar_L0bar_phi1_vs_phi2_cuts_hist = new TH2F("L0bar_L0bar_phi1_vs_phi2_cuts_hist", "L0bar_L0bar_phi1_vs_phi2_cuts_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0bar_L0bar_pT1_vs_pT2_cuts_hist = new TH2F("L0bar_L0bar_pT1_vs_pT2_cuts_hist", "L0bar_L0bar_pT1_vs_pT2_cuts_hist", 20, 0, 5, 20, 0, 5);
  
  
  //--------------------------------------------
  
  //pi kinematics
  TH2F *L0_L0bar_pi_pT1_vs_pi_pT2_cuts_hist = new TH2F("L0_L0bar_pi_pT1_vs_pi_pT2_cuts_hist", "L0_L0bar_pi_pT1_vs_pi_pT2_cuts_hist", 100, 0, 2, 100, 0, 2);


  
  //--------------------------------------------------------------------------------------------------------------------------------------
  
   //mother vs. daughter
  //from L-Lbar pairs
  
  //mother y vs. daughter eta
  TH2F *L0_y_vs_p_eta_cuts[nPtBins+1];
  TH2F *L0_y_vs_pi_eta_cuts[nPtBins+1];
  
  TH2F *L0bar_y_vs_p_eta_cuts[nPtBins+1];
  TH2F *L0bar_y_vs_pi_eta_cuts[nPtBins+1];
  
  //mother y vs. daughter pT
  TH2F *L0_y_vs_p_pT_cuts[nPtBins+1];
  TH2F *L0_y_vs_pi_pT_cuts[nPtBins+1];
  
  TH2F *L0bar_y_vs_p_pT_cuts[nPtBins+1];
  TH2F *L0bar_y_vs_pi_pT_cuts[nPtBins+1];
  
  //----------------------------------------------------------------------------------------------------------------------------------------
  
 
  //mixed event histograms after re-weighing
  //weight from distributions after cuts
  TH1F *L0_L0bar_cosThetaProdPlane_ME_cuts = new TH1F("L0_L0bar_cosThetaProdPlane_ME_cuts", "L0_L0bar_cosThetaProdPlane_ME_cuts", 10, -1, 1);  
  TH1F *L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_ME_cuts = new TH1F("L0_L0_cosThetaProdPlane_ME_cuts", "L0_L0_cosThetaProdPlane_ME_cuts", 10, -1, 1);  
  TH1F *L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_cuts = new TH1F("L0bar_L0bar_cosThetaProdPlane_ME_cuts", "L0bar_L0bar_cosThetaProdPlane_ME_cuts", 10, -1, 1);  
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr]; 
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  //-----------------------------------
  
  //for ME reweight
  TH2F *L0_L0bar_eta1_vs_eta2_ME_cuts_hist = new TH2F("L0_L0bar_eta1_vs_eta2_ME_cuts_hist", "L0_L0bar_eta1_vs_eta2_ME_cuts_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_phi1_vs_phi2_ME_cuts_hist = new TH2F("L0_L0bar_phi1_vs_phi2_ME_cuts_hist", "L0_L0bar_phi1_vs_phi2_ME_cuts_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pT1_vs_pT2_ME_cuts_hist = new TH2F("L0_L0bar_pT1_vs_pT2_ME_cuts_hist", "L0_L0bar_pT1_vs_pT2_ME_cuts_hist", 20, 0, 5, 20, 0, 5);
  
  //----------------------------------------------
  
  TH2F *L0_L0_eta1_vs_eta2_ME_cuts_hist = new TH2F("L0_L0_eta1_vs_eta2_ME_cuts_hist", "L0_L0_eta1_vs_eta2_ME_cuts_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0_phi1_vs_phi2_ME_cuts_hist = new TH2F("L0_L0_phi1_vs_phi2_ME_cuts_hist", "L0_L0_phi1_vs_phi2_ME_cuts_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0_pT1_vs_pT2_ME_cuts_hist = new TH2F("L0_L0_pT1_vs_pT2_ME_cuts_hist", "L0_L0_pT1_vs_pT2_ME_cuts_hist", 20, 0, 5, 20, 0, 5);
  
  //----------------------------------------------
  
  TH2F *L0bar_L0bar_eta1_vs_eta2_ME_cuts_hist = new TH2F("L0bar_L0bar_eta1_vs_eta2_ME_cuts_hist", "L0bar_L0bar_eta1_vs_eta2_ME_cuts_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0bar_L0bar_phi1_vs_phi2_ME_cuts_hist = new TH2F("L0bar_L0bar_phi1_vs_phi2_ME_cuts_hist", "L0bar_L0bar_phi1_vs_phi2_ME_cuts_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0bar_L0bar_pT1_vs_pT2_ME_cuts_hist = new TH2F("L0bar_L0bar_pT1_vs_pT2_ME_cuts_hist", "L0bar_L0bar_pT1_vs_pT2_ME_cuts_hist", 20, 0, 5, 20, 0, 5);
  
  
  //--------------------------------------------
  
  //pi kinematics - with ME reweight
  //first is with weight from distributions after cuts
  //second is with weight from distributions before cuts
  TH2F *L0_L0bar_pi_pT1_vs_pi_pT2_ME_cuts_hist = new TH2F("L0_L0bar_pi_pT1_vs_pi_pT2_ME_cuts_hist", "L0_L0bar_pi_pT1_vs_pi_pT2_ME_cuts_hist", 100, 0, 2, 100, 0, 2);
  TH2F *L0_L0bar_pi_pT1_vs_pi_pT2_ME_cuts_2_hist = new TH2F("L0_L0bar_pi_pT1_vs_pi_pT2_ME_cuts_2_hist", "L0_L0bar_pi_pT1_vs_pi_pT2_ME_cuts_2_hist", 100, 0, 2, 100, 0, 2);

  
  
  //_____________________________________________________________________________________________________________________________________________________________________

    
  
  //QA histograms - mother vs. daughter kinematics in mother pT bins
  
  for(unsigned int pT_bin = 0; pT_bin < nPtBins; pT_bin++)
  {
    //mother vs. daughter
    //from L-Lbar pairs
    
    //mother y vs. daughter eta
    L0_y_vs_p_eta[pT_bin] = new TH2F(Form("L0_y_vs_p_eta_pT_bin_%i", pT_bin), Form("L0_y_vs_p_eta_pT_bin_%i", pT_bin), 100, -1, 1, 100, -5, 5);
    L0_y_vs_pi_eta[pT_bin] = new TH2F(Form("L0_y_vs_pi_eta_pT_bin_%i", pT_bin), Form("L0_y_vs_pi_eta_pT_bin_%i", pT_bin), 100, -1, 1, 100, -5, 5);
    
    L0bar_y_vs_p_eta[pT_bin] = new TH2F(Form("L0bar_y_vs_p_eta_pT_bin_%i", pT_bin), Form("L0bar_y_vs_p_eta_pT_bin_%i", pT_bin), 100, -1, 1, 100, -5, 5);
    L0bar_y_vs_pi_eta[pT_bin] = new TH2F(Form("L0bar_y_vs_pi_eta_pT_bin_%i", pT_bin), Form("L0bar_y_vs_pi_eta_pT_bin_%i", pT_bin), 100, -1, 1, 100, -5, 5);
    
    //mother y vs. daughter pT
    L0_y_vs_p_pT[pT_bin] = new TH2F(Form("L0_y_vs_p_pT_pT_bin_%i", pT_bin), Form("L0_y_vs_p_pT_pT_bin_%i", pT_bin), 100, -1, 1, 100, 0, 5);
    L0_y_vs_pi_pT[pT_bin] = new TH2F(Form("L0_y_vs_pi_pT_pT_bin_%i", pT_bin), Form("L0_y_vs_pi_pT_pT_bin_%i", pT_bin), 100, -1, 1, 100, 0, 5);
    
    L0bar_y_vs_p_pT[pT_bin] = new TH2F(Form("L0bar_y_vs_p_pT_pT_bin_%i", pT_bin), Form("L0bar_y_vs_p_pT_pT_bin_%i", pT_bin), 100, -1, 1, 100, 0, 5);
    L0bar_y_vs_pi_pT[pT_bin] = new TH2F(Form("L0bar_y_vs_pi_pT_pT_bin_%i", pT_bin), Form("L0bar_y_vs_pi_pT_pT_bin_%i", pT_bin), 100, -1, 1, 100, 0, 5);
    
    
    //____________________________________________________________________
    
    //mother y vs. daughter eta after cuts
    L0_y_vs_p_eta_cuts[pT_bin] = new TH2F(Form("L0_y_vs_p_eta_cuts_pT_bin_%i", pT_bin), Form("L0_y_vs_p_eta_cuts_pT_bin_%i", pT_bin), 100, -1, 1, 100, -5, 5);
    L0_y_vs_pi_eta_cuts[pT_bin] = new TH2F(Form("L0_y_vs_pi_eta_cuts_pT_bin_%i", pT_bin), Form("L0_y_vs_pi_eta_cuts_pT_bin_%i", pT_bin), 100, -1, 1, 100, -5, 5);
    
    L0bar_y_vs_p_eta_cuts[pT_bin] = new TH2F(Form("L0bar_y_vs_p_eta_cuts_pT_bin_%i", pT_bin), Form("L0bar_y_vs_p_eta_cuts_pT_bin_%i", pT_bin), 100, -1, 1, 100, -5, 5);
    L0bar_y_vs_pi_eta_cuts[pT_bin] = new TH2F(Form("L0bar_y_vs_pi_eta_cuts_pT_bin_%i", pT_bin), Form("L0bar_y_vs_pi_eta_cuts_pT_bin_%i", pT_bin), 100, -1, 1, 100, -5, 5);
    
    //mother y vs. daughter pT
    L0_y_vs_p_pT_cuts[pT_bin] = new TH2F(Form("L0_y_vs_p_pT_cuts_pT_bin_%i", pT_bin), Form("L0_y_vs_p_pT_cuts_pT_bin_%i", pT_bin), 100, -1, 1, 100, 0, 5);
    L0_y_vs_pi_pT_cuts[pT_bin] = new TH2F(Form("L0_y_vs_pi_pT_cuts_pT_bin_%i", pT_bin), Form("L0_y_vs_pi_pT_cuts_pT_bin_%i", pT_bin), 100, -1, 1, 100, 0, 5);
    
    L0bar_y_vs_p_pT_cuts[pT_bin] = new TH2F(Form("L0bar_y_vs_p_pT_cuts_pT_bin_%i", pT_bin), Form("L0bar_y_vs_p_pT_cuts_pT_bin_%i", pT_bin), 100, -1, 1, 100, 0, 5);
    L0bar_y_vs_pi_pT_cuts[pT_bin] = new TH2F(Form("L0bar_y_vs_pi_pT_cuts_pT_bin_%i", pT_bin), Form("L0bar_y_vs_pi_pT_cuts_pT_bin_%i", pT_bin), 100, -1, 1, 100, 0, 5);
  
  }
  
  
  
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
    
      //mixed event after re-weight after
      L0_L0bar_cosThetaProdPlane_ME_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
      L0_L0_cosThetaProdPlane_ME_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
      L0bar_L0bar_cosThetaProdPlane_ME_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
      //---------------------------------------------------------
      
      
      
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
      L0_L0bar_cosThetaProdPlane_ME_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_ME_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_ME_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      
      //---------------------------------------------------------
         
    }
  }
  //_____________________________________________________________________________________________________________________________________________________________________
  
  
  //mixed event vectors
  
  //SE L-Lbar pairs sotred to be used in mixed event
  //only one of the particles will be used in each ME pair
  vector<TVector3> L_Lbar_L_vector_ME_SE;
  vector<TVector3> L_Lbar_L_cuts_vector_ME_SE;
  
  vector<TVector3> L_Lbar_p_star_vector_ME_SE;
  vector<TVector3> L_Lbar_p_star_cuts_vector_ME_SE;
  
  vector<int> L_Lbar_L_pT_bin_vector_ME_SE;
  vector<int> L_Lbar_L_pT_bin_cuts_vector_ME_SE;
  
  vector<int> L_Lbar_L_eta_bin_vector_ME_SE;
  vector<int> L_Lbar_L_eta_bin_cuts_vector_ME_SE;  
  
  vector<float> L_Lbar_L_pi_pT_vector_ME_SE;
  vector<float> L_Lbar_L_pi_pT_cuts_vector_ME_SE;
  
  
  vector<TVector3> L_Lbar_Lbar_vector_ME_SE;
  vector<TVector3> L_Lbar_Lbar_cuts_vector_ME_SE;
  
  vector<TVector3> L_Lbar_pbar_star_vector_ME_SE;
  vector<TVector3> L_Lbar_pbar_star_cuts_vector_ME_SE;
  
  vector<int> L_Lbar_Lbar_pT_bin_vector_ME_SE;
  vector<int> L_Lbar_Lbar_pT_bin_cuts_vector_ME_SE;
  
  vector<int> L_Lbar_Lbar_eta_bin_vector_ME_SE;
  vector<int> L_Lbar_Lbar_eta_bin_cuts_vector_ME_SE;  
  
  vector<float> L_Lbar_Lbar_pi_pT_vector_ME_SE;
  vector<float> L_Lbar_Lbar_pi_pT_cuts_vector_ME_SE;
  
  //-------------------------------------
  
  //SE L-L pairs sotred to be used in mixed event
  //only one of the particles will be used in each ME pair
  vector<TVector3> L_L_L1_vector_ME_SE;
  vector<TVector3> L_L_L1_cuts_vector_ME_SE;
  
  vector<TVector3> L_L_p1_star_vector_ME_SE;
  vector<TVector3> L_L_p1_star_cuts_vector_ME_SE;
  
  vector<int> L_L_L1_pT_bin_vector_ME_SE;
  vector<int> L_L_L1_pT_bin_cuts_vector_ME_SE;
  
  vector<int> L_L_L1_eta_bin_vector_ME_SE;
  vector<int> L_L_L1_eta_bin_cuts_vector_ME_SE;  
  
  vector<float> L_L_L1_pi_pT_vector_ME_SE;
  vector<float> L_L_L1_pi_pT_cuts_vector_ME_SE;
  
  
  vector<TVector3> L_L_L2_vector_ME_SE;
  vector<TVector3> L_L_L2_cuts_vector_ME_SE;
  
  vector<TVector3> L_L_p2_star_vector_ME_SE;
  vector<TVector3> L_L_p2_star_cuts_vector_ME_SE;
  
  vector<int> L_L_L2_pT_bin_vector_ME_SE;
  vector<int> L_L_L2_pT_bin_cuts_vector_ME_SE;
  
  vector<int> L_L_L2_eta_bin_vector_ME_SE;
  vector<int> L_L_L2_eta_bin_cuts_vector_ME_SE;  
  
  vector<float> L_L_L2_pi_pT_vector_ME_SE;
  vector<float> L_L_L2_pi_pT_cuts_vector_ME_SE;
  
  
  //-------------------------------------
  
  //SE Lbar-Lbar pairs sotred to be used in mixed event
  //only one of the particles will be used in each ME pair
  vector<TVector3> Lbar_Lbar_Lbar1_vector_ME_SE;
  vector<TVector3> Lbar_Lbar_Lbar1_cuts_vector_ME_SE;
  
  vector<TVector3> Lbar_Lbar_pbar1_star_vector_ME_SE;
  vector<TVector3> Lbar_Lbar_pbar1_star_cuts_vector_ME_SE;
  
  vector<int> Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE;
  vector<int> Lbar_Lbar_Lbar1_pT_bin_cuts_vector_ME_SE;
  
  vector<int> Lbar_Lbar_Lbar1_eta_bin_vector_ME_SE;
  vector<int> Lbar_Lbar_Lbar1_eta_bin_cuts_vector_ME_SE;  
  
  vector<float> Lbar_Lbar_Lbar1_pi_pT_vector_ME_SE;
  vector<float> Lbar_Lbar_Lbar1_pi_pT_cuts_vector_ME_SE;
  
  
  vector<TVector3> Lbar_Lbar_Lbar2_vector_ME_SE;
  vector<TVector3> Lbar_Lbar_Lbar2_cuts_vector_ME_SE;
  
  vector<TVector3> Lbar_Lbar_pbar2_star_vector_ME_SE;
  vector<TVector3> Lbar_Lbar_pbar2_star_cuts_vector_ME_SE;
  
  vector<int> Lbar_Lbar_Lbar2_pT_bin_vector_ME_SE;
  vector<int> Lbar_Lbar_Lbar2_pT_bin_cuts_vector_ME_SE;
  
  vector<int> Lbar_Lbar_Lbar2_eta_bin_vector_ME_SE;
  vector<int> Lbar_Lbar_Lbar2_eta_bin_cuts_vector_ME_SE;  
  
  vector<float> Lbar_Lbar_Lbar2_pi_pT_vector_ME_SE;
  vector<float> Lbar_Lbar_Lbar2_pi_pT_cuts_vector_ME_SE; 
  
  
  //--------------------------------------------------------------------------------------------------
  
  //vectors to store L/Lbar from events with just one L/Lbar
  //These will be used as second particle in ME pairs
  vector<TVector3> L_vector_ME;
  vector<TVector3> L_cuts_vector_ME;
    
  vector<TVector3> p_star_vector_ME;
  vector<TVector3> p_star_cuts_vector_ME;
  
  vector<int> L_pT_bin_vector_ME;
  vector<int> L_pT_bin_cuts_vector_ME;
  
  vector<int> L_eta_bin_vector_ME;
  vector<int> L_eta_bin_cuts_vector_ME;  
  
  vector<float> L_pi_pT_vector_ME;
  vector<float> L_pi_pT_cuts_vector_ME;

  //-------------------------------------
  
  vector<TVector3> Lbar_vector_ME;
  vector<TVector3> Lbar_cuts_vector_ME;  
   
  vector<TVector3> pBar_star_vector_ME;
  vector<TVector3> pBar_star_cuts_vector_ME;  
 
  vector<int> Lbar_pT_bin_vector_ME;
  vector<int> Lbar_pT_bin_cuts_vector_ME;
  
  vector<int> Lbar_eta_bin_vector_ME;
  vector<int> Lbar_eta_bin_cuts_vector_ME;
  
  vector<float> Lbar_pi_pT_vector_ME;  
  vector<float> Lbar_pi_pT_cuts_vector_ME;
  
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
    
    vector<int> L_pT_bin_QA_vector;
    vector<int> L_pT_bin_QA_cuts_vector;
    
    vector<int> L_eta_bin_vector;
    vector<int> L_eta_bin_cuts_vector;
    
    
    vector<float> L_pi_pT_vector;
    
    
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
    
    vector<int> Lbar_pT_bin_QA_vector;
    vector<int> Lbar_pT_bin_QA_cuts_vector;
    
    vector<int> Lbar_eta_bin_vector;
    vector<int> Lbar_eta_bin_cuts_vector;
    

    vector<float> Lbar_pi_pT_vector;


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
      TVector3 L_mom;
      TVector3 L_mom_MC; //for pointing angle calculation with smearing
      TLorentzVector L_fourmom;
      
      //daughters
      TLorentzVector p_fourmom;
      p_fourmom.SetPtEtaPhiM(Lambda_MC.p_pT_MC[iLambda_MC], Lambda_MC.p_eta_MC[iLambda_MC], Lambda_MC.p_phi_MC[iLambda_MC], p_mass_PDG);
      
      TLorentzVector pi_fourmom;
      pi_fourmom.SetPtEtaPhiM(Lambda_MC.pi_pT_MC[iLambda_MC], Lambda_MC.pi_eta_MC[iLambda_MC], Lambda_MC.pi_phi_MC[iLambda_MC], pi_mass_PDG);
      
      TLorentzVector p_fourmom_smeared, pi_fourmom_smeared;

      //p momentum boosted to mother rest frame
      TVector3 p_mom_star;
    
      //no smearing
      if( daughterSmearFlag == 0 )
      {
        //Lambda momentum
        L_mom.SetXYZ(Lambda_MC.L_px_MC[iLambda_MC], Lambda_MC.L_py_MC[iLambda_MC], Lambda_MC.L_pz_MC[iLambda_MC]);
        L_mom_MC = L_mom; //without smearing, these are the same
        
        L_fourmom.SetVectM(L_mom, L_mass_PDG);
              
        //if(fabs(L_mom.Eta()) > 0.3 ) continue; //for testing of acceptance effect
        
        p_fourmom_smeared = p_fourmom;
        pi_fourmom_smeared = pi_fourmom;

        //p momentum boosted to mother rest frame
        p_mom_star.SetXYZ(Lambda_MC.p_pxStar_MC[iLambda_MC], Lambda_MC.p_pyStar_MC[iLambda_MC], Lambda_MC.p_pzStar_MC[iLambda_MC]);      
      }      
      else
      {     
        
        L_mom_MC.SetXYZ(Lambda_MC.L_px_MC[iLambda_MC], Lambda_MC.L_py_MC[iLambda_MC], Lambda_MC.L_pz_MC[iLambda_MC]);
        
        if( daughterSmearFlag == 1 ) //pT smearing
        {
          if(Lambda_MC.L_charge_MC[iLambda_MC] > 0) //L
          {
            p_fourmom_smeared = smearMom(p_fourmom, KPlusMomResFit);
            pi_fourmom_smeared = smearMom(pi_fourmom, PiMinusMomResFit);
          }
          
          if(Lambda_MC.L_charge_MC[iLambda_MC] < 0) //Lbar
          {
            p_fourmom_smeared = smearMom(p_fourmom, KMinusMomResFit);
            pi_fourmom_smeared = smearMom(pi_fourmom, PiPlusMomResFit);
          }       
        
        }
        
        if( daughterSmearFlag == 2 ) //pT and eta smearing
        {
          if(Lambda_MC.L_charge_MC[iLambda_MC] > 0) //L
          {
            p_fourmom_smeared = smearMomEtaPhi(p_fourmom, KPlusMomResFit, FitEtaRes, FitPhiRes);
            pi_fourmom_smeared = smearMomEtaPhi(pi_fourmom, PiMinusMomResFit, FitEtaRes, FitPhiRes);
          }
          
          if(Lambda_MC.L_charge_MC[iLambda_MC] < 0) //Lbar
          {
            p_fourmom_smeared = smearMomEtaPhi(p_fourmom, KMinusMomResFit, FitEtaRes, FitPhiRes);
            pi_fourmom_smeared = smearMomEtaPhi(pi_fourmom, PiPlusMomResFit, FitEtaRes, FitPhiRes);
          }
        
        }
        
        if( daughterSmearFlag == 3 ) //phi smearing
        {
          if(Lambda_MC.L_charge_MC[iLambda_MC] > 0) //L
          {
            p_fourmom_smeared = smearPhi(p_fourmom, FitPhiRes);
            pi_fourmom_smeared = smearPhi(pi_fourmom, FitPhiRes);
          }
          
          if(Lambda_MC.L_charge_MC[iLambda_MC] < 0) //Lbar
          {
            p_fourmom_smeared = smearPhi(p_fourmom, FitPhiRes);
            pi_fourmom_smeared = smearPhi(pi_fourmom, FitPhiRes);
          }
        
        }
        
        L_fourmom = p_fourmom_smeared + pi_fourmom_smeared;
        //L_fourmom.SetXYZM(Lambda_MC.L_px_MC[iLambda_MC], Lambda_MC.L_py_MC[iLambda_MC], Lambda_MC.L_pz_MC[iLambda_MC],L_mass_PDG); //to test momentum smearing
               
        L_mom = L_fourmom.Vect();
        
        //calculate p* for smeared momenta
        TLorentzVector L_fourmom_reverse;
        L_fourmom_reverse.SetVectM( -L_fourmom.Vect(), L_fourmom.M() );
        
        TLorentzVector p_fourmom_star = p_fourmom_smeared;
        p_fourmom_star.Boost(L_fourmom_reverse.BoostVector());
        
        p_mom_star = p_fourmom_star.Vect();        
      
      }
    
      
           
        
      //find bins               
      
      int pT_bin_corr = findBinPt( L_mom, pT_bins_corr, nPtBins_corr );
      
      if( pT_bin_corr == -1 ) continue;
      
      
      int pT_bin_QA = findBinPt(L_mom, pT_bins, nPtBins );
      
      if(pT_bin_QA == -1) continue;

      
      int eta_bin = findBinEta( L_mom, eta_bins, nEtaBins );

      if( eta_bin == -1 ) continue;


      //Lambda
      if( Lambda_MC.L_charge_MC[iLambda_MC] > 0 )
      {
        //QA histograms
        L_p_pT->Fill(p_fourmom_smeared.Pt());
        L_pi_pT->Fill(pi_fourmom_smeared.Pt());
      
        L_Minv_hist->Fill(L_fourmom.M());
        
        L0_y_vs_p_eta[pT_bin_QA]->Fill(L_fourmom.Rapidity(), p_fourmom_smeared.Eta());
        L0_y_vs_pi_eta[pT_bin_QA]->Fill(L_fourmom.Rapidity(), pi_fourmom_smeared.Eta());
        
        L0_y_vs_p_pT[pT_bin_QA]->Fill(L_fourmom.Rapidity(), p_fourmom_smeared.Pt());
        L0_y_vs_pi_pT[pT_bin_QA]->Fill(L_fourmom.Rapidity(), pi_fourmom_smeared.Pt());
        
        //---------------------------------------
        
        //vectors for pair analysis
        L_vector.push_back(L_mom);
        
        L_pi_pT_vector.push_back(pi_fourmom_smeared.Pt());
        
        L_decayL_MC_vector.push_back(Lambda_MC.L_decayL_MC[iLambda_MC]);
        
        L_pT_bin_vector.push_back(pT_bin_corr);
        L_pT_bin_QA_vector.push_back(pT_bin_QA);
        L_eta_bin_vector.push_back(eta_bin);
      
        p_star_vector.push_back(p_mom_star);
      }

      //Lambda-bar
      if( Lambda_MC.L_charge_MC[iLambda_MC] < 0 )
      { 
        //QA histograms
        Lbar_p_pT->Fill(p_fourmom_smeared.Pt());
        Lbar_pi_pT->Fill(pi_fourmom_smeared.Pt());
        
        Lbar_Minv_hist->Fill(L_fourmom.M());
        
        L0bar_y_vs_p_eta[pT_bin_QA]->Fill(L_fourmom.Rapidity(), p_fourmom_smeared.Eta());
        L0bar_y_vs_pi_eta[pT_bin_QA]->Fill(L_fourmom.Rapidity(), pi_fourmom_smeared.Eta());
        
        L0bar_y_vs_p_pT[pT_bin_QA]->Fill(L_fourmom.Rapidity(), p_fourmom_smeared.Pt());
        L0bar_y_vs_pi_pT[pT_bin_QA]->Fill(L_fourmom.Rapidity(), pi_fourmom_smeared.Pt());
        
        //---------------------------------------
        
        //vectors for pair analysis
        Lbar_vector.push_back(L_mom);
        
        Lbar_pi_pT_vector.push_back(pi_fourmom_smeared.Pt());
        
        Lbar_decayL_MC_vector.push_back(Lambda_MC.L_decayL_MC[iLambda_MC]);
        
        Lbar_pT_bin_vector.push_back(pT_bin_corr);
        Lbar_pT_bin_QA_vector.push_back(pT_bin_QA);
        Lbar_eta_bin_vector.push_back(eta_bin);

      
        pBar_star_vector.push_back(p_mom_star);      
      }
      
      
      //cosine of the pointing angle (for test of smearing)
      float theta = L_mom_MC.Angle(L_mom);
      
      L_cos_theta->Fill(cos(theta));
      
      if( cos(theta) < 0.998 ) continue;
      
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
        //QA histograms
        L_p_pT_cuts->Fill(p_fourmom_smeared.Pt());
        L_pi_pT_cuts->Fill(pi_fourmom_smeared.Pt());
        
        L_Minv_hist_cuts->Fill(L_fourmom.M());
        
        L0_y_vs_p_eta_cuts[pT_bin_QA]->Fill(L_fourmom.Rapidity(), p_fourmom_smeared.Eta());
        L0_y_vs_pi_eta_cuts[pT_bin_QA]->Fill(L_fourmom.Rapidity(), pi_fourmom_smeared.Eta());
        
        L0_y_vs_p_pT_cuts[pT_bin_QA]->Fill(L_fourmom.Rapidity(), p_fourmom_smeared.Pt());
        L0_y_vs_pi_pT_cuts[pT_bin_QA]->Fill(L_fourmom.Rapidity(), pi_fourmom_smeared.Pt());
        
        //---------------------------------------
        
        //vectors for pair analysis      
        L_cuts_vector.push_back(L_mom);
        
        L_decayL_MC_cuts_vector.push_back(Lambda_MC.L_decayL_MC[iLambda_MC]);
        
        L_pT_bin_cuts_vector.push_back(pT_bin_corr);
        L_pT_bin_QA_cuts_vector.push_back(pT_bin_QA);
        L_eta_bin_cuts_vector.push_back(eta_bin);
      
        p_star_cuts_vector.push_back(p_mom_star);
        
        L_p_pT_cuts_vector.push_back(p_fourmom_smeared.Pt()); 
        L_pi_pT_cuts_vector.push_back(pi_fourmom_smeared.Pt());
        
        L_p_eta_cuts_vector.push_back(p_fourmom_smeared.Eta()); 
        L_pi_eta_cuts_vector.push_back(pi_fourmom_smeared.Eta());
        
        L_p_phi_cuts_vector.push_back(p_fourmom_smeared.Phi()); 
        L_pi_phi_cuts_vector.push_back(pi_fourmom_smeared.Phi());
     
      }

      //Lambda-bar
      if( Lambda_MC.L_charge_MC[iLambda_MC] < 0 )
      {
        //QA histograms
        Lbar_p_pT_cuts->Fill(p_fourmom_smeared.Pt());
        Lbar_pi_pT_cuts->Fill(pi_fourmom_smeared.Pt());
        
        Lbar_Minv_hist_cuts->Fill(L_fourmom.M());
        
        L0bar_y_vs_p_eta_cuts[pT_bin_QA]->Fill(L_fourmom.Rapidity(), p_fourmom_smeared.Eta());
        L0bar_y_vs_pi_eta_cuts[pT_bin_QA]->Fill(L_fourmom.Rapidity(), pi_fourmom_smeared.Eta());
        
        L0bar_y_vs_p_pT_cuts[pT_bin_QA]->Fill(L_fourmom.Rapidity(), p_fourmom_smeared.Pt());
        L0bar_y_vs_pi_pT_cuts[pT_bin_QA]->Fill(L_fourmom.Rapidity(), pi_fourmom_smeared.Pt());
        
        //---------------------------------------
        
        //vectors for pair analysis        
        Lbar_cuts_vector.push_back(L_mom);
        
        Lbar_decayL_MC_cuts_vector.push_back(Lambda_MC.L_decayL_MC[iLambda_MC]);
        
        Lbar_pT_bin_cuts_vector.push_back(pT_bin_corr);
        Lbar_pT_bin_QA_cuts_vector.push_back(pT_bin_QA);
        Lbar_eta_bin_cuts_vector.push_back(eta_bin);
      
        pBar_star_cuts_vector.push_back(p_mom_star);
        
        Lbar_p_pT_cuts_vector.push_back(p_fourmom_smeared.Pt()); 
        Lbar_pi_pT_cuts_vector.push_back(pi_fourmom_smeared.Pt());
        
        Lbar_p_eta_cuts_vector.push_back(p_fourmom_smeared.Eta()); 
        Lbar_pi_eta_cuts_vector.push_back(pi_fourmom_smeared.Eta());
        
        Lbar_p_phi_cuts_vector.push_back(p_fourmom_smeared.Phi()); 
        Lbar_pi_phi_cuts_vector.push_back(pi_fourmom_smeared.Phi());
      
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
          
          
          L0_L0bar_pi_pT1_vs_pi_pT2_hist->Fill(L_pi_pT_vector.at(iLambda), Lbar_pi_pT_vector.at(iLambdaBar));
         
         
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
          
          //for ME reweight
          L0_L0_eta1_vs_eta2_hist->Fill(L_vector.at(iLambda1).Eta(), L_vector.at(iLambda2).Eta() );
          L0_L0_phi1_vs_phi2_hist->Fill(L_vector.at(iLambda1).Phi(), L_vector.at(iLambda2).Phi() );
          L0_L0_pT1_vs_pT2_hist->Fill(L_vector.at(iLambda1).Pt(), L_vector.at(iLambda2).Pt() );
                      
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
          
          
          //for ME re-weighing          
          L0bar_L0bar_eta1_vs_eta2_hist->Fill(Lbar_vector.at(iLambdaBar1).Eta(), Lbar_vector.at(iLambdaBar2).Eta() );
          L0bar_L0bar_phi1_vs_phi2_hist->Fill(Lbar_vector.at(iLambdaBar1).Phi(), Lbar_vector.at(iLambdaBar2).Phi() );
          L0bar_L0bar_pT1_vs_pT2_hist->Fill(Lbar_vector.at(iLambdaBar1).Pt(), Lbar_vector.at(iLambdaBar2).Pt() );
        }     
      }     
    }
    //_____________________________________________________________________________________________________________________________________________________________________

    //mixed event before cuts

    //SE L-Lbar for ME
    if( p_star_vector.size() > 0 && pBar_star_vector.size() > 0 && L_Lbar_L_vector_ME_SE.size() < 1e4) //test ME - L is from good L-Lbar pair and will be paired with Lbar from different event with correct kinematics
    {
      L_Lbar_L_vector_ME_SE.push_back(L_vector.at(0));

      L_Lbar_p_star_vector_ME_SE.push_back(p_star_vector.at(0));

      L_Lbar_L_pT_bin_vector_ME_SE.push_back(L_pT_bin_vector.at(0));
      L_Lbar_L_eta_bin_vector_ME_SE.push_back(L_eta_bin_vector.at(0));

      L_Lbar_L_pi_pT_vector_ME_SE.push_back(L_pi_pT_vector.at(0));

      //------------------------------------------------------------------

      L_Lbar_Lbar_vector_ME_SE.push_back(Lbar_vector.at(0));

      L_Lbar_pbar_star_vector_ME_SE.push_back(pBar_star_vector.at(0));

      L_Lbar_Lbar_pT_bin_vector_ME_SE.push_back(Lbar_pT_bin_vector.at(0));
      L_Lbar_Lbar_eta_bin_vector_ME_SE.push_back(Lbar_eta_bin_vector.at(0));

      L_Lbar_Lbar_pi_pT_vector_ME_SE.push_back(Lbar_pi_pT_vector.at(0));

    }

    //SE L-L for ME
    if( p_star_vector.size() > 1  && L_L_L1_vector_ME_SE.size() < 1e4) //test ME - L is from good L-Lbar pair and will be paired with Lbar from different event with correct kinematics
    {
      L_L_L1_vector_ME_SE.push_back(L_vector.at(0));

      L_L_p1_star_vector_ME_SE.push_back(p_star_vector.at(0));

      L_L_L1_pT_bin_vector_ME_SE.push_back(L_pT_bin_vector.at(0));
      L_L_L1_eta_bin_vector_ME_SE.push_back(L_eta_bin_vector.at(0));

      L_L_L1_pi_pT_vector_ME_SE.push_back(L_pi_pT_vector.at(0));

      //------------------------------------------------------------------

      L_L_L2_vector_ME_SE.push_back(L_vector.at(1));

      L_L_p2_star_vector_ME_SE.push_back(p_star_vector.at(1));

      L_L_L2_pT_bin_vector_ME_SE.push_back(L_pT_bin_vector.at(1));
      L_L_L2_eta_bin_vector_ME_SE.push_back(L_eta_bin_vector.at(1));

      L_L_L2_pi_pT_vector_ME_SE.push_back(L_pi_pT_vector.at(1));

    }


    //SE Lbar-Lbar for ME
    if( pBar_star_vector.size() > 1  && Lbar_Lbar_Lbar1_vector_ME_SE.size() < 1e4) //test ME - L is from good L-Lbar pair and will be paired with Lbar from different event with correct kinematics
    {
      Lbar_Lbar_Lbar1_vector_ME_SE.push_back(Lbar_vector.at(0));

      Lbar_Lbar_pbar1_star_vector_ME_SE.push_back(pBar_star_vector.at(0));

      Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE.push_back(Lbar_pT_bin_vector.at(0));
      Lbar_Lbar_Lbar1_eta_bin_vector_ME_SE.push_back(Lbar_eta_bin_vector.at(0));

      Lbar_Lbar_Lbar1_pi_pT_vector_ME_SE.push_back(Lbar_pi_pT_vector.at(0));

      //------------------------------------------------------------------

      Lbar_Lbar_Lbar2_vector_ME_SE.push_back(Lbar_vector.at(1));

      Lbar_Lbar_pbar2_star_vector_ME_SE.push_back(pBar_star_vector.at(1));

      Lbar_Lbar_Lbar2_pT_bin_vector_ME_SE.push_back(Lbar_pT_bin_vector.at(1));
      Lbar_Lbar_Lbar2_eta_bin_vector_ME_SE.push_back(Lbar_eta_bin_vector.at(1));

      Lbar_Lbar_Lbar2_pi_pT_vector_ME_SE.push_back(Lbar_pi_pT_vector.at(1));

    }


    //---------------------------------------------------------------------------------------------------------------


    //signle L after cuts for ME
    if( p_star_vector.size() == 1 && pBar_star_vector.size() == 0 && p_star_vector_ME.size() < 1e5) //this will be the same for both versions of ME
    {
      L_vector_ME.push_back(L_vector.at(0));

      p_star_vector_ME.push_back(p_star_vector.at(0));

      L_pT_bin_vector_ME.push_back(L_pT_bin_vector.at(0));
      L_eta_bin_vector_ME.push_back(L_eta_bin_vector.at(0));

      L_pi_pT_vector_ME.push_back(L_pi_pT_vector.at(0));

    }

    //signle Lbar after cuts for ME
    if( p_star_vector.size() == 0 && pBar_star_vector.size() == 1 && pBar_star_vector_ME.size() < 1e5) //this will be the same for both versions of ME
    {
      Lbar_vector_ME.push_back(Lbar_vector.at(0));

      pBar_star_vector_ME.push_back(pBar_star_vector.at(0));

      Lbar_pT_bin_vector_ME.push_back(Lbar_pT_bin_vector.at(0));
      Lbar_eta_bin_vector_ME.push_back(Lbar_eta_bin_vector.at(0));

      Lbar_pi_pT_vector_ME.push_back(Lbar_pi_pT_vector.at(0));

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
          
          //for ME re-weighing
          L0_L0bar_eta1_vs_eta2_cuts_hist->Fill(L_cuts_vector.at(iLambda).Eta() , Lbar_cuts_vector.at(iLambdaBar).Eta());
          L0_L0bar_phi1_vs_phi2_cuts_hist->Fill(L_cuts_vector.at(iLambda).Phi() , Lbar_cuts_vector.at(iLambdaBar).Phi());
          L0_L0bar_pT1_vs_pT2_cuts_hist->Fill(L_cuts_vector.at(iLambda).Pt() , Lbar_cuts_vector.at(iLambdaBar).Pt());
          
          
          L0_L0bar_pi_pT1_vs_pi_pT2_cuts_hist->Fill(L_pi_pT_cuts_vector.at(iLambda), Lbar_pi_pT_cuts_vector.at(iLambdaBar));
          
        }   
      
      }
    }
    

   
    //L0-L0
    if(p_star_cuts_vector.size() > 1)
    {
      for(unsigned int iLambda1 = 0; iLambda1 < p_star_cuts_vector.size(); iLambda1++)
      {
        for(unsigned int iLambda2 = iLambda1+1; iLambda2 < p_star_cuts_vector.size(); iLambda2++)
        {
          double theta_star = p_star_cuts_vector.at(iLambda1).Angle(p_star_cuts_vector.at(iLambda2));
          
          L0_L0_cosThetaProdPlane_cuts->Fill(TMath::Cos(theta_star));
          L0_L0_cosThetaProdPlane_pT_cuts_hist[L_pT_bin_cuts_vector.at(iLambda1)][L_pT_bin_cuts_vector.at(iLambda2)]->Fill(TMath::Cos(theta_star));
          L0_L0_cosThetaProdPlane_eta_cuts_hist[L_eta_bin_cuts_vector.at(iLambda1)][L_eta_bin_cuts_vector.at(iLambda2)]->Fill(TMath::Cos(theta_star));           

          //for ME re-weighing
          L0_L0_eta1_vs_eta2_cuts_hist->Fill(L_cuts_vector.at(iLambda1).Eta(), L_cuts_vector.at(iLambda2).Eta() );
          L0_L0_phi1_vs_phi2_cuts_hist->Fill(L_cuts_vector.at(iLambda1).Phi(), L_cuts_vector.at(iLambda2).Phi() );
          L0_L0_pT1_vs_pT2_cuts_hist->Fill(L_cuts_vector.at(iLambda1).Pt(), L_cuts_vector.at(iLambda2).Pt() );
          
        }   
      
      }
    
    }
    
    //L0bar-L0bar
    if(pBar_star_cuts_vector.size() > 1)
    {
      for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pBar_star_cuts_vector.size(); iLambdaBar1++)
      {
        for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < pBar_star_cuts_vector.size(); iLambdaBar2++)
        {
          double theta_star = pBar_star_cuts_vector.at(iLambdaBar1).Angle(pBar_star_cuts_vector.at(iLambdaBar2));
          
          L0bar_L0bar_cosThetaProdPlane_cuts->Fill(TMath::Cos(theta_star));
          L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[Lbar_pT_bin_cuts_vector.at(iLambdaBar1)][Lbar_pT_bin_cuts_vector.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
          L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[Lbar_eta_bin_cuts_vector.at(iLambdaBar1)][Lbar_eta_bin_cuts_vector.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
          
          //for ME re-weighing
          L0bar_L0bar_eta1_vs_eta2_cuts_hist->Fill(Lbar_cuts_vector.at(iLambdaBar1).Eta(), Lbar_cuts_vector.at(iLambdaBar2).Eta() );
          L0bar_L0bar_phi1_vs_phi2_cuts_hist->Fill(Lbar_cuts_vector.at(iLambdaBar1).Phi(), Lbar_cuts_vector.at(iLambdaBar2).Phi() );
          L0bar_L0bar_pT1_vs_pT2_cuts_hist->Fill(Lbar_cuts_vector.at(iLambdaBar1).Pt(), Lbar_cuts_vector.at(iLambdaBar2).Pt() );
        }   
      
      } 
    
    }
    
    //mixed event after cuts
    
    //SE L-Lbar for ME
    if( p_star_cuts_vector.size() > 0 && pBar_star_cuts_vector.size() > 0 && L_Lbar_L_cuts_vector_ME_SE.size() < 1e4) //test ME - L is from good L-Lbar pair and will be paired with Lbar from different event with correct kinematics
    {
      L_Lbar_L_cuts_vector_ME_SE.push_back(L_cuts_vector.at(0));      
      
      L_Lbar_p_star_cuts_vector_ME_SE.push_back(p_star_cuts_vector.at(0));

      L_Lbar_L_pT_bin_cuts_vector_ME_SE.push_back(L_pT_bin_cuts_vector.at(0));
      L_Lbar_L_eta_bin_cuts_vector_ME_SE.push_back(L_eta_bin_cuts_vector.at(0));
      
      L_Lbar_L_pi_pT_cuts_vector_ME_SE.push_back(L_pi_pT_cuts_vector.at(0));
      
      //------------------------------------------------------------------
      
      L_Lbar_Lbar_cuts_vector_ME_SE.push_back(Lbar_cuts_vector.at(0));
      
      L_Lbar_pbar_star_cuts_vector_ME_SE.push_back(pBar_star_cuts_vector.at(0));

      L_Lbar_Lbar_pT_bin_cuts_vector_ME_SE.push_back(Lbar_pT_bin_cuts_vector.at(0));
      L_Lbar_Lbar_eta_bin_cuts_vector_ME_SE.push_back(Lbar_eta_bin_cuts_vector.at(0));
      
      L_Lbar_Lbar_pi_pT_cuts_vector_ME_SE.push_back(Lbar_pi_pT_cuts_vector.at(0));

    }
    
    //SE L-L for ME
    if( p_star_cuts_vector.size() > 1  && L_L_L1_cuts_vector_ME_SE.size() < 1e4) //test ME - L is from good L-Lbar pair and will be paired with Lbar from different event with correct kinematics
    {
      L_L_L1_cuts_vector_ME_SE.push_back(L_cuts_vector.at(0));      
      
      L_L_p1_star_cuts_vector_ME_SE.push_back(p_star_cuts_vector.at(0));

      L_L_L1_pT_bin_cuts_vector_ME_SE.push_back(L_pT_bin_cuts_vector.at(0));
      L_L_L1_eta_bin_cuts_vector_ME_SE.push_back(L_eta_bin_cuts_vector.at(0));
      
      L_L_L1_pi_pT_cuts_vector_ME_SE.push_back(L_pi_pT_cuts_vector.at(0));
      
      //------------------------------------------------------------------
      
      L_L_L2_cuts_vector_ME_SE.push_back(L_cuts_vector.at(1));
      
      L_L_p2_star_cuts_vector_ME_SE.push_back(p_star_cuts_vector.at(1));

      L_L_L2_pT_bin_cuts_vector_ME_SE.push_back(L_pT_bin_cuts_vector.at(1));
      L_L_L2_eta_bin_cuts_vector_ME_SE.push_back(L_eta_bin_cuts_vector.at(1));
      
      L_L_L2_pi_pT_cuts_vector_ME_SE.push_back(L_pi_pT_cuts_vector.at(1));

    }
    
    
    //SE Lbar-Lbar for ME
    if( pBar_star_cuts_vector.size() > 1  && Lbar_Lbar_Lbar1_cuts_vector_ME_SE.size() < 1e4) //test ME - L is from good L-Lbar pair and will be paired with Lbar from different event with correct kinematics
    {
      Lbar_Lbar_Lbar1_cuts_vector_ME_SE.push_back(Lbar_cuts_vector.at(0));      
      
      Lbar_Lbar_pbar1_star_cuts_vector_ME_SE.push_back(pBar_star_cuts_vector.at(0));

      Lbar_Lbar_Lbar1_pT_bin_cuts_vector_ME_SE.push_back(Lbar_pT_bin_cuts_vector.at(0));
      Lbar_Lbar_Lbar1_eta_bin_cuts_vector_ME_SE.push_back(Lbar_eta_bin_cuts_vector.at(0));
      
      Lbar_Lbar_Lbar1_pi_pT_cuts_vector_ME_SE.push_back(Lbar_pi_pT_cuts_vector.at(0));
      
      //------------------------------------------------------------------
      
      Lbar_Lbar_Lbar2_cuts_vector_ME_SE.push_back(Lbar_cuts_vector.at(1));
      
      Lbar_Lbar_pbar2_star_cuts_vector_ME_SE.push_back(pBar_star_cuts_vector.at(1));

      Lbar_Lbar_Lbar2_pT_bin_cuts_vector_ME_SE.push_back(Lbar_pT_bin_cuts_vector.at(1));
      Lbar_Lbar_Lbar2_eta_bin_cuts_vector_ME_SE.push_back(Lbar_eta_bin_cuts_vector.at(1));
      
      Lbar_Lbar_Lbar2_pi_pT_cuts_vector_ME_SE.push_back(Lbar_pi_pT_cuts_vector.at(1));

    }
    
    
    //---------------------------------------------------------------------------------------------------------------


    //signle L after cuts for ME
    if( p_star_cuts_vector.size() == 1 && pBar_star_cuts_vector.size() == 0 && p_star_cuts_vector_ME.size() < 1e4) //this will be the same for both versions of ME
    {    
      L_cuts_vector_ME.push_back(L_cuts_vector.at(0));
            
      p_star_cuts_vector_ME.push_back(p_star_cuts_vector.at(0));

      L_pT_bin_cuts_vector_ME.push_back(L_pT_bin_cuts_vector.at(0));
      L_eta_bin_cuts_vector_ME.push_back(L_eta_bin_cuts_vector.at(0));

      L_pi_pT_cuts_vector_ME.push_back(L_pi_pT_cuts_vector.at(0));

    }

    //signle Lbar after cuts for ME
    if( p_star_cuts_vector.size() == 0 && pBar_star_cuts_vector.size() == 1 && pBar_star_cuts_vector_ME.size() < 1e4) //this will be the same for both versions of ME
    {    
      Lbar_cuts_vector_ME.push_back(Lbar_cuts_vector.at(0));
            
      pBar_star_cuts_vector_ME.push_back(pBar_star_cuts_vector.at(0));

      Lbar_pT_bin_cuts_vector_ME.push_back(Lbar_pT_bin_cuts_vector.at(0));
      Lbar_eta_bin_cuts_vector_ME.push_back(Lbar_eta_bin_cuts_vector.at(0));

      Lbar_pi_pT_cuts_vector_ME.push_back(Lbar_pi_pT_cuts_vector.at(0));

    }
    
    //---------------------------------------------------------------------------------------------------------------------
    
    L_vector.clear();
    L_cuts_vector.clear();
    
    L_pi_pT_vector.clear();
    
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
    
    Lbar_pi_pT_vector.clear();
    
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


  //_____________________________________________________________________________________________________________________________________________________________________
  
  //cout<<p_star_vector_ME.size()<<endl;
  //cout<<pBar_star_vector_ME.size()<<endl;
  
  cout<<"Start MC ME"<<endl;  
    
  //mixed event before cuts
  //L-Lbar
  //L from SE mixed with single Lbar
  for(unsigned int iLambda = 0; iLambda < L_Lbar_p_star_vector_ME_SE.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_ME.size(); iLambdaBar++)
    {
      //limit kinematics of ME Lbar based on kinematics of same event Lbar - double check precision
      if( fabs( L_Lbar_Lbar_vector_ME_SE.at(iLambda).Eta() - Lbar_vector_ME.at(iLambdaBar).Eta()) > 0.1 ) continue;
      if( fabs( L_Lbar_Lbar_vector_ME_SE.at(iLambda).Phi() - Lbar_vector_ME.at(iLambdaBar).Phi()) > 0.1 ) continue;
      if( fabs( L_Lbar_Lbar_vector_ME_SE.at(iLambda).Pt() - Lbar_vector_ME.at(iLambdaBar).Pt()) > 0.1 ) continue;


      double theta_star = L_Lbar_p_star_vector_ME_SE.at(iLambda).Angle(pBar_star_vector_ME.at(iLambdaBar));


      L0_L0bar_cosThetaProdPlane_ME->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_ME_pT_hist[L_Lbar_L_pT_bin_vector_ME_SE.at(iLambda)][Lbar_pT_bin_vector_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_ME_eta_hist[L_Lbar_L_eta_bin_vector_ME_SE.at(iLambda)][Lbar_eta_bin_vector_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));


      //QA histograms with L-Lbar kinematics
      L0_L0bar_eta1_vs_eta2_ME_hist->Fill(L_Lbar_L_vector_ME_SE.at(iLambda).Eta() , Lbar_vector_ME.at(iLambdaBar).Eta());
      L0_L0bar_phi1_vs_phi2_ME_hist->Fill(L_Lbar_L_vector_ME_SE.at(iLambda).Phi() , Lbar_vector_ME.at(iLambdaBar).Phi());
      L0_L0bar_pT1_vs_pT2_ME_hist->Fill(L_Lbar_L_vector_ME_SE.at(iLambda).Pt() , Lbar_vector_ME.at(iLambdaBar).Pt());

      //decay pi kinematics
      L0_L0bar_pi_pT1_vs_pi_pT2_ME_hist->Fill(L_Lbar_L_pi_pT_vector_ME_SE.at(iLambda), Lbar_pi_pT_vector_ME.at(iLambdaBar));

    }

  }

  //Lbar from SE mixed with single L
  for(unsigned int iLambda = 0; iLambda < p_star_vector_ME.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < L_Lbar_pbar_star_vector_ME_SE.size(); iLambdaBar++)
    {
      //limit kinematics of ME L based on kinematics of same event L - double check precision
      if( fabs( L_vector_ME.at(iLambda).Eta() - L_Lbar_L_vector_ME_SE.at(iLambdaBar).Eta() ) > 0.1 ) continue;
      if( fabs( L_vector_ME.at(iLambda).Phi() - L_Lbar_L_vector_ME_SE.at(iLambdaBar).Phi() ) > 0.1 ) continue;
      if( fabs( L_vector_ME.at(iLambda).Pt() -  L_Lbar_L_vector_ME_SE.at(iLambdaBar).Pt() ) > 0.1 ) continue;


      double theta_star = p_star_vector_ME.at(iLambda).Angle(L_Lbar_pbar_star_vector_ME_SE.at(iLambdaBar));


      L0_L0bar_cosThetaProdPlane_ME->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_ME_pT_hist[L_pT_bin_vector_ME.at(iLambda)][L_Lbar_Lbar_pT_bin_vector_ME_SE.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_ME_eta_hist[L_eta_bin_vector_ME.at(iLambda)][L_Lbar_Lbar_eta_bin_vector_ME_SE.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));


      //QA histograms with L-Lbar kinematics
      L0_L0bar_eta1_vs_eta2_ME_hist->Fill(L_vector_ME.at(iLambda).Eta() , L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar).Eta());
      L0_L0bar_phi1_vs_phi2_ME_hist->Fill(L_vector_ME.at(iLambda).Phi() , L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar).Phi());
      L0_L0bar_pT1_vs_pT2_ME_hist->Fill(L_vector_ME.at(iLambda).Pt() , L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar).Pt());

      //decay pi kinematics
      L0_L0bar_pi_pT1_vs_pi_pT2_ME_hist->Fill(L_pi_pT_vector_ME.at(iLambda), L_Lbar_Lbar_pi_pT_vector_ME_SE.at(iLambdaBar));

    }

  }

  //---------------------------------------------------------------------------------------------------------------------------

  //L0-L0
  //L1 is from SE, L2 is from ME
  for(unsigned int iLambda1 = 0; iLambda1 < L_L_p1_star_vector_ME_SE.size(); iLambda1++)
  {
    for(unsigned int iLambda2 = 0; iLambda2 < p_star_vector_ME.size(); iLambda2++)
    {
      //limit kinematics of ME L2 based on kinematics of same event L2 - double check precision
      if( fabs( L_L_L2_vector_ME_SE.at(iLambda1).Eta() - L_vector_ME.at(iLambda2).Eta()) > 0.1 ) continue;
      if( fabs( L_L_L2_vector_ME_SE.at(iLambda1).Phi() - L_vector_ME.at(iLambda2).Phi()) > 0.1 ) continue;
      if( fabs( L_L_L2_vector_ME_SE.at(iLambda1).Pt() - L_vector_ME.at(iLambda2).Pt()) > 0.1 ) continue;

      double theta_star = L_L_p1_star_vector_ME_SE.at(iLambda1).Angle(p_star_vector_ME.at(iLambda2));

      L0_L0_cosThetaProdPlane_ME->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_ME_pT_hist[L_L_L1_pT_bin_vector_ME_SE.at(iLambda1)][L_pT_bin_vector_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_ME_eta_hist[L_L_L1_eta_bin_vector_ME_SE.at(iLambda1)][L_eta_bin_vector_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star));
      
      //QA histograms with L-Lbar kinematics
      L0_L0_eta1_vs_eta2_ME_hist->Fill(L_L_L1_vector_ME_SE.at(iLambda1).Eta() , L_vector_ME.at(iLambda2).Eta());
      L0_L0_phi1_vs_phi2_ME_hist->Fill(L_L_L1_vector_ME_SE.at(iLambda1).Phi() , L_vector_ME.at(iLambda2).Phi());
      L0_L0_pT1_vs_pT2_ME_hist->Fill(L_L_L1_vector_ME_SE.at(iLambda1).Pt() , L_vector_ME.at(iLambda2).Pt());

    }

  }

  //L1 is from ME, L2 is from SE
  for(unsigned int iLambda1 = 0; iLambda1 < p_star_vector_ME.size(); iLambda1++)
  {
    for(unsigned int iLambda2 = 0; iLambda2 < L_L_p2_star_vector_ME_SE.size(); iLambda2++)
    {
      //limit kinematics of ME L2 based on kinematics of same event L2 - double check precision
      if( fabs( L_vector_ME.at(iLambda1).Eta() - L_L_L1_vector_ME_SE.at(iLambda2).Eta()) > 0.1 ) continue;
      if( fabs( L_vector_ME.at(iLambda1).Phi() - L_L_L1_vector_ME_SE.at(iLambda2).Phi()) > 0.1 ) continue;
      if( fabs( L_vector_ME.at(iLambda1).Pt() - L_L_L1_vector_ME_SE.at(iLambda2).Pt()) > 0.1 ) continue;

      double theta_star = p_star_vector_ME.at(iLambda1).Angle(L_L_p2_star_vector_ME_SE.at(iLambda2));

      L0_L0_cosThetaProdPlane_ME->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_ME_pT_hist[L_pT_bin_vector_ME.at(iLambda1)][L_L_L2_pT_bin_vector_ME_SE.at(iLambda2)]->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_ME_eta_hist[L_eta_bin_vector_ME.at(iLambda1)][L_L_L2_eta_bin_vector_ME_SE.at(iLambda2)]->Fill(TMath::Cos(theta_star));
      
      //QA histograms with L-Lbar kinematics
      L0_L0_eta1_vs_eta2_ME_hist->Fill(L_vector_ME.at(iLambda1).Eta() , L_L_L2_vector_ME_SE.at(iLambda2).Eta());
      L0_L0_phi1_vs_phi2_ME_hist->Fill(L_vector_ME.at(iLambda1).Phi() , L_L_L2_vector_ME_SE.at(iLambda2).Phi());
      L0_L0_pT1_vs_pT2_ME_hist->Fill(L_vector_ME.at(iLambda1).Pt() , L_L_L2_vector_ME_SE.at(iLambda2).Pt());

    }

  }

  //---------------------------------------------------------------------------------------------------------------------------

  //L0bar-L0bar
  //Lbar1 is from SE, Lbar2 is from ME
  for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < Lbar_Lbar_pbar1_star_vector_ME_SE.size(); iLambdaBar1++)
  {
    for(unsigned int iLambdaBar2 = 0; iLambdaBar2 < pBar_star_vector_ME.size(); iLambdaBar2++)
    {
      //limit kinematics of ME L2 based on kinematics of same event L2 - double check precision
      if( fabs( Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar1).Eta() - Lbar_vector_ME.at(iLambdaBar2).Eta()) > 0.1 ) continue;
      if( fabs( Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar1).Phi() - Lbar_vector_ME.at(iLambdaBar2).Phi()) > 0.1 ) continue;
      if( fabs( Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar1).Pt() - Lbar_vector_ME.at(iLambdaBar2).Pt()) > 0.1 ) continue;

      double theta_star = Lbar_Lbar_pbar1_star_vector_ME_SE.at(iLambdaBar1).Angle(pBar_star_vector_ME.at(iLambdaBar2));

      L0bar_L0bar_cosThetaProdPlane_ME->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_ME_pT_hist[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE.at(iLambdaBar1)][Lbar_pT_bin_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_ME_eta_hist[Lbar_Lbar_Lbar1_eta_bin_vector_ME_SE.at(iLambdaBar1)][Lbar_eta_bin_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
      
      
      //QA histograms with L-Lbar kinematics
      L0bar_L0bar_eta1_vs_eta2_ME_hist->Fill(Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar1).Eta() , Lbar_vector_ME.at(iLambdaBar2).Eta());
      L0bar_L0bar_phi1_vs_phi2_ME_hist->Fill(Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar1).Phi() , Lbar_vector_ME.at(iLambdaBar2).Phi());
      L0bar_L0bar_pT1_vs_pT2_ME_hist->Fill(Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar1).Pt() , Lbar_vector_ME.at(iLambdaBar2).Pt());

    }

  }

  //Lbar1 is from ME, Lbar2 is from SE
  for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pBar_star_vector_ME.size(); iLambdaBar1++)
  {
    for(unsigned int iLambdaBar2 = 0; iLambdaBar2 < Lbar_Lbar_pbar2_star_vector_ME_SE.size(); iLambdaBar2++)
    {
      //limit kinematics of ME L2 based on kinematics of same event L2 - double check precision
      if( fabs( Lbar_vector_ME.at(iLambdaBar1).Eta() - Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar2).Eta()) > 0.1 ) continue;
      if( fabs( Lbar_vector_ME.at(iLambdaBar1).Phi() - Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar2).Phi()) > 0.1 ) continue;
      if( fabs( Lbar_vector_ME.at(iLambdaBar1).Pt() - Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar2).Pt()) > 0.1 ) continue;

      double theta_star = pBar_star_vector_ME.at(iLambdaBar1).Angle(Lbar_Lbar_pbar2_star_vector_ME_SE.at(iLambdaBar2));

      L0_L0_cosThetaProdPlane_ME->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_ME_pT_hist[Lbar_pT_bin_vector_ME.at(iLambdaBar1)][Lbar_Lbar_Lbar2_pT_bin_vector_ME_SE.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_ME_eta_hist[Lbar_eta_bin_vector_ME.at(iLambdaBar1)][Lbar_Lbar_Lbar2_eta_bin_vector_ME_SE.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
      
      //QA histograms with L-Lbar kinematics
      L0bar_L0bar_eta1_vs_eta2_ME_hist->Fill(Lbar_vector_ME.at(iLambdaBar1).Eta() , Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar2).Eta());
      L0bar_L0bar_phi1_vs_phi2_ME_hist->Fill(Lbar_vector_ME.at(iLambdaBar1).Phi() , Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar2).Phi());
      L0bar_L0bar_pT1_vs_pT2_ME_hist->Fill(Lbar_vector_ME.at(iLambdaBar1).Pt() , Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar2).Pt());

    }

  }
  //_____________________________________________________________________________________________________________________________________________________________________
  

  //mixed event after cuts
  
  //L-Lbar
  //L from SE mixed with single Lbar 
  for(unsigned int iLambda = 0; iLambda < L_Lbar_p_star_cuts_vector_ME_SE.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_cuts_vector_ME.size(); iLambdaBar++)
    {
      //limit kinematics of ME Lbar based on kinematics of same event Lbar - double check precision      
      if( fabs( L_Lbar_Lbar_cuts_vector_ME_SE.at(iLambda).Eta() - Lbar_cuts_vector_ME.at(iLambdaBar).Eta()) > 0.1 ) continue;
      if( fabs( L_Lbar_Lbar_cuts_vector_ME_SE.at(iLambda).Phi() - Lbar_cuts_vector_ME.at(iLambdaBar).Phi()) > 0.1 ) continue;
      if( fabs( L_Lbar_Lbar_cuts_vector_ME_SE.at(iLambda).Pt() - Lbar_cuts_vector_ME.at(iLambdaBar).Pt()) > 0.1 ) continue;
    
            
      double theta_star = L_Lbar_p_star_cuts_vector_ME_SE.at(iLambda).Angle(pBar_star_cuts_vector_ME.at(iLambdaBar));
           
      
      L0_L0bar_cosThetaProdPlane_ME_cuts->Fill(TMath::Cos(theta_star));          
      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[L_Lbar_L_pT_bin_cuts_vector_ME_SE.at(iLambda)][Lbar_pT_bin_cuts_vector_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[L_Lbar_L_eta_bin_cuts_vector_ME_SE.at(iLambda)][Lbar_eta_bin_cuts_vector_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      
      
      //QA histograms with L-Lbar kinematics
      L0_L0bar_eta1_vs_eta2_ME_cuts_hist->Fill(L_Lbar_L_cuts_vector_ME_SE.at(iLambda).Eta() , Lbar_cuts_vector_ME.at(iLambdaBar).Eta());
      L0_L0bar_phi1_vs_phi2_ME_cuts_hist->Fill(L_Lbar_L_cuts_vector_ME_SE.at(iLambda).Phi() , Lbar_cuts_vector_ME.at(iLambdaBar).Phi());
      L0_L0bar_pT1_vs_pT2_ME_cuts_hist->Fill(L_Lbar_L_cuts_vector_ME_SE.at(iLambda).Pt() , Lbar_cuts_vector_ME.at(iLambdaBar).Pt());
      
      //decay pi kinematics
      L0_L0bar_pi_pT1_vs_pi_pT2_ME_cuts_hist->Fill(L_Lbar_L_pi_pT_cuts_vector_ME_SE.at(iLambda), Lbar_pi_pT_cuts_vector_ME.at(iLambdaBar));

    }  
  
  }
  
  //Lbar from SE mixed with single L
  for(unsigned int iLambda = 0; iLambda < p_star_cuts_vector_ME.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < L_Lbar_pbar_star_cuts_vector_ME_SE.size(); iLambdaBar++)
    {
      //limit kinematics of ME L based on kinematics of same event L - double check precision      
      if( fabs( L_cuts_vector_ME.at(iLambda).Eta() - L_Lbar_L_cuts_vector_ME_SE.at(iLambdaBar).Eta() ) > 0.1 ) continue;
      if( fabs( L_cuts_vector_ME.at(iLambda).Phi() - L_Lbar_L_cuts_vector_ME_SE.at(iLambdaBar).Phi() ) > 0.1 ) continue;
      if( fabs( L_cuts_vector_ME.at(iLambda).Pt() -  L_Lbar_L_cuts_vector_ME_SE.at(iLambdaBar).Pt() ) > 0.1 ) continue;
    
            
      double theta_star = p_star_cuts_vector_ME.at(iLambda).Angle(L_Lbar_pbar_star_cuts_vector_ME_SE.at(iLambdaBar));
           
      
      L0_L0bar_cosThetaProdPlane_ME_cuts->Fill(TMath::Cos(theta_star));          
      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[L_pT_bin_cuts_vector_ME.at(iLambda)][L_Lbar_Lbar_pT_bin_cuts_vector_ME_SE.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[L_eta_bin_cuts_vector_ME.at(iLambda)][L_Lbar_Lbar_eta_bin_cuts_vector_ME_SE.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      
      
      //QA histograms with L-Lbar kinematics
      L0_L0bar_eta1_vs_eta2_ME_cuts_hist->Fill(L_cuts_vector_ME.at(iLambda).Eta() , L_Lbar_Lbar_cuts_vector_ME_SE.at(iLambdaBar).Eta());
      L0_L0bar_phi1_vs_phi2_ME_cuts_hist->Fill(L_cuts_vector_ME.at(iLambda).Phi() , L_Lbar_Lbar_cuts_vector_ME_SE.at(iLambdaBar).Phi());
      L0_L0bar_pT1_vs_pT2_ME_cuts_hist->Fill(L_cuts_vector_ME.at(iLambda).Pt() , L_Lbar_Lbar_cuts_vector_ME_SE.at(iLambdaBar).Pt());
      
      //decay pi kinematics
      L0_L0bar_pi_pT1_vs_pi_pT2_ME_cuts_hist->Fill(L_pi_pT_cuts_vector_ME.at(iLambda), L_Lbar_Lbar_pi_pT_cuts_vector_ME_SE.at(iLambdaBar));

    }  
  
  }

  //---------------------------------------------------------------------------------------------------------------------------

  //L0-L0
  //L1 is from SE, L2 is from ME
  for(unsigned int iLambda1 = 0; iLambda1 < L_L_p1_star_cuts_vector_ME_SE.size(); iLambda1++)
  {
    for(unsigned int iLambda2 = 0; iLambda2 < p_star_cuts_vector_ME.size(); iLambda2++)
    {
      //limit kinematics of ME L2 based on kinematics of same event L2 - double check precision      
      if( fabs( L_L_L2_cuts_vector_ME_SE.at(iLambda1).Eta() - L_cuts_vector_ME.at(iLambda2).Eta()) > 0.1 ) continue;
      if( fabs( L_L_L2_cuts_vector_ME_SE.at(iLambda1).Phi() - L_cuts_vector_ME.at(iLambda2).Phi()) > 0.1 ) continue;
      if( fabs( L_L_L2_cuts_vector_ME_SE.at(iLambda1).Pt() - L_cuts_vector_ME.at(iLambda2).Pt()) > 0.1 ) continue;
    
      double theta_star = L_L_p1_star_cuts_vector_ME_SE.at(iLambda1).Angle(p_star_cuts_vector_ME.at(iLambda2)); 
      
      L0_L0_cosThetaProdPlane_ME_cuts->Fill(TMath::Cos(theta_star));          
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[L_L_L1_pT_bin_cuts_vector_ME_SE.at(iLambda1)][L_pT_bin_cuts_vector_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[L_L_L1_eta_bin_cuts_vector_ME_SE.at(iLambda1)][L_eta_bin_cuts_vector_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star));   
      
      //QA histograms with L-Lbar kinematics
      L0_L0_eta1_vs_eta2_ME_cuts_hist->Fill(L_L_L1_cuts_vector_ME_SE.at(iLambda1).Eta() , L_cuts_vector_ME.at(iLambda2).Eta());
      L0_L0_phi1_vs_phi2_ME_cuts_hist->Fill(L_L_L1_cuts_vector_ME_SE.at(iLambda1).Phi() , L_cuts_vector_ME.at(iLambda2).Phi());
      L0_L0_pT1_vs_pT2_ME_cuts_hist->Fill(L_L_L1_cuts_vector_ME_SE.at(iLambda1).Pt() , L_cuts_vector_ME.at(iLambda2).Pt());   
      
    }   
  
  }    
  
  //L1 is from ME, L2 is from SE
  for(unsigned int iLambda1 = 0; iLambda1 < p_star_cuts_vector_ME.size(); iLambda1++)
  {
    for(unsigned int iLambda2 = 0; iLambda2 < L_L_p2_star_cuts_vector_ME_SE.size(); iLambda2++)
    {
      //limit kinematics of ME L2 based on kinematics of same event L2 - double check precision      
      if( fabs( L_cuts_vector_ME.at(iLambda1).Eta() - L_L_L1_cuts_vector_ME_SE.at(iLambda2).Eta()) > 0.1 ) continue;
      if( fabs( L_cuts_vector_ME.at(iLambda1).Phi() - L_L_L1_cuts_vector_ME_SE.at(iLambda2).Phi()) > 0.1 ) continue;
      if( fabs( L_cuts_vector_ME.at(iLambda1).Pt() - L_L_L1_cuts_vector_ME_SE.at(iLambda2).Pt()) > 0.1 ) continue;
    
      double theta_star = p_star_cuts_vector_ME.at(iLambda1).Angle(L_L_p2_star_cuts_vector_ME_SE.at(iLambda2)); 
      
      L0_L0_cosThetaProdPlane_ME_cuts->Fill(TMath::Cos(theta_star));          
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[L_pT_bin_cuts_vector_ME.at(iLambda1)][L_L_L2_pT_bin_cuts_vector_ME_SE.at(iLambda2)]->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[L_eta_bin_cuts_vector_ME.at(iLambda1)][L_L_L2_eta_bin_cuts_vector_ME_SE.at(iLambda2)]->Fill(TMath::Cos(theta_star));  
      
      //QA histograms with L-Lbar kinematics
      L0_L0_eta1_vs_eta2_ME_cuts_hist->Fill(L_cuts_vector_ME.at(iLambda1).Eta() , L_L_L2_cuts_vector_ME_SE.at(iLambda2).Eta());
      L0_L0_phi1_vs_phi2_ME_cuts_hist->Fill(L_cuts_vector_ME.at(iLambda1).Phi() , L_L_L2_cuts_vector_ME_SE.at(iLambda2).Phi());
      L0_L0_pT1_vs_pT2_ME_cuts_hist->Fill(L_cuts_vector_ME.at(iLambda1).Pt() , L_L_L2_cuts_vector_ME_SE.at(iLambda2).Pt());    
      
    }   
  
  }

  //---------------------------------------------------------------------------------------------------------------------------

  //L0bar-L0bar
  //L1 is from SE, L2 is from ME
  for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < Lbar_Lbar_pbar1_star_cuts_vector_ME_SE.size(); iLambdaBar1++)
  {
    for(unsigned int iLambdaBar2 = 0; iLambdaBar2 < pBar_star_cuts_vector_ME.size(); iLambdaBar2++)
    {
      //limit kinematics of ME L2 based on kinematics of same event L2 - double check precision      
      if( fabs( Lbar_Lbar_Lbar2_cuts_vector_ME_SE.at(iLambdaBar1).Eta() - Lbar_cuts_vector_ME.at(iLambdaBar2).Eta()) > 0.1 ) continue;
      if( fabs( Lbar_Lbar_Lbar2_cuts_vector_ME_SE.at(iLambdaBar1).Phi() - Lbar_cuts_vector_ME.at(iLambdaBar2).Phi()) > 0.1 ) continue;
      if( fabs( Lbar_Lbar_Lbar2_cuts_vector_ME_SE.at(iLambdaBar1).Pt() - Lbar_cuts_vector_ME.at(iLambdaBar2).Pt()) > 0.1 ) continue;
    
      double theta_star = Lbar_Lbar_pbar1_star_cuts_vector_ME_SE.at(iLambdaBar1).Angle(pBar_star_cuts_vector_ME.at(iLambdaBar2)); 
      
      L0bar_L0bar_cosThetaProdPlane_ME_cuts->Fill(TMath::Cos(theta_star));          
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[Lbar_Lbar_Lbar1_pT_bin_cuts_vector_ME_SE.at(iLambdaBar1)][Lbar_pT_bin_cuts_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[Lbar_Lbar_Lbar1_eta_bin_cuts_vector_ME_SE.at(iLambdaBar1)][Lbar_eta_bin_cuts_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));      
      
      //QA histograms with L-Lbar kinematics
      L0bar_L0bar_eta1_vs_eta2_ME_cuts_hist->Fill(Lbar_Lbar_Lbar1_cuts_vector_ME_SE.at(iLambdaBar1).Eta() , Lbar_cuts_vector_ME.at(iLambdaBar2).Eta());
      L0bar_L0bar_phi1_vs_phi2_ME_cuts_hist->Fill(Lbar_Lbar_Lbar1_cuts_vector_ME_SE.at(iLambdaBar1).Phi() , Lbar_cuts_vector_ME.at(iLambdaBar2).Phi());
      L0bar_L0bar_pT1_vs_pT2_ME_cuts_hist->Fill(Lbar_Lbar_Lbar1_cuts_vector_ME_SE.at(iLambdaBar1).Pt() , Lbar_cuts_vector_ME.at(iLambdaBar2).Pt());
      
    }   
  
  }    
  
  //L1 is from ME, L2 is from SE
  for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pBar_star_cuts_vector_ME.size(); iLambdaBar1++)
  {
    for(unsigned int iLambdaBar2 = 0; iLambdaBar2 < Lbar_Lbar_pbar2_star_cuts_vector_ME_SE.size(); iLambdaBar2++)
    {
      //limit kinematics of ME L2 based on kinematics of same event L2 - double check precision      
      if( fabs( Lbar_cuts_vector_ME.at(iLambdaBar1).Eta() - Lbar_Lbar_Lbar1_cuts_vector_ME_SE.at(iLambdaBar2).Eta()) > 0.1 ) continue;
      if( fabs( Lbar_cuts_vector_ME.at(iLambdaBar1).Phi() - Lbar_Lbar_Lbar1_cuts_vector_ME_SE.at(iLambdaBar2).Phi()) > 0.1 ) continue;
      if( fabs( Lbar_cuts_vector_ME.at(iLambdaBar1).Pt() - Lbar_Lbar_Lbar1_cuts_vector_ME_SE.at(iLambdaBar2).Pt()) > 0.1 ) continue;
    
      double theta_star = pBar_star_cuts_vector_ME.at(iLambdaBar1).Angle(Lbar_Lbar_pbar2_star_cuts_vector_ME_SE.at(iLambdaBar2)); 
      
      L0bar_L0bar_cosThetaProdPlane_ME_cuts->Fill(TMath::Cos(theta_star));          
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[Lbar_pT_bin_cuts_vector_ME.at(iLambdaBar1)][Lbar_Lbar_Lbar2_pT_bin_cuts_vector_ME_SE.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[Lbar_eta_bin_cuts_vector_ME.at(iLambdaBar1)][Lbar_Lbar_Lbar2_eta_bin_cuts_vector_ME_SE.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));  
      
      //QA histograms with L-Lbar kinematics
      L0bar_L0bar_eta1_vs_eta2_ME_cuts_hist->Fill(Lbar_cuts_vector_ME.at(iLambdaBar1).Eta() , Lbar_Lbar_Lbar2_cuts_vector_ME_SE.at(iLambdaBar2).Eta());
      L0bar_L0bar_phi1_vs_phi2_ME_cuts_hist->Fill(Lbar_cuts_vector_ME.at(iLambdaBar1).Phi() , Lbar_Lbar_Lbar2_cuts_vector_ME_SE.at(iLambdaBar2).Phi());
      L0bar_L0bar_pT1_vs_pT2_ME_cuts_hist->Fill(Lbar_cuts_vector_ME.at(iLambdaBar1).Pt() , Lbar_Lbar_Lbar2_cuts_vector_ME_SE.at(iLambdaBar2).Pt());    
      
    }   
  
  }
   
  //_____________________________________________________________________________________________________________________________________________________________________
  
  
  outFile->cd();
  outFile->Write();
  outFile->Close();

  
  //outFile->Close();

  // Done.
  return;
}
