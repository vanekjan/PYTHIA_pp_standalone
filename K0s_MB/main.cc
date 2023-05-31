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

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

// ROOT, for histogramming.
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TTree.h"
//#include "TChain.h"
#include "TLorentzVector.h"
#include "TF1.h"
//#include "TRandom3.h"
//#include "TString.h"
//#include "TMath.h"


// ROOT, for saving file.
//#include "TFile.h"

using namespace Pythia8;

using namespace std;

int findBinPt( TLorentzVector FourMom, float const *bins, const int nBins )
{
  int bin = -1;

  //find pT bin of Lambda
  for(int j = 0; j < nBins; j++) //loop over pT bins
  {
    if(FourMom.Pt() > bins[j] && FourMom.Pt() <= bins[j+1])
    {
      bin = j;
      break; //stop after bin is found
    }
  }
  
  return bin;
}

int findBinEta( TLorentzVector FourMom, float const *bins, const int nBins )
{
  int bin = -1;

  //find pT bin of Lambda
  for(int j = 0; j < nBins; j++) //loop over pT bins
  {
    if(FourMom.Eta() > bins[j] && FourMom.Eta() <= bins[j+1])
    {
      bin = j;
      break; //stop after bin is found
    }
  }
  
  return bin;
}


int main(int argc, char *argv[])
{

  istringstream nEventsStream(argv[1]);
  int nEvents;
  if (!(nEventsStream >> nEvents))
  {
    cout<<"Invalid first argument!"<<endl;
    return 0;
  }

  istringstream mEnergyStream(argv[2]);
  int mEnergy;
  if (!(mEnergyStream >> mEnergy))
  {
    cout<<"Invalid second argument!"<<endl;
    return 0;
  }
  
  

  cout<<"Number of events: "<<nEvents<<endl;
  cout<<"Energy: "<<mEnergy<<endl;
  cout<<"Output file: "<<argv[3]<<endl;

  
  // Create Pythia instance and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
  //pythia.readString("Beams:eCM = 510"); //beam energy in GeV
  if( mEnergy == 510) pythia.readString("Beams:eCM = 510"); //beam energy in GeV
  else if( mEnergy == 200) pythia.readString("Beams:eCM = 200");
  else
  {
    cout<<"Invalid input energy!"<<endl;
    return 0;
  
  }
  pythia.readString("SoftQCD:nonDiffractive = on"); //equivalent to MB
  pythia.readString("310:onMode=0");
  pythia.readString("310:OnIfMatch=211 -211");
  pythia.readString("-310:onMode=0");
  pythia.readString("-310:OnIfMatch=-211 211");
  
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0"); //should set random seed at start, after calling init()
  
  pythia.init();

  // Create file on which histogram(s) can be saved.
  TFile* outFile = new TFile(argv[3], "RECREATE");

  // Book histogram.
  //TH1F *mult = new TH1F("mult","charged multiplicity", 100, -0.5, 799.5);
  
  const int nPtBins = 8;
  float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5.};
  
  const int nPtBins_corr = 2;
  float const pT_bins_corr[nPtBins_corr+1] = { 0.5, 1.5, 5.};

  const int nEtaBins = 3;
  float const eta_bins[nEtaBins+1] = { -1, -0.2, 0.2, 1 };

  const int K0sPDGid = 310;
  const int K0sbarPDGid = -310;


  //variables for event and particle loop
  TLorentzVector K0_fourmom;
  TLorentzVector K0_fourmom_reverse; //for proton boost
  
  TVector3 K0s_production_vertex;
  TVector3 K0s_decay_vertex;

  TLorentzVector pi1_fourmom;
  TLorentzVector pi2_fourmom;


  //histograms
  
  //true MC
  
  //before cuts
  TH1D *K0s_K0s_cosThetaProdPlane = new TH1D("K0s_K0s_cosThetaProdPlane", "K0s_K0s_cosThetaProdPlane", 10, -1, 1);
  
  TH1D *K0s_K0s_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];
  
  
  //after cuts  
  TH1D *K0s_K0s_cosThetaProdPlane_cuts = new TH1D("K0s_K0s_cosThetaProdPlane_cuts", "K0s_K0s_cosThetaProdPlane_cuts", 10, -1, 1);
  
  TH1D *K0s_K0s_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];
  
  //--------------------------------------------------------------------------------------------------------------------------------------
  
  //true MC mixed event
  //before cuts
  TH1D *K0s_K0s_cosThetaProdPlane_ME = new TH1D("K0s_K0s_cosThetaProdPlane_ME", "K0s_K0s_cosThetaProdPlane_ME", 10, -1, 1);
  
  TH1D *K0s_K0s_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];
  
  
  //after cuts  
  TH1D *K0s_K0s_cosThetaProdPlane_ME_cuts = new TH1D("K0s_K0s_cosThetaProdPlane_ME_cuts", "K0s_K0s_cosThetaProdPlane_ME_cuts", 10, -1, 1);
  
  TH1D *K0s_K0s_cosThetaProdPlane_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_ME_eta_cuts_hist[nEtaBins][nEtaBins];
  
  //_____________________________________________________________________________________________________________________________________________________________________ 
  
  //pi pi pairs to simulate combinatorial background
  
  //US matched to US before cuts
  TH1D *K0s_K0s_cosThetaProdPlane_US = new TH1D("K0s_K0s_cosThetaProdPlane_US", "K0s_K0s_cosThetaProdPlane_US", 10, -1, 1);
  
  TH1D *K0s_K0s_cosThetaProdPlane_US_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_US_eta_hist[nEtaBins][nEtaBins];
  
  
  //after cuts  
  TH1D *K0s_K0s_cosThetaProdPlane_US_cuts = new TH1D("K0s_K0s_cosThetaProdPlane_US_cuts", "K0s_K0s_cosThetaProdPlane_US_cuts", 10, -1, 1);
  
  TH1D *K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[nEtaBins][nEtaBins];
  
  //--------------------------------------------------------------------------------------------------------------------------------------
  
  //mixed event
  //US matched to US before cuts
  TH1D *K0s_K0s_cosThetaProdPlane_US_ME = new TH1D("K0s_K0s_cosThetaProdPlane_US_ME", "K0s_K0s_cosThetaProdPlane_US_ME", 10, -1, 1);
  
  TH1D *K0s_K0s_cosThetaProdPlane_US_ME_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_US_ME_eta_hist[nEtaBins][nEtaBins];
  
  
  //after cuts  
  TH1D *K0s_K0s_cosThetaProdPlane_US_ME_cuts = new TH1D("K0s_K0s_cosThetaProdPlane_US_ME_cuts", "K0s_K0s_cosThetaProdPlane_US_ME_cuts", 10, -1, 1);
  
  TH1D *K0s_K0s_cosThetaProdPlane_US_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_US_ME_eta_cuts_hist[nEtaBins][nEtaBins];
  
  //-------------------------------------------------------------------------------------------------------------------------------------
  
  //US matched to LS before cuts
  TH1D *K0s_K0s_cosThetaProdPlane_US_LS = new TH1D("K0s_K0s_cosThetaProdPlane_US_LS", "K0s_K0s_cosThetaProdPlane_US_LS", 10, -1, 1);
  
  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_pT_hist_tight_eta[nPtBins_corr][nPtBins_corr];

  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_eta_hist[nEtaBins][nEtaBins];
  
  
  //after cuts  
  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_cuts = new TH1D("K0s_K0s_cosThetaProdPlane_US_LS_cuts", "K0s_K0s_cosThetaProdPlane_US_LS_cuts", 10, -1, 1);
  
  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_pT_cuts_hist_tight_eta[nPtBins_corr][nPtBins_corr];

  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_eta_cuts_hist[nEtaBins][nEtaBins];
  
  //--------------------------------------------------------------------------------------------------------------------------------------
  
  //mixed event
  //US matched to LS before cuts
  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_ME = new TH1D("K0s_K0s_cosThetaProdPlane_US_LS_ME", "K0s_K0s_cosThetaProdPlane_US_LS_ME", 10, -1, 1);
  
  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_ME_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_ME_eta_hist[nEtaBins][nEtaBins];
  
  
  //after cuts  
  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_ME_cuts = new TH1D("K0s_K0s_cosThetaProdPlane_US_LS_ME_cuts", "K0s_K0s_cosThetaProdPlane_US_LS_ME_cuts", 10, -1, 1);
  
  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[nEtaBins][nEtaBins];
  
  //_____________________________________________________________________________________________________________________________________________________________________ 
  
  
  
  //L-L correlation histograms in bins
  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      //true MC
      //before cuts
      K0s_K0s_cosThetaProdPlane_pT_hist[pTbin1][pTbin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
      //after cuts
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
      //---------------------------------------------------------------------------------
      
      //true MC mixed event
      //before cuts
      K0s_K0s_cosThetaProdPlane_ME_pT_hist[pTbin1][pTbin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
      //after cuts
      K0s_K0s_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
      //_________________________________________________________________________________
      
      //pi pi pairs to simulate combinatorial background
      //US matched to US
      //before cuts
      K0s_K0s_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
      //after cuts
      K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_US_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_US_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
      //---------------------------------------------------------------------------------
      
      //mixed event
      //before cuts
      K0s_K0s_cosThetaProdPlane_US_ME_pT_hist[pTbin1][pTbin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_US_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_US_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
      //after cuts
      K0s_K0s_cosThetaProdPlane_US_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_US_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_US_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
      //---------------------------------------------------------------------------------
      
      //US matched to LS
      //before cuts
      K0s_K0s_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_US_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_US_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
      //after cuts
      K0s_K0s_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_US_LS_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_US_LS_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
      //---------------------------------------------------------------------------------
      
      //mixed event
      //before cuts
      K0s_K0s_cosThetaProdPlane_US_LS_ME_pT_hist[pTbin1][pTbin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_US_LS_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_US_LS_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
      //after cuts
      K0s_K0s_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_US_LS_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_US_LS_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);      
      

    }
  }

  for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
  {
    for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
    {
      //true MC
      //before cuts
      K0s_K0s_cosThetaProdPlane_eta_hist[etaBin1][etaBin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      
      //after cuts
      K0s_K0s_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      
      //---------------------------------------------------------------------------------
      
      //true MC mixed event
      //before cuts
      K0s_K0s_cosThetaProdPlane_ME_eta_hist[etaBin1][etaBin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      
      //after cuts
      K0s_K0s_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      
      //_________________________________________________________________________________
      
      //pi pi pairs to simulate combinatorial background
      
      //US matched to US
      //before cuts
      K0s_K0s_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      
      //after cuts
      K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_US_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_US_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      
      //---------------------------------------------------------------------------------
      
      //mixed event
      //before cuts
      K0s_K0s_cosThetaProdPlane_US_ME_eta_hist[etaBin1][etaBin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_US_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_US_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      
      //after cuts
      K0s_K0s_cosThetaProdPlane_US_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_US_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_US_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      
      //---------------------------------------------------------------------------------
      
      //US matched to LS
      //before cuts
      K0s_K0s_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_US_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_US_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      
      //after cuts
      K0s_K0s_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_US_LS_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_US_LS_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      
      //---------------------------------------------------------------------------------
      
      //mixed event
      //before cuts
      K0s_K0s_cosThetaProdPlane_US_LS_ME_eta_hist[etaBin1][etaBin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_US_LS_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_US_LS_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      
      //after cuts
      K0s_K0s_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_US_LS_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_US_LS_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

    }
  }

  //mixed event true K0s
  
  vector<TLorentzVector> pi_star_vector_ME;
  
  vector<int> K0s_pT_bin_vector_ME;
  vector<int> K0s_eta_bin_vector_ME;
  
  
  vector<TLorentzVector> pi_star_cuts_vector_ME;
  
  vector<int> K0s_pT_bin_cuts_vector_ME;
  vector<int> K0s_eta_bin_cuts_vector_ME;
  
  //-----------------------------------------
  
  //vectors for creating US and LS pi pairs
  //pi origins to calculate SV of pairs
  
  vector<TLorentzVector> piPlus_vector_event;
  vector<TLorentzVector> piMinus_vector_event;
  
  vector<TVector3> piPlus_origin_vector_event;
  vector<TVector3> piMinus_origin_vector_event;
  
  //mixed event of pi pairs
  //US
  vector<TLorentzVector> pi_star_vector_US_ME;
  vector<int> K0s_cuts_flag_US_ME;
  
  vector<int> K0s_pT_bin_vector_US_ME;  
  vector<int> K0s_eta_bin_vector_US_ME;
  
  
  //US matched to LS
  vector<TLorentzVector> pi_star_vector_LS_ME;
  vector<int> K0s_cuts_flag_LS_ME;
  
  vector<int> K0s_pT_bin_vector_LS_ME;  
  vector<int> K0s_eta_bin_vector_LS_ME;  
  
  
  // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < nEvents; ++iEvent)
  {
    if (!pythia.next()) continue;
    
    vector<float> K0s_y_vector;
    vector<float> K0s_y_cuts_vector;
    
    vector<int> K0s_pT_bin_vector;
    vector<int> K0s_pT_bin_cuts_vector;
    
    vector<int> K0s_eta_bin_vector;
    vector<int> K0s_eta_bin_cuts_vector;
    
    vector<TLorentzVector> pi_star_vector;
    vector<TLorentzVector> pi_star_cuts_vector;

    //loop over particles in event
    for (int i = 0; i < pythia.event.size(); ++i)
    {
      if( fabs(pythia.event[i].id()) == 211 )
      {
        TLorentzVector pi_fourmom;
        pi_fourmom.SetPxPyPzE(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
        
        TVector3 pi_origin = TVector3(pythia.event[i].xProd(), pythia.event[i].yProd(), pythia.event[i].zProd()); //proton production vertex
        
        if( pythia.event[i].id() > 0 ) 
        {
          piPlus_vector_event.push_back(pi_fourmom);
          piPlus_origin_vector_event.push_back(pi_origin);
        }
        
        if( pythia.event[i].id() < 0 )
        {
          piMinus_vector_event.push_back(pi_fourmom);
          piMinus_origin_vector_event.push_back(pi_origin);
        }      
      }    
    
      if( fabs(pythia.event[i].id()) == K0sPDGid)
      {
        //cout<<"Lambda found!"<<endl;
        
        //find index of decay daughters
        const int daughter1_Id = pythia.event[i].daughter1();
        const int daughter2_Id = pythia.event[i].daughter2();
        
        
        //Lambda fourmomentum
        K0_fourmom.SetPxPyPzE(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
        
        //Lambda MC decay veretex
        K0s_decay_vertex.SetXYZ(pythia.event[i].xDec(), pythia.event[i].yDec(), pythia.event[i].zDec());
        
        //Lambda MC production veretex (PV)
        K0s_production_vertex.SetXYZ(pythia.event[i].xProd(), pythia.event[i].yProd(), pythia.event[i].zProd());
        
        TVector3 K0s_decayK0s_vect = K0s_decay_vertex - K0s_production_vertex;
        float K0s_decayL = K0s_decayK0s_vect.Mag();
        
        if( fabs(K0_fourmom.Rapidity()) >= 1.) continue;
        
        K0_fourmom_reverse.SetPxPyPzE(-pythia.event[i].px(), -pythia.event[i].py(), -pythia.event[i].pz(), pythia.event[i].e());

        pi1_fourmom.SetPxPyPzE(pythia.event[daughter1_Id].px(), pythia.event[daughter1_Id].py(), pythia.event[daughter1_Id].pz(), pythia.event[daughter1_Id].e());
        pi2_fourmom.SetPxPyPzE(pythia.event[daughter2_Id].px(), pythia.event[daughter2_Id].py(), pythia.event[daughter2_Id].pz(), pythia.event[daughter2_Id].e());      

        TLorentzVector pi1_fourmom_star = pi1_fourmom;
        pi1_fourmom_star.Boost(K0_fourmom_reverse.BoostVector());  

        TVector3 beamVector(0.,0.,1.); //unity vector along the beam axis
        TVector3 mProdPlane = beamVector.Cross(K0_fourmom.Vect());
        mProdPlane = ( mProdPlane )*(1./mProdPlane.Mag() );

        float mThetaProdPlane = mProdPlane.Angle(pi1_fourmom_star.Vect());
        
        //fill all histograms for all pT and centrality bins
        int pT_bin = -1;

        //find pT bin of Lambda
        for(int j = 0; j < nPtBins; j++) //loop over pT bins
        {
          if(K0_fourmom.Pt() > pT_bins[j] && K0_fourmom.Pt() <= pT_bins[j+1])
          {
            pT_bin = j;
            break; //stop after pT bin is found
          }
        }

        if( pT_bin == -1 ) continue;
        
        
        int pT_bin_corr = -1;

        //find pT bin of Lambda
        for(int j = 0; j < nPtBins_corr; j++) //loop over pT bins
        {
          if(K0_fourmom.Pt() > pT_bins_corr[j] && K0_fourmom.Pt() <= pT_bins_corr[j+1])
          {
            pT_bin_corr = j;
            break; //stop after pT bin is found
          }
        }
        
        if( pT_bin_corr == -1 ) continue;


        //fill all histograms for all eta and centrality bins
        int eta_bin = -1;

        //find eta bin of Lambda
        for(int j = 0; j < nEtaBins; j++) //loop over eta bins
        {
          if(K0_fourmom.Eta() > eta_bins[j] && K0_fourmom.Eta() <= eta_bins[j+1])
          {
            eta_bin = j;
            break; //stop after eta bin is found
          }
        }

        if( eta_bin == -1 ) continue;
        //_____________________________________________________________________________________________
        
        K0s_y_vector.push_back(K0_fourmom.Rapidity());
        
        K0s_pT_bin_vector.push_back(pT_bin_corr);
        K0s_eta_bin_vector.push_back(eta_bin);
        
        pi_star_vector.push_back(pi1_fourmom_star);
       
        //_____________________________________________________________________________________________
        
        
        //fill histograms after cuts
        
        //K0s cuts
        //decay length cuts in mm (default PYTHIA units)
        if(K0s_decayL < 5 || K0s_decayL > 250) continue;
        
        
        //daughter cuts
        if( fabs(pi1_fourmom.Eta()) >= 1. || fabs(pi2_fourmom.Eta()) >= 1. ) continue;       
        if( pi1_fourmom.Pt() < 0.15 || pi1_fourmom.Pt() > 20. ) continue;
        if( pi2_fourmom.Pt() < 0.15 || pi2_fourmom.Pt() > 20. ) continue;
        
        //add decay length cut on K0s
        
        K0s_y_cuts_vector.push_back(K0_fourmom.Rapidity());
        
        K0s_pT_bin_cuts_vector.push_back(pT_bin_corr);
        K0s_eta_bin_cuts_vector.push_back(eta_bin);
        
        pi_star_cuts_vector.push_back(pi1_fourmom_star);     
        

      }//end if K0s
    
    }//end particle loop
    
    //fill K0s-K0s correlation before cuts
    if(pi_star_vector.size() > 1)
    {
      for(unsigned int iK0s1 = 0; iK0s1 < pi_star_vector.size(); iK0s1++)
      {
        for(unsigned int iK0s2 = iK0s1+1; iK0s2 < pi_star_vector.size(); iK0s2++)
        {
          double theta_star = pi_star_vector.at(iK0s1).Angle(pi_star_vector.at(iK0s2).Vect());
          
          K0s_K0s_cosThetaProdPlane->Fill(TMath::Cos(theta_star));
          K0s_K0s_cosThetaProdPlane_pT_hist[K0s_pT_bin_vector.at(iK0s1)][K0s_pT_bin_vector.at(iK0s2)]->Fill(TMath::Cos(theta_star));
          K0s_K0s_cosThetaProdPlane_eta_hist[K0s_eta_bin_vector.at(iK0s1)][K0s_eta_bin_vector.at(iK0s2)]->Fill(TMath::Cos(theta_star));        
        }
      
      }  
    
    }
    
    //-------------------------------------------------------------------
    
    //mixed event before cuts
    if( pi_star_vector.size() == 1 && pi_star_vector_ME.size() < 1e4)
    {
      pi_star_vector_ME.push_back(pi_star_vector.at(0));

      K0s_pT_bin_vector_ME.push_back(K0s_pT_bin_vector.at(0));
      K0s_eta_bin_vector_ME.push_back(K0s_eta_bin_vector.at(0));
    }

    //_____________________________________________________________________________________________________________________________________________________________________
    
    //fill K0s-K0s correlation after cuts
    if(pi_star_cuts_vector.size() > 1)
    {
      for(unsigned int iK0s1 = 0; iK0s1 < pi_star_cuts_vector.size(); iK0s1++)
      {
        for(unsigned int iK0s2 = iK0s1+1; iK0s2 < pi_star_cuts_vector.size(); iK0s2++)
        {
          double theta_star = pi_star_cuts_vector.at(iK0s1).Angle(pi_star_cuts_vector.at(iK0s2).Vect());
          
          K0s_K0s_cosThetaProdPlane_cuts->Fill(TMath::Cos(theta_star));
          K0s_K0s_cosThetaProdPlane_pT_cuts_hist[K0s_pT_bin_cuts_vector.at(iK0s1)][K0s_pT_bin_cuts_vector.at(iK0s2)]->Fill(TMath::Cos(theta_star));
          K0s_K0s_cosThetaProdPlane_eta_cuts_hist[K0s_eta_bin_cuts_vector.at(iK0s1)][K0s_eta_bin_cuts_vector.at(iK0s2)]->Fill(TMath::Cos(theta_star));        
        }
      
      }  
    
    }
    
    //-------------------------------------------------------------------
    
    //mixed event after cuts
    if( pi_star_cuts_vector.size() == 1 && pi_star_cuts_vector_ME.size() < 1e4)
    {
      pi_star_cuts_vector_ME.push_back(pi_star_cuts_vector.at(0));

      K0s_pT_bin_cuts_vector_ME.push_back(K0s_pT_bin_cuts_vector.at(0));
      K0s_eta_bin_cuts_vector_ME.push_back(K0s_eta_bin_cuts_vector.at(0));
    }
    
    //_____________________________________________________________________________________________________________________________________________________________________
    
    //clear vectors
    
    K0s_y_vector.clear();
        
    K0s_pT_bin_vector.clear();
    K0s_eta_bin_vector.clear();
    
    pi_star_vector.clear();
    
    
    K0s_y_cuts_vector.clear();
        
    K0s_pT_bin_cuts_vector.clear();
    K0s_eta_bin_cuts_vector.clear();
    
    pi_star_cuts_vector.clear();
    
    //_____________________________________________________________________________________________________________________________________________________________________
    
    //analyze pi pairs
    
    //L and Lbar "candidates" from pi p pairs
    vector<TLorentzVector> K0s_vector_US;
    vector<int> K0s_pT_bin_vector_US;
    vector<int> K0s_eta_bin_vector_US;
    vector<int> K0s_cuts_flag_US; //flag if pair passed analysis cuts
    vector<TLorentzVector> pi_star_vector_US; //to store p fourmomentum in L (pi p pair) rest frame
    vector<int> pi1_index_vector_US; //to store pi1 index for auto-correlation check for same-sign pairs
    vector<int> pi2_index_vector_US; //to store pi2 index for auto-correlation check for same-sign pairs
    
    
    vector<TLorentzVector> K0s_vector_LS;
    vector<int> K0s_pT_bin_vector_LS;
    vector<int> K0s_eta_bin_vector_LS;
    vector<int> K0s_cuts_flag_LS; //flag if pair passed analysis cuts
    vector<TLorentzVector> pi_star_vector_LS;
    vector<int> pi1_index_vector_LS; //to store pi1 index for auto-correlation check for same-sign pairs
    vector<int> pi2_index_vector_LS; //to store pi2 index for auto-correlation check for same-sign pairs 
    
    
    //pair pions
    
    for( unsigned int pi1_index = 0; pi1_index < piPlus_vector_event.size(); pi1_index++ )
    {
      //US
      for( unsigned int pi2_index = 0; pi2_index < piMinus_vector_event.size(); pi2_index++ )
      {
        TLorentzVector K0s_fourmom_US_event = piPlus_vector_event.at(pi1_index) + piMinus_vector_event.at(pi2_index);
        
        if( fabs( K0s_fourmom_US_event.Rapidity() ) > 1 ) continue;
        
        int pT_bin_K0s_US = findBinPt( K0s_fourmom_US_event, pT_bins_corr, nPtBins_corr );
        int eta_bin_K0s_US = findBinEta( K0s_fourmom_US_event, eta_bins, nEtaBins );
        
        if( pT_bin_K0s_US < 0 || eta_bin_K0s_US < 0 ) continue;        
        
        TLorentzVector K0s_fourmom_reverse_US_event;
        K0s_fourmom_reverse_US_event.SetPxPyPzE(-K0s_fourmom_US_event.Px(), -K0s_fourmom_US_event.Py(), -K0s_fourmom_US_event.Pz(), K0s_fourmom_US_event.E());
        
        TLorentzVector pi_star_US = piPlus_vector_event.at(pi1_index);
        pi_star_US.Boost(K0s_fourmom_reverse_US_event.BoostVector()); //boost p 4-momntum to L rest frame
        
        
        //save info (before cuts)
        K0s_vector_US.push_back(K0s_fourmom_US_event);
        K0s_pT_bin_vector_US.push_back(pT_bin_K0s_US);
        K0s_eta_bin_vector_US.push_back(eta_bin_K0s_US);
        pi_star_vector_US.push_back(pi_star_US);
        
        pi1_index_vector_US.push_back(pi1_index);
        pi2_index_vector_US.push_back(pi2_index);
        
        
        //cuts        
        int cuts_flag = 1;
        //decay length cuts in mm (default PYTHIA units)
        TVector3 K0s_pairVertex_US = ( piPlus_origin_vector_event.at(pi1_index) + piMinus_origin_vector_event.at(pi2_index) )*0.5; //secondary vertex of pi p pair
        float K0s_decayK0s_US = K0s_pairVertex_US.Mag(); //decay length of pi p pair assuming PV at (0,0,0)
        
        if(K0s_decayK0s_US < 5 || K0s_decayK0s_US > 250) cuts_flag = 0;
        //if(K0s_fourmom.Pt() < 0.5) continue; //added pT cut - pT integrated in data starts at 0.5 GeV/c
        //if(K0s_fourmom.Rapidity() >= 1) continue;
        
        //daughter cuts
        if( fabs(piPlus_vector_event.at(pi1_index).Eta()) >= 1. || fabs(piMinus_vector_event.at(pi2_index).Eta()) >= 1. )  cuts_flag = 0;       
        if( piPlus_vector_event.at(pi1_index).Pt() < 0.15 || piPlus_vector_event.at(pi1_index).Pt() > 20. )  cuts_flag = 0;
        if( piMinus_vector_event.at(pi2_index).Pt() < 0.15 || piMinus_vector_event.at(pi2_index).Pt() > 20. )  cuts_flag = 0;
        
        if( K0s_fourmom_US_event.M() < 0.488 || K0s_fourmom_US_event.M() > 0.51 ) cuts_flag = 0; //approximate Minv window from data
        
        K0s_cuts_flag_US.push_back(cuts_flag);
      
      }
      
      //LS with pi+ only, pi- is lower
      for( unsigned int pi2_index = pi1_index+1; pi2_index < piPlus_vector_event.size(); pi2_index++ )
      {
        TLorentzVector K0s_fourmom_LS_event = piPlus_vector_event.at(pi1_index) + piPlus_vector_event.at(pi2_index);
        
        if( fabs( K0s_fourmom_LS_event.Rapidity() ) > 1 ) continue;
        
        int pT_bin_K0s_LS = findBinPt( K0s_fourmom_LS_event, pT_bins_corr, nPtBins_corr );
        int eta_bin_K0s_LS = findBinEta( K0s_fourmom_LS_event, eta_bins, nEtaBins );
        
        if( pT_bin_K0s_LS < 0 || eta_bin_K0s_LS < 0 ) continue;        
        
        TLorentzVector K0s_fourmom_reverse_LS_event;
        K0s_fourmom_reverse_LS_event.SetPxPyPzE(-K0s_fourmom_LS_event.Px(), -K0s_fourmom_LS_event.Py(), -K0s_fourmom_LS_event.Pz(), K0s_fourmom_LS_event.E());
        
        TLorentzVector pi_star_LS = piPlus_vector_event.at(pi1_index);
        pi_star_LS.Boost(K0s_fourmom_reverse_LS_event.BoostVector()); //boost p 4-momntum to L rest frame
        
        
        //save info (before cuts)
        K0s_vector_LS.push_back(K0s_fourmom_LS_event);
        K0s_pT_bin_vector_LS.push_back(pT_bin_K0s_LS);
        K0s_eta_bin_vector_LS.push_back(eta_bin_K0s_LS);
        pi_star_vector_LS.push_back(pi_star_LS);
        
        pi1_index_vector_LS.push_back(pi1_index);
        pi2_index_vector_LS.push_back(pi2_index);
        
        
        //cuts        
        int cuts_flag = 1;
        //decay length cuts in mm (default PYTHIA units)
        TVector3 K0s_pairVertex_LS = ( piPlus_origin_vector_event.at(pi1_index) + piPlus_origin_vector_event.at(pi2_index) )*0.5; //secondary vertex of pi p pair
        float K0s_decayK0s_LS = K0s_pairVertex_LS.Mag(); //decay length of pi p pair assuming PV at (0,0,0)
        
        if(K0s_decayK0s_LS < 5 || K0s_decayK0s_LS > 250) cuts_flag = 0;
        //if(K0s_fourmom.Pt() < 0.5) continue; //added pT cut - pT integrated in data starts at 0.5 GeV/c
        //if(K0s_fourmom.Rapidity() >= 1) continue;
        
        //daughter cuts
        if( fabs(piPlus_vector_event.at(pi1_index).Eta()) >= 1. || fabs(piPlus_vector_event.at(pi2_index).Eta()) >= 1. )  cuts_flag = 0;       
        if( piPlus_vector_event.at(pi1_index).Pt() < 0.15 || piPlus_vector_event.at(pi1_index).Pt() > 20. )  cuts_flag = 0;
        if( piPlus_vector_event.at(pi2_index).Pt() < 0.15 || piPlus_vector_event.at(pi2_index).Pt() > 20. )  cuts_flag = 0;
        
        if( K0s_fourmom_LS_event.M() < 0.488 || K0s_fourmom_LS_event.M() > 0.51 ) cuts_flag = 0; //approximate Minv window from data
        
        K0s_cuts_flag_LS.push_back(cuts_flag);
      
      
      }
    
    }
    
    //LS, with pi-
    for( unsigned int pi1_index = 0; pi1_index < piMinus_vector_event.size(); pi1_index++ )
    {
      for( unsigned int pi2_index = pi1_index+1; pi2_index < piMinus_vector_event.size(); pi2_index++ )
      {
        TLorentzVector K0s_fourmom_LS_event = piMinus_vector_event.at(pi1_index) + piMinus_vector_event.at(pi2_index);
        
        if( fabs( K0s_fourmom_LS_event.Rapidity() ) > 1 ) continue;
        
        int pT_bin_K0s_LS = findBinPt( K0s_fourmom_LS_event, pT_bins_corr, nPtBins_corr );
        int eta_bin_K0s_LS = findBinEta( K0s_fourmom_LS_event, eta_bins, nEtaBins );
        
        if( pT_bin_K0s_LS < 0 || eta_bin_K0s_LS < 0 ) continue;        
        
        TLorentzVector K0s_fourmom_reverse_LS_event;
        K0s_fourmom_reverse_LS_event.SetPxPyPzE(-K0s_fourmom_LS_event.Px(), -K0s_fourmom_LS_event.Py(), -K0s_fourmom_LS_event.Pz(), K0s_fourmom_LS_event.E());
        
        TLorentzVector pi_star_LS = piMinus_vector_event.at(pi1_index);
        pi_star_LS.Boost(K0s_fourmom_reverse_LS_event.BoostVector()); //boost p 4-momntum to L rest frame
        
        
        //save info (before cuts)
        K0s_vector_LS.push_back(K0s_fourmom_LS_event);
        K0s_pT_bin_vector_LS.push_back(pT_bin_K0s_LS);
        K0s_eta_bin_vector_LS.push_back(eta_bin_K0s_LS);
        pi_star_vector_LS.push_back(pi_star_LS);
        
        pi1_index_vector_LS.push_back(pi1_index);
        pi2_index_vector_LS.push_back(pi2_index);
        
        
        //cuts        
        int cuts_flag = 1;
        //decay length cuts in mm (default PYTHIA units)
        TVector3 K0s_pairVertex_LS = ( piMinus_origin_vector_event.at(pi1_index) + piMinus_origin_vector_event.at(pi2_index) )*0.5; //secondary vertex of pi p pair
        float K0s_decayK0s_LS = K0s_pairVertex_LS.Mag(); //decay length of pi p pair assuming PV at (0,0,0)
        
        if(K0s_decayK0s_LS < 5 || K0s_decayK0s_LS > 250) cuts_flag = 0;
        //if(K0s_fourmom.Pt() < 0.5) continue; //added pT cut - pT integrated in data starts at 0.5 GeV/c
        //if(K0s_fourmom.Rapidity() >= 1) continue;
        
        //daughter cuts
        if( fabs(piMinus_vector_event.at(pi1_index).Eta()) >= 1. || fabs(piMinus_vector_event.at(pi2_index).Eta()) >= 1. )  cuts_flag = 0;       
        if( piMinus_vector_event.at(pi1_index).Pt() < 0.15 || piMinus_vector_event.at(pi1_index).Pt() > 20. )  cuts_flag = 0;
        if( piMinus_vector_event.at(pi2_index).Pt() < 0.15 || piMinus_vector_event.at(pi2_index).Pt() > 20. )  cuts_flag = 0;
        
        if( K0s_fourmom_LS_event.M() < 0.488 || K0s_fourmom_LS_event.M() > 0.51 ) cuts_flag = 0; //approximate Minv window from data
        
        K0s_cuts_flag_LS.push_back(cuts_flag);
      
      }
      
    }
    
    
    //clear pi vectors for next event
    
    piPlus_vector_event.clear();
    piPlus_origin_vector_event.clear();
    
    piMinus_vector_event.clear();
    piMinus_origin_vector_event.clear();
    
    //_____________________________________________________________________________________________________________________________________________________________________
    
    //pair K0s
    
    //US paired with US
    if( K0s_vector_US.size() > 1 )
    {
      for( unsigned int iK0s1 = 0; iK0s1 < K0s_vector_US.size(); iK0s1++ )
      {
        for( unsigned int iK0s2 = iK0s1+1; iK0s2 < K0s_vector_US.size(); iK0s2++ )
        {
          //check auto-correlation
          if( pi1_index_vector_US.at(iK0s1) == pi1_index_vector_US.at(iK0s2) ) continue;
          if( pi2_index_vector_US.at(iK0s1) == pi2_index_vector_US.at(iK0s2) ) continue;
          
          float theta_star = pi_star_vector_US.at(iK0s1).Angle(pi_star_vector_US.at(iK0s2).Vect());
          
          K0s_K0s_cosThetaProdPlane_US->Fill(TMath::Cos(theta_star));
          K0s_K0s_cosThetaProdPlane_US_pT_hist[K0s_pT_bin_vector_US.at(iK0s1)][K0s_pT_bin_vector_US.at(iK0s2)]->Fill(TMath::Cos(theta_star));
          K0s_K0s_cosThetaProdPlane_US_eta_hist[K0s_eta_bin_vector_US.at(iK0s1)][K0s_eta_bin_vector_US.at(iK0s2)]->Fill(TMath::Cos(theta_star));
          
          if( K0s_cuts_flag_US.at(iK0s1) == 1 && K0s_cuts_flag_US.at(iK0s2) == 1) //both K0s in pair passed cuts
          {
            K0s_K0s_cosThetaProdPlane_US_cuts->Fill(TMath::Cos(theta_star));
            K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[K0s_pT_bin_vector_US.at(iK0s1)][K0s_pT_bin_vector_US.at(iK0s2)]->Fill(TMath::Cos(theta_star));
            K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[K0s_eta_bin_vector_US.at(iK0s1)][K0s_eta_bin_vector_US.at(iK0s2)]->Fill(TMath::Cos(theta_star));          
          }
        
        }
        
      }
    
    }
        
    //-------------------------------------------------
    
    //US paired with LS
    int nFilK0s_K0s_K0s_US_LS = 0;
    
    if( K0s_vector_US.size() > 0 && K0s_vector_LS.size() )
    {
      for( unsigned int iK0s1 = 0; iK0s1 < K0s_vector_US.size(); iK0s1++ )
      {
        for( unsigned int iK0s2 = 0; iK0s2 < K0s_vector_LS.size(); iK0s2++ )
        {
          //check auto-correlation
          if( pi1_index_vector_US.at(iK0s1) == pi1_index_vector_LS.at(iK0s2) ) continue;
          if( pi2_index_vector_US.at(iK0s1) == pi2_index_vector_LS.at(iK0s2) ) continue;
          if( pi1_index_vector_US.at(iK0s1) == pi2_index_vector_LS.at(iK0s2) ) continue;
          if( pi2_index_vector_US.at(iK0s1) == pi1_index_vector_LS.at(iK0s2) ) continue;
          
          if( nFilK0s_K0s_K0s_US_LS % 2 == 0 )
          {
            float theta_star = pi_star_vector_US.at(iK0s1).Angle(pi_star_vector_LS.at(iK0s2).Vect());
          
            K0s_K0s_cosThetaProdPlane_US_LS->Fill(TMath::Cos(theta_star));
            K0s_K0s_cosThetaProdPlane_US_LS_pT_hist[K0s_pT_bin_vector_US.at(iK0s1)][K0s_pT_bin_vector_LS.at(iK0s2)]->Fill(TMath::Cos(theta_star));
            K0s_K0s_cosThetaProdPlane_US_LS_eta_hist[K0s_eta_bin_vector_US.at(iK0s1)][K0s_eta_bin_vector_LS.at(iK0s2)]->Fill(TMath::Cos(theta_star));
            
            if( K0s_cuts_flag_US.at(iK0s1) == 1 && K0s_cuts_flag_LS.at(iK0s2) == 1) //both K0s in pair passed cuts
            {
              K0s_K0s_cosThetaProdPlane_US_LS_cuts->Fill(TMath::Cos(theta_star));
              K0s_K0s_cosThetaProdPlane_US_LS_pT_cuts_hist[K0s_pT_bin_vector_US.at(iK0s1)][K0s_pT_bin_vector_LS.at(iK0s2)]->Fill(TMath::Cos(theta_star));
              K0s_K0s_cosThetaProdPlane_US_LS_eta_cuts_hist[K0s_eta_bin_vector_US.at(iK0s1)][K0s_eta_bin_vector_LS.at(iK0s2)]->Fill(TMath::Cos(theta_star));          
            }
          
          }
          else
          {
            float theta_star = pi_star_vector_LS.at(iK0s2).Angle(pi_star_vector_US.at(iK0s1).Vect());
          
            K0s_K0s_cosThetaProdPlane_US_LS->Fill(TMath::Cos(theta_star));
            K0s_K0s_cosThetaProdPlane_US_LS_pT_hist[K0s_pT_bin_vector_LS.at(iK0s2)][K0s_pT_bin_vector_US.at(iK0s1)]->Fill(TMath::Cos(theta_star));
            K0s_K0s_cosThetaProdPlane_US_LS_eta_hist[K0s_eta_bin_vector_LS.at(iK0s2)][K0s_eta_bin_vector_US.at(iK0s1)]->Fill(TMath::Cos(theta_star));
            
            if( K0s_cuts_flag_US.at(iK0s1) == 1 && K0s_cuts_flag_LS.at(iK0s2) == 1) //both K0s in pair passed cuts
            {
              K0s_K0s_cosThetaProdPlane_US_LS_cuts->Fill(TMath::Cos(theta_star));
              K0s_K0s_cosThetaProdPlane_US_LS_pT_cuts_hist[K0s_pT_bin_vector_LS.at(iK0s2)][K0s_pT_bin_vector_US.at(iK0s1)]->Fill(TMath::Cos(theta_star));
              K0s_K0s_cosThetaProdPlane_US_LS_eta_cuts_hist[K0s_eta_bin_vector_LS.at(iK0s2)][K0s_eta_bin_vector_US.at(iK0s1)]->Fill(TMath::Cos(theta_star));          
            }  
          
          }
          
          nFilK0s_K0s_K0s_US_LS++;
          
        }
        
      }    
    
    }
    
    //_____________________________________________________________________________________________________________________________________________________________________
    
    //select K0s for mixed event with pi pairs
    if( K0s_vector_US.size() == 1 && pi_star_vector_US_ME.size() < 1e4)
    {
      pi_star_vector_US_ME.push_back(pi_star_vector_US.at(0));
      K0s_cuts_flag_US_ME.push_back(K0s_cuts_flag_US.at(0));

      K0s_pT_bin_vector_US_ME.push_back(K0s_pT_bin_vector_US.at(0));
      K0s_eta_bin_vector_US_ME.push_back(K0s_eta_bin_vector_US.at(0));     
    }

    //_________________________________________________________________________________________________________
    
    if( K0s_vector_LS.size() == 1 && pi_star_vector_LS_ME.size() < 1e4)
    {
      pi_star_vector_LS_ME.push_back(pi_star_vector_LS.at(0));
      K0s_cuts_flag_LS_ME.push_back(K0s_cuts_flag_LS.at(0));

      K0s_pT_bin_vector_LS_ME.push_back(K0s_pT_bin_vector_LS.at(0));
      K0s_eta_bin_vector_LS_ME.push_back(K0s_eta_bin_vector_LS.at(0));
    }
    

  }//end event loop
  
  //__________________________________________________________________________________________________________________________________________________________________________________
  
  
  //true MC mixed event
  
  //before cuts
  for(unsigned int iK0s1 = 0; iK0s1 < pi_star_vector_ME.size(); iK0s1++)
  {
    for(unsigned int iK0s2 = iK0s1+1; iK0s2 < pi_star_vector_ME.size(); iK0s2++)
    {
      double theta_star = pi_star_vector_ME.at(iK0s1).Angle(pi_star_vector_ME.at(iK0s2).Vect());
      
      K0s_K0s_cosThetaProdPlane_ME->Fill(TMath::Cos(theta_star));
      K0s_K0s_cosThetaProdPlane_ME_pT_hist[K0s_pT_bin_vector_ME.at(iK0s1)][K0s_pT_bin_vector_ME.at(iK0s2)]->Fill(TMath::Cos(theta_star));
      K0s_K0s_cosThetaProdPlane_ME_eta_hist[K0s_eta_bin_vector_ME.at(iK0s1)][K0s_eta_bin_vector_ME.at(iK0s2)]->Fill(TMath::Cos(theta_star));
    }   
  
  }
  
  //-----------------------------------------------------------------------------------------------------------------
  
  //after cuts
  for(unsigned int iK0s1 = 0; iK0s1 < pi_star_cuts_vector_ME.size(); iK0s1++)
  {
    for(unsigned int iK0s2 = iK0s1+1; iK0s2 < pi_star_cuts_vector_ME.size(); iK0s2++)
    {
      double theta_star = pi_star_cuts_vector_ME.at(iK0s1).Angle(pi_star_cuts_vector_ME.at(iK0s2).Vect());
      
      K0s_K0s_cosThetaProdPlane_ME_cuts->Fill(TMath::Cos(theta_star));
      K0s_K0s_cosThetaProdPlane_ME_pT_cuts_hist[K0s_pT_bin_cuts_vector_ME.at(iK0s1)][K0s_pT_bin_cuts_vector_ME.at(iK0s2)]->Fill(TMath::Cos(theta_star));
      K0s_K0s_cosThetaProdPlane_ME_eta_cuts_hist[K0s_eta_bin_cuts_vector_ME.at(iK0s1)][K0s_eta_bin_cuts_vector_ME.at(iK0s2)]->Fill(TMath::Cos(theta_star));
    }   
  
  }
  
  //_________________________________________________________________________________
  
  //mixed event with pi pairs
  
  //before cuts
  //US matched to US
  //K0s-K0s
  for(unsigned int iK0s1 = 0; iK0s1 < pi_star_vector_US_ME.size(); iK0s1++)
  {
    for(unsigned int iK0s2 = iK0s1+1; iK0s2 < pi_star_vector_US_ME.size(); iK0s2++)
    {
      double theta_star = pi_star_vector_US_ME.at(iK0s1).Angle(pi_star_vector_US_ME.at(iK0s2).Vect());
      
      K0s_K0s_cosThetaProdPlane_US_ME->Fill(TMath::Cos(theta_star));
      K0s_K0s_cosThetaProdPlane_US_ME_pT_hist[K0s_pT_bin_vector_US_ME.at(iK0s1)][K0s_pT_bin_vector_US_ME.at(iK0s2)]->Fill(TMath::Cos(theta_star));
      K0s_K0s_cosThetaProdPlane_US_ME_eta_hist[K0s_eta_bin_vector_US_ME.at(iK0s1)][K0s_eta_bin_vector_US_ME.at(iK0s2)]->Fill(TMath::Cos(theta_star));
      
      if( K0s_cuts_flag_US_ME.at(iK0s1) == 1 && K0s_cuts_flag_US_ME.at(iK0s2) == 1) //both K0s in pair passed cuts
      {
        K0s_K0s_cosThetaProdPlane_US_ME_cuts->Fill(TMath::Cos(theta_star));
        K0s_K0s_cosThetaProdPlane_US_ME_pT_cuts_hist[K0s_pT_bin_vector_US_ME.at(iK0s1)][K0s_pT_bin_vector_US_ME.at(iK0s2)]->Fill(TMath::Cos(theta_star));
        K0s_K0s_cosThetaProdPlane_US_ME_eta_cuts_hist[K0s_eta_bin_vector_US_ME.at(iK0s1)][K0s_eta_bin_vector_LS_ME.at(iK0s2)]->Fill(TMath::Cos(theta_star));          
      }
      
    }   
  
  }
  
  //US matched to LS
  
  int nFillK0s_K0s_K0s_ME;
  
  for( unsigned int iK0s1 = 0; iK0s1 < pi_star_vector_US_ME.size(); iK0s1++ )
  {
    for( unsigned int iK0s2 = 0; iK0s2 < pi_star_vector_LS_ME.size(); iK0s2++ )
    {
      
      if( nFillK0s_K0s_K0s_ME % 2 == 0 )
      {
        float theta_star = pi_star_vector_US_ME.at(iK0s1).Angle(pi_star_vector_LS_ME.at(iK0s2).Vect());
      
        K0s_K0s_cosThetaProdPlane_US_LS_ME->Fill(TMath::Cos(theta_star));
        K0s_K0s_cosThetaProdPlane_US_LS_ME_pT_hist[K0s_pT_bin_vector_US_ME.at(iK0s1)][K0s_pT_bin_vector_LS_ME.at(iK0s2)]->Fill(TMath::Cos(theta_star));
        K0s_K0s_cosThetaProdPlane_US_LS_ME_eta_hist[K0s_eta_bin_vector_US_ME.at(iK0s1)][K0s_eta_bin_vector_LS_ME.at(iK0s2)]->Fill(TMath::Cos(theta_star));
        
        if( K0s_cuts_flag_US_ME.at(iK0s1) == 1 && K0s_cuts_flag_LS_ME.at(iK0s2) == 1) //both K0s in pair passed cuts
        {
          K0s_K0s_cosThetaProdPlane_US_LS_ME_cuts->Fill(TMath::Cos(theta_star));
          K0s_K0s_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[K0s_pT_bin_vector_US_ME.at(iK0s1)][K0s_pT_bin_vector_LS_ME.at(iK0s2)]->Fill(TMath::Cos(theta_star));
          K0s_K0s_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[K0s_eta_bin_vector_US_ME.at(iK0s1)][K0s_eta_bin_vector_LS_ME.at(iK0s2)]->Fill(TMath::Cos(theta_star));          
        }
      
      }
      else
      {
        float theta_star = pi_star_vector_LS_ME.at(iK0s2).Angle(pi_star_vector_US_ME.at(iK0s1).Vect());
      
        K0s_K0s_cosThetaProdPlane_US_LS_ME->Fill(TMath::Cos(theta_star));
        K0s_K0s_cosThetaProdPlane_US_LS_ME_pT_hist[K0s_pT_bin_vector_LS_ME.at(iK0s2)][K0s_pT_bin_vector_US_ME.at(iK0s1)]->Fill(TMath::Cos(theta_star));
        K0s_K0s_cosThetaProdPlane_US_LS_ME_eta_hist[K0s_eta_bin_vector_LS_ME.at(iK0s2)][K0s_eta_bin_vector_US_ME.at(iK0s1)]->Fill(TMath::Cos(theta_star));
        
        if( K0s_cuts_flag_US_ME.at(iK0s1) == 1 && K0s_cuts_flag_LS_ME.at(iK0s2) == 1) //both K0s in pair passed cuts
        {
          K0s_K0s_cosThetaProdPlane_US_LS_ME_cuts->Fill(TMath::Cos(theta_star));
          K0s_K0s_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[K0s_pT_bin_vector_LS_ME.at(iK0s2)][K0s_pT_bin_vector_US_ME.at(iK0s1)]->Fill(TMath::Cos(theta_star));
          K0s_K0s_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[K0s_eta_bin_vector_LS_ME.at(iK0s2)][K0s_eta_bin_vector_US_ME.at(iK0s1)]->Fill(TMath::Cos(theta_star));          
        }  
      
      }
      
      nFillK0s_K0s_K0s_ME++;
      
    }
    
  }
  
  //_____________________________________________________________________________________________________________________________________________________________________
  
 
  
  outFile->cd();
  outFile->Write();
  outFile->Close();

  // Statistics on event generation.
  pythia.stat();
  
  //outFile->Close();

  // Done.
  return 0;
}
