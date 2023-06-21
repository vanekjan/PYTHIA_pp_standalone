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
#include <sstream>

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


//int nEvents = 1e2, int mEnergy = 510, Char_t *outputFile="./output/output_Lambda_pp_510_MB_1M_events_work.root")
int main( int argc, char *argv[])
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
  
  //return 0;
  
  
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
  pythia.readString("3122:onMode=0");
  pythia.readString("3122:OnIfMatch=2212 -211");
  pythia.readString("-3122:onMode=0");
  pythia.readString("-3122:OnIfMatch=-2212 211");
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
  float const eta_bins[nEtaBins+1] = { -1, -0.4, 0.4, 1 };

  const int L0PDGid = 3122;
  const int L0barPDGid = -3122;


  //variables for event and particle loop
  TLorentzVector L_fourmom;
  TLorentzVector L_fourmom_reverse; //for proton boost
  
  TVector3 L_production_vertex;
  TVector3 L_decay_vertex;

  TLorentzVector p_fourmom;
  TLorentzVector pi_fourmom;


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
  
  //mixed event histograms
  TH1F *L0_L0bar_cosThetaProdPlane_ME = new TH1F("L0_L0bar_cosThetaProdPlane_ME", "L0_L0bar_cosThetaProdPlane_ME", 10, -1, 1);
  
  TH1F *L0_L0bar_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0bar_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_ME = new TH1F("L0_L0_cosThetaProdPlane_ME", "L0_L0_cosThetaProdPlane_ME", 10, -1, 1);
  
  TH1F *L0_L0_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME = new TH1F("L0bar_L0bar_cosThetaProdPlane_ME", "L0bar_L0bar_cosThetaProdPlane_ME", 10, -1, 1);
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr];
 
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];
  
  
  //--------------------------------------------------------------------------
  
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
  
  //----------------------------------------------------------------------------------------
  
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
    
 
  //_____________________________________________________________________________________________________________________________________________________________________

  //histograms after cuts	  
  TH1F *L0_L0bar_cosThetaProdPlane_cuts = new TH1F("L0_L0bar_cosThetaProdPlane_cuts", "L0_L0bar_cosThetaProdPlane_cuts", 10, -1, 1);
  
  TH1F *L0_L0bar_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0bar_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0_L0_cosThetaProdPlane_cuts = new TH1F("L0_L0_cosThetaProdPlane_cuts", "L0_L0_cosThetaProdPlane_cuts", 10, -1, 1);
  
  TH1F *L0_L0_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_cuts = new TH1F("L0bar_L0bar_cosThetaProdPlane_cuts", "L0bar_L0bar_cosThetaProdPlane_cuts", 10, -1, 1);  
  
  TH1F *L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];
  
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
  
  //--------------------------------------------------------------------------
  
  //p pi pairs - to simulate combinatorial background
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
  
  //----------------------------------------------------------------------------------------
  
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
      
      //mixed event
      L0_L0bar_cosThetaProdPlane_ME_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_ME_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_ME_pT_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
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
      
      //mixed event
      L0_L0bar_cosThetaProdPlane_ME_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_ME_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_ME_eta_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      
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
      
    }
  }
  //_____________________________________________________________________________________________________________________________________________________________________
  
  
  //mixed event vectors, true L and Lbar
 
  vector<TLorentzVector> p_star_vector_ME;
  vector<TLorentzVector> p_star_cuts_vector_ME;
  
  vector<int> L_pT_bin_vector_ME;
  vector<int> L_pT_bin_cuts_vector_ME;
  
  vector<int> L_eta_bin_vector_ME;
  vector<int> L_eta_bin_cuts_vector_ME;

  
  vector<TLorentzVector> pBar_star_vector_ME;
  vector<TLorentzVector> pBar_star_cuts_vector_ME;  
 
  vector<int> Lbar_pT_bin_vector_ME;
  vector<int> Lbar_pT_bin_cuts_vector_ME;
  
  vector<int> Lbar_eta_bin_vector_ME;
  vector<int> Lbar_eta_bin_cuts_vector_ME;
  
  //---------------------------------------
  //vectors for creating US and LS pi p pairs
  //pi and p origins to calculate SV of pairs
  
  vector<TLorentzVector> p_vector_event;
  vector<TVector3> p_origin_vector_event;
  
  vector<TLorentzVector> pi_vector_event;
  vector<TVector3> pi_origin_vector_event;
  
  
  vector<TLorentzVector> pBar_vector_event;
  vector<TVector3> pBar_origin_vector_event;
  
  vector<TLorentzVector> piBar_vector_event;
  vector<TVector3> piBar_origin_vector_event;
  
  //mixed event of pi p
  //US
  vector<TLorentzVector> p_star_vector_US_ME;
  
  vector<int> L_pT_bin_vector_US_ME;  
  vector<int> L_eta_bin_vector_US_ME;

  
  vector<TLorentzVector> pBar_star_vector_US_ME;
 
  vector<int> Lbar_pT_bin_vector_US_ME;  
  vector<int> Lbar_eta_bin_vector_US_ME;
  
  //after cuts
  vector<TLorentzVector> p_star_vector_US_ME_cuts;
  
  vector<int> L_pT_bin_vector_US_ME_cuts;  
  vector<int> L_eta_bin_vector_US_ME_cuts;

  
  vector<TLorentzVector> pBar_star_vector_US_ME_cuts;
 
  vector<int> Lbar_pT_bin_vector_US_ME_cuts;  
  vector<int> Lbar_eta_bin_vector_US_ME_cuts;

  
  
  //US matched to LS
  vector<TLorentzVector> p_star_vector_LS_ME;
  
  vector<int> L_pT_bin_vector_LS_ME;  
  vector<int> L_eta_bin_vector_LS_ME;

  
  vector<TLorentzVector> pBar_star_vector_LS_ME;
 
  vector<int> Lbar_pT_bin_vector_LS_ME;  
  vector<int> Lbar_eta_bin_vector_LS_ME;
  
  //after cuts
  vector<TLorentzVector> p_star_vector_LS_ME_cuts;
  
  vector<int> L_pT_bin_vector_LS_ME_cuts;  
  vector<int> L_eta_bin_vector_LS_ME_cuts;

  
  vector<TLorentzVector> pBar_star_vector_LS_ME_cuts;
 
  vector<int> Lbar_pT_bin_vector_LS_ME_cuts;  
  vector<int> Lbar_eta_bin_vector_LS_ME_cuts;
  

  // Begin event loop. Generate event; skip if generation aborted.
  //for (int iEvent = 0; iEvent < 2e8; ++iEvent)
  for (int iEvent = 0; iEvent < nEvents; ++iEvent)
  //for (int iEvent = 0; iEvent < 1e2; ++iEvent)
  {
    if (!pythia.next()) continue;
    
    //vectors of (anti-)proton 4-momenta boosted into mother rest frame
    vector<TLorentzVector> L_vector;
    vector<TLorentzVector> L_cuts_vector;
    
    vector<int> L_pT_bin_vector;
    vector<int> L_pT_bin_cuts_vector;
    
    vector<int> L_eta_bin_vector;
    vector<int> L_eta_bin_cuts_vector;
    
    vector<TLorentzVector> p_vector;
    vector<TLorentzVector> p_cuts_vector;
    
    vector<TLorentzVector> pi_vector;
    vector<TLorentzVector> pi_cuts_vector;
    
    vector<TLorentzVector> p_star_vector;
    vector<TLorentzVector> p_star_cuts_vector;
    
    
    vector<TLorentzVector> Lbar_vector;
    vector<TLorentzVector> Lbar_cuts_vector;
    
    vector<int> Lbar_pT_bin_vector;
    vector<int> Lbar_pT_bin_cuts_vector;
    
    vector<int> Lbar_eta_bin_vector;
    vector<int> Lbar_eta_bin_cuts_vector;
    
    vector<TLorentzVector> pBar_vector;
    vector<TLorentzVector> pBar_cuts_vector;
    
    vector<TLorentzVector> piBar_vector;
    vector<TLorentzVector> piBar_cuts_vector;
    
    vector<TLorentzVector> pBar_star_vector;
    vector<TLorentzVector> pBar_star_cuts_vector;    
    
    

    //loop over aprticles in event
    for (int i = 0; i < pythia.event.size(); ++i)
    {   
      //store all p and pi from event      
      if( fabs(pythia.event[i].id()) == 2212)
      {

        TLorentzVector p_fourmom;
        p_fourmom.SetPxPyPzE(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
        
        TVector3 p_origin = TVector3(pythia.event[i].xProd(), pythia.event[i].yProd(), pythia.event[i].zProd()); //proton production vertex
        
        if( pythia.event[i].id() > 0 ) 
        {
          p_vector_event.push_back(p_fourmom);
          p_origin_vector_event.push_back(p_origin);        
        }
        
        if( pythia.event[i].id() < 0 )
        {
          pBar_vector_event.push_back(p_fourmom);   
          pBar_origin_vector_event.push_back(p_origin);        
        }
            
      }
      
      if( fabs(pythia.event[i].id()) == 211)
      {

        TLorentzVector pi_fourmom;
        pi_fourmom.SetPxPyPzE(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
        
        TVector3 pi_origin = TVector3(pythia.event[i].xProd(), pythia.event[i].yProd(), pythia.event[i].zProd()); //proton production vertex
        
        if( pythia.event[i].id() > 0 ) 
        {
          pi_vector_event.push_back(pi_fourmom);
          pi_origin_vector_event.push_back(pi_origin);
        }
        
        if( pythia.event[i].id() < 0 )
        {
          piBar_vector_event.push_back(pi_fourmom);
          piBar_origin_vector_event.push_back(pi_origin);
        }
          
      } 
      
      //----------------------------------------------------------------------------    
      
      //true L and Lbar
      if( fabs(pythia.event[i].id()) == 3122)
      {
        //cout<<"Lambda found!"<<endl;
        
        //find index of decay daughters
        const int daughter1_Id = pythia.event[i].daughter1();
        const int daughter2_Id = pythia.event[i].daughter2();
        
        
        //Lambda fourmomentum
        L_fourmom.SetPxPyPzE(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
        
        if(L_fourmom.Rapidity() >= 1) continue; //want only L within |y| < 1
        
        //Lambda MC decay veretex
        L_decay_vertex.SetXYZ(pythia.event[i].xDec(), pythia.event[i].yDec(), pythia.event[i].zDec());
        
        //Lambda MC production veretex (PV)
        //L_production_vertex.SetXYZ(pythia.event[i].xProd(), pythia.event[i].yProd(), pythia.event[i].zProd()); //true MC production vertex - can be different from PV due to secondary decays
        L_production_vertex.SetXYZ(0, 0, 0); //fix production vertex to MC PV, i.e. (0,0,0)
        
        TVector3 L_decayL_vect = L_decay_vertex - L_production_vertex;
        float L_decayL = L_decayL_vect.Mag();
        
        
        float L_xF = fabs(L_fourmom.Pz())/pythia.info.eCM()*2.;
                       

        p_fourmom.SetPxPyPzE(pythia.event[daughter1_Id].px(), pythia.event[daughter1_Id].py(), pythia.event[daughter1_Id].pz(), pythia.event[daughter1_Id].e());
        pi_fourmom.SetPxPyPzE(pythia.event[daughter2_Id].px(), pythia.event[daughter2_Id].py(), pythia.event[daughter2_Id].pz(), pythia.event[daughter2_Id].e());
        
        
        L_fourmom_reverse.SetPxPyPzE(-pythia.event[i].px(), -pythia.event[i].py(), -pythia.event[i].pz(), pythia.event[i].e());

        TLorentzVector p_fourmom_star = p_fourmom;
        p_fourmom_star.Boost(L_fourmom_reverse.BoostVector());       
    
        
        //find bins               
        /*
        int pT_bin_corr = -1;

        //find pT bin of Lambda
        for(int j = 0; j < nPtBins_corr; j++) //loop over pT bins
        {
          if(L_fourmom.Pt() > pT_bins_corr[j] && L_fourmom.Pt() <= pT_bins_corr[j+1])
          {
            pT_bin_corr = j;
            break; //stop after pT bin is found
          }
        }
        */
        
        int pT_bin_corr = findBinPt( L_fourmom, pT_bins_corr, nPtBins_corr );
        
        if( pT_bin_corr == -1 ) continue;


        /*
        int eta_bin = -1;

        //find eta bin of Lambda
        for(int j = 0; j < nEtaBins; j++) //loop over eta bins
        {
          if(L_fourmom.Eta() > eta_bins[j] && L_fourmom.Eta() <= eta_bins[j+1])
          {
            eta_bin = j;
            break; //stop after eta bin is found
          }
        }
        */
        
        int eta_bin = findBinEta( L_fourmom, eta_bins, nEtaBins );

        if( eta_bin == -1 ) continue;


        //Lambda
        if( pythia.event[i].id() > 0 )
        {
          L_vector.push_back(L_fourmom);
          L_pT_bin_vector.push_back(pT_bin_corr);
          L_eta_bin_vector.push_back(eta_bin);
          
          p_vector.push_back(p_fourmom);
          pi_vector.push_back(pi_fourmom);
        
          p_star_vector.push_back(p_fourmom_star);
        }

        //Lambda-bar
        if( pythia.event[i].id() < 0 )
        {                  
          Lbar_vector.push_back(L_fourmom);
          Lbar_pT_bin_vector.push_back(pT_bin_corr);
          Lbar_eta_bin_vector.push_back(eta_bin);
          
          pBar_vector.push_back(p_fourmom);
          piBar_vector.push_back(pi_fourmom);
        
          pBar_star_vector.push_back(p_fourmom_star);      
        }
        
        
        //Lambda cuts
        //decay length cuts in mm (default PYTHIA units)
        if(L_decayL < 20 || L_decayL > 250) continue;
        //if(L_fourmom.Pt() < 0.5) continue; //added pT cut - pT integrated in data starts at 0.5 GeV/c
        //if(L_fourmom.Rapidity() >= 1) continue;
        
        //daughter cuts
        if( fabs(p_fourmom.Eta()) >= 1. || fabs(pi_fourmom.Eta()) >= 1. ) continue;       
        if( p_fourmom.Pt() < 0.15 || p_fourmom.Pt() > 20. ) continue;
        if( pi_fourmom.Pt() < 0.15 || pi_fourmom.Pt() > 20. ) continue;  
        

        //Lambda
        if( pythia.event[i].id() > 0 )
        {
          L_cuts_vector.push_back(L_fourmom);
          L_pT_bin_cuts_vector.push_back(pT_bin_corr);
          L_eta_bin_cuts_vector.push_back(eta_bin);
          
          p_cuts_vector.push_back(p_fourmom);
          pi_cuts_vector.push_back(pi_fourmom);
        
          p_star_cuts_vector.push_back(p_fourmom_star);        
        }

        //Lambda-bar
        if( pythia.event[i].id() < 0 )
        {
          Lbar_cuts_vector.push_back(L_fourmom);
          Lbar_pT_bin_cuts_vector.push_back(pT_bin_corr);
          Lbar_eta_bin_cuts_vector.push_back(eta_bin);
          
          pBar_cuts_vector.push_back(p_fourmom);
          piBar_cuts_vector.push_back(pi_fourmom);
        
          pBar_star_cuts_vector.push_back(p_fourmom_star);        
        }             
        
      } //end PID if for Lambdas
      
         

    } //end loop over particles in event
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
          double theta_star = p_star_vector.at(iLambda).Angle(pBar_star_vector.at(iLambdaBar).Vect());
         
          L0_L0bar_cosThetaProdPlane->Fill(TMath::Cos(theta_star));          
          L0_L0bar_cosThetaProdPlane_pT_hist[L_pT_bin_vector.at(iLambda)][Lbar_pT_bin_vector.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
          L0_L0bar_cosThetaProdPlane_eta_hist[L_eta_bin_vector.at(iLambda)][Lbar_eta_bin_vector.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));       
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
          double theta_star = p_star_vector.at(iLambda1).Angle(p_star_vector.at(iLambda2).Vect());
          
          L0_L0_cosThetaProdPlane->Fill(TMath::Cos(theta_star));
          L0_L0_cosThetaProdPlane_pT_hist[L_pT_bin_vector.at(iLambda1)][L_pT_bin_vector.at(iLambda2)]->Fill(TMath::Cos(theta_star));
          L0_L0_cosThetaProdPlane_eta_hist[L_eta_bin_vector.at(iLambda1)][L_eta_bin_vector.at(iLambda2)]->Fill(TMath::Cos(theta_star));            
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
          double theta_star = pBar_star_vector.at(iLambdaBar1).Angle(pBar_star_vector.at(iLambdaBar2).Vect());
          
          L0bar_L0bar_cosThetaProdPlane->Fill(TMath::Cos(theta_star));
          L0bar_L0bar_cosThetaProdPlane_pT_hist[Lbar_pT_bin_vector.at(iLambdaBar1)][Lbar_pT_bin_vector.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
          L0bar_L0bar_cosThetaProdPlane_eta_hist[Lbar_eta_bin_vector.at(iLambdaBar1)][Lbar_eta_bin_vector.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));                
        }     
      }     
    }
    //_____________________________________________________________________________________________________________________________________________________________________
    
    //mixed event before cuts
    if( p_star_vector.size() == 1 && pBar_star_vector.size() == 0 && p_star_vector_ME.size() < 1e4)
    {
      p_star_vector_ME.push_back(p_star_vector.at(0));

      L_pT_bin_vector_ME.push_back(L_pT_bin_vector.at(0));
      L_eta_bin_vector_ME.push_back(L_eta_bin_vector.at(0));
    }

    if( p_star_vector.size() == 0 && pBar_star_vector.size() == 1 && pBar_star_vector_ME.size() < 1e4)
    {    
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
          double theta_star = p_star_cuts_vector.at(iLambda).Angle(pBar_star_cuts_vector.at(iLambdaBar).Vect());
          
          L0_L0bar_cosThetaProdPlane_cuts->Fill(TMath::Cos(theta_star));
          L0_L0bar_cosThetaProdPlane_pT_cuts_hist[L_pT_bin_cuts_vector.at(iLambda)][Lbar_pT_bin_cuts_vector.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
          L0_L0bar_cosThetaProdPlane_eta_cuts_hist[L_eta_bin_cuts_vector.at(iLambda)][Lbar_eta_bin_cuts_vector.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));       
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
          double theta_star = p_star_cuts_vector.at(iLambda1).Angle(p_star_cuts_vector.at(iLambda2).Vect());
          
          L0_L0_cosThetaProdPlane_cuts->Fill(TMath::Cos(theta_star));
          L0_L0_cosThetaProdPlane_pT_cuts_hist[L_pT_bin_cuts_vector.at(iLambda1)][L_pT_bin_cuts_vector.at(iLambda2)]->Fill(TMath::Cos(theta_star));
          L0_L0_cosThetaProdPlane_eta_cuts_hist[L_eta_bin_cuts_vector.at(iLambda1)][L_eta_bin_cuts_vector.at(iLambda2)]->Fill(TMath::Cos(theta_star));           
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
          double theta_star = pBar_star_cuts_vector.at(iLambdaBar1).Angle(pBar_star_cuts_vector.at(iLambdaBar2).Vect());
          
          L0bar_L0bar_cosThetaProdPlane_cuts->Fill(TMath::Cos(theta_star));
          L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[Lbar_pT_bin_cuts_vector.at(iLambdaBar1)][Lbar_pT_bin_cuts_vector.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
          L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[Lbar_eta_bin_cuts_vector.at(iLambdaBar1)][Lbar_eta_bin_cuts_vector.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
        }   
      
      } 
    
    }
    
    //mixed event after cuts
    if( p_star_cuts_vector.size() == 1 && pBar_star_cuts_vector.size() == 0 && p_star_cuts_vector_ME.size() < 1e4)
    {
   
      p_star_cuts_vector_ME.push_back(p_star_cuts_vector.at(0));

      L_pT_bin_cuts_vector_ME.push_back(L_pT_bin_cuts_vector.at(0));
      L_eta_bin_cuts_vector_ME.push_back(L_eta_bin_cuts_vector.at(0));

    }

    if( p_star_cuts_vector.size() == 0 && pBar_star_cuts_vector.size() == 1 && pBar_star_cuts_vector_ME.size() < 1e4)
    {    
      pBar_star_cuts_vector_ME.push_back(pBar_star_cuts_vector.at(0));

      Lbar_pT_bin_cuts_vector_ME.push_back(Lbar_pT_bin_cuts_vector.at(0));
      Lbar_eta_bin_cuts_vector_ME.push_back(Lbar_eta_bin_cuts_vector.at(0));

    }
    //_________________________________________________________________________________________________________

    
    L_vector.clear();
    L_cuts_vector.clear();
    
    L_pT_bin_vector.clear();
    L_pT_bin_cuts_vector.clear();
    
    L_eta_bin_vector.clear();
    L_eta_bin_cuts_vector.clear();
    
    p_vector.clear();
    p_cuts_vector.clear();
    
    pi_vector.clear();
    pi_cuts_vector.clear();
    
    p_star_vector.clear();
    p_star_cuts_vector.clear();
    
    
    Lbar_vector.clear();
    Lbar_cuts_vector.clear();
    
    Lbar_pT_bin_vector.clear();
    Lbar_pT_bin_cuts_vector.clear();
    
    Lbar_eta_bin_vector.clear();
    Lbar_eta_bin_cuts_vector.clear();
    
    pBar_vector.clear();
    pBar_cuts_vector.clear();
    
    piBar_vector.clear();
    piBar_cuts_vector.clear();
    
    pBar_star_vector.clear();
    pBar_star_cuts_vector.clear();
    
    
    //---------------------------------------------------------------------------------------------------------------------
    
    // analyze pi p pairs
    
    //L and Lbar "candidates" from pi p pairs
    vector<TLorentzVector> L_vector_US;
    vector<int> L_pT_bin_vector_US;
    vector<int> L_eta_bin_vector_US;
    vector<int> L_cuts_flag_US; //flag if pair passed analysis cuts
    vector<int> L_Minv_flag_US; //invariant mass window flag
    
    vector<TLorentzVector> p_star_vector_US; //to store p fourmomentum in L (pi p pair) rest frame
    vector<int> p_index_vector_US; //to store p index for auto-correlation check for same-sign pairs
    vector<int> pi_index_vector_US; //to store pi index for auto-correlation check for same-sign pairs
    
    
    vector<TLorentzVector> Lbar_vector_US;
    vector<int> Lbar_pT_bin_vector_US;
    vector<int> Lbar_eta_bin_vector_US;
    vector<int> Lbar_cuts_flag_US; //flag if pair passed analysis cuts
    vector<int> Lbar_Minv_flag_US; //invariant mass window flag
    
    vector<TLorentzVector> pBar_star_vector_US; //to store p fourmomentum in Lbar (pi pBar pair) rest frame
    vector<int> pBar_index_vector_US; //to store p index for auto-correlation check for same-sign pairs
    vector<int> piBar_index_vector_US; //to store pi index for auto-correlation check for same-sign pairs
 
    
    
    vector<TLorentzVector> L_vector_LS;
    vector<int> L_pT_bin_vector_LS;
    vector<int> L_eta_bin_vector_LS;
    vector<int> L_cuts_flag_LS; //flag if pair passed analysis cuts
    vector<int> L_Minv_flag_LS; //invariant mass window flag    
    
    vector<TLorentzVector> p_star_vector_LS;
    vector<int> p_index_vector_LS; //to store p index for auto-correlation check for same-sign pairs
    vector<int> pi_index_vector_LS; //to store pi index for auto-correlation check for same-sign pairs 
    
    
    vector<TLorentzVector> Lbar_vector_LS;
    vector<int> Lbar_pT_bin_vector_LS;
    vector<int> Lbar_eta_bin_vector_LS;
    vector<int> Lbar_cuts_flag_LS; //flag if pair passed analysis cuts
    vector<int> Lbar_Minv_flag_LS; //invariant mass window flag
    
    vector<TLorentzVector> pBar_star_vector_LS;
    vector<int> pBar_index_vector_LS; //to store p index for auto-correlation check for same-sign pairs
    vector<int> piBar_index_vector_LS; //to store pi index for auto-correlation check for same-sign pairs 
    
 /*   
    cout<<p_vector_event.size()<<endl;
    cout<<pi_vector_event.size()<<endl;
    
    cout<<pBar_vector_event.size()<<endl;
    cout<<piBar_vector_event.size()<<endl;
 */   
    //pair pi and p
    //L
    for(unsigned int p_index = 0; p_index < p_vector_event.size(); p_index++)
    {
      //US
      for(unsigned int piBar_index = 0; piBar_index < piBar_vector_event.size(); piBar_index++)
      {
        TLorentzVector L_fourmom_US_event = p_vector_event.at(p_index) + piBar_vector_event.at(piBar_index);
        
        if( fabs( L_fourmom_US_event.Rapidity() ) > 1 ) continue;
        
        int pT_bin_L_US = findBinPt( L_fourmom_US_event, pT_bins_corr, nPtBins_corr );
        int eta_bin_L_US = findBinEta( L_fourmom_US_event, eta_bins, nEtaBins );
        
        if( pT_bin_L_US < 0 || eta_bin_L_US < 0 ) continue;        
        
        TLorentzVector L_fourmom_reverse_US_event;
        L_fourmom_reverse_US_event.SetPxPyPzE(-L_fourmom_US_event.Px(), -L_fourmom_US_event.Py(), -L_fourmom_US_event.Pz(), L_fourmom_US_event.E());
        
        TLorentzVector p_star_US = p_vector_event.at(p_index);
        p_star_US.Boost(L_fourmom_reverse_US_event.BoostVector()); //boost p 4-momntum to L rest frame
        
        
        //save info (before cuts)
        L_vector_US.push_back(L_fourmom_US_event);
        L_pT_bin_vector_US.push_back(pT_bin_L_US);
        L_eta_bin_vector_US.push_back(eta_bin_L_US);
        p_star_vector_US.push_back(p_star_US);
        
        p_index_vector_US.push_back(p_index);
        piBar_index_vector_US.push_back(piBar_index);
        
        
        //cuts        
        int cuts_flag = 1;
        int Minv_flag = 1;
        
        //decay length cuts in mm (default PYTHIA units)
        TVector3 L_pairVertex_US = ( p_origin_vector_event.at(p_index) + piBar_origin_vector_event.at(piBar_index) )*0.5; //secondary vertex of pi p pair
        float L_decayL_US = L_pairVertex_US.Mag(); //decay length of pi p pair assuming PV at (0,0,0)
        
        if(L_decayL_US < 20 || L_decayL_US > 250) cuts_flag = 0;
        
        //cos of pointing angle
        float L_point_angle_US = L_pairVertex_US.Angle(L_fourmom_US_event.Vect());
        if(cos(L_point_angle_US) < 0.996) cuts_flag = 0;
                
        //if(L_fourmom.Pt() < 0.5) continue; //added pT cut - pT integrated in data starts at 0.5 GeV/c
        //if(L_fourmom.Rapidity() >= 1) continue;
        
        //daughter cuts
        if( fabs(p_vector_event.at(p_index).Eta()) >= 1. || fabs(piBar_vector_event.at(piBar_index).Eta()) >= 1. )  cuts_flag = 0;       
        if( p_vector_event.at(p_index).Pt() < 0.15 || p_vector_event.at(p_index).Pt() > 20. )  cuts_flag = 0;
        if( piBar_vector_event.at(piBar_index).Pt() < 0.15 || piBar_vector_event.at(piBar_index).Pt() > 20. )  cuts_flag = 0;
                
        L_cuts_flag_US.push_back(cuts_flag);
        
        
        if( L_fourmom_US_event.M() < 1.115 || L_fourmom_US_event.M() > 1.12 ) Minv_flag = 0; //approximate Minv window from data
        
        L_Minv_flag_US.push_back(Minv_flag);
        
      }
      
      //LS
      for(unsigned int pi_index = 0; pi_index < pi_vector_event.size(); pi_index++)
      {
        TLorentzVector L_fourmom_LS_event = p_vector_event.at(p_index) + pi_vector_event.at(pi_index);
        
        if( fabs( L_fourmom_LS_event.Rapidity() ) > 1 ) continue;
        
        
        int pT_bin_L_LS = findBinPt( L_fourmom_LS_event, pT_bins_corr, nPtBins_corr );
        int eta_bin_L_LS = findBinEta( L_fourmom_LS_event, eta_bins, nEtaBins );
        
        if( pT_bin_L_LS < 0 || eta_bin_L_LS < 0 ) continue;
        
        
        TLorentzVector L_fourmom_reverse_LS_event;
        L_fourmom_reverse_LS_event.SetPxPyPzE(-L_fourmom_LS_event.Px(), -L_fourmom_LS_event.Py(), -L_fourmom_LS_event.Pz(), L_fourmom_LS_event.E());
        
        TLorentzVector p_star_LS = p_vector_event.at(p_index);
        p_star_LS.Boost(L_fourmom_reverse_LS_event.BoostVector()); //boost p 4-momntum to L rest frame
        
        //save info (before cuts)
        L_vector_LS.push_back(L_fourmom_LS_event);
        L_pT_bin_vector_LS.push_back(pT_bin_L_LS);
        L_eta_bin_vector_LS.push_back(eta_bin_L_LS);
        p_star_vector_LS.push_back(p_star_LS);
        
        p_index_vector_LS.push_back(p_index);
        pi_index_vector_LS.push_back(pi_index);
        
        //cuts        
        int cuts_flag = 1;
        int Minv_flag = 1;
        
        //decay length cuts in mm (default PYTHIA units)
        TVector3 L_pairVertex_LS = ( p_origin_vector_event.at(p_index) + pi_origin_vector_event.at(pi_index) )*0.5; //secondary vertex of pi p pair
        float L_decayL_LS = L_pairVertex_LS.Mag(); //decay length of pi p pair assuming PV at (0,0,0)
        
        if(L_decayL_LS < 20 || L_decayL_LS > 250) cuts_flag = 0;
        
        //cos of pointing angle
        float L_point_angle_LS = L_pairVertex_LS.Angle(L_fourmom_LS_event.Vect());
        if(cos(L_point_angle_LS) < 0.996) cuts_flag = 0;
        
        //if(L_fourmom.Pt() < 0.5) continue; //added pT cut - pT integrated in data starts at 0.5 GeV/c
        //if(L_fourmom.Rapidity() >= 1) continue;
        
        //daughter cuts
        if( fabs(p_vector_event.at(p_index).Eta()) >= 1. || fabs(pi_vector_event.at(pi_index).Eta()) >= 1. )  cuts_flag = 0;       
        if( p_vector_event.at(p_index).Pt() < 0.15 || p_vector_event.at(p_index).Pt() > 20. )  cuts_flag = 0;
        if( pi_vector_event.at(pi_index).Pt() < 0.15 || pi_vector_event.at(pi_index).Pt() > 20. )  cuts_flag = 0;       
        
        L_cuts_flag_LS.push_back(cuts_flag);
        
        
        if( L_fourmom_LS_event.M() < 1.115 || L_fourmom_LS_event.M() > 1.12 ) Minv_flag = 0; //approximate Minv window from data
        
        L_Minv_flag_LS.push_back(Minv_flag);
        
      }    
    }
    
    //Lbar
    for(unsigned int pBar_index = 0; pBar_index < pBar_vector_event.size(); pBar_index++)
    {
      //US
      for(unsigned int pi_index = 0; pi_index < pi_vector_event.size(); pi_index++)
      {
        TLorentzVector Lbar_fourmom_US_event = pBar_vector_event.at(pBar_index) + pi_vector_event.at(pi_index);
        
        if( fabs( Lbar_fourmom_US_event.Rapidity() ) > 1 ) continue;
        
        int pT_bin_Lbar_US = findBinPt( Lbar_fourmom_US_event, pT_bins_corr, nPtBins_corr );
        int eta_bin_Lbar_US = findBinEta( Lbar_fourmom_US_event, eta_bins, nEtaBins );
        
        if( pT_bin_Lbar_US < 0 || eta_bin_Lbar_US < 0 ) continue;         
        
        TLorentzVector Lbar_fourmom_reverse_US_event;
        Lbar_fourmom_reverse_US_event.SetPxPyPzE(-Lbar_fourmom_US_event.Px(), -Lbar_fourmom_US_event.Py(), -Lbar_fourmom_US_event.Pz(), Lbar_fourmom_US_event.E());
        
        TLorentzVector pBar_star_US = pBar_vector_event.at(pBar_index);
        pBar_star_US.Boost(Lbar_fourmom_reverse_US_event.BoostVector()); //boost p 4-momntum to Lbar rest frame
        
        //save info (before cuts)
        Lbar_vector_US.push_back(Lbar_fourmom_US_event);
        Lbar_pT_bin_vector_US.push_back(pT_bin_Lbar_US);
        Lbar_eta_bin_vector_US.push_back(eta_bin_Lbar_US);
        pBar_star_vector_US.push_back(pBar_star_US);
        
        pBar_index_vector_US.push_back(pBar_index);
        pi_index_vector_US.push_back(pi_index);        
        
        //cuts        
        int cuts_flag = 1;
        int Minv_flag = 1;
        
        //decay length cuts in mm (default PYTHIA units)
        TVector3 Lbar_pairVertex_US = ( pBar_origin_vector_event.at(pBar_index) + pi_origin_vector_event.at(pi_index) )*0.5; //secondary vertex of pi p pair
        float Lbar_decayL_US = Lbar_pairVertex_US.Mag(); //decay length of pi p pair assuming PV at (0,0,0)
        
        if(Lbar_decayL_US < 20 || Lbar_decayL_US > 250) cuts_flag = 0;
        
        //cos of pointing angle
        float Lbar_point_angle_US = Lbar_pairVertex_US.Angle(Lbar_fourmom_US_event.Vect());
        if(cos(Lbar_point_angle_US) < 0.996) cuts_flag = 0;
        
        //if(L_fourmom.Pt() < 0.5) continue; //added pT cut - pT integrated in data starts at 0.5 GeV/c
        //if(L_fourmom.Rapidity() >= 1) continue;
        
        //daughter cuts
        if( fabs(pBar_vector_event.at(pBar_index).Eta()) >= 1. || fabs(pi_vector_event.at(pi_index).Eta()) >= 1. )  cuts_flag = 0;       
        if( pBar_vector_event.at(pBar_index).Pt() < 0.15 || pBar_vector_event.at(pBar_index).Pt() > 20. )  cuts_flag = 0;
        if( pi_vector_event.at(pi_index).Pt() < 0.15 || pi_vector_event.at(pi_index).Pt() > 20. )  cuts_flag = 0;        
        
        Lbar_cuts_flag_US.push_back(cuts_flag);        
        
        
        if( Lbar_fourmom_US_event.M() < 1.115 || Lbar_fourmom_US_event.M() > 1.12 ) Minv_flag = 0; //approximate Minv window from data
        
        Lbar_Minv_flag_US.push_back(Minv_flag);
        
      }
      
      //LS
      for(unsigned int piBar_index = 0; piBar_index < piBar_vector_event.size(); piBar_index++)
      {
        TLorentzVector Lbar_fourmom_LS_event = pBar_vector_event.at(pBar_index) + piBar_vector_event.at(piBar_index);
                
        if( fabs( Lbar_fourmom_LS_event.Rapidity() ) > 1 ) continue;
        
        int pT_bin_Lbar_LS = findBinPt( Lbar_fourmom_LS_event, pT_bins_corr, nPtBins_corr );
        int eta_bin_Lbar_LS = findBinEta( Lbar_fourmom_LS_event, eta_bins, nEtaBins );
        
        if( pT_bin_Lbar_LS < 0 || eta_bin_Lbar_LS < 0 ) continue;   
        
        TLorentzVector Lbar_fourmom_reverse_LS_event;
        Lbar_fourmom_reverse_LS_event.SetPxPyPzE(-Lbar_fourmom_LS_event.Px(), -Lbar_fourmom_LS_event.Py(), -Lbar_fourmom_LS_event.Pz(), Lbar_fourmom_LS_event.E());
        
        TLorentzVector pBar_star_LS = pBar_vector_event.at(pBar_index);
        pBar_star_LS.Boost(Lbar_fourmom_reverse_LS_event.BoostVector()); //boost p 4-momntum to Lbar rest frame
        
        //save info (before cuts)
        Lbar_vector_LS.push_back(Lbar_fourmom_LS_event);
        Lbar_pT_bin_vector_LS.push_back(pT_bin_Lbar_LS);
        Lbar_eta_bin_vector_LS.push_back(eta_bin_Lbar_LS);
        pBar_star_vector_LS.push_back(pBar_star_LS);
        
        pBar_index_vector_LS.push_back(pBar_index);
        piBar_index_vector_LS.push_back(piBar_index);
        
        //cuts        
        int cuts_flag = 1;
        int Minv_flag = 1;
        
        //decay length cuts in mm (default PYTHIA units)
        TVector3 Lbar_pairVertex_LS = ( pBar_origin_vector_event.at(pBar_index) + piBar_origin_vector_event.at(piBar_index) )*0.5; //secondary vertex of pi p pair
        float Lbar_decayL_LS = Lbar_pairVertex_LS.Mag(); //decay length of pi p pair assuming PV at (0,0,0)
        
        if(Lbar_decayL_LS < 20 || Lbar_decayL_LS > 250) cuts_flag = 0;
        
        //cos of pointing angle
        float Lbar_point_angle_LS = Lbar_pairVertex_LS.Angle(Lbar_fourmom_LS_event.Vect());
        if(cos(Lbar_point_angle_LS) < 0.996) cuts_flag = 0;
        
        //if(L_fourmom.Pt() < 0.5) continue; //added pT cut - pT integrated in data starts at 0.5 GeV/c
        //if(L_fourmom.Rapidity() >= 1) continue;
        
        //daughter cuts
        if( fabs(pBar_vector_event.at(pBar_index).Eta()) >= 1. || fabs(piBar_vector_event.at(piBar_index).Eta()) >= 1. )  cuts_flag = 0;       
        if( pBar_vector_event.at(pBar_index).Pt() < 0.15 || pBar_vector_event.at(pBar_index).Pt() > 20. )  cuts_flag = 0;
        if( piBar_vector_event.at(piBar_index).Pt() < 0.15 || piBar_vector_event.at(piBar_index).Pt() > 20. )  cuts_flag = 0;
                
        Lbar_cuts_flag_LS.push_back(cuts_flag);
        
                
        if( Lbar_fourmom_LS_event.M() < 1.115 || Lbar_fourmom_LS_event.M() > 1.12 ) Minv_flag = 0; //approximate Minv window from data
        
        Lbar_Minv_flag_LS.push_back(Minv_flag);
          
      }     
    }
    
    //clear p and pi vectors for next events
       
    
    p_vector_event.clear();
    p_origin_vector_event.clear();
    
    pi_vector_event.clear();
    pi_origin_vector_event.clear();
    
    
    pBar_vector_event.clear();
    pBar_origin_vector_event.clear();
    
    piBar_vector_event.clear();
    piBar_origin_vector_event.clear();
    
   //---------------------------------------------------------------------------------------------    
    
    //analyze L and Lbar candidates created from pi p pairs
    //cout<<"test 1"<<endl;
   /* 
    cout<<L_vector_US.size()<<endl;
    cout<<Lbar_vector_US.size()<<endl;
    
    cout<<L_vector_LS.size()<<endl;
    cout<<Lbar_vector_LS.size()<<endl;
    
   */ 
    //US paired with US
    //L-Lbar
    if( L_vector_US.size() > 0 && Lbar_vector_US.size() > 0 )
    {
      for( unsigned int iLambda = 0; iLambda < L_vector_US.size(); iLambda++)
      {
        for( unsigned int iLambdaBar = 0; iLambdaBar < Lbar_vector_US.size(); iLambdaBar++)
        {
          L0_inv_mass_vs_L0bar_inv_mass_US[L_pT_bin_vector_US.at(iLambda)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill( L_vector_US.at(iLambda).M(), Lbar_vector_US.at(iLambdaBar).M() );
          
          float theta_star = p_star_vector_US.at(iLambda).Angle(pBar_star_vector_US.at(iLambdaBar).Vect());
          
          if( L_Minv_flag_US.at(iLambda) == 1 && Lbar_Minv_flag_US.at(iLambdaBar) == 1) //fill cos theta* for pairs in Minv window
          {
            L0_L0bar_cosThetaProdPlane_US->Fill(TMath::Cos(theta_star));
            L0_L0bar_cosThetaProdPlane_US_pT_hist[L_pT_bin_vector_US.at(iLambda)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
            L0_L0bar_cosThetaProdPlane_US_eta_hist[L_eta_bin_vector_US.at(iLambda)][Lbar_eta_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));          
          }        
          
          
          if( L_cuts_flag_US.at(iLambda) == 1 && Lbar_cuts_flag_US.at(iLambdaBar) == 1) //both L in pair passed cuts
          {
            L0_inv_mass_vs_L0bar_inv_mass_US_cuts[L_pT_bin_vector_US.at(iLambda)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill( L_vector_US.at(iLambda).M(), Lbar_vector_US.at(iLambdaBar).M() );
            
            if( L_Minv_flag_US.at(iLambda) == 1 && Lbar_Minv_flag_US.at(iLambdaBar) == 1) //fill cos theta* for pairs in Minv window
            {            
              L0_L0bar_cosThetaProdPlane_US_cuts->Fill(TMath::Cos(theta_star));
              L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[L_pT_bin_vector_US.at(iLambda)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
              L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[L_eta_bin_vector_US.at(iLambda)][Lbar_eta_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
            }
            
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
          if( p_index_vector_US.at(iLambda1) == p_index_vector_US.at(iLambda2) ) continue;
          if( piBar_index_vector_US.at(iLambda1) == piBar_index_vector_US.at(iLambda2) ) continue;
          
          L0_inv_mass_vs_L0_inv_mass_US[L_pT_bin_vector_US.at(iLambda1)][L_pT_bin_vector_US.at(iLambda2)]->Fill( L_vector_US.at(iLambda1).M(), L_vector_US.at(iLambda2).M() );         
          
          float theta_star = p_star_vector_US.at(iLambda1).Angle(p_star_vector_US.at(iLambda2).Vect());
          
          if( L_Minv_flag_US.at(iLambda1) == 1 && L_Minv_flag_US.at(iLambda2) == 1) //fill cos theta* for pairs in Minv window
          {
            L0_L0_cosThetaProdPlane_US->Fill(TMath::Cos(theta_star));
            L0_L0_cosThetaProdPlane_US_pT_hist[L_pT_bin_vector_US.at(iLambda1)][L_pT_bin_vector_US.at(iLambda2)]->Fill(TMath::Cos(theta_star));
            L0_L0_cosThetaProdPlane_US_eta_hist[L_eta_bin_vector_US.at(iLambda1)][L_eta_bin_vector_US.at(iLambda2)]->Fill(TMath::Cos(theta_star));
          }
          
          if( L_cuts_flag_US.at(iLambda1) == 1 && L_cuts_flag_US.at(iLambda2) == 1) //both L in pair passed cuts
          {
            L0_inv_mass_vs_L0_inv_mass_US_cuts[L_pT_bin_vector_US.at(iLambda1)][L_pT_bin_vector_US.at(iLambda2)]->Fill( L_vector_US.at(iLambda1).M(), L_vector_US.at(iLambda2).M() );
            
            if( L_Minv_flag_US.at(iLambda1) == 1 && L_Minv_flag_US.at(iLambda2) == 1) //fill cos theta* for pairs in Minv window
            {  
              L0_L0_cosThetaProdPlane_US_cuts->Fill(TMath::Cos(theta_star));
              L0_L0_cosThetaProdPlane_US_pT_cuts_hist[L_pT_bin_vector_US.at(iLambda1)][L_pT_bin_vector_US.at(iLambda2)]->Fill(TMath::Cos(theta_star));
              L0_L0_cosThetaProdPlane_US_eta_cuts_hist[L_eta_bin_vector_US.at(iLambda1)][L_eta_bin_vector_US.at(iLambda2)]->Fill(TMath::Cos(theta_star));
            }
            
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
          if( pBar_index_vector_US.at(iLambdaBar1) == pBar_index_vector_US.at(iLambdaBar2) ) continue;
          if( pi_index_vector_US.at(iLambdaBar1) == pi_index_vector_US.at(iLambdaBar2) ) continue;
          
          L0bar_inv_mass_vs_L0bar_inv_mass_US[Lbar_pT_bin_vector_US.at(iLambdaBar1)][Lbar_pT_bin_vector_US.at(iLambdaBar2)]->Fill( Lbar_vector_US.at(iLambdaBar1).M(), Lbar_vector_US.at(iLambdaBar2).M() );
          
          float theta_star = pBar_star_vector_US.at(iLambdaBar1).Angle(pBar_star_vector_US.at(iLambdaBar2).Vect());
          
          if( Lbar_Minv_flag_US.at(iLambdaBar1) == 1 && Lbar_Minv_flag_US.at(iLambdaBar2) == 1) //fill cos theta* for pairs in Minv window
          {
            L0bar_L0bar_cosThetaProdPlane_US->Fill(TMath::Cos(theta_star));
            L0bar_L0bar_cosThetaProdPlane_US_pT_hist[Lbar_pT_bin_vector_US.at(iLambdaBar1)][Lbar_pT_bin_vector_US.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
            L0bar_L0bar_cosThetaProdPlane_US_eta_hist[Lbar_eta_bin_vector_US.at(iLambdaBar1)][Lbar_eta_bin_vector_US.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
          }
          
          if( Lbar_cuts_flag_US.at(iLambdaBar1) == 1 && Lbar_cuts_flag_US.at(iLambdaBar2) == 1) //both L in pair passed cuts
          {
            L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts[Lbar_pT_bin_vector_US.at(iLambdaBar1)][Lbar_pT_bin_vector_US.at(iLambdaBar2)]->Fill( Lbar_vector_US.at(iLambdaBar1).M(), Lbar_vector_US.at(iLambdaBar2).M() );
            
            if( Lbar_Minv_flag_US.at(iLambdaBar1) == 1 && Lbar_Minv_flag_US.at(iLambdaBar2) == 1) //fill cos theta* for pairs in Minv window
            {
              L0bar_L0bar_cosThetaProdPlane_US_cuts->Fill(TMath::Cos(theta_star));
              L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[Lbar_pT_bin_vector_US.at(iLambdaBar1)][Lbar_pT_bin_vector_US.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
              L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[Lbar_eta_bin_vector_US.at(iLambdaBar1)][Lbar_eta_bin_vector_US.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
            }
            
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
          if( piBar_index_vector_US.at(iLambda) == piBar_index_vector_LS.at(iLambdaBar) ) continue;
          
          L0_inv_mass_vs_L0bar_inv_mass_US_LS[L_pT_bin_vector_US.at(iLambda)][Lbar_pT_bin_vector_LS.at(iLambdaBar)]->Fill( L_vector_US.at(iLambda).M(), Lbar_vector_LS.at(iLambdaBar).M() );          
        
          float theta_star = p_star_vector_US.at(iLambda).Angle(pBar_star_vector_LS.at(iLambdaBar).Vect());
          
          if( L_Minv_flag_US.at(iLambda) == 1 && Lbar_Minv_flag_LS.at(iLambdaBar) == 1) //fill cos theta* for pairs in Minv window
          {
            L0_L0bar_cosThetaProdPlane_US_LS->Fill(TMath::Cos(theta_star));
            L0_L0bar_cosThetaProdPlane_US_LS_pT_hist[L_pT_bin_vector_US.at(iLambda)][Lbar_pT_bin_vector_LS.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
            L0_L0bar_cosThetaProdPlane_US_LS_eta_hist[L_eta_bin_vector_US.at(iLambda)][Lbar_eta_bin_vector_LS.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
          }
          
          if( L_cuts_flag_US.at(iLambda) == 1 && Lbar_cuts_flag_LS.at(iLambdaBar) == 1) //both L in pair passed cuts
          {
            L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[L_pT_bin_vector_US.at(iLambda)][Lbar_pT_bin_vector_LS.at(iLambdaBar)]->Fill( L_vector_US.at(iLambda).M(), Lbar_vector_LS.at(iLambdaBar).M() );
          
            if( L_Minv_flag_US.at(iLambda) == 1 && Lbar_Minv_flag_LS.at(iLambdaBar) == 1) //fill cos theta* for pairs in Minv window
            {
              L0_L0bar_cosThetaProdPlane_US_LS_cuts->Fill(TMath::Cos(theta_star));
              L0_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[L_pT_bin_vector_US.at(iLambda)][Lbar_pT_bin_vector_LS.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
              L0_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[L_eta_bin_vector_US.at(iLambda)][Lbar_eta_bin_vector_LS.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));          
            }
            
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
          if( pi_index_vector_LS.at(iLambda) == pi_index_vector_US.at(iLambdaBar) ) continue;
          
          L0_inv_mass_vs_L0bar_inv_mass_US_LS[L_pT_bin_vector_LS.at(iLambda)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill( L_vector_LS.at(iLambda).M(), Lbar_vector_US.at(iLambdaBar).M() );
        
          float theta_star = p_star_vector_LS.at(iLambda).Angle(pBar_star_vector_US.at(iLambdaBar).Vect());
          
          if( L_Minv_flag_LS.at(iLambda) == 1 && Lbar_Minv_flag_US.at(iLambdaBar) == 1) //fill cos theta* for pairs in Minv window
          {
            L0_L0bar_cosThetaProdPlane_US_LS->Fill(TMath::Cos(theta_star));
            L0_L0bar_cosThetaProdPlane_US_LS_pT_hist[L_pT_bin_vector_LS.at(iLambda)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
            L0_L0bar_cosThetaProdPlane_US_LS_eta_hist[L_eta_bin_vector_LS.at(iLambda)][Lbar_eta_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
          }
          
          if( L_cuts_flag_LS.at(iLambda) == 1 && Lbar_cuts_flag_US.at(iLambdaBar) == 1) //both L in pair passed cuts
          {
            L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[L_pT_bin_vector_LS.at(iLambda)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill( L_vector_LS.at(iLambda).M(), Lbar_vector_US.at(iLambdaBar).M() );
            
            if( L_Minv_flag_LS.at(iLambda) == 1 && Lbar_Minv_flag_US.at(iLambdaBar) == 1) //fill cos theta* for pairs in Minv window
            {
              L0_L0bar_cosThetaProdPlane_US_LS_cuts->Fill(TMath::Cos(theta_star));
              L0_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[L_pT_bin_vector_LS.at(iLambda)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
              L0_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[L_eta_bin_vector_LS.at(iLambda)][Lbar_eta_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));          
            }
            
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
          if( p_index_vector_US.at(iLambda) == p_index_vector_LS.at(iLambdaBck) ) continue;
          
          //interate order of US and LS for L
          if( nFills_L_L_US_LS % 2 == 0 )
          {
            L0_inv_mass_vs_L0_inv_mass_US_LS[L_pT_bin_vector_US.at(iLambda)][L_pT_bin_vector_LS.at(iLambdaBck)]->Fill( L_vector_US.at(iLambda).M(), L_vector_LS.at(iLambdaBck).M() );
          
            float theta_star = p_star_vector_US.at(iLambda).Angle(p_star_vector_LS.at(iLambdaBck).Vect());
            
            if( L_Minv_flag_US.at(iLambda) == 1 && L_Minv_flag_LS.at(iLambdaBck) == 1) //fill cos theta* for pairs in Minv window
            {
              L0_L0_cosThetaProdPlane_US_LS->Fill(TMath::Cos(theta_star));
              L0_L0_cosThetaProdPlane_US_LS_pT_hist[L_pT_bin_vector_US.at(iLambda)][L_pT_bin_vector_LS.at(iLambdaBck)]->Fill(TMath::Cos(theta_star));
              L0_L0_cosThetaProdPlane_US_LS_eta_hist[L_eta_bin_vector_US.at(iLambda)][L_eta_bin_vector_LS.at(iLambdaBck)]->Fill(TMath::Cos(theta_star));
            }
            
            if( L_cuts_flag_US.at(iLambda) == 1 && L_cuts_flag_LS.at(iLambdaBck) == 1) //both L in pair passed cuts
            {
              L0_inv_mass_vs_L0_inv_mass_US_LS_cuts[L_pT_bin_vector_US.at(iLambda)][L_pT_bin_vector_LS.at(iLambdaBck)]->Fill( L_vector_US.at(iLambda).M(), L_vector_LS.at(iLambdaBck).M() );
              
              if( L_Minv_flag_US.at(iLambda) == 1 && L_Minv_flag_LS.at(iLambdaBck) == 1) //fill cos theta* for pairs in Minv window
              {
                L0_L0_cosThetaProdPlane_US_LS_cuts->Fill(TMath::Cos(theta_star));
                L0_L0_cosThetaProdPlane_US_LS_pT_cuts_hist[L_pT_bin_vector_US.at(iLambda)][L_pT_bin_vector_LS.at(iLambdaBck)]->Fill(TMath::Cos(theta_star));
                L0_L0_cosThetaProdPlane_US_LS_eta_cuts_hist[L_eta_bin_vector_US.at(iLambda)][L_eta_bin_vector_LS.at(iLambdaBck)]->Fill(TMath::Cos(theta_star));          
              }
              
            }          
          
          }
          else
          {
            L0_inv_mass_vs_L0_inv_mass_US_LS[L_pT_bin_vector_LS.at(iLambdaBck)][L_pT_bin_vector_US.at(iLambda)]->Fill( L_vector_LS.at(iLambdaBck).M(), L_vector_US.at(iLambda).M() );
          
            float theta_star = p_star_vector_LS.at(iLambdaBck).Angle(p_star_vector_US.at(iLambda).Vect());
            
            if( L_Minv_flag_US.at(iLambda) == 1 && L_Minv_flag_LS.at(iLambdaBck) == 1) //fill cos theta* for pairs in Minv window
            {    
              L0_L0_cosThetaProdPlane_US_LS->Fill(TMath::Cos(theta_star));
              L0_L0_cosThetaProdPlane_US_LS_pT_hist[L_pT_bin_vector_LS.at(iLambdaBck)][L_pT_bin_vector_US.at(iLambda)]->Fill(TMath::Cos(theta_star));
              L0_L0_cosThetaProdPlane_US_LS_eta_hist[L_eta_bin_vector_LS.at(iLambdaBck)][L_eta_bin_vector_US.at(iLambda)]->Fill(TMath::Cos(theta_star));
            }
            
            if( L_cuts_flag_LS.at(iLambdaBck) == 1 && L_cuts_flag_US.at(iLambda) == 1) //both L in pair passed cuts
            {
              L0_inv_mass_vs_L0_inv_mass_US_LS_cuts[L_pT_bin_vector_LS.at(iLambdaBck)][L_pT_bin_vector_US.at(iLambda)]->Fill( L_vector_LS.at(iLambdaBck).M(), L_vector_US.at(iLambda).M() );
              
              if( L_Minv_flag_US.at(iLambda) == 1 && L_Minv_flag_LS.at(iLambdaBck) == 1) //fill cos theta* for pairs in Minv window
              {
                L0_L0_cosThetaProdPlane_US_LS_cuts->Fill(TMath::Cos(theta_star));
                L0_L0_cosThetaProdPlane_US_LS_pT_cuts_hist[L_pT_bin_vector_LS.at(iLambdaBck)][L_pT_bin_vector_US.at(iLambda)]->Fill(TMath::Cos(theta_star));
                L0_L0_cosThetaProdPlane_US_LS_eta_cuts_hist[L_eta_bin_vector_LS.at(iLambdaBck)][L_eta_bin_vector_US.at(iLambda)]->Fill(TMath::Cos(theta_star));          
              }
              
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
          if( pBar_index_vector_US.at(iLambdaBar) == pBar_index_vector_LS.at(iLambdaBarBck) ) continue;
          
          //interate order of US and LS for L
          if( nFills_Lbar_Lbar_US_LS % 2 == 0 )
          {
            L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[Lbar_pT_bin_vector_US.at(iLambdaBar)][Lbar_pT_bin_vector_LS.at(iLambdaBarBck)]->Fill( Lbar_vector_US.at(iLambdaBar).M(), Lbar_vector_LS.at(iLambdaBarBck).M() );
          
            float theta_star = pBar_star_vector_US.at(iLambdaBar).Angle(pBar_star_vector_LS.at(iLambdaBarBck).Vect());
            
            if( Lbar_Minv_flag_US.at(iLambdaBar) == 1 && Lbar_Minv_flag_LS.at(iLambdaBarBck) == 1) //fill cos theta* for pairs in Minv window
            {       
              L0bar_L0bar_cosThetaProdPlane_US_LS->Fill(TMath::Cos(theta_star));
              L0bar_L0bar_cosThetaProdPlane_US_LS_pT_hist[Lbar_pT_bin_vector_US.at(iLambdaBar)][Lbar_pT_bin_vector_LS.at(iLambdaBarBck)]->Fill(TMath::Cos(theta_star));
              L0bar_L0bar_cosThetaProdPlane_US_LS_eta_hist[Lbar_eta_bin_vector_US.at(iLambdaBar)][Lbar_eta_bin_vector_LS.at(iLambdaBarBck)]->Fill(TMath::Cos(theta_star));
            }
            
            if( Lbar_cuts_flag_US.at(iLambdaBar) == 1 && Lbar_cuts_flag_LS.at(iLambdaBarBck) == 1) //both L in pair passed cuts
            {
              L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[Lbar_pT_bin_vector_US.at(iLambdaBar)][Lbar_pT_bin_vector_LS.at(iLambdaBarBck)]->Fill( Lbar_vector_US.at(iLambdaBar).M(), Lbar_vector_LS.at(iLambdaBarBck).M() );
              
              if( Lbar_Minv_flag_US.at(iLambdaBar) == 1 && Lbar_Minv_flag_LS.at(iLambdaBarBck) == 1) //fill cos theta* for pairs in Minv window
              { 
                L0bar_L0bar_cosThetaProdPlane_US_LS_cuts->Fill(TMath::Cos(theta_star));
                L0bar_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[Lbar_pT_bin_vector_US.at(iLambdaBar)][Lbar_pT_bin_vector_LS.at(iLambdaBarBck)]->Fill(TMath::Cos(theta_star));
                L0bar_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[Lbar_eta_bin_vector_US.at(iLambdaBar)][Lbar_eta_bin_vector_LS.at(iLambdaBarBck)]->Fill(TMath::Cos(theta_star));          
              }
              
            }          
          
          }
          else
          {
            L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[Lbar_pT_bin_vector_LS.at(iLambdaBarBck)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill( Lbar_vector_LS.at(iLambdaBarBck).M(), Lbar_vector_US.at(iLambdaBar).M() );
          
            float theta_star = pBar_star_vector_LS.at(iLambdaBarBck).Angle(pBar_star_vector_US.at(iLambdaBar).Vect());
                   
            if( Lbar_Minv_flag_US.at(iLambdaBar) == 1 && Lbar_Minv_flag_LS.at(iLambdaBarBck) == 1) //fill cos theta* for pairs in Minv window
            {    
              L0bar_L0bar_cosThetaProdPlane_US_LS->Fill(TMath::Cos(theta_star));
              L0bar_L0bar_cosThetaProdPlane_US_LS_pT_hist[Lbar_pT_bin_vector_LS.at(iLambdaBarBck)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
              L0bar_L0bar_cosThetaProdPlane_US_LS_eta_hist[Lbar_eta_bin_vector_LS.at(iLambdaBarBck)][Lbar_eta_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
            }
            
            if( Lbar_cuts_flag_LS.at(iLambdaBarBck) == 1 && Lbar_cuts_flag_US.at(iLambdaBar) == 1) //both L in pair passed cuts
            {
              L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[Lbar_pT_bin_vector_LS.at(iLambdaBarBck)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill( Lbar_vector_LS.at(iLambdaBarBck).M(), Lbar_vector_US.at(iLambdaBar).M() );
              
              if( Lbar_Minv_flag_US.at(iLambdaBar) == 1 && Lbar_Minv_flag_LS.at(iLambdaBarBck) == 1) //fill cos theta* for pairs in Minv window
              { 
                L0bar_L0bar_cosThetaProdPlane_US_LS_cuts->Fill(TMath::Cos(theta_star));
                L0bar_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[Lbar_pT_bin_vector_LS.at(iLambdaBarBck)][Lbar_pT_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
                L0bar_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[Lbar_eta_bin_vector_LS.at(iLambdaBarBck)][Lbar_eta_bin_vector_US.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));          
              }
              
            }
            
          }//end else
          
          nFills_Lbar_Lbar_US_LS++;
        
        }
      
      }
    
    }
    
    //----------------------------------------------------------------
    
    //select L and Lbar for mixed event from pi p pairs
    //before cuts
    if( L_vector_US.size() == 1 && Lbar_vector_US.size() == 0 && p_star_vector_US_ME.size() < 1e4)
    {
      if( L_Minv_flag_US.at(0) == 0 ) continue;
    
      p_star_vector_US_ME.push_back(p_star_vector_US.at(0));

      L_pT_bin_vector_US_ME.push_back(L_pT_bin_vector_US.at(0));
      L_eta_bin_vector_US_ME.push_back(L_eta_bin_vector_US.at(0));     
    }

    if( L_vector_US.size() == 0 && Lbar_vector_US.size() == 1 && pBar_star_vector_US_ME.size() < 1e4)
    {
      if( Lbar_Minv_flag_US.at(0) == 0 ) continue;
      
      pBar_star_vector_US_ME.push_back(pBar_star_vector_US.at(0));

      Lbar_pT_bin_vector_US_ME.push_back(Lbar_pT_bin_vector_US.at(0));
      Lbar_eta_bin_vector_US_ME.push_back(Lbar_eta_bin_vector_US.at(0));

    }
    
    //after cuts
    if( L_vector_US.size() == 1 && Lbar_vector_US.size() == 0 && p_star_vector_US_ME_cuts.size() < 1e4)
    {
      if( L_cuts_flag_US.at(0) == 0 ) continue;
      if( L_Minv_flag_US.at(0) == 0 ) continue;
      
      p_star_vector_US_ME_cuts.push_back(p_star_vector_US.at(0));

      L_pT_bin_vector_US_ME_cuts.push_back(L_pT_bin_vector_US.at(0));
      L_eta_bin_vector_US_ME_cuts.push_back(L_eta_bin_vector_US.at(0));     
    }

    if( L_vector_US.size() == 0 && Lbar_vector_US.size() == 1 && pBar_star_vector_US_ME_cuts.size() < 1e4)
    {
      //cout<<"fill Lbar ME"<<endl;
      if( Lbar_cuts_flag_US.at(0) == 0 ) continue;
      if( Lbar_Minv_flag_US.at(0) == 0 ) continue;
      
      pBar_star_vector_US_ME_cuts.push_back(pBar_star_vector_US.at(0));

      Lbar_pT_bin_vector_US_ME_cuts.push_back(Lbar_pT_bin_vector_US.at(0));
      Lbar_eta_bin_vector_US_ME_cuts.push_back(Lbar_eta_bin_vector_US.at(0));

    }

    //_________________________________________________________________________________________________________
    
    //before cuts
    if( L_vector_LS.size() == 1 && Lbar_vector_LS.size() == 0 && p_star_vector_LS_ME.size() < 1e4)
    {
      if( L_Minv_flag_LS.at(0) == 0 ) continue;
      
      p_star_vector_LS_ME.push_back(p_star_vector_LS.at(0));

      L_pT_bin_vector_LS_ME.push_back(L_pT_bin_vector_LS.at(0));
      L_eta_bin_vector_LS_ME.push_back(L_eta_bin_vector_LS.at(0));
    }

    if( L_vector_LS.size() == 0 && Lbar_vector_LS.size() == 1 && pBar_star_vector_LS_ME.size() < 1e4)
    {
      if( Lbar_Minv_flag_LS.at(0) == 0 ) continue;
    
      pBar_star_vector_LS_ME.push_back(pBar_star_vector_LS.at(0));

      Lbar_pT_bin_vector_LS_ME.push_back(Lbar_pT_bin_vector_LS.at(0));
      Lbar_eta_bin_vector_LS_ME.push_back(Lbar_eta_bin_vector_LS.at(0));
    }
    
    //after cuts
    if( L_vector_LS.size() == 1 && Lbar_vector_LS.size() == 0 && p_star_vector_LS_ME_cuts.size() < 1e4)
    {
      if( L_cuts_flag_LS.at(0) == 0 ) continue;
      if( L_Minv_flag_LS.at(0) == 0 ) continue;
      
      p_star_vector_LS_ME_cuts.push_back(p_star_vector_LS.at(0));

      L_pT_bin_vector_LS_ME_cuts.push_back(L_pT_bin_vector_LS.at(0));
      L_eta_bin_vector_LS_ME_cuts.push_back(L_eta_bin_vector_LS.at(0));     
    }

    if( L_vector_LS.size() == 0 && Lbar_vector_LS.size() == 1 && pBar_star_vector_LS_ME_cuts.size() < 1e4)
    {
      //cout<<"fill Lbar ME"<<endl;
      if( Lbar_cuts_flag_LS.at(0) == 0 ) continue;
      if( Lbar_Minv_flag_LS.at(0) == 0 ) continue;
      
      pBar_star_vector_LS_ME_cuts.push_back(pBar_star_vector_LS.at(0));

      Lbar_pT_bin_vector_LS_ME_cuts.push_back(Lbar_pT_bin_vector_LS.at(0));
      Lbar_eta_bin_vector_LS_ME_cuts.push_back(Lbar_eta_bin_vector_LS.at(0));
    }

  }//end event loop
  //_____________________________________________________________________________________________________________________________________________________________________
  
  //cout<<p_star_vector_ME.size()<<endl;
  //cout<<pBar_star_vector_ME.size()<<endl;
  
  //mixed event before cuts
  for(unsigned int iLambda = 0; iLambda < p_star_vector_ME.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_ME.size(); iLambdaBar++)
    {
      double theta_star = p_star_vector_ME.at(iLambda).Angle(pBar_star_vector_ME.at(iLambdaBar).Vect());
           
      L0_L0bar_cosThetaProdPlane_ME->Fill(TMath::Cos(theta_star));          
      L0_L0bar_cosThetaProdPlane_ME_pT_hist[L_pT_bin_vector_ME.at(iLambda)][Lbar_pT_bin_vector_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_ME_eta_hist[L_eta_bin_vector_ME.at(iLambda)][Lbar_eta_bin_vector_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));  
    }   
  
  }
  
  //L0-L0
  for(unsigned int iLambda1 = 0; iLambda1 < p_star_vector_ME.size(); iLambda1++)
  {
    for(unsigned int iLambda2 = iLambda1+1; iLambda2 < p_star_vector_ME.size(); iLambda2++)
    {
      double theta_star = p_star_vector_ME.at(iLambda1).Angle(p_star_vector_ME.at(iLambda2).Vect());
      
      L0_L0_cosThetaProdPlane_ME->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_ME_pT_hist[L_pT_bin_vector_ME.at(iLambda1)][L_pT_bin_vector_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_ME_eta_hist[L_eta_bin_vector_ME.at(iLambda1)][L_eta_bin_vector_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star));
    }   
  
  }
  
  
  //L0bar-L0bar
  for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pBar_star_vector_ME.size(); iLambdaBar1++)
  {
    for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < pBar_star_vector_ME.size(); iLambdaBar2++)
    {
      double theta_star = pBar_star_vector_ME.at(iLambdaBar1).Angle(pBar_star_vector_ME.at(iLambdaBar2).Vect());
      
      L0bar_L0bar_cosThetaProdPlane_ME->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_ME_pT_hist[Lbar_pT_bin_vector_ME.at(iLambdaBar1)][Lbar_pT_bin_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_ME_eta_hist[Lbar_eta_bin_vector_ME.at(iLambdaBar1)][Lbar_eta_bin_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));            
    }   
  
  }
  //_____________________________________________________________________________________________________________________________________________________________________
  
  
  //mixed event after cuts
  for(unsigned int iLambda = 0; iLambda < p_star_cuts_vector_ME.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_cuts_vector_ME.size(); iLambdaBar++)
    {
      double theta_star = p_star_cuts_vector_ME.at(iLambda).Angle(pBar_star_cuts_vector_ME.at(iLambdaBar).Vect());
      
      //cout<<L_pT_bin_vector.at(iLambda)<<" "<<Lbar_pT_bin_vector.at(iLambdaBar)<<endl;
      
      L0_L0bar_cosThetaProdPlane_ME_cuts->Fill(TMath::Cos(theta_star));          
      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[L_pT_bin_cuts_vector_ME.at(iLambda)][Lbar_pT_bin_cuts_vector_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[L_eta_bin_cuts_vector_ME.at(iLambda)][Lbar_eta_bin_cuts_vector_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));   
    }   
  
  }
  
  //L0-L0
  for(unsigned int iLambda1 = 0; iLambda1 < p_star_cuts_vector_ME.size(); iLambda1++)
  {
    for(unsigned int iLambda2 = iLambda1+1; iLambda2 < p_star_cuts_vector_ME.size(); iLambda2++)
    {
      double theta_star = p_star_cuts_vector_ME.at(iLambda1).Angle(p_star_cuts_vector_ME.at(iLambda2).Vect());
      
      L0_L0_cosThetaProdPlane_ME_cuts->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[L_pT_bin_cuts_vector_ME.at(iLambda1)][L_pT_bin_cuts_vector_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[L_eta_bin_cuts_vector_ME.at(iLambda1)][L_eta_bin_cuts_vector_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star));
    }   
  
  }    
  
  //L0bar-L0bar
  for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pBar_star_cuts_vector_ME.size(); iLambdaBar1++)
  {
    for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < pBar_star_cuts_vector_ME.size(); iLambdaBar2++)
    {
      double theta_star = pBar_star_cuts_vector_ME.at(iLambdaBar1).Angle(pBar_star_cuts_vector_ME.at(iLambdaBar2).Vect());
      
      L0bar_L0bar_cosThetaProdPlane_ME_cuts->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[Lbar_pT_bin_cuts_vector_ME.at(iLambdaBar1)][Lbar_pT_bin_cuts_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[Lbar_eta_bin_cuts_vector_ME.at(iLambdaBar1)][Lbar_eta_bin_cuts_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
    }   
  
  }  
  //_____________________________________________________________________________________________________________________________________________________________________
  
  //mixed event with pi p pairs
  //US paired with US
  
  //L-Lbar before cuts
  for(unsigned int iLambda = 0; iLambda < p_star_vector_US_ME.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_US_ME.size(); iLambdaBar++)
    {
      double theta_star = p_star_vector_US_ME.at(iLambda).Angle(pBar_star_vector_US_ME.at(iLambdaBar).Vect());
           
      L0_L0bar_cosThetaProdPlane_US_ME->Fill(TMath::Cos(theta_star));          
      L0_L0bar_cosThetaProdPlane_US_ME_pT_hist[L_pT_bin_vector_US_ME.at(iLambda)][Lbar_pT_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_US_ME_eta_hist[L_eta_bin_vector_US_ME.at(iLambda)][Lbar_eta_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      
    }   
  
  }
  
  //L-Lbar after cuts
  for(unsigned int iLambda = 0; iLambda < p_star_vector_US_ME_cuts.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_US_ME_cuts.size(); iLambdaBar++)
    {
      double theta_star = p_star_vector_US_ME_cuts.at(iLambda).Angle(pBar_star_vector_US_ME_cuts.at(iLambdaBar).Vect());

      L0_L0bar_cosThetaProdPlane_US_ME_cuts->Fill(TMath::Cos(theta_star));          
      L0_L0bar_cosThetaProdPlane_US_ME_pT_cuts_hist[L_pT_bin_vector_US_ME.at(iLambda)][Lbar_pT_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_US_ME_eta_cuts_hist[L_eta_bin_vector_US_ME.at(iLambda)][Lbar_eta_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));      
    
    }   
  
  }
  
  //------------------------------------------------------------------------------------------------------------
  
  //L0-L0 before cuts
  for(unsigned int iLambda1 = 0; iLambda1 < p_star_vector_US_ME.size(); iLambda1++)
  {
    for(unsigned int iLambda2 = iLambda1+1; iLambda2 < p_star_vector_US_ME.size(); iLambda2++)
    {
      double theta_star = p_star_vector_US_ME.at(iLambda1).Angle(p_star_vector_US_ME.at(iLambda2).Vect());
      
      L0_L0_cosThetaProdPlane_US_ME->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_US_ME_pT_hist[L_pT_bin_vector_US_ME.at(iLambda1)][L_pT_bin_vector_US_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_US_ME_eta_hist[L_eta_bin_vector_US_ME.at(iLambda1)][L_eta_bin_vector_US_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star));

    }   
  
  }
    
  //L0-L0 after cuts
  for(unsigned int iLambda1 = 0; iLambda1 < p_star_vector_US_ME_cuts.size(); iLambda1++)
  {
    for(unsigned int iLambda2 = iLambda1+1; iLambda2 < p_star_vector_US_ME_cuts.size(); iLambda2++)
    {
      double theta_star = p_star_vector_US_ME_cuts.at(iLambda1).Angle(p_star_vector_US_ME_cuts.at(iLambda2).Vect());
           
      L0_L0_cosThetaProdPlane_US_ME_cuts->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_US_ME_pT_cuts_hist[L_pT_bin_vector_US_ME.at(iLambda1)][L_pT_bin_vector_US_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star));
      L0_L0_cosThetaProdPlane_US_ME_eta_cuts_hist[L_eta_bin_vector_US_ME.at(iLambda1)][L_eta_bin_vector_US_ME.at(iLambda2)]->Fill(TMath::Cos(theta_star));
      
    }   
  
  }
  
  //------------------------------------------------------------------------------------------------------------
  
  //L0bar-L0bar before cuts
  for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pBar_star_vector_US_ME.size(); iLambdaBar1++)
  {
    for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < pBar_star_vector_US_ME.size(); iLambdaBar2++)
    {
      double theta_star = pBar_star_vector_US_ME.at(iLambdaBar1).Angle(pBar_star_vector_US_ME.at(iLambdaBar2).Vect());
      
      L0bar_L0bar_cosThetaProdPlane_US_ME->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_US_ME_pT_hist[Lbar_pT_bin_vector_US_ME.at(iLambdaBar1)][Lbar_pT_bin_vector_US_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_US_ME_eta_hist[Lbar_eta_bin_vector_US_ME.at(iLambdaBar1)][Lbar_eta_bin_vector_US_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
      
    }   
  
  }
  
  //L0bar-L0bar after cuts
  for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pBar_star_vector_US_ME_cuts.size(); iLambdaBar1++)
  {
    for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < pBar_star_vector_US_ME_cuts.size(); iLambdaBar2++)
    {
      double theta_star = pBar_star_vector_US_ME_cuts.at(iLambdaBar1).Angle(pBar_star_vector_US_ME_cuts.at(iLambdaBar2).Vect());

      L0bar_L0bar_cosThetaProdPlane_US_ME_cuts->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_US_ME_pT_cuts_hist[Lbar_pT_bin_vector_US_ME.at(iLambdaBar1)][Lbar_pT_bin_vector_US_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_US_ME_eta_cuts_hist[Lbar_eta_bin_vector_US_ME.at(iLambdaBar1)][Lbar_eta_bin_vector_US_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));   
      
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
      float theta_star = p_star_vector_US_ME.at(iLambda).Angle(pBar_star_vector_LS_ME.at(iLambdaBar).Vect());
      
      L0_L0bar_cosThetaProdPlane_US_LS_ME->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_US_LS_ME_pT_hist[L_pT_bin_vector_US_ME.at(iLambda)][Lbar_pT_bin_vector_LS_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_US_LS_ME_eta_hist[L_eta_bin_vector_US_ME.at(iLambda)][Lbar_eta_bin_vector_LS_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
    
    }
  
  } 
  
  //L - LS, Lbar - US  
  for( unsigned int iLambda = 0; iLambda < p_star_vector_LS_ME.size(); iLambda++)
  {
    for( unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_US_ME.size(); iLambdaBar++)
    {    
      float theta_star = p_star_vector_LS_ME.at(iLambda).Angle(pBar_star_vector_US_ME.at(iLambdaBar).Vect());
      
      L0_L0bar_cosThetaProdPlane_US_LS_ME->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_US_LS_ME_pT_hist[L_pT_bin_vector_LS_ME.at(iLambda)][Lbar_pT_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_US_LS_ME_eta_hist[L_eta_bin_vector_LS_ME.at(iLambda)][Lbar_eta_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
    
    }
  
  }
    
  //L-Lbar after cuts
  //L - US, Lbar - LS 
  for( unsigned int iLambda = 0; iLambda < p_star_vector_US_ME_cuts.size(); iLambda++)
  {
    for( unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_LS_ME_cuts.size(); iLambdaBar++)
    {     
      float theta_star = p_star_vector_US_ME_cuts.at(iLambda).Angle(pBar_star_vector_LS_ME_cuts.at(iLambdaBar).Vect());

      L0_L0bar_cosThetaProdPlane_US_LS_ME_cuts->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[L_pT_bin_vector_US_ME.at(iLambda)][Lbar_pT_bin_vector_LS_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[L_eta_bin_vector_US_ME.at(iLambda)][Lbar_eta_bin_vector_LS_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));          
 
    }
  
  } 
  
  //L - LS, Lbar - US  
  for( unsigned int iLambda = 0; iLambda < p_star_vector_LS_ME_cuts.size(); iLambda++)
  {
    for( unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_US_ME_cuts.size(); iLambdaBar++)
    {    
      float theta_star = p_star_vector_LS_ME_cuts.at(iLambda).Angle(pBar_star_vector_US_ME_cuts.at(iLambdaBar).Vect());

      L0_L0bar_cosThetaProdPlane_US_LS_ME_cuts->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[L_pT_bin_vector_LS_ME.at(iLambda)][Lbar_pT_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
      L0_L0bar_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[L_eta_bin_vector_LS_ME.at(iLambda)][Lbar_eta_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));      
    
    }
  
  } 
  
  //------------------------------------------------------------------------------------------------------------------------
  
  //L-L before cuts
  int nFills_L_L_US_LS_ME = 0;  
  
  for( unsigned int iLambda = 0; iLambda < p_star_vector_US_ME.size(); iLambda++ )
  {
    for( unsigned int iLambdaBck = 0; iLambdaBck < p_star_vector_LS_ME.size(); iLambdaBck++ )
    {
            
      //interate order of US and LS for L
      if( nFills_L_L_US_LS_ME % 2 == 0 )
      {
        float theta_star = p_star_vector_US_ME.at(iLambda).Angle(p_star_vector_LS_ME.at(iLambdaBck).Vect());
                  
        L0_L0_cosThetaProdPlane_US_LS_ME->Fill(TMath::Cos(theta_star));
        L0_L0_cosThetaProdPlane_US_LS_ME_pT_hist[L_pT_bin_vector_US_ME.at(iLambda)][L_pT_bin_vector_LS_ME.at(iLambdaBck)]->Fill(TMath::Cos(theta_star));
        L0_L0_cosThetaProdPlane_US_LS_ME_eta_hist[L_eta_bin_vector_US_ME.at(iLambda)][L_eta_bin_vector_LS_ME.at(iLambdaBck)]->Fill(TMath::Cos(theta_star));
     
      
      }
      else
      {
        float theta_star = p_star_vector_LS_ME.at(iLambdaBck).Angle(p_star_vector_US_ME.at(iLambda).Vect());
                  
        L0_L0_cosThetaProdPlane_US_LS_ME->Fill(TMath::Cos(theta_star));
        L0_L0_cosThetaProdPlane_US_LS_ME_pT_hist[L_pT_bin_vector_LS_ME.at(iLambdaBck)][L_pT_bin_vector_US_ME.at(iLambda)]->Fill(TMath::Cos(theta_star));
        L0_L0_cosThetaProdPlane_US_LS_ME_eta_hist[L_eta_bin_vector_LS_ME.at(iLambdaBck)][L_eta_bin_vector_US_ME.at(iLambda)]->Fill(TMath::Cos(theta_star));
        
      }//end else
      
      nFills_L_L_US_LS_ME++;
    
    }
  
  }
  
  //L-L before cuts
  int nFills_L_L_US_LS_ME_cuts = 0;  
  
  for( unsigned int iLambda = 0; iLambda < p_star_vector_US_ME_cuts.size(); iLambda++ )
  {
    for( unsigned int iLambdaBck = 0; iLambdaBck < p_star_vector_LS_ME_cuts.size(); iLambdaBck++ )
    {
            
      //interate order of US and LS for L
      if( nFills_L_L_US_LS_ME_cuts % 2 == 0 )
      {
        float theta_star = p_star_vector_US_ME_cuts.at(iLambda).Angle(p_star_vector_LS_ME_cuts.at(iLambdaBck).Vect());

        L0_L0_cosThetaProdPlane_US_LS_ME_cuts->Fill(TMath::Cos(theta_star));
        L0_L0_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[L_pT_bin_vector_US_ME.at(iLambda)][L_pT_bin_vector_LS_ME.at(iLambdaBck)]->Fill(TMath::Cos(theta_star));
        L0_L0_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[L_eta_bin_vector_US_ME.at(iLambda)][L_eta_bin_vector_LS_ME.at(iLambdaBck)]->Fill(TMath::Cos(theta_star));          
               
      
      }
      else
      {
        float theta_star = p_star_vector_LS_ME_cuts.at(iLambdaBck).Angle(p_star_vector_US_ME_cuts.at(iLambda).Vect());

        L0_L0_cosThetaProdPlane_US_LS_ME_cuts->Fill(TMath::Cos(theta_star));
        L0_L0_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[L_pT_bin_vector_LS_ME.at(iLambdaBck)][L_pT_bin_vector_US_ME.at(iLambda)]->Fill(TMath::Cos(theta_star));
        L0_L0_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[L_eta_bin_vector_LS_ME.at(iLambdaBck)][L_eta_bin_vector_US_ME.at(iLambda)]->Fill(TMath::Cos(theta_star));         
        
      }//end else
      
      nFills_L_L_US_LS_ME_cuts++;
    
    }
  
  }
  
  //-----------------------------------------------------------------------------------------------------------
  
  //Lbar-Lbar
  int nFills_Lbar_Lbar_US_LS_ME = 0;
  
  for( unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_US_ME.size(); iLambdaBar++ )
  {
    for( unsigned int iLambdaBarBck = 0; iLambdaBarBck < pBar_star_vector_LS_ME.size(); iLambdaBarBck++ )
    {     
      //interate order of US and LS for L
      if( nFills_Lbar_Lbar_US_LS_ME % 2 == 0 )
      {
        float theta_star = pBar_star_vector_US_ME.at(iLambdaBar).Angle(pBar_star_vector_LS_ME.at(iLambdaBarBck).Vect());
                  
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME->Fill(TMath::Cos(theta_star));
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_pT_hist[Lbar_pT_bin_vector_US_ME.at(iLambdaBar)][Lbar_pT_bin_vector_LS_ME.at(iLambdaBarBck)]->Fill(TMath::Cos(theta_star));
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_eta_hist[Lbar_eta_bin_vector_US_ME.at(iLambdaBar)][Lbar_eta_bin_vector_LS_ME.at(iLambdaBarBck)]->Fill(TMath::Cos(theta_star));
        
      
      }
      else
      {
        float theta_star = pBar_star_vector_LS_ME.at(iLambdaBarBck).Angle(pBar_star_vector_US_ME.at(iLambdaBar).Vect());
                  
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME->Fill(TMath::Cos(theta_star));
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_pT_hist[Lbar_pT_bin_vector_LS_ME.at(iLambdaBarBck)][Lbar_pT_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_eta_hist[Lbar_eta_bin_vector_LS_ME.at(iLambdaBarBck)][Lbar_eta_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
        
      }//end else
      
      nFills_Lbar_Lbar_US_LS_ME++;
    
    }
  
  }
  
  
  //Lbar-Lbar
  int nFills_Lbar_Lbar_US_LS_ME_cuts = 0;
  
  for( unsigned int iLambdaBar = 0; iLambdaBar < pBar_star_vector_US_ME_cuts.size(); iLambdaBar++ )
  {
    for( unsigned int iLambdaBarBck = 0; iLambdaBarBck < pBar_star_vector_LS_ME_cuts.size(); iLambdaBarBck++ )
    {     
      //interate order of US and LS for L
      if( nFills_Lbar_Lbar_US_LS_ME_cuts % 2 == 0 )
      {
        float theta_star = pBar_star_vector_US_ME_cuts.at(iLambdaBar).Angle(pBar_star_vector_LS_ME_cuts.at(iLambdaBarBck).Vect());        

        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_cuts->Fill(TMath::Cos(theta_star));
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[Lbar_pT_bin_vector_US_ME.at(iLambdaBar)][Lbar_pT_bin_vector_LS_ME.at(iLambdaBarBck)]->Fill(TMath::Cos(theta_star));
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[Lbar_eta_bin_vector_US_ME.at(iLambdaBar)][Lbar_eta_bin_vector_LS_ME.at(iLambdaBarBck)]->Fill(TMath::Cos(theta_star));              
      
      }
      else
      {
        float theta_star = pBar_star_vector_LS_ME_cuts.at(iLambdaBarBck).Angle(pBar_star_vector_US_ME_cuts.at(iLambdaBar).Vect());
                  
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_cuts->Fill(TMath::Cos(theta_star));
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[Lbar_pT_bin_vector_LS_ME.at(iLambdaBarBck)][Lbar_pT_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
        L0bar_L0bar_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[Lbar_eta_bin_vector_LS_ME.at(iLambdaBarBck)][Lbar_eta_bin_vector_US_ME.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));          
        
      }//end else
      
      nFills_Lbar_Lbar_US_LS_ME_cuts++;
    
    }
  
  }
  
  //_________________________________________________________________________________________________________

  
  
  outFile->cd();
  outFile->Write();
  outFile->Close();

  // Statistics on event generation.
  pythia.stat();
  
  //outFile->Close();

  // Done.
  return 0;
}
