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
  TH1D *L_decayL_hist = new TH1D("L_decayL_hist", "L_decayL_hist", 100, 0, 250);
  
  TH1D* L0_pt_hist = new TH1D("L0_pt_hist","L0_pt_hist",100,0,10);
	TH1D* L0_mass_hist = new TH1D("L0_mass_hist","L0_mass_hist",100,1.,1.2);
	
	
	
  TH1D *L0_thetaProdPlane[nPtBins+1][nEtaBins+1];
  TH1D *L0_cosThetaProdPlane[nPtBins+1][nEtaBins+1];
  
  TH2D *L0_y_vs_p_eta[nPtBins+1];
  TH2D *L0_y_vs_pi_eta[nPtBins+1];
  
  TH1D *L0_pz[nPtBins+1][nEtaBins+1];
  TH1D *L0_xF[nPtBins+1][nEtaBins+1];
  
  TH2D *L0_pT_vs_L0_pz[nEtaBins+1];  
  
  TH2D *L0_p_eta_vs_pi_eta[nPtBins+1];  
  

  TH1D* L0bar_pt_hist = new TH1D("L0bar_pt_hist","L0bar_pt_hist",100,0,10);
	TH1D* L0bar_mass_hist = new TH1D("L0bar_mass_hist","L0bar_mass_hist",100,1.,1.2);

  TH1D *L0bar_thetaProdPlane[nPtBins+1][nEtaBins+1];
  TH1D *L0bar_cosThetaProdPlane[nPtBins+1][nEtaBins+1];
  
  TH2D *L0bar_y_vs_p_eta[nPtBins+1];
  TH2D *L0bar_y_vs_pi_eta[nPtBins+1];  
  
  TH1D *L0bar_pz[nPtBins+1][nEtaBins+1];
  TH1D *L0bar_xF[nPtBins+1][nEtaBins+1];
  
  TH2D *L0bar_pT_vs_L0bar_pz[nEtaBins+1];
  
  TH2D *L0bar_p_eta_vs_pi_eta[nPtBins+1];
  
  
  
  TH1D *L0_L0bar_cosThetaProdPlane = new TH1D("L0_L0bar_cosThetaProdPlane", "L0_L0bar_cosThetaProdPlane", 10, -1, 1);
  
  TH1D *L0_L0bar_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1D *L0_L0bar_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1D *L0_L0_cosThetaProdPlane = new TH1D("L0_L0_cosThetaProdPlane", "L0_L0_cosThetaProdPlane", 10, -1, 1);
  
  TH1D *L0_L0_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1D *L0_L0_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1D *L0bar_L0bar_cosThetaProdPlane = new TH1D("L0bar_L0bar_cosThetaProdPlane", "L0bar_L0bar_cosThetaProdPlane", 10, -1, 1);
  
  TH1D *L0bar_L0bar_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];
 
  TH1D *L0bar_L0bar_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];
  
  //mixed event histograms
  TH1D *L0_L0bar_cosThetaProdPlane_ME = new TH1D("L0_L0bar_cosThetaProdPlane_ME", "L0_L0bar_cosThetaProdPlane_ME", 10, -1, 1);
  
  TH1D *L0_L0bar_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1D *L0_L0bar_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1D *L0_L0_cosThetaProdPlane_ME = new TH1D("L0_L0_cosThetaProdPlane_ME", "L0_L0_cosThetaProdPlane_ME", 10, -1, 1);
  
  TH1D *L0_L0_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1D *L0_L0_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];
  
  
  TH1D *L0bar_L0bar_cosThetaProdPlane_ME = new TH1D("L0bar_L0bar_cosThetaProdPlane_ME", "L0bar_L0bar_cosThetaProdPlane_ME", 10, -1, 1);
  
  TH1D *L0bar_L0bar_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr];
 
  TH1D *L0bar_L0bar_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];
  
 
  //L-L correlation QA histograms
  TH2D *L0_L0bar_delta_eta_vs_delta_phi_hist = new TH2D("L0_L0bar_delta_eta_vs_delta_phi_hist", "L0_L0bar_delta_eta_vs_delta_phi_hist", 100, 0, 2, 100, 0, TMath::TwoPi());

  TH2D *L0_L0bar_y1_vs_y2_hist = new TH2D("L0_L0bar_y1_vs_y2_hist", "L0_L0bar_y1_vs_y2_hist", 100, -1, 1, 100, -1, 1);

  TH2D *L0_L0bar_pT1_vs_pT2_hist = new TH2D("L0_L0bar_pT1_vs_pT2_hist", "L0_L0bar_pT1_vs_pT2_hist", 50, 0, 5, 50, 0, 5);

  TH2D *L0_L0bar_phi1_vs_phi2_hist = new TH2D("L0_L0bar_phi1_vs_phi2_hist", "L0_L0bar_phi1_vs_phi2_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());
  
  
  TH2D *L0_L0_delta_eta_vs_delta_phi_hist = new TH2D("L0_L0_delta_eta_vs_delta_phi_hist", "L0_L0_delta_eta_vs_delta_phi_hist", 100, 0, 2, 100, 0, TMath::TwoPi());

  TH2D *L0_L0_delta_eta_vs_delta_phi_zoom_hist = new TH2D("L0_L0_delta_eta_vs_delta_phi_zoom_hist", "L0_L0_delta_eta_vs_delta_phi_zoom_hist", 100, 0, 0.2, 100, 0, 0.2);

  TH2D *L0_L0_y1_vs_y2_hist = new TH2D("L0_L0_y1_vs_y2_hist", "L0_L0_y1_vs_y2_hist", 100, -1, 1, 100, -1, 1);

  TH2D *L0_p_L0_p_y1_vs_y2_hist = new TH2D("L0_p_L0_p_y1_vs_y2_hist", "L0_p_L0_p_y1_vs_y2_hist", 100, -1, 1, 100, -1, 1);

  TH2D *L0_pi_L0_pi_y1_vs_y2_hist = new TH2D("L0_pi_L0_pi_y1_vs_y2_hist", "L0_pi_L0_pi_y1_vs_y2_hist", 100, -1, 1, 100, -1, 1);

  TH2D *L0_L0_pT1_vs_pT2_hist = new TH2D("L0_L0_pT1_vs_pT2_hist", "L0_L0_pT1_vs_pT2_hist", 50, 0, 5, 50, 0, 5);

  TH2D *L0_L0_phi1_vs_phi2_hist = new TH2D("L0_L0_phi1_vs_phi2_hist", "L0_L0_phi1_vs_phi2_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());

  TH2D *L0_p_L0_p_phi1_vs_phi2_hist = new TH2D("L0_p_L0_p_phi1_vs_phi2_hist", "L0_p_L0_p_phi1_vs_phi2_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());

  TH2D *L0_pi_L0_pi_phi1_vs_phi2_hist = new TH2D("L0_pi_L0_pi_phi1_vs_phi2_hist", "L0_pi_L0_pi_phi1_vs_phi2_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());
  
  
  TH2D *L0bar_L0bar_delta_eta_vs_delta_phi_hist = new TH2D("L0bar_L0bar_delta_eta_vs_delta_phi_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_hist", 100, 0, 2, 100, 0, TMath::TwoPi());

  TH2D *L0bar_L0bar_delta_eta_vs_delta_phi_zoom_hist = new TH2D("L0bar_L0bar_delta_eta_vs_delta_phi_zoom_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_zoom_hist", 100, 0, 0.2, 100, 0, 0.2);

  TH2D *L0bar_L0bar_y1_vs_y2_hist = new TH2D("L0bar_L0bar_y1_vs_y2_hist", "L0bar_L0bar_y1_vs_y2_hist", 100, -1, 1, 100, -1, 1);

  TH2D *L0bar_p_L0bar_p_y1_vs_y2_hist = new TH2D("L0bar_p_L0bar_p_y1_vs_y2_hist", "L0bar_p_L0bar_p_y1_vs_y2_hist", 100, -1, 1, 100, -1, 1);

  TH2D *L0bar_pi_L0bar_pi_y1_vs_y2_hist = new TH2D("L0bar_pi_L0bar_pi_y1_vs_y2_hist", "L0bar_pi_L0bar_pi_y1_vs_y2_hist", 100, -1, 1, 100, -1, 1);

  TH2D *L0bar_L0bar_pT1_vs_pT2_hist = new TH2D("L0bar_L0bar_pT1_vs_pT2_hist", "L0bar_L0bar_pT1_vs_pT2_hist", 50, 0, 5, 50, 0, 5);

  TH2D *L0bar_L0bar_phi1_vs_phi2_hist = new TH2D("L0bar_L0bar_phi1_vs_phi2_hist", "L0bar_L0bar_phi1_vs_phi2_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());

  TH2D *L0bar_p_L0bar_p_phi1_vs_phi2_hist = new TH2D("L0bar_p_L0bar_p_phi1_vs_phi2_hist", "L0bar_p_L0bar_p_phi1_vs_phi2_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());

  TH2D *L0bar_pi_L0bar_pi_phi1_vs_phi2_hist = new TH2D("L0bar_pi_L0bar_pi_phi1_vs_phi2_hist", "L0bar_pi_L0bar_pi_phi1_vs_phi2_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());
  //_____________________________

  //histograms after cuts	
  TH1D *L0_thetaProdPlane_cuts[nPtBins+1][nEtaBins+1];
  TH1D *L0_cosThetaProdPlane_cuts[nPtBins+1][nEtaBins+1];
  
  TH2D *L0_y_vs_p_eta_cuts[nPtBins+1];
  TH2D *L0_y_vs_pi_eta_cuts[nPtBins+1];
  
  TH1D *L0_pz_cuts[nPtBins+1][nEtaBins+1];
  TH1D *L0_xF_cuts[nPtBins+1][nEtaBins+1];
  
  TH2D *L0_pT_vs_L0_pz_cuts[nEtaBins+1];  
  
  TH2D *L0_p_eta_vs_pi_eta_cuts[nPtBins+1];
  
  

  TH1D *L0bar_thetaProdPlane_cuts[nPtBins+1][nEtaBins+1];
  TH1D *L0bar_cosThetaProdPlane_cuts[nPtBins+1][nEtaBins+1];
  
  TH2D *L0bar_y_vs_p_eta_cuts[nPtBins+1];
  TH2D *L0bar_y_vs_pi_eta_cuts[nPtBins+1];
  
  TH1D *L0bar_pz_cuts[nPtBins+1][nEtaBins+1];
  TH1D *L0bar_xF_cuts[nPtBins+1][nEtaBins+1];
  
  TH2D *L0bar_pT_vs_L0bar_pz_cuts[nEtaBins+1];
  
  TH2D *L0bar_p_eta_vs_pi_eta_cuts[nPtBins+1];
  
  
  TH1D *L0_L0bar_cosThetaProdPlane_cuts = new TH1D("L0_L0bar_cosThetaProdPlane_cuts", "L0_L0bar_cosThetaProdPlane_cuts", 10, -1, 1);
  
  TH1D *L0_L0bar_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1D *L0_L0bar_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1D *L0_L0_cosThetaProdPlane_cuts = new TH1D("L0_L0_cosThetaProdPlane_cuts", "L0_L0_cosThetaProdPlane_cuts", 10, -1, 1);
  
  TH1D *L0_L0_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1D *L0_L0_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1D *L0bar_L0bar_cosThetaProdPlane_cuts = new TH1D("L0bar_L0bar_cosThetaProdPlane_cuts", "L0bar_L0bar_cosThetaProdPlane_cuts", 10, -1, 1);  
  
  TH1D *L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0bar_L0bar_cosThetaProdPlane_pT_cuts_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1D *L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];
  
  //mixed event histograms
  TH1D *L0_L0bar_cosThetaProdPlane_ME_cuts = new TH1D("L0_L0bar_cosThetaProdPlane_ME_cuts", "L0_L0bar_cosThetaProdPlane_ME_cuts", 10, -1, 1);
  
  TH1D *L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1D *L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1D *L0_L0_cosThetaProdPlane_ME_cuts = new TH1D("L0_L0_cosThetaProdPlane_ME_cuts", "L0_L0_cosThetaProdPlane_ME_cuts", 10, -1, 1);
  
  TH1D *L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1D *L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1D *L0bar_L0bar_cosThetaProdPlane_ME_cuts = new TH1D("L0bar_L0bar_cosThetaProdPlane_ME_cuts", "L0bar_L0bar_cosThetaProdPlane_ME_cuts", 10, -1, 1);
  
  TH1D *L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
 
  TH1D *L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[nEtaBins][nEtaBins];
 
  
 
  //QA histograms
  TH2D *L0_L0bar_delta_eta_vs_delta_phi_cuts_hist = new TH2D("L0_L0bar_delta_eta_vs_delta_phi_cuts_hist", "L0_L0bar_delta_eta_vs_delta_phi_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());

  TH2D *L0_L0bar_y1_vs_y2_cuts_hist = new TH2D("L0_L0bar_y1_vs_y2_cuts_hist", "L0_L0bar_y1_vs_y2_cuts_hist", 100, -1, 1, 100, -1, 1);

  TH2D *L0_L0bar_pT1_vs_pT2_cuts_hist = new TH2D("L0_L0bar_pT1_vs_pT2_cuts_hist", "L0_L0bar_pT1_vs_pT2_cuts_hist", 50, 0, 5, 50, 0, 5);

  TH2D *L0_L0bar_phi1_vs_phi2_cuts_hist = new TH2D("L0_L0bar_phi1_vs_phi2_cuts_hist", "L0_L0bar_phi1_vs_phi2_cuts_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());
  
  
  TH2D *L0_L0_delta_eta_vs_delta_phi_cuts_hist = new TH2D("L0_L0_delta_eta_vs_delta_phi_cuts_hist", "L0_L0_delta_eta_vs_delta_phi_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());

  TH2D *L0_L0_delta_eta_vs_delta_phi_zoom_cuts_hist = new TH2D("L0_L0_delta_eta_vs_delta_phi_zoom_cuts_hist", "L0_L0_delta_eta_vs_delta_phi_zoom_cuts_hist", 100, 0, 0.2, 100, 0, 0.2);

  TH2D *L0_L0_y1_vs_y2_cuts_hist = new TH2D("L0_L0_y1_vs_y2_cuts_hist", "L0_L0_y1_vs_y2_cuts_hist", 100, -1, 1, 100, -1, 1);

  TH2D *L0_p_L0_p_y1_vs_y2_cuts_hist = new TH2D("L0_p_L0_p_y1_vs_y2_cuts_hist", "L0_p_L0_p_y1_vs_y2_cuts_hist", 100, -1, 1, 100, -1, 1);

  TH2D *L0_pi_L0_pi_y1_vs_y2_cuts_hist = new TH2D("L0_pi_L0_pi_y1_vs_y2_cuts_hist", "L0_pi_L0_pi_y1_vs_y2_cuts_hist", 100, -1, 1, 100, -1, 1);

  TH2D *L0_L0_pT1_vs_pT2_cuts_hist = new TH2D("L0_L0_pT1_vs_pT2_cuts_hist", "L0_L0_pT1_vs_pT2_cuts_hist", 50, 0, 5, 50, 0, 5);

  TH2D *L0_L0_phi1_vs_phi2_cuts_hist = new TH2D("L0_L0_phi1_vs_phi2_cuts_hist", "L0_L0_phi1_vs_phi2_cuts_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());

  TH2D *L0_p_L0_p_phi1_vs_phi2_cuts_hist = new TH2D("L0_p_L0_p_phi1_vs_phi2_cuts_hist", "L0_p_L0_p_phi1_vs_phi2_cuts_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());

  TH2D *L0_pi_L0_pi_phi1_vs_phi2_cuts_hist = new TH2D("L0_pi_L0_pi_phi1_vs_phi2_cuts_hist", "L0_pi_L0_pi_phi1_vs_phi2_cuts_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());
  
  
  TH2D *L0bar_L0bar_delta_eta_vs_delta_phi_cuts_hist = new TH2D("L0bar_L0bar_delta_eta_vs_delta_phi_cuts_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_cuts_hist", 100, 0, 2, 100, 0, TMath::TwoPi());

  TH2D *L0bar_L0bar_delta_eta_vs_delta_phi_zoom_cuts_hist = new TH2D("L0bar_L0bar_delta_eta_vs_delta_phi_zoom_cuts_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_zoom_cuts_hist", 100, 0, 0.2, 100, 0, 0.2);

  TH2D *L0bar_L0bar_y1_vs_y2_cuts_hist = new TH2D("L0bar_L0bar_y1_vs_y2_cuts_hist", "L0bar_L0bar_y1_vs_y2_cuts_hist", 100, -1, 1, 100, -1, 1);

  TH2D *L0bar_p_L0bar_p_y1_vs_y2_cuts_hist = new TH2D("L0bar_p_L0bar_p_y1_vs_y2_cuts_hist", "L0bar_p_L0bar_p_y1_vs_y2_cuts_hist", 100, -1, 1, 100, -1, 1);

  TH2D *L0bar_pi_L0bar_pi_y1_vs_y2_cuts_hist = new TH2D("L0bar_pi_L0bar_pi_y1_vs_y2_cuts_hist", "L0bar_pi_L0bar_pi_y1_vs_y2_cuts_hist", 100, -1, 1, 100, -1, 1);

  TH2D *L0bar_L0bar_pT1_vs_pT2_cuts_hist = new TH2D("L0bar_L0bar_pT1_vs_pT2_cuts_hist", "L0bar_L0bar_pT1_vs_pT2_cuts_hist", 50, 0, 5, 50, 0, 5);

  TH2D *L0bar_L0bar_phi1_vs_phi2_cuts_hist = new TH2D("L0bar_L0bar_phi1_vs_phi2_cuts_hist", "L0bar_L0bar_phi1_vs_phi2_cuts_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());

  TH2D *L0bar_p_L0bar_p_phi1_vs_phi2_cuts_hist = new TH2D("L0bar_p_L0bar_p_phi1_vs_phi2_cuts_hist", "L0bar_p_L0bar_p_phi1_vs_phi2_cuts_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());

  TH2D *L0bar_pi_L0bar_pi_phi1_vs_phi2_cuts_hist = new TH2D("L0bar_pi_L0bar_pi_phi1_vs_phi2_cuts_hist", "L0bar_pi_L0bar_pi_phi1_vs_phi2_cuts_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());
  //_____________________________
  

  
  for(unsigned int pTbin = 0; pTbin < nPtBins+1; pTbin++)
  {
  
    L0_y_vs_p_eta[pTbin] = new TH2D(Form("L0_y_vs_p_eta_pT_%i", pTbin), Form("L0_y_vs_p_eta_pT_%i", pTbin), 100, -1, 1, 500, -5, 5);
    L0_y_vs_pi_eta[pTbin] = new TH2D(Form("L0_y_vs_pi_eta_pT_%i", pTbin), Form("L0_y_vs_pi_eta_pT_%i", pTbin), 100, -1, 1, 500, -5, 5);
    
    L0_p_eta_vs_pi_eta[pTbin] = new TH2D(Form("L0_p_eta_vs_pi_eta_L0_pT_%i", pTbin), Form("L0_p_eta_vs_pi_eta_L0_pT_%i", pTbin), 500, -5, 5, 500, -5, 5);
    
    L0bar_y_vs_p_eta[pTbin] = new TH2D(Form("L0bar_y_vs_p_eta_pT_%i", pTbin), Form("L0bar_y_vs_p_eta_pT_%i", pTbin), 100, -1, 1, 500, -5, 5);
    L0bar_y_vs_pi_eta[pTbin] = new TH2D(Form("L0bar_y_vs_pi_eta_pT_%i", pTbin), Form("L0bar_y_vs_pi_eta_pT_%i", pTbin), 100, -1, 1, 500, -5, 5);    
    
    L0bar_p_eta_vs_pi_eta[pTbin] = new TH2D(Form("L0bar_p_eta_vs_pi_eta_L0_pT_%i", pTbin), Form("L0bar_p_eta_vs_pi_eta_L0_pT_%i", pTbin), 500, -5, 5, 500, -5, 5);
    
        
    
    L0_y_vs_p_eta_cuts[pTbin] = new TH2D(Form("L0_y_vs_p_eta_cuts_pT_%i", pTbin), Form("L0_y_vs_p_eta_cuts_pT_%i", pTbin), 100, -1, 1, 500, -5, 5);
    L0_y_vs_pi_eta_cuts[pTbin] = new TH2D(Form("L0_y_vs_pi_eta_cuts_pT_%i", pTbin), Form("L0_y_vs_pi_eta_cuts_pT_%i", pTbin), 100, -1, 1, 500, -5, 5);
    
    L0_p_eta_vs_pi_eta_cuts[pTbin] = new TH2D(Form("L0_p_eta_vs_pi_eta_cuts_L0_pT_%i", pTbin), Form("L0_p_eta_vs_pi_eta_cuts_L0_pT_%i", pTbin), 500, -5, 5, 500, -5, 5);
    
    L0bar_y_vs_p_eta_cuts[pTbin] = new TH2D(Form("L0bar_y_vs_p_eta_cuts_pT_%i", pTbin), Form("L0bar_y_vs_p_eta_cuts_pT_%i", pTbin), 100, -1, 1, 500, -5, 5);
    L0bar_y_vs_pi_eta_cuts[pTbin] = new TH2D(Form("L0bar_y_vs_pi_eta_cuts_pT_%i", pTbin), Form("L0bar_y_vs_pi_eta_cuts_pT_%i", pTbin), 100, -1, 1, 500, -5, 5);
    
    L0bar_p_eta_vs_pi_eta_cuts[pTbin] = new TH2D(Form("L0bar_p_eta_vs_pi_eta_cuts_L0_pT_%i", pTbin), Form("L0bar_p_eta_vs_pi_eta_cuts_L0_pT_%i", pTbin), 500, -5, 5, 500, -5, 5);
 
    for(unsigned int etaBin = 0; etaBin < nEtaBins+1; etaBin++)
    {
      if(pTbin == 0)
      {
      
        L0_pT_vs_L0_pz[etaBin] = new TH2D(Form("L0_pT_vs_L0_pz_eta_%i", etaBin), Form("L0_pT_vs_L0_pz_eta_%i", etaBin), 100, 0, 10, 200, -10, 10);
        
        L0bar_pT_vs_L0bar_pz[etaBin] = new TH2D(Form("L0bar_pT_vs_L0bar_pz_eta_%i", etaBin), Form("L0bar_pT_vs_L0bar_pz_eta_%i", etaBin), 100, 0, 10, 200, -10, 10);
        
        
        L0_pT_vs_L0_pz_cuts[etaBin] = new TH2D(Form("L0_pT_vs_L0_pz_cuts_eta_%i", etaBin), Form("L0_pT_vs_L0_pz_cuts_eta_%i", etaBin), 100, 0, 10, 200, -10, 10);
        
        L0bar_pT_vs_L0bar_pz_cuts[etaBin] = new TH2D(Form("L0bar_pT_vs_L0bar_pz_cuts_eta_%i", etaBin), Form("L0bar_pT_vs_L0bar_pz_cuts_eta_%i", etaBin), 100, 0, 10, 200, -10, 10); 
 
      }

      L0_thetaProdPlane[pTbin][etaBin] = new TH1D(Form("L0_thetaProdPlane_pT_%i_eta_%i", pTbin, etaBin), Form("L0_thetaProdPlane_pT_%i_eta_%i", pTbin, etaBin), 20, 0, TMath::Pi());
      L0_cosThetaProdPlane[pTbin][etaBin] = new TH1D(Form("L0_cosThetaProdPlane_pT_%i_eta_%i", pTbin, etaBin), Form("L0_cosThetaProdPlane_pT_%i_eta_%i", pTbin, etaBin), 20, -1, 1);
      
      L0_pz[pTbin][etaBin] = new TH1D(Form("L0_pz_pT_%i_eta_%i", pTbin, etaBin), Form("L0_pz_pT_%i_eta_%i", pTbin, etaBin), 200,-10, 10);
      L0_xF[pTbin][etaBin] = new TH1D(Form("L0_xF_pT_%i_eta_%i", pTbin, etaBin), Form("L0_xF_pT_%i_eta_%i", pTbin, etaBin), 100, 0, 0.01);

      L0bar_thetaProdPlane[pTbin][etaBin] = new TH1D(Form("L0bar_thetaProdPlane_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_thetaProdPlane_pT_%i_eta_%i", pTbin, etaBin), 20, 0, TMath::Pi());
      L0bar_cosThetaProdPlane[pTbin][etaBin] = new TH1D(Form("L0bar_cosThetaProdPlane_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_cosThetaProdPlane_pT_%i_eta_%i", pTbin, etaBin), 20, -1, 1);
      
      L0bar_pz[pTbin][etaBin] = new TH1D(Form("L0bar_pz_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_pz_pT_%i_eta_%i", pTbin, etaBin), 200, -10, 10);
      L0bar_xF[pTbin][etaBin] = new TH1D(Form("L0bar_xF_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_xF_pT_%i_eta_%i", pTbin, etaBin), 100, 0, 0.01);
      
      
      L0_thetaProdPlane_cuts[pTbin][etaBin] = new TH1D(Form("L0_thetaProdPlane_cuts_pT_%i_eta_%i", pTbin, etaBin), Form("L0_thetaProdPlane_cuts_pT_%i_eta_%i", pTbin, etaBin), 20, 0, TMath::Pi());
      L0_cosThetaProdPlane_cuts[pTbin][etaBin] = new TH1D(Form("L0_cosThetaProdPlane_cuts_pT_%i_eta_%i", pTbin, etaBin), Form("L0_cosThetaProdPlane_cuts_pT_%i_eta_%i", pTbin, etaBin), 20, -1, 1);
      
      L0_pz_cuts[pTbin][etaBin] = new TH1D(Form("L0_pz_cuts_pT_%i_eta_%i", pTbin, etaBin), Form("L0_pz_cuts_pT_%i_eta_%i", pTbin, etaBin), 200,-10, 10);
      L0_xF_cuts[pTbin][etaBin] = new TH1D(Form("L0_xF_cuts_pT_%i_eta_%i", pTbin, etaBin), Form("L0_xF_cuts_pT_%i_eta_%i", pTbin, etaBin), 100, 0, 0.01);

      L0bar_thetaProdPlane_cuts[pTbin][etaBin] = new TH1D(Form("L0bar_thetaProdPlane_cuts_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_thetaProdPlane_cuts_pT_%i_eta_%i", pTbin, etaBin), 20, 0, TMath::Pi());
      L0bar_cosThetaProdPlane_cuts[pTbin][etaBin] = new TH1D(Form("L0bar_cosThetaProdPlane_cuts_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_cosThetaProdPlane_cuts_pT_%i_eta_%i", pTbin, etaBin), 20, -1, 1);
      
      L0bar_pz_cuts[pTbin][etaBin] = new TH1D(Form("L0bar_pz_cuts_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_pz_cuts_pT_%i_eta_%i", pTbin, etaBin), 200, -10, 10);
      L0bar_xF_cuts[pTbin][etaBin] = new TH1D(Form("L0bar_xF_cuts_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_xF_cuts_pT_%i_eta_%i", pTbin, etaBin), 100, 0, 0.01);
    
    }

  }
  
  
  //L-L correlation histograms in bins
  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      L0_L0bar_cosThetaProdPlane_pT_hist[pTbin1][pTbin2] = new TH1D(Form("L0_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2] = new TH1D(Form("L0_L0bar_cosThetaProdPlane_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_pT_hist[pTbin1][pTbin2] = new TH1D(Form("L0_L0_cosThetaProdPlane_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2] = new TH1D(Form("L0_L0_cosThetaProdPlane_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_pT_hist[pTbin1][pTbin2] = new TH1D(Form("L0bar_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2] = new TH1D(Form("L0bar_L0bar_cosThetaProdPlane_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      
      //mixed event
      L0_L0bar_cosThetaProdPlane_ME_pT_hist[pTbin1][pTbin2] = new TH1D(Form("L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1D(Form("L0_L0bar_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_ME_pT_hist[pTbin1][pTbin2] = new TH1D(Form("L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1D(Form("L0_L0_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_ME_pT_hist[pTbin1][pTbin2] = new TH1D(Form("L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2] = new TH1D(Form("L0bar_L0bar_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
    }
  }

  for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
  {
    for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
    {
      L0_L0bar_cosThetaProdPlane_eta_hist[etaBin1][etaBin2] = new TH1D(Form("L0_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2] = new TH1D(Form("L0_L0bar_cosThetaProdPlane_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_eta_hist[etaBin1][etaBin2] = new TH1D(Form("L0_L0_cosThetaProdPlane_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2] = new TH1D(Form("L0_L0_cosThetaProdPlane_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_eta_hist[etaBin1][etaBin2] = new TH1D(Form("L0bar_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2] = new TH1D(Form("L0bar_L0bar_cosThetaProdPlane_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      
      //mixed event
      L0_L0bar_cosThetaProdPlane_ME_eta_hist[etaBin1][etaBin2] = new TH1D(Form("L0_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1D(Form("L0_L0bar_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_ME_eta_hist[etaBin1][etaBin2] = new TH1D(Form("L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1D(Form("L0_L0_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_ME_eta_hist[etaBin1][etaBin2] = new TH1D(Form("L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2] = new TH1D(Form("L0bar_L0bar_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
    }
  }
  //_____________________________________________________________________________________________________________________________________________________________________
  
  
  //mixed event vectors       
  vector<TLorentzVector> p_star_vector_ME;
  vector<TLorentzVector> p_star_cuts_vector_ME;
  
  vector<int> L_pT_bin_vector_ME;
  vector<int> L_pT_bin_cuts_vector_ME;
  
  vector<int> L_eta_bin_vector_ME;
  vector<int> L_eta_bin_cuts_vector_ME;

  
  
  vector<TLorentzVector> pbar_star_vector_ME;
  vector<TLorentzVector> pbar_star_cuts_vector_ME;  
 
  vector<int> Lbar_pT_bin_vector_ME;
  vector<int> Lbar_pT_bin_cuts_vector_ME;
  
  vector<int> Lbar_eta_bin_vector_ME;
  vector<int> Lbar_eta_bin_cuts_vector_ME;  

  

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
    
    vector<TLorentzVector> pbar_vector;
    vector<TLorentzVector> pbar_cuts_vector;
    
    vector<TLorentzVector> pibar_vector;
    vector<TLorentzVector> pibar_cuts_vector;
    
    vector<TLorentzVector> pbar_star_vector;
    vector<TLorentzVector> pbar_star_cuts_vector;
    
    
    

    //loop over aprticles in event
    for (int i = 0; i < pythia.event.size(); ++i)
    {
    
    
      if( fabs(pythia.event[i].id()) == 3122)
      {
        //cout<<"Lambda found!"<<endl;
        
        //find index of decay daughters
        const int daughter1_Id = pythia.event[i].daughter1();
        const int daughter2_Id = pythia.event[i].daughter2();
        
        
        //Lambda fourmomentum
        L_fourmom.SetPxPyPzE(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
        
        //Lambda MC decay veretex
        L_decay_vertex.SetXYZ(pythia.event[i].xDec(), pythia.event[i].yDec(), pythia.event[i].zDec());
        
        //Lambda MC production veretex (PV)
        L_production_vertex.SetXYZ(pythia.event[i].xProd(), pythia.event[i].yProd(), pythia.event[i].zProd());
        
        TVector3 L_decayL_vect = L_decay_vertex - L_production_vertex;
        float L_decayL = L_decayL_vect.Mag();
        
        L_decayL_hist->Fill(L_decayL);
        
        float L_xF = fabs(L_fourmom.Pz())/pythia.info.eCM()*2.;
        
        if( fabs(L_fourmom.Eta()) >= 1. ) continue;
        
        L_fourmom_reverse.SetPxPyPzE(-pythia.event[i].px(), -pythia.event[i].py(), -pythia.event[i].pz(), pythia.event[i].e());
        

        p_fourmom.SetPxPyPzE(pythia.event[daughter1_Id].px(), pythia.event[daughter1_Id].py(), pythia.event[daughter1_Id].pz(), pythia.event[daughter1_Id].e());
        pi_fourmom.SetPxPyPzE(pythia.event[daughter2_Id].px(), pythia.event[daughter2_Id].py(), pythia.event[daughter2_Id].pz(), pythia.event[daughter2_Id].e());
        

        TLorentzVector p_fourmom_star = p_fourmom;
        p_fourmom_star.Boost(L_fourmom_reverse.BoostVector());       
    

        TVector3 beamVector(0.,0.,1.); //unity vector along the beam axis
        TVector3 mProdPlane = beamVector.Cross(L_fourmom.Vect());
        mProdPlane = ( mProdPlane )*(1./mProdPlane.Mag() );

        float mThetaProdPlane = mProdPlane.Angle(p_fourmom_star.Vect());


        
        //fill all histograms for all pT and centrality bins
        int pT_bin = -1;

        //find pT bin of Lambda
        for(int j = 0; j < nPtBins; j++) //loop over pT bins
        {
          if(L_fourmom.Pt() > pT_bins[j] && L_fourmom.Pt() <= pT_bins[j+1])
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
          if(L_fourmom.Pt() > pT_bins_corr[j] && L_fourmom.Pt() <= pT_bins_corr[j+1])
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
          if(L_fourmom.Eta() > eta_bins[j] && L_fourmom.Eta() <= eta_bins[j+1])
          {
            eta_bin = j;
            break; //stop after eta bin is found
          }
        }

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
        
        
          L0_mass_hist->Fill(pythia.event[i].m());
          L0_pt_hist->Fill(L_fourmom.Pt());

          L0_thetaProdPlane[pT_bin][eta_bin]->Fill(mThetaProdPlane);
          L0_thetaProdPlane[nPtBins][eta_bin]->Fill(mThetaProdPlane);  //pT integrated, eta bins
          L0_thetaProdPlane[pT_bin][nEtaBins]->Fill(mThetaProdPlane);  //pT bins, -1 < eta < 1
          L0_thetaProdPlane[nPtBins][nEtaBins]->Fill(mThetaProdPlane); //pT integrated and -1 < eta < 1


          L0_cosThetaProdPlane[pT_bin][eta_bin]->Fill(cos(mThetaProdPlane));
          L0_cosThetaProdPlane[nPtBins][eta_bin]->Fill(cos(mThetaProdPlane));
          L0_cosThetaProdPlane[pT_bin][nEtaBins]->Fill(cos(mThetaProdPlane));
          L0_cosThetaProdPlane[nPtBins][nEtaBins]->Fill(cos(mThetaProdPlane));
          
          
          L0_pz[pT_bin][eta_bin]->Fill(L_fourmom.Pz());
          L0_pz[nPtBins][eta_bin]->Fill(L_fourmom.Pz());
          L0_pz[pT_bin][nEtaBins]->Fill(L_fourmom.Pz());
          L0_pz[nPtBins][nEtaBins]->Fill(L_fourmom.Pz());
          
          L0_xF[pT_bin][eta_bin]->Fill(L_xF);
          L0_xF[nPtBins][eta_bin]->Fill(L_xF);
          L0_xF[pT_bin][nEtaBins]->Fill(L_xF);
          L0_xF[nPtBins][nEtaBins]->Fill(L_xF);          
    
          L0_pT_vs_L0_pz[eta_bin]->Fill(L_fourmom.Pt(), L_fourmom.Pz());
          L0_pT_vs_L0_pz[nEtaBins]->Fill(L_fourmom.Pt(), L_fourmom.Pz());
          
          L0_y_vs_p_eta[pT_bin]->Fill(L_fourmom.Rapidity(), p_fourmom.Eta());
          L0_y_vs_pi_eta[pT_bin]->Fill(L_fourmom.Rapidity(), pi_fourmom.Eta());
          
          L0_y_vs_p_eta[nPtBins]->Fill(L_fourmom.Rapidity(), p_fourmom.Eta());
          L0_y_vs_pi_eta[nPtBins]->Fill(L_fourmom.Rapidity(), pi_fourmom.Eta());
          
          L0_p_eta_vs_pi_eta[pT_bin]->Fill(p_fourmom.Eta(), pi_fourmom.Eta());
          L0_p_eta_vs_pi_eta[nPtBins]->Fill(p_fourmom.Eta(), pi_fourmom.Eta());

        }

        //Lambda-bar
        if( pythia.event[i].id() < 0 )
        {
          Lbar_vector.push_back(L_fourmom);
          Lbar_pT_bin_vector.push_back(pT_bin_corr);
          Lbar_eta_bin_vector.push_back(eta_bin);
          
          pbar_vector.push_back(p_fourmom);
          pibar_vector.push_back(pi_fourmom);
          
          pbar_star_vector.push_back(p_fourmom_star);
          
          L0bar_mass_hist->Fill(pythia.event[i].m());
          L0bar_pt_hist->Fill(L_fourmom.Pt());

          L0bar_thetaProdPlane[pT_bin][eta_bin]->Fill(mThetaProdPlane);
          L0bar_thetaProdPlane[nPtBins][eta_bin]->Fill(mThetaProdPlane);
          L0bar_thetaProdPlane[pT_bin][nEtaBins]->Fill(mThetaProdPlane);
          L0bar_thetaProdPlane[nPtBins][nEtaBins]->Fill(mThetaProdPlane);

          L0bar_cosThetaProdPlane[pT_bin][eta_bin]->Fill(cos(mThetaProdPlane));
          L0bar_cosThetaProdPlane[nPtBins][eta_bin]->Fill(cos(mThetaProdPlane));
          L0bar_cosThetaProdPlane[pT_bin][nEtaBins]->Fill(cos(mThetaProdPlane));
          L0bar_cosThetaProdPlane[nPtBins][nEtaBins]->Fill(cos(mThetaProdPlane));
          
          
          L0bar_pz[pT_bin][eta_bin]->Fill(L_fourmom.Pz());
          L0bar_pz[nPtBins][eta_bin]->Fill(L_fourmom.Pz());
          L0bar_pz[pT_bin][nEtaBins]->Fill(L_fourmom.Pz());
          L0bar_pz[nPtBins][nEtaBins]->Fill(L_fourmom.Pz());
          
          L0bar_xF[pT_bin][eta_bin]->Fill(L_xF);
          L0bar_xF[nPtBins][eta_bin]->Fill(L_xF);
          L0bar_xF[pT_bin][nEtaBins]->Fill(L_xF);
          L0bar_xF[nPtBins][nEtaBins]->Fill(L_xF);
        
          L0bar_pT_vs_L0bar_pz[eta_bin]->Fill(L_fourmom.Pt(), L_fourmom.Pz());
          L0bar_pT_vs_L0bar_pz[nEtaBins]->Fill(L_fourmom.Pt(), L_fourmom.Pz());
          
          L0bar_y_vs_p_eta[pT_bin]->Fill(L_fourmom.Rapidity(), p_fourmom.Eta());
          L0bar_y_vs_pi_eta[pT_bin]->Fill(L_fourmom.Rapidity(), pi_fourmom.Eta());
          
          L0bar_y_vs_p_eta[nPtBins]->Fill(L_fourmom.Rapidity(), p_fourmom.Eta());
          L0bar_y_vs_pi_eta[nPtBins]->Fill(L_fourmom.Rapidity(), pi_fourmom.Eta());
          
          L0bar_p_eta_vs_pi_eta[pT_bin]->Fill(p_fourmom.Eta(), pi_fourmom.Eta());
          L0bar_p_eta_vs_pi_eta[nPtBins]->Fill(p_fourmom.Eta(), pi_fourmom.Eta());
       
          
        }
        
        
        //Lambda cuts
        //decay length cuts in mm (default PYTHIA units)
        if(L_decayL < 20 || L_decayL > 250) continue;
        //if(L_fourmom.Pt() < 0.5) continue; //added pT cut - pT integrated in data starts at 0.5 GeV/c
        if(L_fourmom.Rapidity() >= 1) continue;
        
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

          L0_thetaProdPlane_cuts[pT_bin][eta_bin]->Fill(mThetaProdPlane);
          L0_thetaProdPlane_cuts[nPtBins][eta_bin]->Fill(mThetaProdPlane);  //pT integrated, eta bins
          L0_thetaProdPlane_cuts[pT_bin][nEtaBins]->Fill(mThetaProdPlane);  //pT bins, -1 < eta < 1
          L0_thetaProdPlane_cuts[nPtBins][nEtaBins]->Fill(mThetaProdPlane); //pT integrated and -1 < eta < 1


          L0_cosThetaProdPlane_cuts[pT_bin][eta_bin]->Fill(cos(mThetaProdPlane));
          L0_cosThetaProdPlane_cuts[nPtBins][eta_bin]->Fill(cos(mThetaProdPlane));
          L0_cosThetaProdPlane_cuts[pT_bin][nEtaBins]->Fill(cos(mThetaProdPlane));
          L0_cosThetaProdPlane_cuts[nPtBins][nEtaBins]->Fill(cos(mThetaProdPlane));
          
          
          L0_pz_cuts[pT_bin][eta_bin]->Fill(L_fourmom.Pz());
          L0_pz_cuts[nPtBins][eta_bin]->Fill(L_fourmom.Pz());
          L0_pz_cuts[pT_bin][nEtaBins]->Fill(L_fourmom.Pz());
          L0_pz_cuts[nPtBins][nEtaBins]->Fill(L_fourmom.Pz());
          
          L0_xF_cuts[pT_bin][eta_bin]->Fill(L_xF);
          L0_xF_cuts[nPtBins][eta_bin]->Fill(L_xF);
          L0_xF_cuts[pT_bin][nEtaBins]->Fill(L_xF);
          L0_xF_cuts[nPtBins][nEtaBins]->Fill(L_xF);
       
          L0_pT_vs_L0_pz_cuts[eta_bin]->Fill(L_fourmom.Pt(), L_fourmom.Pz());
          L0_pT_vs_L0_pz_cuts[nEtaBins]->Fill(L_fourmom.Pt(), L_fourmom.Pz());
          
          L0_y_vs_p_eta_cuts[pT_bin]->Fill(L_fourmom.Rapidity(), p_fourmom.Eta());
          L0_y_vs_pi_eta_cuts[pT_bin]->Fill(L_fourmom.Rapidity(), pi_fourmom.Eta());
          
          L0_y_vs_p_eta_cuts[nPtBins]->Fill(L_fourmom.Rapidity(), p_fourmom.Eta());
          L0_y_vs_pi_eta_cuts[nPtBins]->Fill(L_fourmom.Rapidity(), pi_fourmom.Eta());
          
          L0_p_eta_vs_pi_eta_cuts[pT_bin]->Fill(p_fourmom.Eta(), pi_fourmom.Eta());
          L0_p_eta_vs_pi_eta_cuts[nPtBins]->Fill(p_fourmom.Eta(), pi_fourmom.Eta());

        }

        //Lambda-bar
        if( pythia.event[i].id() < 0 )
        {
          Lbar_cuts_vector.push_back(L_fourmom);
          Lbar_pT_bin_cuts_vector.push_back(pT_bin_corr);
          Lbar_eta_bin_cuts_vector.push_back(eta_bin);
          
          pbar_cuts_vector.push_back(p_fourmom);
          pibar_cuts_vector.push_back(pi_fourmom);
        
          pbar_star_cuts_vector.push_back(p_fourmom_star);

          L0bar_thetaProdPlane_cuts[pT_bin][eta_bin]->Fill(mThetaProdPlane);
          L0bar_thetaProdPlane_cuts[nPtBins][eta_bin]->Fill(mThetaProdPlane);
          L0bar_thetaProdPlane_cuts[pT_bin][nEtaBins]->Fill(mThetaProdPlane);
          L0bar_thetaProdPlane_cuts[nPtBins][nEtaBins]->Fill(mThetaProdPlane);

          L0bar_cosThetaProdPlane_cuts[pT_bin][eta_bin]->Fill(cos(mThetaProdPlane));
          L0bar_cosThetaProdPlane_cuts[nPtBins][eta_bin]->Fill(cos(mThetaProdPlane));
          L0bar_cosThetaProdPlane_cuts[pT_bin][nEtaBins]->Fill(cos(mThetaProdPlane));
          L0bar_cosThetaProdPlane_cuts[nPtBins][nEtaBins]->Fill(cos(mThetaProdPlane));
          
          
          L0bar_pz_cuts[pT_bin][eta_bin]->Fill(L_fourmom.Pz());
          L0bar_pz_cuts[nPtBins][eta_bin]->Fill(L_fourmom.Pz());
          L0bar_pz_cuts[pT_bin][nEtaBins]->Fill(L_fourmom.Pz());
          L0bar_pz_cuts[nPtBins][nEtaBins]->Fill(L_fourmom.Pz());
          
          L0bar_xF_cuts[pT_bin][eta_bin]->Fill(L_xF);
          L0bar_xF_cuts[nPtBins][eta_bin]->Fill(L_xF);
          L0bar_xF_cuts[pT_bin][nEtaBins]->Fill(L_xF);
          L0bar_xF_cuts[nPtBins][nEtaBins]->Fill(L_xF);
       
          L0bar_pT_vs_L0bar_pz_cuts[eta_bin]->Fill(L_fourmom.Pt(), L_fourmom.Pz());
          L0bar_pT_vs_L0bar_pz_cuts[nEtaBins]->Fill(L_fourmom.Pt(), L_fourmom.Pz());
          
          L0bar_y_vs_p_eta_cuts[pT_bin]->Fill(L_fourmom.Rapidity(), p_fourmom.Eta());
          L0bar_y_vs_pi_eta_cuts[pT_bin]->Fill(L_fourmom.Rapidity(), pi_fourmom.Eta());
          
          L0bar_y_vs_p_eta_cuts[nPtBins]->Fill(L_fourmom.Rapidity(), p_fourmom.Eta());
          L0bar_y_vs_pi_eta_cuts[nPtBins]->Fill(L_fourmom.Rapidity(), pi_fourmom.Eta());
          
          L0bar_p_eta_vs_pi_eta_cuts[pT_bin]->Fill(p_fourmom.Eta(), pi_fourmom.Eta());
          L0bar_p_eta_vs_pi_eta_cuts[nPtBins]->Fill(p_fourmom.Eta(), pi_fourmom.Eta());

        }
        
        
        
      } //end PID if for Lambdas
    
    }//end loop over particles in event
    //______________________
    
     
    
    //fill L-L correlation histograms before cuts
    
    //L0-L0bar
    if(p_star_vector.size() > 0 && pbar_star_vector.size() > 0)
    {
      for(unsigned int iLambda = 0; iLambda < p_star_vector.size(); iLambda++)
      {
        for(unsigned int iLambdaBar = 0; iLambdaBar < pbar_star_vector.size(); iLambdaBar++)
        {
          double theta_star = p_star_vector.at(iLambda).Angle(pbar_star_vector.at(iLambdaBar).Vect());
          
          //cout<<L_pT_bin_vector.at(iLambda)<<" "<<Lbar_pT_bin_vector.at(iLambdaBar)<<endl;
          
          L0_L0bar_cosThetaProdPlane->Fill(TMath::Cos(theta_star));          
          L0_L0bar_cosThetaProdPlane_pT_hist[L_pT_bin_vector.at(iLambda)][Lbar_pT_bin_vector.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
          L0_L0bar_cosThetaProdPlane_eta_hist[L_eta_bin_vector.at(iLambda)][Lbar_eta_bin_vector.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
          
          
          //____________________

          double delta_eta = fabs( L_vector.at(iLambda).Eta() - Lbar_vector.at(iLambdaBar).Eta() );
          double delta_phi = fabs( L_vector.at(iLambda).Phi() - Lbar_vector.at(iLambdaBar).Phi() );

          L0_L0bar_delta_eta_vs_delta_phi_hist->Fill(delta_eta, delta_phi);
          //_____________________

          L0_L0bar_y1_vs_y2_hist->Fill( L_vector.at(iLambda).Rapidity(), Lbar_vector.at(iLambdaBar).Rapidity() );

          L0_L0bar_pT1_vs_pT2_hist->Fill( L_vector.at(iLambda).Pt(), Lbar_vector.at(iLambdaBar).Pt() );

          L0_L0bar_phi1_vs_phi2_hist->Fill( L_vector.at(iLambda).Phi(), Lbar_vector.at(iLambdaBar).Phi() );
       
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
          //____________________

          double delta_eta = fabs( L_vector.at(iLambda1).Eta() - L_vector.at(iLambda2).Eta() );
          double delta_phi = fabs( L_vector.at(iLambda1).Phi() - L_vector.at(iLambda2).Phi() );

          L0_L0_delta_eta_vs_delta_phi_hist->Fill(delta_eta, delta_phi);
          L0_L0_delta_eta_vs_delta_phi_zoom_hist->Fill(delta_eta, delta_phi);
          //_____________________

          L0_L0_y1_vs_y2_hist->Fill( L_vector.at(iLambda1).Rapidity(), L_vector.at(iLambda2).Rapidity() );
          L0_p_L0_p_y1_vs_y2_hist->Fill( p_vector.at(iLambda1).Rapidity(), p_vector.at(iLambda2).Rapidity() );
          L0_pi_L0_pi_y1_vs_y2_hist->Fill( pi_vector.at(iLambda1).Rapidity(), pi_vector.at(iLambda2).Rapidity() );

          L0_L0_pT1_vs_pT2_hist->Fill( L_vector.at(iLambda1).Pt(), L_vector.at(iLambda2).Pt() );

          L0_L0_phi1_vs_phi2_hist->Fill( L_vector.at(iLambda1).Phi(), L_vector.at(iLambda2).Phi() );
          L0_p_L0_p_phi1_vs_phi2_hist->Fill( p_vector.at(iLambda1).Phi(), p_vector.at(iLambda2).Phi() );
          L0_pi_L0_pi_phi1_vs_phi2_hist->Fill( pi_vector.at(iLambda1).Phi(), pi_vector.at(iLambda2).Phi() );
             
        }   
      
      }
    
    }
    
    //L0bar-L0bar
    if(pbar_star_vector.size() > 1)
    {
      for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pbar_star_vector.size(); iLambdaBar1++)
      {
        for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < pbar_star_vector.size(); iLambdaBar2++)
        {
          double theta_star = pbar_star_vector.at(iLambdaBar1).Angle(pbar_star_vector.at(iLambdaBar2).Vect());
          
          L0bar_L0bar_cosThetaProdPlane->Fill(TMath::Cos(theta_star));
          L0bar_L0bar_cosThetaProdPlane_pT_hist[Lbar_pT_bin_vector.at(iLambdaBar1)][Lbar_pT_bin_vector.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
          L0bar_L0bar_cosThetaProdPlane_eta_hist[Lbar_eta_bin_vector.at(iLambdaBar1)][Lbar_eta_bin_vector.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
          //____________________

          double delta_eta = fabs( Lbar_vector.at(iLambdaBar1).Eta() - Lbar_vector.at(iLambdaBar2).Eta() );
          double delta_phi = fabs( Lbar_vector.at(iLambdaBar1).Phi() - Lbar_vector.at(iLambdaBar2).Phi() );

          L0bar_L0bar_delta_eta_vs_delta_phi_hist->Fill(delta_eta, delta_phi);
          L0bar_L0bar_delta_eta_vs_delta_phi_zoom_hist->Fill(delta_eta, delta_phi);
          //_____________________

          L0bar_L0bar_y1_vs_y2_hist->Fill( Lbar_vector.at(iLambdaBar1).Rapidity(), Lbar_vector.at(iLambdaBar2).Rapidity() );
          L0bar_p_L0bar_p_y1_vs_y2_hist->Fill( pbar_vector.at(iLambdaBar1).Rapidity(), pbar_vector.at(iLambdaBar2).Rapidity() );
          L0bar_pi_L0bar_pi_y1_vs_y2_hist->Fill( pibar_vector.at(iLambdaBar1).Rapidity(), pibar_vector.at(iLambdaBar2).Rapidity() );

          L0bar_L0bar_pT1_vs_pT2_hist->Fill( Lbar_vector.at(iLambdaBar1).Pt(), Lbar_vector.at(iLambdaBar2).Pt() );

          L0bar_L0bar_phi1_vs_phi2_hist->Fill( Lbar_vector.at(iLambdaBar1).Phi(), Lbar_vector.at(iLambdaBar2).Phi() );
          L0bar_p_L0bar_p_phi1_vs_phi2_hist->Fill( pbar_vector.at(iLambdaBar1).Phi(), pbar_vector.at(iLambdaBar2).Phi() );
          L0bar_pi_L0bar_pi_phi1_vs_phi2_hist->Fill( pibar_vector.at(iLambdaBar1).Phi(), pibar_vector.at(iLambdaBar2).Phi() );
                
        }   
      
      } 
    
    }
    //_____________________________________________________________________________________________________________________________________________________________________
    
    //mixed event before cuts
    if( p_star_vector.size() == 1 && pbar_star_vector.size() == 0 && p_star_vector_ME.size() < 1e4)
    {
      p_star_vector_ME.push_back(p_star_vector.at(0));

      L_pT_bin_vector_ME.push_back(L_pT_bin_vector.at(0));
      L_eta_bin_vector_ME.push_back(L_eta_bin_vector.at(0));

    }

    if( p_star_vector.size() == 0 && pbar_star_vector.size() == 1 && pbar_star_vector_ME.size() < 1e4)
    {
      pbar_star_vector_ME.push_back(pbar_star_vector.at(0));

      Lbar_pT_bin_vector_ME.push_back(Lbar_pT_bin_vector.at(0));
      Lbar_eta_bin_vector_ME.push_back(Lbar_eta_bin_vector.at(0));

    }
    //_________________________________________________________________________________________________________
    
     
    
    //fill L-L correlation histograms afrer cuts
    
    //L0-L0bar
    if(p_star_cuts_vector.size() > 0 && pbar_star_cuts_vector.size() > 0)
    {
      for(unsigned int iLambda = 0; iLambda < p_star_cuts_vector.size(); iLambda++)
      {
        for(unsigned int iLambdaBar = 0; iLambdaBar < pbar_star_cuts_vector.size(); iLambdaBar++)
        {
          double theta_star = p_star_cuts_vector.at(iLambda).Angle(pbar_star_cuts_vector.at(iLambdaBar).Vect());
          
          L0_L0bar_cosThetaProdPlane_cuts->Fill(TMath::Cos(theta_star));
          L0_L0bar_cosThetaProdPlane_pT_cuts_hist[L_pT_bin_cuts_vector.at(iLambda)][Lbar_pT_bin_cuts_vector.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
          L0_L0bar_cosThetaProdPlane_eta_cuts_hist[L_eta_bin_cuts_vector.at(iLambda)][Lbar_eta_bin_cuts_vector.at(iLambdaBar)]->Fill(TMath::Cos(theta_star));
          //____________________

          double delta_eta = fabs( L_cuts_vector.at(iLambda).Eta() - Lbar_cuts_vector.at(iLambdaBar).Eta() );
          double delta_phi = fabs( L_cuts_vector.at(iLambda).Phi() - Lbar_cuts_vector.at(iLambdaBar).Phi() );

          L0_L0bar_delta_eta_vs_delta_phi_cuts_hist->Fill(delta_eta, delta_phi);
          //_____________________

          L0_L0bar_y1_vs_y2_cuts_hist->Fill( L_cuts_vector.at(iLambda).Rapidity(), Lbar_cuts_vector.at(iLambdaBar).Rapidity() );

          L0_L0bar_pT1_vs_pT2_cuts_hist->Fill( L_cuts_vector.at(iLambda).Pt(), Lbar_cuts_vector.at(iLambdaBar).Pt() );

          L0_L0bar_phi1_vs_phi2_cuts_hist->Fill( L_cuts_vector.at(iLambda).Phi(), Lbar_cuts_vector.at(iLambdaBar).Phi() );
         
        
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
          //____________________

          double delta_eta = fabs( L_cuts_vector.at(iLambda1).Eta() - L_cuts_vector.at(iLambda2).Eta() );
          double delta_phi = fabs( L_cuts_vector.at(iLambda1).Phi() - L_cuts_vector.at(iLambda2).Phi() );

          L0_L0_delta_eta_vs_delta_phi_cuts_hist->Fill(delta_eta, delta_phi);
          L0_L0_delta_eta_vs_delta_phi_zoom_cuts_hist->Fill(delta_eta, delta_phi);
          //_____________________

          L0_L0_y1_vs_y2_cuts_hist->Fill( L_cuts_vector.at(iLambda1).Rapidity(), L_cuts_vector.at(iLambda2).Rapidity() );
          L0_p_L0_p_y1_vs_y2_cuts_hist->Fill( p_cuts_vector.at(iLambda1).Rapidity(), p_cuts_vector.at(iLambda2).Rapidity() );
          L0_pi_L0_pi_y1_vs_y2_cuts_hist->Fill( pi_cuts_vector.at(iLambda1).Rapidity(), pi_cuts_vector.at(iLambda2).Rapidity() );

          L0_L0_pT1_vs_pT2_cuts_hist->Fill( L_cuts_vector.at(iLambda1).Pt(), L_cuts_vector.at(iLambda2).Pt() );

          L0_L0_phi1_vs_phi2_cuts_hist->Fill( L_cuts_vector.at(iLambda1).Phi(), L_cuts_vector.at(iLambda2).Phi() );
          L0_p_L0_p_phi1_vs_phi2_cuts_hist->Fill( p_cuts_vector.at(iLambda1).Phi(), p_cuts_vector.at(iLambda2).Phi() );
          L0_pi_L0_pi_phi1_vs_phi2_cuts_hist->Fill( pi_cuts_vector.at(iLambda1).Phi(), pi_cuts_vector.at(iLambda2).Phi() );
              
        }   
      
      }
    
    }
    
    //L0bar-L0bar
    if(pbar_star_cuts_vector.size() > 0)
    {
      for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pbar_star_cuts_vector.size(); iLambdaBar1++)
      {
        for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < pbar_star_cuts_vector.size(); iLambdaBar2++)
        {
          double theta_star = pbar_star_cuts_vector.at(iLambdaBar1).Angle(pbar_star_cuts_vector.at(iLambdaBar2).Vect());
          
          L0bar_L0bar_cosThetaProdPlane_cuts->Fill(TMath::Cos(theta_star));
          L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[Lbar_pT_bin_cuts_vector.at(iLambdaBar1)][Lbar_pT_bin_cuts_vector.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
          L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[Lbar_eta_bin_cuts_vector.at(iLambdaBar1)][Lbar_eta_bin_cuts_vector.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
          //____________________

          double delta_eta = fabs( Lbar_cuts_vector.at(iLambdaBar1).Eta() - Lbar_cuts_vector.at(iLambdaBar2).Eta() );
          double delta_phi = fabs( Lbar_cuts_vector.at(iLambdaBar1).Phi() - Lbar_cuts_vector.at(iLambdaBar2).Phi() );

          L0bar_L0bar_delta_eta_vs_delta_phi_cuts_hist->Fill(delta_eta, delta_phi);
          L0bar_L0bar_delta_eta_vs_delta_phi_zoom_cuts_hist->Fill(delta_eta, delta_phi);
          //_____________________

          L0bar_L0bar_y1_vs_y2_cuts_hist->Fill( Lbar_cuts_vector.at(iLambdaBar1).Rapidity(), Lbar_cuts_vector.at(iLambdaBar2).Rapidity() );
          L0bar_p_L0bar_p_y1_vs_y2_cuts_hist->Fill( pbar_cuts_vector.at(iLambdaBar1).Rapidity(), pbar_cuts_vector.at(iLambdaBar2).Rapidity() );
          L0bar_pi_L0bar_pi_y1_vs_y2_cuts_hist->Fill( pibar_cuts_vector.at(iLambdaBar1).Rapidity(), pibar_cuts_vector.at(iLambdaBar2).Rapidity() );

          L0bar_L0bar_pT1_vs_pT2_cuts_hist->Fill( Lbar_cuts_vector.at(iLambdaBar1).Pt(), Lbar_cuts_vector.at(iLambdaBar2).Pt() );

          L0bar_L0bar_phi1_vs_phi2_cuts_hist->Fill( Lbar_cuts_vector.at(iLambdaBar1).Phi(), Lbar_cuts_vector.at(iLambdaBar2).Phi() );
          L0bar_p_L0bar_p_phi1_vs_phi2_cuts_hist->Fill( pbar_cuts_vector.at(iLambdaBar1).Phi(), pbar_cuts_vector.at(iLambdaBar2).Phi() );
          L0bar_pi_L0bar_pi_phi1_vs_phi2_cuts_hist->Fill( pibar_cuts_vector.at(iLambdaBar1).Phi(), pibar_cuts_vector.at(iLambdaBar2).Phi() );
                
        }   
      
      } 
    
    }
    
    //mixed event after cuts
    if( p_star_cuts_vector.size() == 1 && pbar_star_cuts_vector.size() == 0 && p_star_cuts_vector_ME.size() < 1e4)
    {
      p_star_cuts_vector_ME.push_back(p_star_cuts_vector.at(0));

      L_pT_bin_cuts_vector_ME.push_back(L_pT_bin_cuts_vector.at(0));
      L_eta_bin_cuts_vector_ME.push_back(L_eta_bin_cuts_vector.at(0));

    }

    if( p_star_cuts_vector.size() == 0 && pbar_star_cuts_vector.size() == 1 && pbar_star_cuts_vector_ME.size() < 1e4)
    {
      pbar_star_cuts_vector_ME.push_back(pbar_star_cuts_vector.at(0));

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
    
    pbar_vector.clear();
    pbar_cuts_vector.clear();
    
    pibar_vector.clear();
    pibar_cuts_vector.clear();
    
    pbar_star_vector.clear();
    pbar_star_cuts_vector.clear();      

  }//end event loop
  //____________________________________________________________________________________________________________
  
  //cout<<p_star_vector_ME.size()<<endl;
  //cout<<pbar_star_vector_ME.size()<<endl;
  
  //mixed event before cuts
  for(unsigned int iLambda = 0; iLambda < p_star_vector_ME.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < pbar_star_vector_ME.size(); iLambdaBar++)
    {
      double theta_star = p_star_vector_ME.at(iLambda).Angle(pbar_star_vector_ME.at(iLambdaBar).Vect());
      
      //cout<<L_pT_bin_vector.at(iLambda)<<" "<<Lbar_pT_bin_vector.at(iLambdaBar)<<endl;
      
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
      //____________________
        
    }   
  
  }
  
  
  //L0bar-L0bar
  for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pbar_star_vector_ME.size(); iLambdaBar1++)
  {
    for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < pbar_star_vector_ME.size(); iLambdaBar2++)
    {
      double theta_star = pbar_star_vector_ME.at(iLambdaBar1).Angle(pbar_star_vector_ME.at(iLambdaBar2).Vect());
      
      L0bar_L0bar_cosThetaProdPlane_ME->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_ME_pT_hist[Lbar_pT_bin_vector_ME.at(iLambdaBar1)][Lbar_pT_bin_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_ME_eta_hist[Lbar_eta_bin_vector_ME.at(iLambdaBar1)][Lbar_eta_bin_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
      //____________________
            
    }   
  
  }
  //_____________________________________________________________________________________________________________________________________________________________________
  
  
  //mixed event after cuts
  for(unsigned int iLambda = 0; iLambda < p_star_cuts_vector_ME.size(); iLambda++)
  {
    for(unsigned int iLambdaBar = 0; iLambdaBar < pbar_star_cuts_vector_ME.size(); iLambdaBar++)
    {
      double theta_star = p_star_cuts_vector_ME.at(iLambda).Angle(pbar_star_cuts_vector_ME.at(iLambdaBar).Vect());
      
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
      //____________________
        
    }   
  
  }    
  
  //L0bar-L0bar
  for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pbar_star_cuts_vector_ME.size(); iLambdaBar1++)
  {
    for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < pbar_star_cuts_vector_ME.size(); iLambdaBar2++)
    {
      double theta_star = pbar_star_cuts_vector_ME.at(iLambdaBar1).Angle(pbar_star_cuts_vector_ME.at(iLambdaBar2).Vect());
      
      L0bar_L0bar_cosThetaProdPlane_ME_cuts->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[Lbar_pT_bin_cuts_vector_ME.at(iLambdaBar1)][Lbar_pT_bin_cuts_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[Lbar_eta_bin_cuts_vector_ME.at(iLambdaBar1)][Lbar_eta_bin_cuts_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
      //____________________
            
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
