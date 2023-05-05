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
 
  //_____________________________

  //histograms after cuts	  
  TH1D *L0_L0bar_cosThetaProdPlane_cuts = new TH1D("L0_L0bar_cosThetaProdPlane_cuts", "L0_L0bar_cosThetaProdPlane_cuts", 10, -1, 1);
  
  TH1D *L0_L0bar_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1D *L0_L0bar_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1D *L0_L0_cosThetaProdPlane_cuts = new TH1D("L0_L0_cosThetaProdPlane_cuts", "L0_L0_cosThetaProdPlane_cuts", 10, -1, 1);
  
  TH1D *L0_L0_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1D *L0_L0_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];
  
  
  TH1D *L0bar_L0bar_cosThetaProdPlane_cuts = new TH1D("L0bar_L0bar_cosThetaProdPlane_cuts", "L0bar_L0bar_cosThetaProdPlane_cuts", 10, -1, 1);  
  
  TH1D *L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

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

  //_____________________________
    
  
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
  vector<float> L_y_vector_ME;
  vector<float> L_y_cuts_vector_ME;
  
  vector<TLorentzVector> p_star_vector_ME;
  vector<TLorentzVector> p_star_cuts_vector_ME;
  
  vector<int> L_pT_bin_vector_ME;
  vector<int> L_pT_bin_cuts_vector_ME;
  
  vector<int> L_eta_bin_vector_ME;
  vector<int> L_eta_bin_cuts_vector_ME;

  
  vector<float> Lbar_y_vector_ME;
  vector<float> Lbar_y_cuts_vector_ME;
  
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
        
        if(L_fourmom.Rapidity() >= 1) continue; //want only L within |y| < 1
        
        //Lambda MC decay veretex
        L_decay_vertex.SetXYZ(pythia.event[i].xDec(), pythia.event[i].yDec(), pythia.event[i].zDec());
        
        //Lambda MC production veretex (PV)
        L_production_vertex.SetXYZ(pythia.event[i].xProd(), pythia.event[i].yProd(), pythia.event[i].zProd());
        
        TVector3 L_decayL_vect = L_decay_vertex - L_production_vertex;
        float L_decayL = L_decayL_vect.Mag();
        
        L_decayL_hist->Fill(L_decayL);
        
        float L_xF = fabs(L_fourmom.Pz())/pythia.info.eCM()*2.;
                       

        p_fourmom.SetPxPyPzE(pythia.event[daughter1_Id].px(), pythia.event[daughter1_Id].py(), pythia.event[daughter1_Id].pz(), pythia.event[daughter1_Id].e());
        pi_fourmom.SetPxPyPzE(pythia.event[daughter2_Id].px(), pythia.event[daughter2_Id].py(), pythia.event[daughter2_Id].pz(), pythia.event[daughter2_Id].e());
        
        
        L_fourmom_reverse.SetPxPyPzE(-pythia.event[i].px(), -pythia.event[i].py(), -pythia.event[i].pz(), pythia.event[i].e());

        TLorentzVector p_fourmom_star = p_fourmom;
        p_fourmom_star.Boost(L_fourmom_reverse.BoostVector());       
    
        
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
          
          pbar_cuts_vector.push_back(p_fourmom);
          pibar_cuts_vector.push_back(pi_fourmom);
        
          pbar_star_cuts_vector.push_back(p_fourmom_star);        
        }             
        
      } //end PID if for Lambdas
    
    } //end loop over particles in event
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
        }   
      
      } 
    
    }
    //_____________________________________________________________________________________________________________________________________________________________________
    
    //mixed event before cuts
    if( p_star_vector.size() == 1 && pbar_star_vector.size() == 0 && p_star_vector_ME.size() < 1e4)
    {
      L_y_vector_ME.push_back(L_vector.at(0).Rapidity());
    
      p_star_vector_ME.push_back(p_star_vector.at(0));

      L_pT_bin_vector_ME.push_back(L_pT_bin_vector.at(0));
      L_eta_bin_vector_ME.push_back(L_eta_bin_vector.at(0));
    }

    if( p_star_vector.size() == 0 && pbar_star_vector.size() == 1 && pbar_star_vector_ME.size() < 1e4)
    {
      Lbar_y_vector_ME.push_back(Lbar_vector.at(0).Rapidity());
    
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
        }   
      
      } 
    
    }
    
    //mixed event after cuts
    if( p_star_cuts_vector.size() == 1 && pbar_star_cuts_vector.size() == 0 && p_star_cuts_vector_ME.size() < 1e4)
    {
      L_y_cuts_vector_ME.push_back(L_vector.at(0).Rapidity());
    
      p_star_cuts_vector_ME.push_back(p_star_cuts_vector.at(0));

      L_pT_bin_cuts_vector_ME.push_back(L_pT_bin_cuts_vector.at(0));
      L_eta_bin_cuts_vector_ME.push_back(L_eta_bin_cuts_vector.at(0));

    }

    if( p_star_cuts_vector.size() == 0 && pbar_star_cuts_vector.size() == 1 && pbar_star_cuts_vector_ME.size() < 1e4)
    {
      Lbar_y_cuts_vector_ME.push_back(Lbar_vector.at(0).Rapidity());
    
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
  for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < pbar_star_vector_ME.size(); iLambdaBar1++)
  {
    for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < pbar_star_vector_ME.size(); iLambdaBar2++)
    {
      double theta_star = pbar_star_vector_ME.at(iLambdaBar1).Angle(pbar_star_vector_ME.at(iLambdaBar2).Vect());
      
      L0bar_L0bar_cosThetaProdPlane_ME->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_ME_pT_hist[Lbar_pT_bin_vector_ME.at(iLambdaBar1)][Lbar_pT_bin_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));
      L0bar_L0bar_cosThetaProdPlane_ME_eta_hist[Lbar_eta_bin_vector_ME.at(iLambdaBar1)][Lbar_eta_bin_vector_ME.at(iLambdaBar2)]->Fill(TMath::Cos(theta_star));            
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
