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


struct K0s_MC {
  // general information
  // lab frame
  enum {
    nK0s_MAX_MC=10 //check that thi is enough
  };

  //number of L
  Int_t nK0s_MC;

  //L momentum for eta and pT bins and for Delta eta and Delta phi
  Float_t K0s_px_MC[nK0s_MAX_MC];
  Float_t K0s_py_MC[nK0s_MAX_MC];
  Float_t K0s_pz_MC[nK0s_MAX_MC];

  //proton momntum in L rest frame
  Float_t pi_pxStar_MC[nK0s_MAX_MC];
  Float_t pi_pyStar_MC[nK0s_MAX_MC];
  Float_t pi_pzStar_MC[nK0s_MAX_MC];

  //for cuts
  Float_t K0s_decayL_MC[nK0s_MAX_MC];
  
  Float_t pi1_pT_MC[nK0s_MAX_MC];
  Float_t pi1_eta_MC[nK0s_MAX_MC];
  
  Float_t pi2_pT_MC[nK0s_MAX_MC];
  Float_t pi2_eta_MC[nK0s_MAX_MC];

};


struct K0s_from_pairs {
  // general information
  // lab frame
  enum {
    nK0s_MAX=100 //check that thi is enough
  };

  //number of L
  Int_t nK0s;

  //L momentum for eta and pT bins and for Delta eta and Delta phi
  Float_t K0s_px[nK0s_MAX];
  Float_t K0s_py[nK0s_MAX];
  Float_t K0s_pz[nK0s_MAX];
  Float_t K0s_Minv[nK0s_MAX];

  //proton momntum in L rest frame
  Float_t pi_pxStar[nK0s_MAX];
  Float_t pi_pyStar[nK0s_MAX];
  Float_t pi_pzStar[nK0s_MAX];

  //for cuts
  Float_t K0s_decayL[nK0s_MAX];
  Float_t K0s_theta[nK0s_MAX]; //pointing angle
  
  Float_t pi1_pT[nK0s_MAX];
  Float_t pi1_eta[nK0s_MAX];
  
  Float_t pi2_pT[nK0s_MAX];
  Float_t pi2_eta[nK0s_MAX];
  
  //pipi charge combinations
  Int_t K0s_charge[nK0s_MAX]; //1 - pi+pi+, 0 - pi+pi-, -1 - pi-pi-
};


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
  
  //-----------------------------------------
  
  //vectors for creating US and LS pi pairs
  //pi origins to calculate SV of pairs
  
  vector<TLorentzVector> piPlus_vector_event;
  vector<TLorentzVector> piMinus_vector_event;
  
  vector<TVector3> piPlus_origin_vector_event;
  vector<TVector3> piMinus_origin_vector_event;
  
  //-----------------------------------------
  
  //trees to save L info
  
  TTree *K0s_MC_tree = new TTree("K0s_MC_tree", "K0s_MC_tree");
  K0s_MC Kaon_MC;
  
  //number of L
  K0s_MC_tree->Branch( "nK0s_MC", &Kaon_MC.nK0s_MC, "nK0s_MC/I" );

  //L momentum for eta and pT bins and for Delta eta and Delta phi
  K0s_MC_tree->Branch( "K0s_px_MC", &Kaon_MC.K0s_px_MC, "K0s_px_MC[nK0s_MC]/F" );
  K0s_MC_tree->Branch( "K0s_py_MC", &Kaon_MC.K0s_py_MC, "K0s_py_MC[nK0s_MC]/F" );
  K0s_MC_tree->Branch( "K0s_pz_MC", &Kaon_MC.K0s_pz_MC, "K0s_pz_MC[nK0s_MC]/F" );

  //proton momntum in L rest frame
  K0s_MC_tree->Branch( "pi_pxStar_MC", &Kaon_MC.pi_pxStar_MC, "pi_pxStar_MC[nK0s_MC]/F" );
  K0s_MC_tree->Branch( "pi_pyStar_MC", &Kaon_MC.pi_pyStar_MC, "pi_pyStar_MC[nK0s_MC]/F" );
  K0s_MC_tree->Branch( "pi_pzStar_MC", &Kaon_MC.pi_pzStar_MC, "pi_pzStar_MC[nK0s_MC]/F" );

  //for cuts
  K0s_MC_tree->Branch( "K0s_decayL_MC", &Kaon_MC.K0s_decayL_MC, "K0s_decayL_MC[nK0s_MC]/F" );
  
  K0s_MC_tree->Branch( "pi1_pT_MC", &Kaon_MC.pi1_pT_MC, "pi1_pT_MC[nK0s_MC]/F" );
  K0s_MC_tree->Branch( "pi1_eta_MC", &Kaon_MC.pi1_eta_MC, "pi1_eta_MC[nK0s_MC]/F" );
  
  K0s_MC_tree->Branch( "pi2_pT_MC", &Kaon_MC.pi2_pT_MC, "pi2_pT_MC[nK0s_MC]/F" );
  K0s_MC_tree->Branch( "pi2_eta_MC", &Kaon_MC.pi2_eta_MC, "pi2_eta_MC[nK0s_MC]/F" );
  
  
  //-------------------------------------------------------------------------------
  
  TTree *K0s_tree = new TTree("K0s_tree", "K0s_tree");
  K0s_from_pairs Kaon_from_pair;
  
  
  //number of L
  K0s_tree->Branch( "nK0s", &Kaon_from_pair.nK0s, "nK0s/I" );

  //L momentum for eta and pT bins and for Delta eta and Delta phi
  K0s_tree->Branch( "K0s_px", &Kaon_from_pair.K0s_px, "K0s_px[nK0s]/F" );
  K0s_tree->Branch( "K0s_py", &Kaon_from_pair.K0s_py, "K0s_py[nK0s]/F" );
  K0s_tree->Branch( "K0s_pz", &Kaon_from_pair.K0s_pz, "K0s_pz[nK0s]/F" );
  K0s_tree->Branch( "K0s_Minv", &Kaon_from_pair.K0s_Minv, "K0s_Minv[nK0s]/F" );

  //proton momntum in L rest frame
  K0s_tree->Branch( "pi_pxStar", &Kaon_from_pair.pi_pxStar, "pi_pxStar[nK0s]/F" );
  K0s_tree->Branch( "pi_pyStar", &Kaon_from_pair.pi_pyStar, "pi_pyStar[nK0s]/F" );
  K0s_tree->Branch( "pi_pzStar", &Kaon_from_pair.pi_pzStar, "pi_pzStar[nK0s]/F" );

  //for cuts
  K0s_tree->Branch( "K0s_decayL", &Kaon_from_pair.K0s_decayL, "K0s_decayL[nK0s]/F" );
  K0s_tree->Branch( "K0s_theta", &Kaon_from_pair.K0s_theta, "K0s_theta[nK0s]/F" );
  
  K0s_tree->Branch( "pi1_pT", &Kaon_from_pair.pi1_pT, "pi1_pT[nK0s]/F" );
  K0s_tree->Branch( "pi1_eta", &Kaon_from_pair.pi1_eta, "pi1_eta[nK0s]/F" );
  
  K0s_tree->Branch( "pi2_pT", &Kaon_from_pair.pi2_pT, "pi2_pT[nK0s]/F" );
  K0s_tree->Branch( "pi2_eta", &Kaon_from_pair.pi2_eta, "pi2_eta[nK0s]/F" );
  
  //L charge
  K0s_tree->Branch( "K0s_charge", &Kaon_from_pair.K0s_charge , "K0s_charge[nK0s]/I" ); //1 - L, -1 - Lba
 
  
  
  // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < nEvents; ++iEvent)
  {
    if (!pythia.next()) continue;
    
    int iK0s_MC = 0;
    
    int iK0s_from_pairs = 0;

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
    
      if( fabs(pythia.event[i].id()) == K0sPDGid && iK0s_MC < 10)
      {
        //cout<<"Lambda found!"<<endl;
        
        //find index of decay daughters
        const int daughter1_Id = pythia.event[i].daughter1();
        const int daughter2_Id = pythia.event[i].daughter2();
        
        
        //Lambda fourmomentum
        K0_fourmom.SetPxPyPzE(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
        
        if( fabs(K0_fourmom.Rapidity()) >= 1. ) continue;
        
        //save K0s momentum
        Kaon_MC.K0s_px_MC[iK0s_MC] = K0_fourmom.Px();
        Kaon_MC.K0s_py_MC[iK0s_MC] = K0_fourmom.Py();
        Kaon_MC.K0s_pz_MC[iK0s_MC] = K0_fourmom.Pz();
        
        //daughter kinematics
        pi1_fourmom.SetPxPyPzE(pythia.event[daughter1_Id].px(), pythia.event[daughter1_Id].py(), pythia.event[daughter1_Id].pz(), pythia.event[daughter1_Id].e());
        
        Kaon_MC.pi1_pT_MC[iK0s_MC] = pi1_fourmom.Pt();
        Kaon_MC.pi1_eta_MC[iK0s_MC] = pi1_fourmom.Eta();
        
        
        pi2_fourmom.SetPxPyPzE(pythia.event[daughter2_Id].px(), pythia.event[daughter2_Id].py(), pythia.event[daughter2_Id].pz(), pythia.event[daughter2_Id].e());
        
        Kaon_MC.pi2_pT_MC[iK0s_MC] = pi2_fourmom.Pt();
        Kaon_MC.pi2_eta_MC[iK0s_MC] = pi2_fourmom.Eta();
        
        
        //Lambda MC decay veretex
        K0s_decay_vertex.SetXYZ(pythia.event[i].xDec(), pythia.event[i].yDec(), pythia.event[i].zDec());
        
        //Lambda MC production veretex (PV)
        K0s_production_vertex.SetXYZ(pythia.event[i].xProd(), pythia.event[i].yProd(), pythia.event[i].zProd());
        
        TVector3 K0s_decayK0s_vect = K0s_decay_vertex - K0s_production_vertex;
        Kaon_MC.K0s_decayL_MC[iK0s_MC] = K0s_decayK0s_vect.Mag();        
        
        
        K0_fourmom_reverse.SetPxPyPzE(-pythia.event[i].px(), -pythia.event[i].py(), -pythia.event[i].pz(), pythia.event[i].e());              

        TLorentzVector pi1_fourmom_star = pi1_fourmom;
        pi1_fourmom_star.Boost(K0_fourmom_reverse.BoostVector());  
        
        Kaon_MC.pi_pxStar_MC[iK0s_MC] = pi1_fourmom_star.Px();
        Kaon_MC.pi_pyStar_MC[iK0s_MC] = pi1_fourmom_star.Py();
        Kaon_MC.pi_pzStar_MC[iK0s_MC] = pi1_fourmom_star.Pz();
        
        
        iK0s_MC++;

      }//end if K0s
    
    }//end particle loop
    
    Kaon_MC.nK0s_MC = iK0s_MC;
    
    if(iK0s_MC > 0) K0s_MC_tree->Fill(); //fill tree only for events with more than 1 K0s
    
    
    //_____________________________________________________________________________________________

    
    //pair pions
    
    for( unsigned int pi1_index = 0; pi1_index < piPlus_vector_event.size(); pi1_index++ )
    {
      //US
      for( unsigned int pi2_index = 0; pi2_index < piMinus_vector_event.size(); pi2_index++ )
      {
        TLorentzVector K0s_fourmom_US_event = piPlus_vector_event.at(pi1_index) + piMinus_vector_event.at(pi2_index);
        
        if( iK0s_from_pairs >= 100  ) continue;
        if( fabs( K0s_fourmom_US_event.Rapidity() ) > 1 ) continue;
        if( K0s_fourmom_US_event.M() < 0.488 || K0s_fourmom_US_event.M() > 0.51 ) continue;
        
        Kaon_from_pair.K0s_px[iK0s_from_pairs] = K0s_fourmom_US_event.Px();
        Kaon_from_pair.K0s_py[iK0s_from_pairs] = K0s_fourmom_US_event.Py();
        Kaon_from_pair.K0s_pz[iK0s_from_pairs] = K0s_fourmom_US_event.Pz();
        Kaon_from_pair.K0s_Minv[iK0s_from_pairs] = K0s_fourmom_US_event.M();
        
        
        //daughter kinematics
        Kaon_from_pair.pi1_pT[iK0s_from_pairs] = piPlus_vector_event.at(pi1_index).Pt();
        Kaon_from_pair.pi1_eta[iK0s_from_pairs] = piPlus_vector_event.at(pi1_index).Eta();
        
        Kaon_from_pair.pi2_pT[iK0s_from_pairs] = piMinus_vector_event.at(pi2_index).Pt();
        Kaon_from_pair.pi2_eta[iK0s_from_pairs] = piMinus_vector_event.at(pi2_index).Eta();
        
        
        
        TLorentzVector K0s_fourmom_reverse_US_event;
        K0s_fourmom_reverse_US_event.SetPxPyPzE(-K0s_fourmom_US_event.Px(), -K0s_fourmom_US_event.Py(), -K0s_fourmom_US_event.Pz(), K0s_fourmom_US_event.E());
        
        TLorentzVector pi_star_US = piPlus_vector_event.at(pi1_index);
        pi_star_US.Boost(K0s_fourmom_reverse_US_event.BoostVector()); //boost p 4-momntum to L rest frame
        
        Kaon_from_pair.pi_pxStar[iK0s_from_pairs] = pi_star_US.Px();
        Kaon_from_pair.pi_pyStar[iK0s_from_pairs] = pi_star_US.Py();
        Kaon_from_pair.pi_pzStar[iK0s_from_pairs] = pi_star_US.Pz();
        
        
        //decay length cuts in mm (default PYTHIA units)
        TVector3 K0s_pairVertex_US = ( piPlus_origin_vector_event.at(pi1_index) + piMinus_origin_vector_event.at(pi2_index) )*0.5; //secondary vertex of pi p pair
        Kaon_from_pair.K0s_decayL[iK0s_from_pairs] = K0s_pairVertex_US.Mag(); //decay length of pi p pair assuming PV at (0,0,0)
        
        
        
        float K0s_point_angle_US = K0s_pairVertex_US.Angle(K0s_fourmom_US_event.Vect());
        Kaon_from_pair.K0s_theta[iK0s_from_pairs] = K0s_point_angle_US;
        
        
        Kaon_from_pair.K0s_charge[iK0s_from_pairs] = 0; //US pi pair
        
        
        iK0s_from_pairs++;
     
      }
      
      //LS with pi+ only, with pi- is lower
      for( unsigned int pi2_index = pi1_index+1; pi2_index < piPlus_vector_event.size(); pi2_index++ )
      {
        TLorentzVector K0s_fourmom_LS_event = piPlus_vector_event.at(pi1_index) + piPlus_vector_event.at(pi2_index);
        
        if( iK0s_from_pairs >= 100  ) continue;
        if( fabs( K0s_fourmom_LS_event.Rapidity() ) > 1 ) continue;
        if( K0s_fourmom_LS_event.M() < 0.488 || K0s_fourmom_LS_event.M() > 0.51 ) continue;
        
        Kaon_from_pair.K0s_px[iK0s_from_pairs] = K0s_fourmom_LS_event.Px();
        Kaon_from_pair.K0s_py[iK0s_from_pairs] = K0s_fourmom_LS_event.Py();
        Kaon_from_pair.K0s_pz[iK0s_from_pairs] = K0s_fourmom_LS_event.Pz();
        Kaon_from_pair.K0s_Minv[iK0s_from_pairs] = K0s_fourmom_LS_event.M();
        
        
        //daughter kinematics
        Kaon_from_pair.pi1_pT[iK0s_from_pairs] = piPlus_vector_event.at(pi1_index).Pt();
        Kaon_from_pair.pi1_eta[iK0s_from_pairs] = piPlus_vector_event.at(pi1_index).Eta();
        
        Kaon_from_pair.pi2_pT[iK0s_from_pairs] = piPlus_vector_event.at(pi2_index).Pt();
        Kaon_from_pair.pi2_eta[iK0s_from_pairs] = piPlus_vector_event.at(pi2_index).Eta();
        
            
        
        TLorentzVector K0s_fourmom_reverse_LS_event;
        K0s_fourmom_reverse_LS_event.SetPxPyPzE(-K0s_fourmom_LS_event.Px(), -K0s_fourmom_LS_event.Py(), -K0s_fourmom_LS_event.Pz(), K0s_fourmom_LS_event.E());
        
        TLorentzVector pi_star_LS = piPlus_vector_event.at(pi1_index);
        pi_star_LS.Boost(K0s_fourmom_reverse_LS_event.BoostVector()); //boost p 4-momntum to L rest frame
        
        Kaon_from_pair.pi_pxStar[iK0s_from_pairs] = pi_star_LS.Px();
        Kaon_from_pair.pi_pyStar[iK0s_from_pairs] = pi_star_LS.Py();
        Kaon_from_pair.pi_pzStar[iK0s_from_pairs] = pi_star_LS.Pz();

        
        //decay length cuts in mm (default PYTHIA units)
        TVector3 K0s_pairVertex_LS = ( piPlus_origin_vector_event.at(pi1_index) + piPlus_origin_vector_event.at(pi2_index) )*0.5; //secondary vertex of pi p pair
        Kaon_from_pair.K0s_decayL[iK0s_from_pairs] = K0s_pairVertex_LS.Mag(); //decay length of pi p pair assuming PV at (0,0,0)
   
        float K0s_point_angle_LS = K0s_pairVertex_LS.Angle(K0s_fourmom_LS_event.Vect());
        Kaon_from_pair.K0s_theta[iK0s_from_pairs] = K0s_point_angle_LS;
        
        
        Kaon_from_pair.K0s_charge[iK0s_from_pairs] = 1; //LS pi+pi+ pair
        
        
        iK0s_from_pairs++;
      
      
      }
    
    }
    
    //LS, with pi-
    for( unsigned int pi1_index = 0; pi1_index < piMinus_vector_event.size(); pi1_index++ )
    {
      for( unsigned int pi2_index = pi1_index+1; pi2_index < piMinus_vector_event.size(); pi2_index++ )
      {
        TLorentzVector K0s_fourmom_LS_event = piMinus_vector_event.at(pi1_index) + piMinus_vector_event.at(pi2_index);
        
        if( iK0s_from_pairs >= 100  ) continue;
        if( fabs( K0s_fourmom_LS_event.Rapidity() ) > 1 ) continue;
        if( K0s_fourmom_LS_event.M() < 0.488 || K0s_fourmom_LS_event.M() > 0.51 ) continue;
        
        
        Kaon_from_pair.K0s_px[iK0s_from_pairs] = K0s_fourmom_LS_event.Px();
        Kaon_from_pair.K0s_py[iK0s_from_pairs] = K0s_fourmom_LS_event.Py();
        Kaon_from_pair.K0s_pz[iK0s_from_pairs] = K0s_fourmom_LS_event.Pz();
        Kaon_from_pair.K0s_Minv[iK0s_from_pairs] = K0s_fourmom_LS_event.M();
        
        
        //daughter kinematics
        Kaon_from_pair.pi1_pT[iK0s_from_pairs] = piMinus_vector_event.at(pi1_index).Pt();
        Kaon_from_pair.pi1_eta[iK0s_from_pairs] = piMinus_vector_event.at(pi1_index).Eta();
        
        Kaon_from_pair.pi2_pT[iK0s_from_pairs] = piMinus_vector_event.at(pi2_index).Pt();
        Kaon_from_pair.pi2_eta[iK0s_from_pairs] = piMinus_vector_event.at(pi2_index).Eta();
        
        
        TLorentzVector K0s_fourmom_reverse_LS_event;
        K0s_fourmom_reverse_LS_event.SetPxPyPzE(-K0s_fourmom_LS_event.Px(), -K0s_fourmom_LS_event.Py(), -K0s_fourmom_LS_event.Pz(), K0s_fourmom_LS_event.E());
        
        TLorentzVector pi_star_LS = piMinus_vector_event.at(pi1_index);
        pi_star_LS.Boost(K0s_fourmom_reverse_LS_event.BoostVector()); //boost p 4-momntum to L rest frame
        
        Kaon_from_pair.pi_pxStar[iK0s_from_pairs] = pi_star_LS.Px();
        Kaon_from_pair.pi_pyStar[iK0s_from_pairs] = pi_star_LS.Py();
        Kaon_from_pair.pi_pzStar[iK0s_from_pairs] = pi_star_LS.Pz();
        
        
        //decay length cuts in mm (default PYTHIA units)
        TVector3 K0s_pairVertex_LS = ( piMinus_origin_vector_event.at(pi1_index) + piMinus_origin_vector_event.at(pi2_index) )*0.5; //secondary vertex of pi p pair
        Kaon_from_pair.K0s_decayL[iK0s_from_pairs] = K0s_pairVertex_LS.Mag(); //decay length of pi p pair assuming PV at (0,0,0)
              
        
        float K0s_point_angle_LS = K0s_pairVertex_LS.Angle(K0s_fourmom_LS_event.Vect());
        Kaon_from_pair.K0s_theta[iK0s_from_pairs] = K0s_point_angle_LS;
        
        
        Kaon_from_pair.K0s_charge[iK0s_from_pairs] = -1; //LS pi-pi- pair
        
        
        iK0s_from_pairs++;
        
      
      }
      
    }
    
    
    Kaon_from_pair.nK0s = iK0s_from_pairs;
    
    if(iK0s_from_pairs > 0) K0s_tree->Fill(); //fill tree only for events with more than 1 K0s
    
    
    //clear pi vectors for next event
    
    piPlus_vector_event.clear();
    piPlus_origin_vector_event.clear();
    
    piMinus_vector_event.clear();
    piMinus_origin_vector_event.clear();
    
    //_____________________________________________________________________________________________________________________________________________________________________

  }//end event loop
  
  
  outFile->cd();
  outFile->Write();
  outFile->Close();

  // Statistics on event generation.
  pythia.stat();
  
  //outFile->Close();

  // Done.
  return 0;
}
