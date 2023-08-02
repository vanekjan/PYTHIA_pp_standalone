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
  
  istringstream mySeedStream(argv[4]);
  int mySeed;
  if (!(mySeedStream >> mySeed))
  {
    cout<<"Invalid fourth argument!"<<endl;
    return 0;
  }
  
  mySeed += 1; //need to sart at 1, argument starts at 0
  
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
  //pythia.readString("Random:seed = 0"); //should set random seed at start, after calling init() - based on time - sets same seed for several jobs
  pythia.readString(Form("Random:seed = %i", mySeed)); //set own seed - now set based on process ID of submitted job
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

  vector<TLorentzVector> p_vector_event;
  vector<TVector3> p_origin_vector_event;
  
  vector<TLorentzVector> pi_vector_event;
  vector<TVector3> pi_origin_vector_event;
  
  
  vector<TLorentzVector> pBar_vector_event;
  vector<TVector3> pBar_origin_vector_event;
  
  vector<TLorentzVector> piBar_vector_event;
  vector<TVector3> piBar_origin_vector_event;

  //variables for event and particle loop
  TLorentzVector L_fourmom;
  TLorentzVector L_fourmom_reverse; //for proton boost
  
  TVector3 L_production_vertex;
  TVector3 L_decay_vertex;

  TLorentzVector p_fourmom;
  TLorentzVector pi_fourmom;
  
  //_____________________________________________________________
  
  //trees to save L info
  
  TTree *L_MC_tree = new TTree("L_MC_tree", "L_MC_tree");
  L_MC Lambda_MC;
  
  //number of L
  L_MC_tree->Branch( "nL_MC", &Lambda_MC.nL_MC, "nL_MC/I" );

  //L momentum for eta and pT bins and for Delta eta and Delta phi
  L_MC_tree->Branch( "L_px_MC", &Lambda_MC.L_px_MC, "L_px_MC[nL_MC]/F" );
  L_MC_tree->Branch( "L_py_MC", &Lambda_MC.L_py_MC, "L_py_MC[nL_MC]/F" );
  L_MC_tree->Branch( "L_pz_MC", &Lambda_MC.L_pz_MC, "L_pz_MC[nL_MC]/F" );

  //proton momntum in L rest frame
  L_MC_tree->Branch( "p_pxStar_MC", &Lambda_MC.p_pxStar_MC, "p_pxStar_MC[nL_MC]/F" );
  L_MC_tree->Branch( "p_pyStar_MC", &Lambda_MC.p_pyStar_MC, "p_pyStar_MC[nL_MC]/F" );
  L_MC_tree->Branch( "p_pzStar_MC", &Lambda_MC.p_pzStar_MC, "p_pzStar_MC[nL_MC]/F" );

  //for cuts
  L_MC_tree->Branch( "L_decayL_MC", &Lambda_MC.L_decayL_MC, "L_decayL_MC[nL_MC]/F" );
  
  L_MC_tree->Branch( "p_pT_MC", &Lambda_MC.p_pT_MC, "p_pT_MC[nL_MC]/F" );
  L_MC_tree->Branch( "p_eta_MC", &Lambda_MC.p_eta_MC, "p_eta_MC[nL_MC]/F" );
  L_MC_tree->Branch( "p_phi_MC", &Lambda_MC.p_phi_MC, "p_phi_MC[nL_MC]/F" );
  
  L_MC_tree->Branch( "pi_pT_MC", &Lambda_MC.pi_pT_MC, "pi_pT_MC[nL_MC]/F" );
  L_MC_tree->Branch( "pi_eta_MC", &Lambda_MC.pi_eta_MC, "pi_eta_MC[nL_MC]/F" );
  L_MC_tree->Branch( "pi_phi_MC", &Lambda_MC.pi_phi_MC, "pi_phi_MC[nL_MC]/F" );
  
  //L charge
  L_MC_tree->Branch( "L_charge_MC", &Lambda_MC.L_charge_MC , "L_charge_MC[nL_MC]/I" ); //1 - L, -1 - Lbar
  
  //-------------------------------------------------------------------------------
  
  TTree *L_tree = new TTree("L_tree", "L_tree");
  L_from_pairs Lambda_from_pair;
  
  
  //number of L
  L_tree->Branch( "nL", &Lambda_from_pair.nL, "nL/I" );

  //L momentum for eta and pT bins and for Delta eta and Delta phi
  L_tree->Branch( "L_px", &Lambda_from_pair.L_px, "L_px[nL]/F" );
  L_tree->Branch( "L_py", &Lambda_from_pair.L_py, "L_py[nL]/F" );
  L_tree->Branch( "L_pz", &Lambda_from_pair.L_pz, "L_pz[nL]/F" );
  L_tree->Branch( "L_Minv", &Lambda_from_pair.L_Minv, "L_Minv[nL]/F" );

  //proton momntum in L rest frame
  L_tree->Branch( "p_pxStar", &Lambda_from_pair.p_pxStar, "p_pxStar[nL]/F" );
  L_tree->Branch( "p_pyStar", &Lambda_from_pair.p_pyStar, "p_pyStar[nL]/F" );
  L_tree->Branch( "p_pzStar", &Lambda_from_pair.p_pzStar, "p_pzStar[nL]/F" );

  //for cuts
  L_tree->Branch( "L_decayL", &Lambda_from_pair.L_decayL, "L_decayL[nL]/F" );
  L_tree->Branch( "L_theta", &Lambda_from_pair.L_theta, "L_theta[nL]/F" );
  
  L_tree->Branch( "p_pT", &Lambda_from_pair.p_pT, "p_pT[nL]/F" );
  L_tree->Branch( "p_eta", &Lambda_from_pair.p_eta, "p_eta[nL]/F" );
  L_tree->Branch( "p_phi", &Lambda_from_pair.p_phi, "p_phi[nL]/F" );
  
  L_tree->Branch( "pi_pT", &Lambda_from_pair.pi_pT, "pi_pT[nL]/F" );
  L_tree->Branch( "pi_eta", &Lambda_from_pair.pi_eta, "pi_eta[nL]/F" );
  L_tree->Branch( "pi_phi", &Lambda_from_pair.pi_phi, "pi_phi[nL]/F" );
  
  //L charge
  L_tree->Branch( "L_charge", &Lambda_from_pair.L_charge , "L_charge[nL]/I" ); //1 - L, -1 - Lbar
  L_tree->Branch( "L_US_LS_flag", &Lambda_from_pair.L_US_LS_flag , "L_US_LS_flag[nL]/I" ); // 1 - US, 0 - LS
   

  // Begin event loop. Generate event; skip if generation aborted.
  //for (int iEvent = 0; iEvent < 2e8; ++iEvent)
  for (int iEvent = 0; iEvent < nEvents; ++iEvent)
  //for (int iEvent = 0; iEvent < 1e2; ++iEvent)
  {
    if (!pythia.next()) continue;
    
    
    int iLambda = 0; //index of MC Lambdas
    
    int iLambda_from_pairs = 0; //index of p-pi pairs (total)

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
      if( fabs(pythia.event[i].id()) == 3122 && iLambda < 10 )
      {
        //cout<<"Lambda found!"<<endl;
        
        //find index of decay daughters
        const int daughter1_Id = pythia.event[i].daughter1();
        const int daughter2_Id = pythia.event[i].daughter2();
        
        
        //Lambda fourmomentum
        L_fourmom.SetPxPyPzE(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
        
        if(L_fourmom.Rapidity() >= 1) continue; //want only L within |y| < 1
        
        //save L momentum
        Lambda_MC.L_px_MC[iLambda] = L_fourmom.Px();
        Lambda_MC.L_py_MC[iLambda] = L_fourmom.Py();
        Lambda_MC.L_pz_MC[iLambda] = L_fourmom.Pz();        
               
        //daughter kinematics
        p_fourmom.SetPxPyPzE(pythia.event[daughter1_Id].px(), pythia.event[daughter1_Id].py(), pythia.event[daughter1_Id].pz(), pythia.event[daughter1_Id].e());
        
        Lambda_MC.p_pT_MC[iLambda] = p_fourmom.Pt();
        Lambda_MC.p_eta_MC[iLambda] = p_fourmom.Eta();
        Lambda_MC.p_phi_MC[iLambda] = p_fourmom.Phi();
        
        pi_fourmom.SetPxPyPzE(pythia.event[daughter2_Id].px(), pythia.event[daughter2_Id].py(), pythia.event[daughter2_Id].pz(), pythia.event[daughter2_Id].e());
        
        Lambda_MC.pi_pT_MC[iLambda] = pi_fourmom.Pt();
        Lambda_MC.pi_eta_MC[iLambda] = pi_fourmom.Eta();
        Lambda_MC.pi_phi_MC[iLambda] = pi_fourmom.Phi();
        
        
        L_fourmom_reverse.SetPxPyPzE(-pythia.event[i].px(), -pythia.event[i].py(), -pythia.event[i].pz(), pythia.event[i].e());

        TLorentzVector p_fourmom_star = p_fourmom;
        p_fourmom_star.Boost(L_fourmom_reverse.BoostVector());       
    
        
        Lambda_MC.p_pxStar_MC[iLambda] = p_fourmom_star.Px();
        Lambda_MC.p_pyStar_MC[iLambda] = p_fourmom_star.Py();
        Lambda_MC.p_pzStar_MC[iLambda] = p_fourmom_star.Pz();
        
        
        //Lambda MC decay veretex
        L_decay_vertex.SetXYZ(pythia.event[i].xDec(), pythia.event[i].yDec(), pythia.event[i].zDec());
        
        //Lambda MC production veretex (PV)
        //L_production_vertex.SetXYZ(pythia.event[i].xProd(), pythia.event[i].yProd(), pythia.event[i].zProd()); //true MC production vertex - can be different from PV due to secondary decays
        L_production_vertex.SetXYZ(0, 0, 0); //fix production vertex to MC PV, i.e. (0,0,0)
        
        TVector3 L_decayL_vect = L_decay_vertex - L_production_vertex;
        Lambda_MC.L_decayL_MC[iLambda] = L_decayL_vect.Mag();


        //Lambda
        if( pythia.event[i].id() > 0 )
        {
          Lambda_MC.L_charge_MC[iLambda] = 1;
        }

        //Lambda-bar
        if( pythia.event[i].id() < 0 )
        {
          Lambda_MC.L_charge_MC[iLambda] = -1;        
        }
               
        iLambda++;
        
      } //end PID if for Lambdas
         

    } //end loop over particles in event
    
    
    Lambda_MC.nL_MC = iLambda; //nLambas = index of last Lambda + 1 (index starts at 0)
      
    if(iLambda > 0) L_MC_tree->Fill(); //fill tree only if there is at least 1 L (Lbar) in the event
    
    
    //_____________________________________________________________________________________________________________________________________________________________________
    
    //pair pi and p
    //L
    for(unsigned int p_index = 0; p_index < p_vector_event.size(); p_index++)
    {
      //US
      for(unsigned int piBar_index = 0; piBar_index < piBar_vector_event.size(); piBar_index++)
      {
        TLorentzVector L_fourmom_US_event = p_vector_event.at(p_index) + piBar_vector_event.at(piBar_index);
        
        if( iLambda_from_pairs >= 100 ) //overfill protection
        if( fabs( L_fourmom_US_event.Rapidity() ) > 1 ) continue;
        if( L_fourmom_US_event.M() < 1.115 || L_fourmom_US_event.M() > 1.12 ) continue;
        
        Lambda_from_pair.L_px[iLambda_from_pairs] = L_fourmom_US_event.Px();
        Lambda_from_pair.L_py[iLambda_from_pairs] = L_fourmom_US_event.Py();
        Lambda_from_pair.L_pz[iLambda_from_pairs] = L_fourmom_US_event.Pz();
        Lambda_from_pair.L_Minv[iLambda_from_pairs] = L_fourmom_US_event.M();     
        
        
        //daughter kinematics
        Lambda_from_pair.p_pT[iLambda_from_pairs] = p_vector_event.at(p_index).Pt();
        Lambda_from_pair.p_eta[iLambda_from_pairs] = p_vector_event.at(p_index).Eta();
        Lambda_from_pair.p_phi[iLambda_from_pairs] = p_vector_event.at(p_index).Phi();
        
        Lambda_from_pair.pi_pT[iLambda_from_pairs] = piBar_vector_event.at(piBar_index).Pt();
        Lambda_from_pair.pi_eta[iLambda_from_pairs] = piBar_vector_event.at(piBar_index).Eta();   
        Lambda_from_pair.pi_phi[iLambda_from_pairs] = piBar_vector_event.at(piBar_index).Phi();   
        
        
        TLorentzVector L_fourmom_reverse_US_event;
        L_fourmom_reverse_US_event.SetPxPyPzE(-L_fourmom_US_event.Px(), -L_fourmom_US_event.Py(), -L_fourmom_US_event.Pz(), L_fourmom_US_event.E());      
 
        TLorentzVector p_star_US = p_vector_event.at(p_index);
        p_star_US.Boost(L_fourmom_reverse_US_event.BoostVector()); //boost p 4-momntum to L rest frame
        
        Lambda_from_pair.p_pxStar[iLambda_from_pairs] = p_star_US.Px();
        Lambda_from_pair.p_pyStar[iLambda_from_pairs] = p_star_US.Py();
        Lambda_from_pair.p_pzStar[iLambda_from_pairs] = p_star_US.Pz(); 
        
                
        //decay length cuts in mm (default PYTHIA units)
        TVector3 L_pairVertex_US = ( p_origin_vector_event.at(p_index) + piBar_origin_vector_event.at(piBar_index) )*0.5; //secondary vertex of pi p pair
        Lambda_from_pair.L_decayL[iLambda_from_pairs] = L_pairVertex_US.Mag(); //decay length of pi p pair assuming PV at (0,0,0)  
        
        //cos of pointing angle
        float L_point_angle_US = L_pairVertex_US.Angle(L_fourmom_US_event.Vect());
        Lambda_from_pair.L_theta[iLambda_from_pairs] = L_point_angle_US;
                
        Lambda_from_pair.L_charge[iLambda_from_pairs] = 1; //L
        Lambda_from_pair.L_US_LS_flag[iLambda_from_pairs] = 1; //US pair
        
        
        iLambda_from_pairs++;        
      }
      
      //LS
      for(unsigned int pi_index = 0; pi_index < pi_vector_event.size(); pi_index++)
      {
        TLorentzVector L_fourmom_LS_event = p_vector_event.at(p_index) + pi_vector_event.at(pi_index);
        
        if( iLambda_from_pairs >= 100 ) //overfill protection
        if( fabs( L_fourmom_LS_event.Rapidity() ) > 1 ) continue;
        if( L_fourmom_LS_event.M() < 1.115 || L_fourmom_LS_event.M() > 1.12 ) continue;
        
        Lambda_from_pair.L_px[iLambda_from_pairs] = L_fourmom_LS_event.Px();
        Lambda_from_pair.L_py[iLambda_from_pairs] = L_fourmom_LS_event.Py();
        Lambda_from_pair.L_pz[iLambda_from_pairs] = L_fourmom_LS_event.Pz();
        Lambda_from_pair.L_Minv[iLambda_from_pairs] = L_fourmom_LS_event.M(); 

        
        //daughter kinematics
        Lambda_from_pair.p_pT[iLambda_from_pairs] = p_vector_event.at(p_index).Pt();
        Lambda_from_pair.p_eta[iLambda_from_pairs] = p_vector_event.at(p_index).Eta();
        Lambda_from_pair.p_phi[iLambda_from_pairs] = p_vector_event.at(p_index).Phi();
        
        Lambda_from_pair.pi_pT[iLambda_from_pairs] = pi_vector_event.at(pi_index).Pt();
        Lambda_from_pair.pi_eta[iLambda_from_pairs] = pi_vector_event.at(pi_index).Eta();
        Lambda_from_pair.pi_phi[iLambda_from_pairs] = pi_vector_event.at(pi_index).Phi();
        
        
        TLorentzVector L_fourmom_reverse_LS_event;
        L_fourmom_reverse_LS_event.SetPxPyPzE(-L_fourmom_LS_event.Px(), -L_fourmom_LS_event.Py(), -L_fourmom_LS_event.Pz(), L_fourmom_LS_event.E());
        
        TLorentzVector p_star_LS = p_vector_event.at(p_index);
        p_star_LS.Boost(L_fourmom_reverse_LS_event.BoostVector()); //boost p 4-momntum to L rest frame
        
        Lambda_from_pair.p_pxStar[iLambda_from_pairs] = p_star_LS.Px();
        Lambda_from_pair.p_pyStar[iLambda_from_pairs] = p_star_LS.Py();
        Lambda_from_pair.p_pzStar[iLambda_from_pairs] = p_star_LS.Pz(); 

        
        //decay length cuts in mm (default PYTHIA units)
        TVector3 L_pairVertex_LS = ( p_origin_vector_event.at(p_index) + pi_origin_vector_event.at(pi_index) )*0.5; //secondary vertex of pi p pair
        Lambda_from_pair.L_decayL[iLambda_from_pairs] = L_pairVertex_LS.Mag();
        
        
        //cos of pointing angle
        float L_point_angle_LS = L_pairVertex_LS.Angle(L_fourmom_LS_event.Vect());
        Lambda_from_pair.L_theta[iLambda_from_pairs] = L_point_angle_LS;
        
        
        Lambda_from_pair.L_charge[iLambda_from_pairs] = 1; //L
        Lambda_from_pair.L_US_LS_flag[iLambda_from_pairs] = 0; //LS pair
        
        
        iLambda_from_pairs++;      

      }    
    }
    
    //Lbar
    for(unsigned int pBar_index = 0; pBar_index < pBar_vector_event.size(); pBar_index++)
    {
      //US
      for(unsigned int pi_index = 0; pi_index < pi_vector_event.size(); pi_index++)
      {
        TLorentzVector Lbar_fourmom_US_event = pBar_vector_event.at(pBar_index) + pi_vector_event.at(pi_index);
        
        if( iLambda_from_pairs >= 100 ) //overfill protection
        if( fabs( Lbar_fourmom_US_event.Rapidity() ) > 1 ) continue;
        if( Lbar_fourmom_US_event.M() < 1.115 || Lbar_fourmom_US_event.M() > 1.12 ) continue;
        
        Lambda_from_pair.L_px[iLambda_from_pairs] = Lbar_fourmom_US_event.Px();
        Lambda_from_pair.L_py[iLambda_from_pairs] = Lbar_fourmom_US_event.Py();
        Lambda_from_pair.L_pz[iLambda_from_pairs] = Lbar_fourmom_US_event.Pz();
        Lambda_from_pair.L_Minv[iLambda_from_pairs] = Lbar_fourmom_US_event.M();
        
        
        //daughter kinematics
        Lambda_from_pair.p_pT[iLambda_from_pairs] = pBar_vector_event.at(pBar_index).Pt();
        Lambda_from_pair.p_eta[iLambda_from_pairs] = pBar_vector_event.at(pBar_index).Eta();
        Lambda_from_pair.p_phi[iLambda_from_pairs] = pBar_vector_event.at(pBar_index).Phi();
        
        Lambda_from_pair.pi_pT[iLambda_from_pairs] = pi_vector_event.at(pi_index).Pt();
        Lambda_from_pair.pi_eta[iLambda_from_pairs] = pi_vector_event.at(pi_index).Eta();   
        Lambda_from_pair.pi_phi[iLambda_from_pairs] = pi_vector_event.at(pi_index).Phi();
        
        
        TLorentzVector Lbar_fourmom_reverse_US_event;
        Lbar_fourmom_reverse_US_event.SetPxPyPzE(-Lbar_fourmom_US_event.Px(), -Lbar_fourmom_US_event.Py(), -Lbar_fourmom_US_event.Pz(), Lbar_fourmom_US_event.E());
        
        TLorentzVector pBar_star_US = pBar_vector_event.at(pBar_index);
        pBar_star_US.Boost(Lbar_fourmom_reverse_US_event.BoostVector()); //boost p 4-momntum to Lbar rest frame

        Lambda_from_pair.p_pxStar[iLambda_from_pairs] = pBar_star_US.Px();
        Lambda_from_pair.p_pyStar[iLambda_from_pairs] = pBar_star_US.Py();
        Lambda_from_pair.p_pzStar[iLambda_from_pairs] = pBar_star_US.Pz();         

        
        //decay length cuts in mm (default PYTHIA units)
        TVector3 Lbar_pairVertex_US = ( pBar_origin_vector_event.at(pBar_index) + pi_origin_vector_event.at(pi_index) )*0.5; //secondary vertex of pi p pair
        Lambda_from_pair.L_decayL[iLambda_from_pairs] = Lbar_pairVertex_US.Mag(); //decay length of pi p pair assuming PV at (0,0,0)
        
        
        //cos of pointing angle
        float Lbar_point_angle_US = Lbar_pairVertex_US.Angle(Lbar_fourmom_US_event.Vect());
        Lambda_from_pair.L_theta[iLambda_from_pairs] = Lbar_point_angle_US;
                
        Lambda_from_pair.L_charge[iLambda_from_pairs] = -1; //Lbar
        Lambda_from_pair.L_US_LS_flag[iLambda_from_pairs] = 1; //US pair
        
        
        iLambda_from_pairs++;
    
      }
      
      //LS
      for(unsigned int piBar_index = 0; piBar_index < piBar_vector_event.size(); piBar_index++)
      {
        TLorentzVector Lbar_fourmom_LS_event = pBar_vector_event.at(pBar_index) + piBar_vector_event.at(piBar_index);
        
        if( iLambda_from_pairs >= 100 ) //overfill protection
        if( fabs( Lbar_fourmom_LS_event.Rapidity() ) > 1 ) continue;
        if( Lbar_fourmom_LS_event.M() < 1.115 || Lbar_fourmom_LS_event.M() > 1.12 ) continue;
        
        Lambda_from_pair.L_px[iLambda_from_pairs] = Lbar_fourmom_LS_event.Px();
        Lambda_from_pair.L_py[iLambda_from_pairs] = Lbar_fourmom_LS_event.Py();
        Lambda_from_pair.L_pz[iLambda_from_pairs] = Lbar_fourmom_LS_event.Pz();
        Lambda_from_pair.L_Minv[iLambda_from_pairs] = Lbar_fourmom_LS_event.M();
        
        
        //daughter kinematics
        Lambda_from_pair.p_pT[iLambda_from_pairs] = pBar_vector_event.at(pBar_index).Pt();
        Lambda_from_pair.p_eta[iLambda_from_pairs] = pBar_vector_event.at(pBar_index).Eta();
        Lambda_from_pair.p_phi[iLambda_from_pairs] = pBar_vector_event.at(pBar_index).Phi();
        
        Lambda_from_pair.pi_pT[iLambda_from_pairs] = piBar_vector_event.at(piBar_index).Pt();
        Lambda_from_pair.pi_eta[iLambda_from_pairs] = piBar_vector_event.at(piBar_index).Eta(); 
        Lambda_from_pair.pi_phi[iLambda_from_pairs] = piBar_vector_event.at(piBar_index).Phi();
  
        
        TLorentzVector Lbar_fourmom_reverse_LS_event;
        Lbar_fourmom_reverse_LS_event.SetPxPyPzE(-Lbar_fourmom_LS_event.Px(), -Lbar_fourmom_LS_event.Py(), -Lbar_fourmom_LS_event.Pz(), Lbar_fourmom_LS_event.E());
        
        TLorentzVector pBar_star_LS = pBar_vector_event.at(pBar_index);
        pBar_star_LS.Boost(Lbar_fourmom_reverse_LS_event.BoostVector()); //boost p 4-momntum to Lbar rest frame
        
     
        
        //decay length cuts in mm (default PYTHIA units)
        TVector3 Lbar_pairVertex_LS = ( pBar_origin_vector_event.at(pBar_index) + piBar_origin_vector_event.at(piBar_index) )*0.5; //secondary vertex of pi p pair
        Lambda_from_pair.L_decayL[iLambda_from_pairs] = Lbar_pairVertex_LS.Mag(); //decay length of pi p pair assuming PV at (0,0,0)
        
        
        //cos of pointing angle
        float Lbar_point_angle_LS = Lbar_pairVertex_LS.Angle(Lbar_fourmom_LS_event.Vect());
        Lambda_from_pair.L_theta[iLambda_from_pairs] = Lbar_point_angle_LS;
                
        Lambda_from_pair.L_charge[iLambda_from_pairs] = -1; //Lbar
        Lambda_from_pair.L_US_LS_flag[iLambda_from_pairs] = 0; //LS pair
        
        
        iLambda_from_pairs++;        
          
      }     
    }
    
    Lambda_from_pair.nL = iLambda_from_pairs;
    
    if( iLambda_from_pairs > 0 ) L_tree->Fill(); //fill tree if there is more than one US or LS p-pi pair in the event
    
    //_________________________________________________________________________________________________________
    
    //clear p and pi vectors
    p_vector_event.clear();
    p_origin_vector_event.clear();
    
    pi_vector_event.clear();
    pi_origin_vector_event.clear();
    
    
    pBar_vector_event.clear();
    pBar_origin_vector_event.clear();
    
    piBar_vector_event.clear();
    piBar_origin_vector_event.clear();
    

  }//end event loop
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
