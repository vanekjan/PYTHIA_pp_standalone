// main92.cc is a part of the PYTHIA event generator.
// Copyright (C) 2022 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: analysis; root;

// This is a simple test program.
// Modified by Rene Brun and Axel Naumann to put the Pythia::event
// into a TTree.

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

// ROOT, for saving Pythia events as trees in a file.
//#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLorentzVector.h"

using namespace Pythia8;
using namespace TMath;


Float_t ptHatMin = 0;
Float_t ptHatMax = 128;

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

  const int nPtBins = 6;
  float const pT_bins[nPtBins+1] = { 0,1,2,3,4,5,6 };

  // Create Pythia instance and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
  //pythia.readString("HardQCD:all = on");
  //pythia.readString("PhaseSpace:pTHatMin = 20.");
  if( mEnergy == 510) pythia.readString("Beams:eCM = 510"); //beam energy in GeV
  else if( mEnergy == 200) pythia.readString("Beams:eCM = 200");
  else
  {
    cout<<"Invalid input energy!"<<endl;
    return 0;
  
  }
  
  pythia.readString("HardQCD:hardccbar = on");
  pythia.readString("HardQCD:gg2ccbar = on");
  pythia.readString("HardQCD:qqbar2ccbar = on");
  
  
  
  pythia.readString("413:onMode = 0");
  pythia.readString("413:onIfAll = 421 211");
  pythia.readString("421:onMode = 0");
  pythia.readString("421:onIfAll = -321 211");

  pythia.readString("-413:onMode = 0");
  pythia.readString("-413:onIfAll = -421 -211");
  pythia.readString("-421:onMode = 0");
  pythia.readString("-421:onIfAll = 321 -211");
  
/*
  //old tune - not all options work locally
  pythia.readString("Tune:pp = 6");

  //http://home.thep.lu.se/~torbjorn/pythia81html/Tunes.html
  //option 6 : "Tune 4Cx", based on tune 4C, but using the x-dependent matter profile,
  //MultipartonInteractions:bProfile = 4 and an increased MultipartonInteractions:pT0Ref [Cor11]. 

  pythia.readString("SigmaProcess:renormScale2 = 3");
  pythia.readString("SigmaProcess:factorScale2 = 3");
  pythia.readString("SigmaProcess:renormMultFac = 2"); //2mT
  pythia.readString("SigmaProcess:factorMultFac = 2");
  pythia.readString(Form("PhaseSpace:pTHatMin = %f",ptHatMin));
  pythia.readString(Form("PhaseSpace:pTHatMax = %f",ptHatMax));

  pythia.readString("PDF:useLHAPDF = on");
  pythia.readString("PDF:LHAPDFset = MRSTMCal.LHgrid");
  pythia.readString("PDF:extrapolateLHAPDF = on");
  pythia.readString("PartonLevel:MI = on");
  pythia.readString("PartonLevel:ISR = on");
  pythia.readString("BeamRemnants:primordialKT = on");
  pythia.readString("PartonLevel:FSR = on");
  pythia.readString("StringFlav:mesonCvector = 1.5");
  pythia.readString("StringFlav:mesonBvector = 3");
  pythia.readString("4:m0 = 1.43");
  pythia.readString("5:m0 = 4.30");
  
  
  //______________________________________________________________
  //new tune - runs slow

  pythia.readString("PDF:pSet = 17");
  pythia.readString("MultipartonInteractions:ecmRef = 200");
  pythia.readString("MultipartonInteractions:bprofile = 2");
  pythia.readString("MultipartonInteractions:pT0Ref = 0.140");
  pythia.readString("MultipartonInteractions:ecmPow  = 0.135");
  pythia.readString("MultipartonInteractions:coreRadius = 0.56");
  pythia.readString("MultipartonInteractions:coreFraction = 0.78");
  pythia.readString("ColourReconnection:range = 5.4");

  pythia.readString(Form("PhaseSpace:pTHatMin = %f",ptHatMin));
  pythia.readString(Form("PhaseSpace:pTHatMax = %f",ptHatMax));


*/
  pythia.init();

  // Set up the ROOT TFile and TTree.
  TFile *file = TFile::Open(argv[3],"recreate");
  
  Event *event = &pythia.event;

  //TTree *T = new TTree("T","ev1 Tree");
  //T->Branch("event",&event);
  
  TH2D *kaon_phi_vs_pion_phi_D0[nPtBins];
  TH2D *kaon_phi_vs_pion_phi_Dstar[nPtBins];
  TH2D *pion_phi_Dstar_vs_pion_phi_D0[nPtBins];
  
  TH2D *kaon_eta_vs_pion_eta_D0[nPtBins];
  TH2D *kaon_eta_vs_pion_eta_Dstar[nPtBins];
  TH2D *pion_eta_Dstar_vs_pion_eta_D0[nPtBins];
  
  TH1D *nCorr_K_pi_D0[nPtBins];
  TH1D *nCorr_K_pi_Dstar[nPtBins];
  TH1D *nCorr_pi_Dstar_pi_D0[nPtBins];

  TH1D *nCorr_clust_K_pi_D0[nPtBins];
  TH1D *nCorr_clust_K_pi_Dstar[nPtBins];
  TH1D *nCorr_clust_pi_Dstar_pi_D0[nPtBins];
  
  for(unsigned int pTbin = 0; pTbin < nPtBins; pTbin++)
  {
    kaon_phi_vs_pion_phi_D0[pTbin] = new TH2D(Form("kaon_phi_vs_pion_phi_D0_pT_%i", pTbin), Form("kaon_phi_vs_pion_phi_D0_pT_%i", pTbin), 120, -Pi(), Pi(), 120, -Pi(), Pi()); //phi of K and pi from decay of D0
    kaon_phi_vs_pion_phi_Dstar[pTbin] = new TH2D(Form("kaon_phi_vs_pion_phi_Dstar_pT_%i", pTbin), Form("kaon_phi_vs_pion_phi_Dstar_pT_%i", pTbin), 120, -Pi(), Pi(), 120, -Pi(), Pi()); //phi of K from decay of D0, and pi from decay of D*
    pion_phi_Dstar_vs_pion_phi_D0[pTbin] = new TH2D(Form("pion_phi_Dstar_vs_pion_phi_D0_pT_%i", pTbin), Form("pion_phi_Dstar_vs_pion_phi_D0_pT_%i", pTbin), 120, -Pi(), Pi(), 120, -Pi(), Pi()); //phi of decay of pi from decay of D0 and pi from decay of D*
    
    kaon_eta_vs_pion_eta_D0[pTbin] = new TH2D(Form("kaon_eta_vs_pion_eta_D0_pT_%i", pTbin), Form("kaon_eta_vs_pion_eta_D0_pT_%i", pTbin), 40, -1, 1, 40, -1, 1); //phi of K and pi from decay of D0
    kaon_eta_vs_pion_eta_Dstar[pTbin] = new TH2D(Form("kaon_eta_vs_pion_eta_Dstar_pT_%i", pTbin), Form("kaon_eta_vs_pion_eta_Dstar_pT_%i", pTbin), 40, -1, 1, 40, -1, 1); //phi of K from decay of D0, and pi from decay of D*
    pion_eta_Dstar_vs_pion_eta_D0[pTbin] = new TH2D(Form("pion_eta_Dstar_vs_pion_eta_D0_pT_%i", pTbin), Form("pion_eta_Dstar_vs_pion_eta_D0_pT_%i", pTbin), 40, -1, 1, 40, -1, 1); //phi of decay of pi from decay of D0 and pi from decay of D*  
    
    nCorr_K_pi_D0[pTbin] = new TH1D(Form("nCorr_K_pi_D0_pT_%i", pTbin), Form("nCorr_K_pi_D0_pT_%i", pTbin), 2, 0, 2);
    nCorr_K_pi_Dstar[pTbin] = new TH1D(Form("nCorr_K_pi_Dstar_pT_%i", pTbin), Form("nCorr_K_pi_Dstar_pT_%i", pTbin), 2, 0, 2);
    nCorr_pi_Dstar_pi_D0[pTbin] = new TH1D(Form("nCorr_pi_Dstar_pi_D0_pT_%i", pTbin), Form("nCorr_pi_Dstar_pi_D0_pT_%i", pTbin), 2, 0, 2);

    nCorr_clust_K_pi_D0[pTbin] = new TH1D(Form("nCorr_clust_K_pi_D0_pT_%i", pTbin), Form("nCorr_clust_K_pi_D0_pT_%i", pTbin), 2, 0, 2);
    nCorr_clust_K_pi_Dstar[pTbin] = new TH1D(Form("nCorr_clust_K_pi_Dstar_pT_%i", pTbin), Form("nCorr_clust_K_pi_Dstar_pT_%i", pTbin), 2, 0, 2);
    nCorr_clust_pi_Dstar_pi_D0[pTbin] = new TH1D(Form("nCorr_clust_pi_Dstar_pi_D0_pT_%i", pTbin), Form("nCorr_clust_pi_Dstar_pi_D0_pT_%i", pTbin), 2, 0, 2);
  
  }

 // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < nEvents; ++iEvent)
  {
    if (!pythia.next()) continue;
    
    //if( (iEvent % 10) == 0) cout<<iEvent<<endl;
    
    for(unsigned int iPart = 0; iPart < pythia.event.size(); iPart++)
    {
      if( fabs( pythia.event[iPart].id() ) != 413  ) continue; //continue only for D*+/-
      
      if( fabs(pythia.event[iPart].eta()) > 1 ) continue;
      
      TLorentzVector D_star_fourmom;
      D_star_fourmom.SetPxPyPzE(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), pythia.event[iPart].e());
      
      //fill all histograms for all pT and centrality bins
      int pT_bin = -1;

      //find pT bin of Lambda
      for(int j = 0; j < nPtBins; j++) //loop over pT bins
      {
        if(D_star_fourmom.Pt() > pT_bins[j] && D_star_fourmom.Pt() <= pT_bins[j+1])
        {
          pT_bin = j;
          break; //stop after pT bin is found
        }
      }

      if( pT_bin == -1 ) continue;
     
      //cout<<"D* found!"<<endl;
      
      int firstDaughter = pythia.event[iPart].daughter1(); //ID of D0 from decay of D*
      int secondDaughter = pythia.event[iPart].daughter2(); // ID of pi from decay of D*
      
      if( fabs( pythia.event[firstDaughter].id() ) == 421 && fabs( pythia.event[secondDaughter].id() ) == 211 )
      {
        //cout<<"First daughter is D0"<<endl;
        
        if( fabs(pythia.event[firstDaughter].eta()) > 1 ) continue;
        if( fabs(pythia.event[secondDaughter].eta()) > 1 ) continue;
        
        int firstDaughter_D0 = pythia.event[firstDaughter].daughter1(); //ID of K from decay of D0
        int secondDaughter_D0 = pythia.event[firstDaughter].daughter2(); //ID of pi from decay of D0
        
        if( fabs(pythia.event[firstDaughter_D0].id()) == 321 && fabs( pythia.event[secondDaughter_D0].id() ) == 211 )
        {
          //cout<<"First daughter is K"<<endl;
          
          if( fabs(pythia.event[firstDaughter_D0].eta()) > 1 ) continue;
          if( fabs(pythia.event[secondDaughter_D0].eta()) > 1 ) continue;
          
          kaon_phi_vs_pion_phi_D0[pT_bin]->Fill(pythia.event[firstDaughter_D0].phi(), pythia.event[secondDaughter_D0].phi());
          kaon_phi_vs_pion_phi_Dstar[pT_bin]->Fill(pythia.event[firstDaughter_D0].phi(), pythia.event[secondDaughter].phi());
          pion_phi_Dstar_vs_pion_phi_D0[pT_bin]->Fill(pythia.event[secondDaughter].phi(), pythia.event[secondDaughter_D0].phi());
          
          kaon_eta_vs_pion_eta_D0[pT_bin]->Fill(pythia.event[firstDaughter_D0].eta(), pythia.event[secondDaughter_D0].eta());
          kaon_eta_vs_pion_eta_Dstar[pT_bin]->Fill(pythia.event[firstDaughter_D0].eta(), pythia.event[secondDaughter].eta());
          pion_eta_Dstar_vs_pion_eta_D0[pT_bin]->Fill(pythia.event[secondDaughter].eta(), pythia.event[secondDaughter_D0].eta());
          
          int bin_phi_pi_D0 = kaon_phi_vs_pion_phi_D0[pT_bin]->GetXaxis()->FindBin(pythia.event[secondDaughter_D0].phi());
          int bin_phi_K_D0 = kaon_phi_vs_pion_phi_D0[pT_bin]->GetXaxis()->FindBin(pythia.event[firstDaughter_D0].phi());
          int bin_phi_pi_D_star = kaon_phi_vs_pion_phi_Dstar[pT_bin]->GetXaxis()->FindBin(pythia.event[secondDaughter].phi());
          
          int bin_eta_pi_D0 = kaon_eta_vs_pion_eta_D0[pT_bin]->GetXaxis()->FindBin(pythia.event[secondDaughter_D0].eta());
          int bin_eta_K_D0 = kaon_eta_vs_pion_eta_D0[pT_bin]->GetXaxis()->FindBin(pythia.event[firstDaughter_D0].eta());
          int bin_eta_pi_D_star = kaon_eta_vs_pion_eta_Dstar[pT_bin]->GetXaxis()->FindBin(pythia.event[secondDaughter].eta());

          // K and pi from decay of D0 in the same BEMC tower
          if( (bin_phi_K_D0 == bin_phi_pi_D0) && ( bin_eta_K_D0 == bin_eta_pi_D0 ) )
          {
            nCorr_K_pi_D0[pT_bin]->Fill(1.5);          
          }
          else
          {
            nCorr_K_pi_D0[pT_bin]->Fill(0.5);          
          }
          
          // K from decay of D0, pi from decay of D* in the same BEMC tower
          if( (bin_phi_K_D0 == bin_phi_pi_D_star) && ( bin_eta_K_D0 == bin_eta_pi_D_star ) )
          {
            nCorr_K_pi_Dstar[pT_bin]->Fill(1.5);          
          }
          else
          {
            nCorr_K_pi_Dstar[pT_bin]->Fill(0.5);          
          }
          
          // pi from decay of D0, pi from decay of D* in the same BEMC tower
          if( (bin_phi_pi_D0 == bin_phi_pi_D_star) && ( bin_eta_pi_D0 == bin_eta_pi_D_star ) )
          {
            nCorr_pi_Dstar_pi_D0[pT_bin]->Fill(1.5);          
          }
          else
          {
            nCorr_pi_Dstar_pi_D0[pT_bin]->Fill(0.5);          
          }
          //__________________________________________________
          
          // K and pi from decay of D0 in the same BEMC cluster
          if( fabs(bin_phi_K_D0 - bin_phi_pi_D0) < 2 && fabs( bin_eta_K_D0 - bin_eta_pi_D0 ) < 2 )
          {
            nCorr_clust_K_pi_D0[pT_bin]->Fill(1.5);          
          }
          else
          {
            nCorr_clust_K_pi_D0[pT_bin]->Fill(0.5);          
          }
          
          // K from decay of D0, pi from decay of D* in the same BEMC cluster
          if( fabs(bin_phi_K_D0 - bin_phi_pi_D_star) < 2 && fabs( bin_eta_K_D0 - bin_eta_pi_D_star ) < 2 )
          {
            nCorr_clust_K_pi_Dstar[pT_bin]->Fill(1.5);          
          }
          else
          {
            nCorr_clust_K_pi_Dstar[pT_bin]->Fill(0.5);          
          }
          
          // pi from decay of D0, pi from decay of D* in the same BEMC cluster
          if( fabs(bin_phi_pi_D0 - bin_phi_pi_D_star) < 2 && fabs( bin_eta_pi_D0 - bin_eta_pi_D_star ) < 2 )
          {
            nCorr_clust_pi_Dstar_pi_D0[pT_bin]->Fill(1.5);          
          }
          else
          {
            nCorr_clust_pi_Dstar_pi_D0[pT_bin]->Fill(0.5);          
          }
          
          
          break; //stop going through event when good D* decay is found
        
        }
        
        break;       
      
      }
      
      
      
      break;  
          
    
    }
    
    pythia.event.free();

  // End event loop.
  }

  // Statistics on event generation.
  pythia.stat();

  //  Write tree.
  //T->Print();
  //T->Write();
  file->Write();
  file->Close();

  // Done.
  return 0;
}
