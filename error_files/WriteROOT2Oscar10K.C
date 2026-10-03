#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TVector.h>

#include "/phenix/plhf/roli/HELIOS/WriteEvent.h"
using namespace std;

void WriteROOT2Oscar10K(TString infile = "single_pi0_HELIOS_1B.root", TString output = "oscar.txt", Int_t nfolder = 0, int cent = 2){
  gSystem->Load("/phenix/plhf/roli/HELIOS/simulation/WriteEvent_C.so");
  //gSystem->Load("/phenix/plhf/tongzhouguo/occu/gen_sim/WriteEvent_C.so");
  TFile* input = new TFile(infile,"READ");
  if(!(input)){
    cout << "no input file" << endl;
    exit(1);
  }

  ofstream file(output);
  int nstart = nfolder * 10000;
  //int nend = nstart + 10000;

  //Read in the TTrees 
  TTree* T = (TTree*)input->Get("T");
  TBranch* br = T->GetBranch("MyEvent");

  WriteEvent* event = 0;
  br->SetAddress(&event);

  file << "# OSC1999A" << endl;
  file << "# final_id_p_x" << endl;
  file << "# SimName 1.0" << endl;
  file << "#" << endl;
  file << "# Some comments..." << endl;
  file << endl;

  TFile* fvz = new TFile("/phenix/plhf/tongzhouguo/error_files/hvz_CuAu.root", "READ");
  TH1D* hvz = (TH1D*) fvz->Get(Form("hvz_cent%d",cent));
  hvz->SetDirectory(0);
int iwritten =0;
int ievt = nstart;
  while (iwritten < 10000){
    double zpos = hvz->GetRandom();
    //event->ClearEvent();
    //int ievt = nstart+iwritten;
    br->GetEntry(ievt);
    ievt++;
    //file << 0 << "\t" << event->GetNStable() << endl;
    bool do_procedure = false;

    for(int i = 0; i < event->GetNEntries(); i++){
	    WriteTrack Track = event->GetWriteTrack(i);
	    if (Track.GetID()==111){
	      if (sqrt(pow(Track.GetPx(),2)+pow(Track.GetPy(),2))<1000000.0) do_procedure = true;
	    }
    }

    if (!do_procedure) continue;
file << 0 << "\t" << event->GetNStable() << endl;
    for(int i = 0; i < event->GetNEntries(); i++){
      
      WriteTrack Track = event->GetWriteTrack(i);
      if(Track.GetFinal() == 1){
        file << Track.GetNum() << "\t" << Track.GetID() << "\t" << 0 << "\t" << Track.GetPx() << "\t" << Track.GetPy() << "\t" << Track.GetPz() << "\t" << Track.GetEnergy() << "\t" << Track.GetMass() << "\t" << "0" << "\t" << "0" << "\t" << zpos*pow(10,13) << "\t" << "0" << endl;
      }
    }
    file << "0" << "\t" << "0" << endl;
  iwritten++;
  }
}
