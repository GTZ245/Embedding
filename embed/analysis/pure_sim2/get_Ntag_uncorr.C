#include <math.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH3D.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TMatrixD.h>
#include <TNtuple.h>
#include <TF1.h>
#include <TRandom3.h>

#include "/direct/phenix+u/tongzhouguo/install/include/DileptonAnalysis/MyEvent.h"
#include "/direct/phenix+u/roli/scratch/HELIOS/WriteEvent.h"
#include "/direct/phenix+u/tongzhouguo/install/include/DileptonAnalysis/Reconstruction.h"
#include "dc_dead_map.h"

using namespace std;
using namespace DileptonAnalysis;

static const double PI = 3.14159;
static const double Me = 0.000510998918;  // Particle Data Group 7-25-2004
static const double Me2 = Me*Me;
static const double c = 299792458;

const double PC1_DPHI_CUT = 0.02, PC1_DZ_CUT = 0.5; // PC1 ghost cut
const int centbin = 1;
const int sectbin = 1;
const int cutbin = 13;
float w_pt = 1;

const char* cut_name[cutbin] = {
  "nominal", "conv_solution", "eid0", "eid1", "eid3",
  "pte02", "pte04", "ecore03", "ecore05", "chi2_4",
  "chi2_5", "escale099", "escale101"
};

// cut0: conversion solution + dzed < 4
// cut1: conversion solution only
const bool APPLY_DZED_CUT[cutbin] = {
  true, false, true, true, true,
  true, true, true, true, true,
  true, true, true
};

const double DZED_CUT[cutbin] = {
  4.0, 4.0, 4.0, 4.0, 4.0,
  4.0, 4.0, 4.0, 4.0, 4.0,
  4.0, 4.0, 4.0
};

const double DCENTER_CUT[cutbin] = {
  10.0, 10.0, 10.0, 10.0, 10.0,
  10.0, 10.0, 10.0, 10.0, 10.0,
  10.0, 10.0, 10.0
};

const double R_LO_CUT[cutbin] = {
  1.0, 1.0, 1.0, 1.0, 1.0,
  1.0, 1.0, 1.0, 1.0, 1.0,
  1.0, 1.0, 1.0
};

const double R_HI_CUT[cutbin] = {
  29.0, 29.0, 29.0, 29.0, 29.0,
  29.0, 29.0, 29.0, 29.0, 29.0,
  29.0, 29.0, 29.0
};

const double DPHIR_CUT[cutbin] = {
  0.005, 0.005, 0.005, 0.005, 0.005,
  0.005, 0.005, 0.005, 0.005, 0.005,
  0.005, 0.005, 0.005
};

const double ECORE_CUT[cutbin] = {
  0.4, 0.4, 0.4, 0.4, 0.4,
  0.4, 0.4, 0.3, 0.5, 0.4,
  0.4, 0.4, 0.4
};

const double CHI2_CUT[cutbin] = {
  3.0, 3.0, 3.0, 3.0, 3.0,
  3.0, 3.0, 3.0, 3.0, 4.0,
  5.0, 3.0, 3.0
};

const double TRK_PT_CUT[cutbin] = {
  0.30, 0.30, 0.30, 0.30, 0.30,
  0.20, 0.40, 0.30, 0.30, 0.30,
  0.30, 0.30, 0.30
};

const double EP_LO_CUT[cutbin] = {
  0.7, 0.7, 0.5, 0.6, 0.8,
  0.7, 0.7, 0.7, 0.7, 0.7,
  0.7, 0.7, 0.7
};

const double EP_HI_CUT[cutbin] = {
  1.0e9, 1.0e9, 1.0e9, 1.0e9, 1.0e9,
  1.0e9, 1.0e9, 1.0e9, 1.0e9, 1.0e9,
  1.0e9, 1.0e9, 1.0e9
};

const double N0_CUT[cutbin] = {
  2.0, 2.0, 2.0, 2.0, 4.0,
  2.0, 2.0, 2.0, 2.0, 2.0,
  2.0, 2.0, 2.0
};

const double EMCDZ_CUT[cutbin] = {
  20.0, 20.0, 20.0, 25.0, 10.0,
  20.0, 20.0, 20.0, 20.0, 20.0,
  20.0, 20.0, 20.0
};

const double EMCDPHI_CUT[cutbin] = {
  0.05, 0.05, 0.05, 0.06, 0.03,
  0.05, 0.05, 0.05, 0.05, 0.05,
  0.05, 0.05, 0.05
};

const double PHOTON_ESCALE[cutbin] = {
  1.00, 1.00, 1.00, 1.00, 1.00,
  1.00, 1.00, 1.00, 1.00, 1.00,
  1.00, 0.99, 1.01
};

float ecore_cut = 0.5;
float chi2_cut = 3.0;
double bbcz_cut_val = 10;
const double bbcq_cut_val = 1800;
float trk_pT_cut = 0.4;
float dz_cut = 10;
float dphi_cut = 0.02;

float EMCMAP[8][48][96];
const double POOLD = 100; // pool depth

// Scale/smear arrays are intentionally left unused in this version.
// The single-electron Ecore treatment follows the embedding code: use GetEcore() directly.
const double scale[8] = {
  0.988, 0.990, 0.985, 0.981,
  0.980, 0.985, 1.025, 1.023
};

const double smear_c1[8] = {
  0.055, 0.066, 0.055, 0.059,
  0.057, 0.056, 0.072, 0.072
};

const double smear_c2[8] = {
  0.011, 0.017, 0.011, 0.013,
  0.012, 0.012, 0.063, 0.062
};

TH1D* h1d_truepiz;
TH1D* h1d_truepiz_unw;
TH2D* eemass[cutbin][centbin][sectbin];
TH2D* eegmass[cutbin][centbin][sectbin];
TH2D* eemass_unw[cutbin][centbin][sectbin];
TH2D* eegmass_unw[cutbin][centbin][sectbin];
TFile *ratio[10];
TGraph *gr[10];
TLorentzVector pion; 

double getDcenter(double phi1, double z1, double phi2, double z2){
  double dcenter_phi_sigma = 0.01, dcenter_phi_offset = 0.;
  double dcenter_z_sigma = 3.6, dcenter_z_offset = 0.;

  double dcenter_z = (z1-z2-dcenter_z_offset)/dcenter_z_sigma;
  double dcenter_phi = (phi1-phi2-dcenter_phi_offset)/dcenter_phi_sigma;

  return sqrt(dcenter_phi*dcenter_phi+dcenter_z*dcenter_z);
}

void read_in_emcmap(){
  const char* directory_file =
    "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/pi02gg/directory.txt";

  ifstream directory_stream(directory_file);
  if (!directory_stream.is_open()){
    cerr << "ERROR: cannot open EMCal dead-map directory file: "
         << directory_file << endl;
    exit(1);
  }

  string emcmap_path;
  getline(directory_stream, emcmap_path);
  directory_stream.close();

  while (!emcmap_path.empty() &&
         (emcmap_path[emcmap_path.size()-1] == '\r' ||
          emcmap_path[emcmap_path.size()-1] == '\n' ||
          emcmap_path[emcmap_path.size()-1] == ' ' ||
          emcmap_path[emcmap_path.size()-1] == '\t')){
    emcmap_path.erase(emcmap_path.size()-1);
  }

  if (emcmap_path.empty()){
    cerr << "ERROR: no EMCal dead-map path found in: "
         << directory_file << endl;
    exit(1);
  }

  ifstream readmap(emcmap_path.c_str());
  if (!readmap.is_open()){
    cerr << "ERROR: cannot open EMCal dead map: "
         << emcmap_path << endl;
    exit(1);
  }

  for (int i = 0; i < 8; ++i){
    for (int j = 0; j < 48; ++j){
      for (int k = 0; k < 96; ++k){
        EMCMAP[i][j][k] = 0;
      }
    }
  }

  int armsect = 0, ypos = 0, zpos = 0, status = 0;
  cout << "reading EMCal dead map: " << emcmap_path << endl;
  while (readmap >> armsect >> ypos >> zpos >> status){
    EMCMAP[armsect][ypos][zpos] = status;
  }
  readmap.close();
}

bool apply_emcmap(int sect, int y, int z){
  if (y==0 || z==0) return false;
  if ( sect < 6 && ( y == 35 || z == 71) ) return false;
  if ( sect > 5 && ( y == 47 || z == 95) ) return false;
  if (EMCMAP[sect][y][z]) return false;
  return true;
}

bool single_cut(MyTrack* mytrk, int icut){
  int _side = mytrk->GetDCside();
  int _arm = mytrk->GetArm();
  float _phi = mytrk->GetPhiDC();
  float _alpha = mytrk->GetAlpha(); 
  float _zed = mytrk->GetZDC();

  float _board = -999;
  if(_arm == 0) _board = (3.72402-_phi+0.008047*cos(_phi+0.87851))/0.01963496;
  else if(_arm == 1) _board = (0.573231+_phi-0.0046*cos(_phi+0.05721))/0.01963496;
  if(_board < -99) return false;

  bool dead_area = DCdeadArea(_alpha, _phi, _board, _zed);
  if(dead_area) return false;

  double _p = sqrt(mytrk->GetPx()*mytrk->GetPx() + mytrk->GetPy()*mytrk->GetPy() + mytrk->GetPz()*mytrk->GetPz());
  double _theta = mytrk->GetThe0();
  double _pt = _p * sin(_theta);
  if (_pt < TRK_PT_CUT[icut]) return false;

  double _e = mytrk->GetEcore();
  double _ep = _e/_p;
  if(_ep < EP_LO_CUT[icut]) return false;
  if(_ep > EP_HI_CUT[icut]) return false;

  double _n0 = mytrk->GetN0();
  if(_n0 < N0_CUT[icut]) return false;

  int _charge = mytrk->GetCharge();
  if(fabs(_charge) != 1) return false;

  double _emcdz = mytrk->GetEmcdz();
  if(fabs(_emcdz) > EMCDZ_CUT[icut]) return false;

  double _emcdphi = mytrk->GetEmcdphi();
  if(fabs(_emcdphi) > EMCDPHI_CUT[icut]) return false;

  return true;
}

float radius_r;
float phi_positron;
float phi_electron;
float theta_positron;  
float theta_electron;  
bool solution(MyTrack* mytrk1, MyTrack* mytrk2, MyPair* mypair, Reconstruction* reco, float zVtx, int icut){
  reco->findIntersection(mytrk1, mytrk2, mypair, zVtx);
  radius_r = mypair->GetRPair();
  phi_electron = mypair->GetPhiElectron();
  phi_positron = mypair->GetPhiPositron();
  theta_electron = mypair->GetThetaElectron();
  theta_positron = mypair->GetThetaPositron();
  if ( (radius_r<=R_LO_CUT[icut]) || (radius_r>=R_HI_CUT[icut]) ) return false;

  float dphi_r = mypair->GetPhiElectron()-mypair->GetPhiPositron();
  if ( dphi_r>TMath::Pi() )  dphi_r = 2*TMath::Pi()-dphi_r;
  if ( dphi_r<-TMath::Pi() ) dphi_r = -2*TMath::Pi()-dphi_r;
  if ( TMath::Abs(dphi_r)>=DPHIR_CUT[icut] ) return false;
  
  return true;
}

//void get_Ntag_uncorr(const char* inFile = "/phenix/plhf/tongzhouguo/test_tree1/trees/sim_trees_00182.root", const char* outFile = "./out.root", const int system = 0, int num = 182, int iter = 0){
void get_Ntag_uncorr(const char* inFile = "/phenix/plhf/tongzhouguo/test_tree1/trees/sim_trees_08475.root", const char* outFile = "./out.root", const int system = 0, int num = 8475, int iter = 0){	
  gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libDileptonAnalysisEvent");
  gSystem->Load("/direct/phenix+u/roli/scratch/HELIOS/WriteEvent_C.so");
  gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libDileptonAnalysisReco");
  Reconstruction reco("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/yuri_embed/embed/analysis/emb/lookup_3D_one_phi.root");

  read_in_emcmap();

  TF1* fHagedorn = new TF1("fHagedorn","6.2831853*x*((1/(1+exp((x-[5])/[6])))*[0]/(1+x/[1])**[2] + (1-(1/(1+exp((x-[5])/[6]))))*[3]/x**[4])",0.0,10.0);
  fHagedorn->SetParameters(1.922170e+02, 2.106180e+00, 1.284030e+01, 1.630000e+01, 8.060480e+00, 4.033000e+00, 6.395340e-02);
  fHagedorn->SetNpx(10000);
  //fHagedorn->SetParameters(2.990550e+01, 1.768200e+00, 1.146540e+01, 3.837460e+00, 8.222340e+00, 4.555530e+00, 5.845620e-02);

  TString inFile_A = inFile;
  TString inFile_B = inFile;

  TFile* input_A = new TFile(inFile_A,"READ");
  if(!(input_A)){
    cout << "no input_A file" << endl;
    exit(1);
  }
  TFile* input_B = new TFile(inFile_B,"READ");
  if(!(input_B)){
    cout << "no input_B file" << endl;
    exit(1);
  }
  
  TTree* T_A = (TTree*)input_A->Get("T");
  TBranch* br_A = T_A->GetBranch("MyEvent");
  MyEvent* event = 0;
  br_A->SetAddress(&event);
  int nevt = T_A->GetEntries();
  //if(nevt < 100) exit(1);
  cout << "Input Tree read in" << endl;
  TTree* T_B = (TTree*)input_B->Get("T");
  TBranch* br_B = T_B->GetBranch("MyEvent");
  MyEvent* event_B = 0;
  br_B->SetAddress(&event_B);
  cout << "Trees read!" << endl; 
  TFile* output = new TFile(outFile,"RECREATE");

  int nfolder = num; //root file number
  int nstart = nfolder * 10000;
  TFile* helios = new TFile("/direct/phenix+u/roli/scratch/HELIOS/simulation/singlePi0.root","READ");
  TTree* helios_T = (TTree*)helios->Get("T");
  TBranch* helios_br = helios_T->GetBranch("MyEvent");
  WriteEvent* helios_event = 0;
  helios_br->SetAddress(&helios_event);
  cout << "HELIOS Tree read in" << endl;

  for(int icent = 0; icent < centbin; ++icent){
    for (int icut = 0; icut < cutbin; ++icut){
      for (int isect = 0; isect < sectbin; ++isect){
        eemass[icut][icent][isect] = new TH2D(Form("eemass_%s_cent%d_sect%d",cut_name[icut],icent,isect), "ee mass versus ee pt simulation", 200, 0.0, 1.0, 400, 0.0, 40.0);
        eegmass[icut][icent][isect] = new TH2D(Form("eegmass_%s_cent%d_sect%d",cut_name[icut],icent,isect), "pi0 mass versus pi0 pt simulation", 200, 0.0, 1.0, 400, 0.0, 40.0);
        eemass_unw[icut][icent][isect] = new TH2D(Form("eemass_unw_%s_cent%d_sect%d",cut_name[icut],icent,isect), "ee mass versus ee pt simulation", 200, 0.0, 1.0, 400, 0.0, 40.0);
        eegmass_unw[icut][icent][isect] = new TH2D(Form("eegmass_unw_%s_cent%d_sect%d",cut_name[icut],icent,isect), "pi0 mass versus pi0 pt simulation", 200, 0.0, 1.0, 400, 0.0, 40.0);
        eemass[icut][icent][isect]->Sumw2();
        eegmass[icut][icent][isect]->Sumw2();
        eemass_unw[icut][icent][isect]->Sumw2();
        eegmass_unw[icut][icent][isect]->Sumw2();
      }
    }
  }

  h1d_truepiz = new TH1D("h1d_truepiz", "number of pizero generated", 400, 0.0, 40.0);
  h1d_truepiz_unw = new TH1D("h1d_truepiz_unw", "number of pizero generated", 400, 0.0, 40.0); 
  h1d_truepiz->Sumw2();
  h1d_truepiz_unw->Sumw2();
  cout << "Histograms defined" << endl;

  for(int ievt = 0; ievt < 5000; ievt++){
    helios_br->GetEntry(abs(nstart+ievt));
    for(int i = 0; i < helios_event->GetNEntries(); i++){
      WriteTrack helios_track = helios_event->GetWriteTrack(i);
      if(helios_track.GetID() == 111){
        pion.SetPxPyPzE(helios_track.GetPx(), helios_track.GetPy(), helios_track.GetPz(), helios_track.GetEnergy());
      }
    }
    w_pt = fHagedorn->Eval(pion.Pt());
    for(int i = 0; i < iter; i++){
      w_pt = w_pt * gr[i]->Eval(pion.Pt());
    }
    h1d_truepiz->Fill(pion.Pt(), w_pt);
    h1d_truepiz_unw->Fill(pion.Pt());
  }

  //==============================
  //         Event A Loop
  //==============================
  for (int ievent_A = 0; ievent_A < nevt; ievent_A++){
    if (ievent_A%1000==0) cout << "Event: " << ievent_A << " / " << nevt << endl;

    event->ClearEvent();
    br_A->GetEntry(ievent_A);
  
    int ntrack_A = event->GetNtrack();
    const int nclust_A = event->GetNcluster();
    double zVtx_A = event->GetVtxZ();
    float bbcq_A = event->GetBBCcharge();

    int id = event->GetEvtNo() - 1;
    helios_br->GetEntry(abs(nstart+id));

    for(int i = 0; i < helios_event->GetNEntries(); i++){
      WriteTrack helios_track = helios_event->GetWriteTrack(i);
      if(helios_track.GetID() == 111){
        pion.SetPxPyPzE(helios_track.GetPx(), helios_track.GetPy(), helios_track.GetPz(), helios_track.GetEnergy());
      }
    }

    w_pt = fHagedorn->Eval(pion.Pt()); 
    for(int i = 0; i < iter; i++){
      w_pt = w_pt * gr[i]->Eval(pion.Pt());
    }
    if(nclust_A == 0) continue;

    if ( ntrack_A == 0 ) continue;
    if ( TMath::Abs(zVtx_A)>bbcz_cut_val ) continue;
    if ( bbcq_A>bbcq_cut_val ) continue;

    int icent = 0;


    for (int icut = 0; icut < cutbin; ++icut){
      //===========================================
      //       Track Selection for Event A
      //===========================================
      int ghost_A[40], real_A[40];
      int all_A[40];
      int nghost_A = 0, nreal_A = 0, nall_A = 0;
      for (int k1 = 0; k1 < ntrack_A; ++k1){
        if ( !single_cut(&(event->GetEntry(k1)), icut)) continue; // eid cuts
        all_A[nall_A] = k1;
        ++nall_A;        
        int ghost_flag = 0;

        for (int k2 = 0; k2 < ntrack_A; ++k2){
          if ( k1==k2 ) continue;
          if ( !single_cut(&(event->GetEntry(k2)), icut) ) continue;
          // found a RICH ghost pair
          if ( getDcenter((event->GetEntry(k1)).GetCrkphi(), (event->GetEntry(k1)).GetCrkz(), (event->GetEntry(k2)).GetCrkphi(), (event->GetEntry(k2)).GetCrkz()) < DCENTER_CUT[icut] ) { ghost_flag = 1; break; }
        }
        if ( ghost_flag ) { ghost_A[nghost_A] = k1; ++nghost_A; }
        else { real_A[nreal_A] = k1; ++nreal_A; }
      }
      if(nreal_A == 0) continue;

      //======================================
      //      Electron/Positron A1 Loop
      //======================================
      for (int ireal_A1 = 0; ireal_A1 < nreal_A; ireal_A1++){
        MyTrack mytrk_A1 = event->GetEntry(real_A[ireal_A1]);
        int id_A1 = mytrk_A1.GetEmcId();
        int arm_A1 = mytrk_A1.GetArm();
        int charge_A1 = mytrk_A1.GetCharge();
        float the_A1 = mytrk_A1.GetThe0();
        float phi_A1 = mytrk_A1.GetPhi0();
        float phidc_A1 = mytrk_A1.GetPhiDC();
        float zed_A1 = mytrk_A1.GetZDC();
        float mom_A1 = sqrt(mytrk_A1.GetPx()*mytrk_A1.GetPx() + mytrk_A1.GetPy()*mytrk_A1.GetPy() + mytrk_A1.GetPz()*mytrk_A1.GetPz());
        float pt_A1 = mom_A1 * sin(the_A1);

        float momx_A1 = pt_A1 * cos(phi_A1);
        float momy_A1 = pt_A1 * sin(phi_A1);
        float momz_A1 = pt_A1/tan(the_A1);

        //===================================
        //     Electron/Positron A2 Loop
        //===================================
        for (int ireal_A2 = ireal_A1+1; ireal_A2 < nreal_A; ireal_A2++){
          MyTrack mytrk_A2 = event->GetEntry(real_A[ireal_A2]);
          int id_A2 = mytrk_A2.GetEmcId();
          int arm_A2 = mytrk_A2.GetArm();
          int charge_A2 = mytrk_A2.GetCharge();
          float the_A2 = mytrk_A2.GetThe0();
          float phi_A2 = mytrk_A2.GetPhi0();
          float phidc_A2 = mytrk_A2.GetPhiDC();
          float zed_A2 = mytrk_A2.GetZDC();
          float mom_A2 = sqrt(mytrk_A2.GetPx()*mytrk_A2.GetPx() + mytrk_A2.GetPy()*mytrk_A2.GetPy() + mytrk_A2.GetPz()*mytrk_A2.GetPz());
          float pt_A2 = mom_A2 * sin(the_A2);

          float momx_A2 = pt_A2 * cos(phi_A2);
          float momy_A2 = pt_A2 * sin(phi_A2);
          float momz_A2 = pt_A2/tan(the_A2);

          if (fabs(phidc_A1-phidc_A2)<PC1_DPHI_CUT && fabs(zed_A1-zed_A2)<PC1_DZ_CUT) continue;
          if (charge_A1 == charge_A2) continue;
          if (arm_A1 != arm_A2) continue; 
          if (APPLY_DZED_CUT[icut] && fabs(zed_A1-zed_A2)>DZED_CUT[icut]) continue;

          MyPair mypair_AA;
          if (!solution(&mytrk_A1, &mytrk_A2, &mypair_AA, &reco, zVtx_A, icut)) continue;
          
          TLorentzVector p_A1_DC;
          p_A1_DC.Clear();
          p_A1_DC.SetPx(momx_A1);
          p_A1_DC.SetPy(momy_A1);
          p_A1_DC.SetPz(momz_A1);
          p_A1_DC.SetE(sqrt(pow(p_A1_DC.P(),2) + Me2));

          TLorentzVector p_A2_DC;
          p_A2_DC.Clear();
          p_A2_DC.SetPx(momx_A2);
          p_A2_DC.SetPy(momy_A2);
          p_A2_DC.SetPz(momz_A2);
          p_A2_DC.SetE(sqrt(pow(p_A2_DC.P(),2) + Me2));
          
          TLorentzVector convPhoton_AA;
          convPhoton_AA.Clear();
          convPhoton_AA = p_A1_DC + p_A2_DC;
          float Pt_AA = convPhoton_AA.Pt();

          TLorentzVector real_convPhoton_AA;
          real_convPhoton_AA.Clear();
          real_convPhoton_AA.SetPx(momx_A1+momx_A2);
          real_convPhoton_AA.SetPy(momy_A1+momy_A2);
          real_convPhoton_AA.SetPz(momz_A1+momz_A2);
          real_convPhoton_AA.SetE(sqrt(pow(convPhoton_AA.P(),2)));
          float Pt_AA_real = real_convPhoton_AA.Pt();

          eemass[icut][icent][0]->Fill(convPhoton_AA.M(), Pt_AA, w_pt);
          eemass_unw[icut][icent][0]->Fill(convPhoton_AA.M(), Pt_AA);
          if((convPhoton_AA.M()<0.01)||(convPhoton_AA.M()>0.15)) continue;

          // combine with photons from same event after distance-based charged veto
          for (int i_A = 0; i_A < nclust_A; ++i_A){
            MyCluster myclust_A = event->GetClusterEntry(i_A);
            int id_cluster =  myclust_A.GetID();
            if((id_cluster==id_A1)||(id_cluster==id_A2)) continue;
            //if(myclust_A.GetMcId()!=1) continue;
            float xclust_A = myclust_A.GetX();
            float yclust_A = myclust_A.GetY();
            float zclust_A = myclust_A.GetZ()-zVtx_A; 
            float ecore_A = myclust_A.GetEcore() * PHOTON_ESCALE[icut];
            float chi2_A   = myclust_A.GetChi2();
            int sector = myclust_A.GetSect();
            //====================================================================
            //                    apply cluster cuts
            //====================================================================
            if(!apply_emcmap(sector,myclust_A.GetIY(),myclust_A.GetIZ())) continue;
            if(ecore_A < ECORE_CUT[icut]) continue;
            if(chi2_A > CHI2_CUT[icut]) continue;

            TLorentzVector emcPhoton_A;
            emcPhoton_A.Clear();
            emcPhoton_A.SetPx( ecore_A*xclust_A/sqrt(xclust_A*xclust_A+yclust_A*yclust_A+zclust_A*zclust_A) );
            emcPhoton_A.SetPy( ecore_A*yclust_A/sqrt(xclust_A*xclust_A+yclust_A*yclust_A+zclust_A*zclust_A) );
            emcPhoton_A.SetPz( ecore_A*zclust_A/sqrt(xclust_A*xclust_A+yclust_A*yclust_A+zclust_A*zclust_A) );
            emcPhoton_A.SetE( ecore_A );

            TLorentzVector pi0_AA;
            pi0_AA.Clear();
            pi0_AA = real_convPhoton_AA + emcPhoton_A;
            float Mass_AA = pi0_AA.M();
            float Pt_pi0 = pi0_AA.Pt();


            eegmass[icut][icent][0]->Fill(Mass_AA, Pt_AA, w_pt);
            eegmass_unw[icut][icent][0]->Fill(Mass_AA, Pt_AA);
          }
        }
      }
    }
  }

  for (int icut = 0; icut < cutbin; ++icut) {
    for (int icent = 0; icent < centbin; ++icent) {
      for (int isect = 0; isect < sectbin; ++isect) {
        int nbinsx_ee = eemass[icut][icent][isect]->GetNbinsX();
        int nbinsy_ee = eemass[icut][icent][isect]->GetNbinsY();

        for (int ibin_x = 1; ibin_x <= nbinsx_ee; ++ibin_x) {
          for (int ibin_y = 1; ibin_y <= nbinsy_ee; ++ibin_y) {
            double n_unw = eemass_unw[icut][icent][isect]->GetBinContent(ibin_x, ibin_y);
            if (n_unw == 0) continue;
            double content = eemass[icut][icent][isect]->GetBinContent(ibin_x, ibin_y);
            double err = content / sqrt(n_unw);
            eemass[icut][icent][isect]->SetBinError(ibin_x, ibin_y, err);
          }
        }

        int nbinsx_eeg = eegmass[icut][icent][isect]->GetNbinsX();
        int nbinsy_eeg = eegmass[icut][icent][isect]->GetNbinsY();

        for (int ibin_x = 1; ibin_x <= nbinsx_eeg; ++ibin_x) {
          for (int ibin_y = 1; ibin_y <= nbinsy_eeg; ++ibin_y) {
            double n_unw = eegmass_unw[icut][icent][isect]->GetBinContent(ibin_x, ibin_y);
            if (n_unw == 0) continue;
            double content = eegmass[icut][icent][isect]->GetBinContent(ibin_x, ibin_y);
            double err = content / sqrt(n_unw);
            eegmass[icut][icent][isect]->SetBinError(ibin_x, ibin_y, err);
          }
        }
      }
    }
  }

  output->cd();
  h1d_truepiz->Write();
  h1d_truepiz_unw->Write();
  for(int icent = 0; icent < centbin; icent++){
    for (int icut = 0; icut < cutbin; ++icut){
      for(int isect = 0; isect < sectbin; isect++){
        eemass[icut][icent][isect]->Write();
        eegmass[icut][icent][isect]->Write();
        eemass_unw[icut][icent][isect]->Write();
        eegmass_unw[icut][icent][isect]->Write();
      }
    }
  }
  input_A->Close();
  input_B->Close();
  output->Close();
} // End make_TTree
