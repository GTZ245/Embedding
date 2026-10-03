#define ReadNtpCntTrk_cxx
#include "ReadNtpCntTrk.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

void ReadNtpCntTrk::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L ReadNtpCntTrk.C
//      Root > ReadNtpCntTrk t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if(jentry%10000==0) cout<<jentry<<endl;

      float pt = mom*sin(the0);
      
      if( (pt>2&&ecore>0.5)       || // high pt track
          (n0>2&&cch2/npe0<10)    || // electron
          (sn0>2&&scch2/snpe0<10)    // electron bg
        )
        {
          m_ntp_cnttrk->Fill();
          m_nfilled++;
        }

      if(jentry==0) cout<<jentry<<" "<<mom<<" "<<n0<<" "<<endl;
   }
}

void ReadNtpCntTrk::LoopCntTrk(TTree *tree)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if(jentry%10000==0) cout<<jentry<<endl;

      float pt = mom*sin(the0);
      
      if( (pt>2&&ecore>0.5)       || // high pt track
          (n0>2&&cch2/npe0<10)    || // electron
          (sn0>2&&scch2/snpe0<10)    // electron bg
        )
        {
          tree->Fill();
          m_nfilled++;
        }

      if(jentry==0) cout<<jentry<<" "<<mom<<" "<<n0<<" "<<endl;
   }
}

void ReadNtpCntTrk::LoopEvt(TTree *tree)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if(jentry%10000==0) cout<<jentry<<endl;

      tree->Fill();
      m_nfilled++;
   }
}


void ReadNtpCntTrk::initOutput(const char* outname){
  m_outfile = TFile::Open(outname, "recreate");

/*
  m_ntpevt = new TTree("ntpevt", "Event Info Tree");
  m_ntpevt->Branch("run",     &m_e_run,     "run/I");
  m_ntpevt->Branch("event",   &m_e_event,   "event/I");
  m_ntpevt->Branch("strig",   &m_e_strig,   "strig/I");
  m_ntpevt->Branch("eseq",    &m_e_eseq,    "eseq/I");
  m_ntpevt->Branch("xctr",    &m_e_xctr,    "xctr/I");
  m_ntpevt->Branch("n0",      &m_e_n0,      "n0/I");
  m_ntpevt->Branch("n1",      &m_e_n1,      "n1/I");
  m_ntpevt->Branch("n2",      &m_e_n2,      "n2/I");
  m_ntpevt->Branch("n3",      &m_e_n3,      "n3/I");
  m_ntpevt->Branch("zvtx",    &m_e_zvtx,    "zvtx/F");
  m_ntpevt->Branch("zbbc",    &m_e_zbbc,    "zbbc/F");
  m_ntpevt->Branch("zzdc",    &m_e_zzdc,    "zzdc/F");
  m_ntpevt->Branch("nseg",    &m_e_nseg,    "nseg/I");
  m_ntpevt->Branch("ntrk",    &m_e_ntrk,    "trk/I");
  m_ntpevt->Branch("bbcq",    &m_e_bbcq,    "bbcq/F");
  m_ntpevt->Branch("zdcq",    &m_e_zdcq,    "zdcq/F");
  m_ntpevt->Branch("xvtxs",   &m_e_xvtxs,   "xvtxs/F");
  m_ntpevt->Branch("yvtxs",   &m_e_yvtxs,   "yvtxs/F");
  m_ntpevt->Branch("zvtxs",   &m_e_zvtxs,   "zvtxs/F");
  m_ntpevt->Branch("zvtxsw",  &m_e_zvtxsw,  "zvtxsw/F");
  m_ntpevt->Branch("zvtxse",  &m_e_zvtxse,  "zvtxse/F");
  m_ntpevt->Branch("xvtxp",   &m_e_xvtxp,   "xvtxp/F");
  m_ntpevt->Branch("yvtxp",   &m_e_yvtxp,   "yvtxp/F");
  m_ntpevt->Branch("zvtxp",   &m_e_zvtxp,   "zvtxp/F");
  m_ntpevt->Branch("xvtxpw",  &m_e_xvtxpw,   "xvtxpw/F");
  m_ntpevt->Branch("yvtxpw",  &m_e_yvtxpw,   "yvtxpw/F");
  m_ntpevt->Branch("zvtxpw",  &m_e_zvtxpw,   "zvtxpw/F");
  m_ntpevt->Branch("xvtxpe",  &m_e_xvtxpe,   "xvtxpe/F");
  m_ntpevt->Branch("yvtxpe",  &m_e_yvtxpe,   "yvtxpe/F");
  m_ntpevt->Branch("zvtxpe",  &m_e_zvtxpe,   "zvtxpe/F");
*/



  m_ntp_cnttrk    = new TTree("ntp_cnttrk",  "CNTtrack with VTX info tree");
  m_ntp_cnttrk_bg = new TTree("ntp_cnttrk_bg",  "CNTtrack with VTX info tree (BG)");


  TTree *t[2] = {m_ntp_cnttrk, m_ntp_cnttrk_bg};
  for(int i=0; i<2; i++){
    t[i]->Branch("bbcz",    &bbcz,    "bbcz/F");
    t[i]->Branch("bbcq",    &bbcq,    "bbcq/F");
    t[i]->Branch("t0",      &t0,      "t0/F");
    t[i]->Branch("xvtx",    &xvtx,    "xvtx/F");
    t[i]->Branch("yvtx",    &yvtx,    "yvtx/F");
    t[i]->Branch("zvtx",    &zvtx,    "zvtx/F");
    t[i]->Branch("zvtxp",   &zvtxp,   "zvtxp/F");
    t[i]->Branch("zvtxs",   &zvtxs,   "zvtxs/F");
    t[i]->Branch("mom",     &mom,     "mom/F");
    t[i]->Branch("the0",    &the0,    "the0/F");
    t[i]->Branch("phi0",    &phi0,    "phi0/F");
    t[i]->Branch("zed",     &zed,     "zed/F");
    t[i]->Branch("c",       &c,       "c/F");
    t[i]->Branch("dcqual",  &dcqual,  "dcqual/I");
    t[i]->Branch("emcdphi", &emcdphi, "emcdphi/F");
    t[i]->Branch("emcdz",   &emcdz,   "emcdz/F");
    t[i]->Branch("pc3dphi", &pc3dphi, "pc3dphi/F");
    t[i]->Branch("pc3dz",   &pc3dz,   "pc3dz/F");
    t[i]->Branch("pc2dphi", &pc2dphi, "pc2dphi/F");
    t[i]->Branch("pc2dz",   &pc2dz,   "pc2dz/F");
    t[i]->Branch("n0",      &n0,      "n0/F");
    t[i]->Branch("cch2",    &cch2,    "cch2/F");
    t[i]->Branch("npe0",    &npe0,    "npe0/F");
    t[i]->Branch("disp",    &disp,    "disp/F");
    t[i]->Branch("sn0",     &sn0,     "sn0/F");
    t[i]->Branch("scch2",   &scch2,   "scch2/F");
    t[i]->Branch("snpe0",   &snpe0,   "snpe0/F");
    t[i]->Branch("sdisp",   &sdisp,   "sdisp/F");
    t[i]->Branch("ecore",   &ecore,   "ecore/F");

    t[i]->Branch("ly3",     &ly3,      "ly3/F");
    t[i]->Branch("ld3",     &ld3,      "ld3/F");
    t[i]->Branch("zproj3",  &zproj3,   "zproj3/F");
    t[i]->Branch("dproj3",  &dproj3,   "dproj3/F");
    t[i]->Branch("bend3",   &bend3,    "bend3/F");
    t[i]->Branch("zv3",     &zv3,      "zv3/F");
    t[i]->Branch("ph3",     &ph3,      "ph3/F");
    t[i]->Branch("r3",      &r3,       "r3/F");
    t[i]->Branch("fitdp3",  &fitdp3,   "fitdp3/F");
    t[i]->Branch("fitdz3",  &fitdz3,   "fitdz3/F");
    t[i]->Branch("ly2",     &ly2,      "ly2/F");
    t[i]->Branch("ld2",     &ld2,      "ld2/F");
    t[i]->Branch("zproj2",  &zproj2,   "zproj2/F");
    t[i]->Branch("dproj2",  &dproj2,   "dproj2/F");
    t[i]->Branch("bend2",   &bend2,    "bend2/F");
    t[i]->Branch("zv2",     &zv2,      "zv2/F");
    t[i]->Branch("ph2",     &ph2,      "ph2/F");
    t[i]->Branch("r2",      &r2,       "r2/F");
    t[i]->Branch("fitdp2",  &fitdp2,   "fitdp2/F");
    t[i]->Branch("fitdz2",  &fitdz2,   "fitdz2/F");
    t[i]->Branch("ly1",     &ly1,      "ly1/F");
    t[i]->Branch("ld1",     &ld1,      "ld1/F");
    t[i]->Branch("zproj1",  &zproj1,   "zproj1/F");
    t[i]->Branch("dproj1",  &dproj1,   "dproj1/F");
    t[i]->Branch("bend1",   &bend1,    "bend1/F");
    t[i]->Branch("zv1",     &zv1,      "zv1/F");
    t[i]->Branch("ph1",     &ph1,      "ph1/F");
    t[i]->Branch("r1",      &r1,       "r1/F");
    t[i]->Branch("fitdp1",  &fitdp1,   "fitdp1/F");
    t[i]->Branch("fitdz1",  &fitdz1,   "fitdz1/F");
    t[i]->Branch("ly0",     &ly0,      "ly0/F");
    t[i]->Branch("ld0",     &ld0,      "ld0/F");
    t[i]->Branch("zproj0",  &zproj0,   "zproj0/F");
    t[i]->Branch("dproj0",  &dproj0,   "dproj0/F");
    t[i]->Branch("bend0",   &bend0,    "bend0/F");
    t[i]->Branch("zv0",     &zv0,      "zv0/F");
    t[i]->Branch("ph0",     &ph0,      "ph0/F");
    t[i]->Branch("r0",      &r0,       "r0/F");
    t[i]->Branch("fitdp0",   &fitdp0,    "fitdp0/F");
    t[i]->Branch("fitdz0",   &fitdz0,    "fitdz0/F");
    t[i]->Branch("chi2",     &chi2,      "chi2/F");
    t[i]->Branch("ndf",      &ndf,       "ndf/I");
    t[i]->Branch("nhit",     &nhit,      "nhit/I");
    t[i]->Branch("chi22",    &chi22,     "chi22/F");
    t[i]->Branch("unique",   &unique,    "unique/I");
    t[i]->Branch("dpchi2",   &dpchi2,    "dpchi2/F");
    t[i]->Branch("dzchi2",   &dzchi2,    "dzchi2/F");
    t[i]->Branch("dpndf",    &dpndf,     "dpndf/I");
    t[i]->Branch("dzndf",    &dzndf,     "dzndf/I");
    t[i]->Branch("dpchi2p",  &dpchi2p,   "dpchi2p/F");
    t[i]->Branch("dzchi2p",  &dzchi2p,   "dpchi2z/F");
    t[i]->Branch("d2dca",    &d2dca,     "d2dca/F");
    t[i]->Branch("d2dca0",   &d2dca0,    "d2dca0/F");
    t[i]->Branch("d2dcab",   &d2dcab,    "d2dcab/F");
    t[i]->Branch("zdca",     &zdca,      "zdca/F");
    t[i]->Branch("dedx1",    &dedx1,     "dedx1/F");
    t[i]->Branch("dedx2",    &dedx2,     "dedx2/F");
    t[i]->Branch("exn3",     &exn3,      "exn3/I");
    t[i]->Branch("exdz3",     exdz3,     "exdz3[exn3]/F");
    t[i]->Branch("exdp3",     exdp3,     "exdp3[exn3]/F");
    t[i]->Branch("exn2",     &exn2,      "exn2/I");
    t[i]->Branch("exdz2",     exdz2,     "exdz2[exn2]/F");
    t[i]->Branch("exdp2",     exdp2,     "exdp2[exn2]/F");
    t[i]->Branch("exn1",    &exn1,       "exn1/I");
    t[i]->Branch("exdz1",    exdz1,      "exdz1[exn1]/F");
    t[i]->Branch("exdp1",    exdp1,      "exdp1[exn1]/F");
    t[i]->Branch("exn0",    &exn0,       "exn0/I");
    t[i]->Branch("exdz0",    exdz0,      "exdz0[exn0]/F");
    t[i]->Branch("exdp0",    exdp0,      "exdp0[exn0]/F");
  }
}


void ReadNtpCntTrk::loop_a_file(const char *listfile){
  ifstream fin(listfile);

  if(fin.bad()){
    cout<<"listfile is bad: "<<listfile<<endl;
  }

  string fname;
  while( getline(fin, fname) ){
    cout<<fname.c_str()<<" "<<fname.size()<<endl;
    if(fname.size()==0) {
      //cout<<"File size is zero : "<<fname.c_str()<<endl;
      continue;
    }

    TFile *fntp = TFile::Open(fname.c_str());
    if(fntp==NULL){
      cout<<"File does not exist : "<<fname.c_str()<<endl;
      continue;
    }

/*
    // evt
    TTree *ntpevt = (TTree*)fntp->Get("ntpevt");
    cout<<"start evt"<<endl;
    InitEvt(ntpevt);
    LoopEvt(m_ntpevt);
    cout<<"Nfilled : "<<m_nfilled<<endl;
*/

    // cnttrk
    TTree *ntp_cnttrk = (TTree*)fntp->Get("ntp_cnttrk");
    cout<<"start cnttrk"<<endl;
    Init(ntp_cnttrk);
    LoopCntTrk(m_ntp_cnttrk);
    cout<<"Nfilled : "<<m_nfilled<<endl;

    // cnttrk_bg
    TTree *ntp_cnttrk_bg = (TTree*)fntp->Get("ntp_cnttrk_bg");
    cout<<"start cnttrk_bg : "<<ntp_cnttrk_bg<<endl;
    Init(ntp_cnttrk_bg);
    LoopCntTrk(m_ntp_cnttrk_bg);
    cout<<"Nfilled : "<<m_nfilled<<endl;

    fntp->Close();
    delete fntp;
    fChain = NULL;
  }
  fin.close();

  cout<<"end of loop_a_file"<<endl;
}
