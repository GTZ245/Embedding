//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 14 13:20:42 2012 by ROOT version 5.30/03
// from TTree ntp_cnttrk/CNTtrack with VTX info tree
// found on file: /phenix/zdata05/phnxreco/run11auau_200GeV_pro90/run_0000349000_0000350000/SvxCntQA/SvxCntQA_run11auau_200GeV_pro90-0000349000-0000.root
//////////////////////////////////////////////////////////

#ifndef ReadNtpCntTrk_h
#define ReadNtpCntTrk_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>

class ReadNtpCntTrk {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   //////////////////////////////
   // event tree
   // Declaration of leaf types
   Int_t          m_e_run;
   Int_t          m_e_event;
   Int_t          m_e_strig;
   Int_t          m_e_eseq;
   Int_t          m_e_xctr;
   Int_t          m_e_n0;
   Int_t          m_e_n1;
   Int_t          m_e_n2;
   Int_t          m_e_n3;
   Float_t        m_e_zvtx;
   Float_t        m_e_zbbc;
   Float_t        m_e_zzdc;
   Int_t          m_e_nseg;
   Int_t          m_e_ntrk;
   Float_t        m_e_bbcq;
   Float_t        m_e_zdcq;
   Float_t        m_e_xvtxs;
   Float_t        m_e_yvtxs;
   Float_t        m_e_zvtxs;
   Float_t        m_e_zvtxsw;
   Float_t        m_e_zvtxse;
   Float_t        m_e_xvtxp;
   Float_t        m_e_yvtxp;
   Float_t        m_e_zvtxp;
   Float_t        m_e_xvtxpw;
   Float_t        m_e_yvtxpw;
   Float_t        m_e_zvtxpw;
   Float_t        m_e_xvtxpe;
   Float_t        m_e_yvtxpe;
   Float_t        m_e_zvtxpe;

// List of branches
   TBranch        *b_e_run;   //!
   TBranch        *b_e_event;   //!
   TBranch        *b_e_strig;   //!
   TBranch        *b_e_eseq;   //!
   TBranch        *b_e_xctr;   //!
   TBranch        *b_e_n0;   //!
   TBranch        *b_e_n1;   //!
   TBranch        *b_e_n2;   //!
   TBranch        *b_e_n3;   //!
   TBranch        *b_e_zvtx;   //!
   TBranch        *b_e_zbbc;   //!
   TBranch        *b_e_zzdc;   //!
   TBranch        *b_e_nseg;   //!
   TBranch        *b_e_trk;   //!
   TBranch        *b_e_bbcq;   //!
   TBranch        *b_e_zdcq;   //!
   TBranch        *b_e_xvtxs;   //!
   TBranch        *b_e_yvtxs;   //!
   TBranch        *b_e_zvtxs;   //!
   TBranch        *b_e_zvtxsw;   //!
   TBranch        *b_e_zvtxse;   //!
   TBranch        *b_e_xvtxp;   //!
   TBranch        *b_e_yvtxp;   //!
   TBranch        *b_e_zvtxp;   //!
   TBranch        *b_e_xvtxpw;   //!
   TBranch        *b_e_yvtxpw;   //!
   TBranch        *b_e_zvtxpw;   //!
   TBranch        *b_e_xvtxpe;   //!
   TBranch        *b_e_yvtxpe;   //!
   TBranch        *b_e_zvtxpe;   //!


   //////////////////////////////
   // track tree
   // Declaration of leaf types
   Float_t         bbcz;
   Float_t         bbcq;
   Float_t         t0;
   Float_t         xvtx;
   Float_t         yvtx;
   Float_t         zvtx;
   Float_t         zvtxp;
   Float_t         zvtxs;
   Float_t         mom;
   Float_t         the0;
   Float_t         phi0;
   Float_t         zed;
   Float_t         c;
   Int_t           dcqual;
   Float_t         emcdphi;
   Float_t         emcdz;
   Float_t         pc3dphi;
   Float_t         pc3dz;
   Float_t         pc2dphi;
   Float_t         pc2dz;
   Float_t         n0;
   Float_t         cch2;
   Float_t         npe0;
   Float_t         disp;
   Float_t         sn0;
   Float_t         scch2;
   Float_t         snpe0;
   Float_t         sdisp;
   Float_t         ecore;
   Float_t         ly3;
   Float_t         ld3;
   Float_t         zproj3;
   Float_t         dproj3;
   Float_t         bend3;
   Float_t         zv3;
   Float_t         ph3;
   Float_t         r3;
   Float_t         fitdp3;
   Float_t         fitdz3;
   Float_t         ly2;
   Float_t         ld2;
   Float_t         zproj2;
   Float_t         dproj2;
   Float_t         bend2;
   Float_t         zv2;
   Float_t         ph2;
   Float_t         r2;
   Float_t         fitdp2;
   Float_t         fitdz2;
   Float_t         ly1;
   Float_t         ld1;
   Float_t         zproj1;
   Float_t         dproj1;
   Float_t         bend1;
   Float_t         zv1;
   Float_t         ph1;
   Float_t         r1;
   Float_t         fitdp1;
   Float_t         fitdz1;
   Float_t         ly0;
   Float_t         ld0;
   Float_t         zproj0;
   Float_t         dproj0;
   Float_t         bend0;
   Float_t         zv0;
   Float_t         ph0;
   Float_t         r0;
   Float_t         fitdp0;
   Float_t         fitdz0;
   Float_t         chi2;
   Int_t           ndf;
   Int_t           nhit;
   Float_t         chi22;
   Int_t           unique;
   Float_t         dpchi2;
   Float_t         dzchi2;
   Int_t           dpndf;
   Int_t           dzndf;
   Float_t         dpchi2p;
   Float_t         dzchi2p;
   Float_t         d2dca;
   Float_t         d2dca0;
   Float_t         d2dcab;
   Float_t         zdca;
   Float_t         dedx1;
   Float_t         dedx2;
   Int_t           exn3;
   Float_t         exdz3[3];   //[exn3]
   Float_t         exdp3[3];   //[exn3]
   Int_t           exn2;
   Float_t         exdz2[3];   //[exn2]
   Float_t         exdp2[3];   //[exn2]
   Int_t           exn1;
   Float_t         exdz1[3];   //[exn1]
   Float_t         exdp1[3];   //[exn1]
   Int_t           exn0;
   Float_t         exdz0[3];   //[exn0]
   Float_t         exdp0[3];   //[exn0]

   // List of branches
   TBranch        *b_bbcz;   //!
   TBranch        *b_bbcq;   //!
   TBranch        *b_t0;   //!
   TBranch        *b_xvtx;   //!
   TBranch        *b_yvtx;   //!
   TBranch        *b_zvtx;   //!
   TBranch        *b_zvtxp;   //!
   TBranch        *b_zvtxs;   //!
   TBranch        *b_mom;   //!
   TBranch        *b_the0;   //!
   TBranch        *b_phi0;   //!
   TBranch        *b_zed;   //!
   TBranch        *b_c;   //!
   TBranch        *b_dcqual;   //!
   TBranch        *b_emcdphi;   //!
   TBranch        *b_emcdz;   //!
   TBranch        *b_pc3dphi;   //!
   TBranch        *b_pc3dz;   //!
   TBranch        *b_pc2dphi;   //!
   TBranch        *b_pc2dz;   //!
   TBranch        *b_n0;   //!
   TBranch        *b_cch2;   //!
   TBranch        *b_npe0;   //!
   TBranch        *b_disp;   //!
   TBranch        *b_sn0;   //!
   TBranch        *b_scch2;   //!
   TBranch        *b_snpe0;   //!
   TBranch        *b_sdisp;   //!
   TBranch        *b_ecore;   //!
   TBranch        *b_ly3;   //!
   TBranch        *b_ld3;   //!
   TBranch        *b_zproj3;   //!
   TBranch        *b_dproj3;   //!
   TBranch        *b_bend3;   //!
   TBranch        *b_zv3;   //!
   TBranch        *b_ph3;   //!
   TBranch        *b_r3;   //!
   TBranch        *b_fitdp3;   //!
   TBranch        *b_fitdz3;   //!
   TBranch        *b_ly2;   //!
   TBranch        *b_ld2;   //!
   TBranch        *b_zproj2;   //!
   TBranch        *b_dproj2;   //!
   TBranch        *b_bend2;   //!
   TBranch        *b_zv2;   //!
   TBranch        *b_ph2;   //!
   TBranch        *b_r2;   //!
   TBranch        *b_fitdp2;   //!
   TBranch        *b_fitdz2;   //!
   TBranch        *b_ly1;   //!
   TBranch        *b_ld1;   //!
   TBranch        *b_zproj1;   //!
   TBranch        *b_dproj1;   //!
   TBranch        *b_bend1;   //!
   TBranch        *b_zv1;   //!
   TBranch        *b_ph1;   //!
   TBranch        *b_r1;   //!
   TBranch        *b_fitdp1;   //!
   TBranch        *b_fitdz1;   //!
   TBranch        *b_ly0;   //!
   TBranch        *b_ld0;   //!
   TBranch        *b_zproj0;   //!
   TBranch        *b_dproj0;   //!
   TBranch        *b_bend0;   //!
   TBranch        *b_zv0;   //!
   TBranch        *b_ph0;   //!
   TBranch        *b_r0;   //!
   TBranch        *b_fitdp0;   //!
   TBranch        *b_fitdz0;   //!
   TBranch        *b_chi2;   //!
   TBranch        *b_ndf;   //!
   TBranch        *b_nhit;   //!
   TBranch        *b_chi22;   //!
   TBranch        *b_unique;   //!
   TBranch        *b_dpchi2;   //!
   TBranch        *b_dzchi2;   //!
   TBranch        *b_dpndf;   //!
   TBranch        *b_dzndf;   //!
   TBranch        *b_dpchi2p;   //!
   TBranch        *b_dpchi2z;   //!
   TBranch        *b_d2dca;   //!
   TBranch        *b_d2dca0;   //!
   TBranch        *b_d2dcab;   //!
   TBranch        *b_zdca;   //!
   TBranch        *b_dedx1;   //!
   TBranch        *b_dedx2;   //!
   TBranch        *b_exn3;   //!
   TBranch        *b_exdz3;   //!
   TBranch        *b_exdp3;   //!
   TBranch        *b_exn2;   //!
   TBranch        *b_exdz2;   //!
   TBranch        *b_exdp2;   //!
   TBranch        *b_exn1;   //!
   TBranch        *b_exdz1;   //!
   TBranch        *b_exdp1;   //!
   TBranch        *b_exn0;   //!
   TBranch        *b_exdz0;   //!
   TBranch        *b_exdp0;   //!

   TTree          *m_ntpevt;
   TTree          *m_ntp_cnttrk;
   TTree          *m_ntp_cnttrk_bg;
   TFile          *m_outfile;
   int             m_nfilled;

   ReadNtpCntTrk(TTree *tree=0);
   virtual ~ReadNtpCntTrk();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     InitEvt(TTree *tree);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual void     LoopEvt(TTree *tree);
   virtual void     LoopCntTrk(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void initOutput(const char* outname="readCntTrk.root");
   void end(){ if(m_outfile!=NULL) { m_outfile->Write();  m_outfile->Close(); } }

   void loop_a_file(const char* listfile="readCntTrk.list");
};

#endif

#ifdef ReadNtpCntTrk_cxx
ReadNtpCntTrk::ReadNtpCntTrk(TTree *tree) : fChain(NULL), m_ntp_cnttrk(NULL), m_ntp_cnttrk_bg(NULL), m_outfile(NULL), m_nfilled(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
/*
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/phenix/zdata05/phnxreco/run11auau_200GeV_pro90/run_0000349000_0000350000/SvxCntQA/SvxCntQA_run11auau_200GeV_pro90-0000349000-0000.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/phenix/zdata05/phnxreco/run11auau_200GeV_pro90/run_0000349000_0000350000/SvxCntQA/SvxCntQA_run11auau_200GeV_pro90-0000349000-0000.root");
      }
      f->GetObject("ntp_cnttrk",tree);

   }
   Init(tree);
*/
}

ReadNtpCntTrk::~ReadNtpCntTrk()
{
   //if(m_outtree!=NULL) delete m_outtree;
   if(m_outfile!=NULL) delete m_outfile;


   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ReadNtpCntTrk::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ReadNtpCntTrk::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ReadNtpCntTrk::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) {
     std::cout<<"No tree exist"<<std::endl;
     return;
   }
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("bbcz", &bbcz, &b_bbcz);
   fChain->SetBranchAddress("bbcq", &bbcq, &b_bbcq);
   fChain->SetBranchAddress("t0", &t0, &b_t0);
   fChain->SetBranchAddress("xvtx", &xvtx, &b_xvtx);
   fChain->SetBranchAddress("yvtx", &yvtx, &b_yvtx);
   fChain->SetBranchAddress("zvtx", &zvtx, &b_zvtx);
   fChain->SetBranchAddress("zvtxp", &zvtxp, &b_zvtxp);
   fChain->SetBranchAddress("zvtxs", &zvtxs, &b_zvtxs);
   fChain->SetBranchAddress("mom", &mom, &b_mom);
   fChain->SetBranchAddress("the0", &the0, &b_the0);
   fChain->SetBranchAddress("phi0", &phi0, &b_phi0);
   fChain->SetBranchAddress("zed", &zed, &b_zed);
   fChain->SetBranchAddress("c", &c, &b_c);
   fChain->SetBranchAddress("dcqual", &dcqual, &b_dcqual);
   fChain->SetBranchAddress("emcdphi", &emcdphi, &b_emcdphi);
   fChain->SetBranchAddress("emcdz", &emcdz, &b_emcdz);
   fChain->SetBranchAddress("pc3dphi", &pc3dphi, &b_pc3dphi);
   fChain->SetBranchAddress("pc3dz", &pc3dz, &b_pc3dz);
   fChain->SetBranchAddress("pc2dphi", &pc2dphi, &b_pc2dphi);
   fChain->SetBranchAddress("pc2dz", &pc2dz, &b_pc2dz);
   fChain->SetBranchAddress("n0", &n0, &b_n0);
   fChain->SetBranchAddress("cch2", &cch2, &b_cch2);
   fChain->SetBranchAddress("npe0", &npe0, &b_npe0);
   fChain->SetBranchAddress("disp", &disp, &b_disp);
   fChain->SetBranchAddress("sn0", &sn0, &b_sn0);
   fChain->SetBranchAddress("scch2", &scch2, &b_scch2);
   fChain->SetBranchAddress("snpe0", &snpe0, &b_snpe0);
   fChain->SetBranchAddress("sdisp", &sdisp, &b_sdisp);
   fChain->SetBranchAddress("ecore", &ecore, &b_ecore);
   fChain->SetBranchAddress("ly3", &ly3, &b_ly3);
   fChain->SetBranchAddress("ld3", &ld3, &b_ld3);
   fChain->SetBranchAddress("zproj3", &zproj3, &b_zproj3);
   fChain->SetBranchAddress("dproj3", &dproj3, &b_dproj3);
   fChain->SetBranchAddress("bend3", &bend3, &b_bend3);
   fChain->SetBranchAddress("zv3", &zv3, &b_zv3);
   fChain->SetBranchAddress("ph3", &ph3, &b_ph3);
   fChain->SetBranchAddress("r3", &r3, &b_r3);
   fChain->SetBranchAddress("fitdp3", &fitdp3, &b_fitdp3);
   fChain->SetBranchAddress("fitdz3", &fitdz3, &b_fitdz3);
   fChain->SetBranchAddress("ly2", &ly2, &b_ly2);
   fChain->SetBranchAddress("ld2", &ld2, &b_ld2);
   fChain->SetBranchAddress("zproj2", &zproj2, &b_zproj2);
   fChain->SetBranchAddress("dproj2", &dproj2, &b_dproj2);
   fChain->SetBranchAddress("bend2", &bend2, &b_bend2);
   fChain->SetBranchAddress("zv2", &zv2, &b_zv2);
   fChain->SetBranchAddress("ph2", &ph2, &b_ph2);
   fChain->SetBranchAddress("r2", &r2, &b_r2);
   fChain->SetBranchAddress("fitdp2", &fitdp2, &b_fitdp2);
   fChain->SetBranchAddress("fitdz2", &fitdz2, &b_fitdz2);
   fChain->SetBranchAddress("ly1", &ly1, &b_ly1);
   fChain->SetBranchAddress("ld1", &ld1, &b_ld1);
   fChain->SetBranchAddress("zproj1", &zproj1, &b_zproj1);
   fChain->SetBranchAddress("dproj1", &dproj1, &b_dproj1);
   fChain->SetBranchAddress("bend1", &bend1, &b_bend1);
   fChain->SetBranchAddress("zv1", &zv1, &b_zv1);
   fChain->SetBranchAddress("ph1", &ph1, &b_ph1);
   fChain->SetBranchAddress("r1", &r1, &b_r1);
   fChain->SetBranchAddress("fitdp1", &fitdp1, &b_fitdp1);
   fChain->SetBranchAddress("fitdz1", &fitdz1, &b_fitdz1);
   fChain->SetBranchAddress("ly0", &ly0, &b_ly0);
   fChain->SetBranchAddress("ld0", &ld0, &b_ld0);
   fChain->SetBranchAddress("zproj0", &zproj0, &b_zproj0);
   fChain->SetBranchAddress("dproj0", &dproj0, &b_dproj0);
   fChain->SetBranchAddress("bend0", &bend0, &b_bend0);
   fChain->SetBranchAddress("zv0", &zv0, &b_zv0);
   fChain->SetBranchAddress("ph0", &ph0, &b_ph0);
   fChain->SetBranchAddress("r0", &r0, &b_r0);
   fChain->SetBranchAddress("fitdp0", &fitdp0, &b_fitdp0);
   fChain->SetBranchAddress("fitdz0", &fitdz0, &b_fitdz0);
   fChain->SetBranchAddress("chi2", &chi2, &b_chi2);
   fChain->SetBranchAddress("ndf", &ndf, &b_ndf);
   fChain->SetBranchAddress("nhit", &nhit, &b_nhit);
   fChain->SetBranchAddress("chi22", &chi22, &b_chi22);
   fChain->SetBranchAddress("unique", &unique, &b_unique);
   fChain->SetBranchAddress("dpchi2", &dpchi2, &b_dpchi2);
   fChain->SetBranchAddress("dzchi2", &dzchi2, &b_dzchi2);
   fChain->SetBranchAddress("dpndf", &dpndf, &b_dpndf);
   fChain->SetBranchAddress("dzndf", &dzndf, &b_dzndf);
   fChain->SetBranchAddress("dpchi2p", &dpchi2p, &b_dpchi2p);
   fChain->SetBranchAddress("dzchi2p", &dzchi2p, &b_dpchi2z);
   fChain->SetBranchAddress("d2dca", &d2dca, &b_d2dca);
   fChain->SetBranchAddress("d2dca0", &d2dca0, &b_d2dca0);
   fChain->SetBranchAddress("d2dcab", &d2dcab, &b_d2dcab);
   fChain->SetBranchAddress("zdca", &zdca, &b_zdca);
   fChain->SetBranchAddress("dedx1", &dedx1, &b_dedx1);
   fChain->SetBranchAddress("dedx2", &dedx2, &b_dedx2);
   fChain->SetBranchAddress("exn3", &exn3, &b_exn3);
   fChain->SetBranchAddress("exdz3", &exdz3, &b_exdz3);
   fChain->SetBranchAddress("exdp3", &exdp3, &b_exdp3);
   fChain->SetBranchAddress("exn2", &exn2, &b_exn2);
   fChain->SetBranchAddress("exdz2", &exdz2, &b_exdz2);
   fChain->SetBranchAddress("exdp2", &exdp2, &b_exdp2);
   fChain->SetBranchAddress("exn1", &exn1, &b_exn1);
   fChain->SetBranchAddress("exdz1", exdz1, &b_exdz1);
   fChain->SetBranchAddress("exdp1", exdp1, &b_exdp1);
   fChain->SetBranchAddress("exn0", &exn0, &b_exn0);
   fChain->SetBranchAddress("exdz0", exdz0, &b_exdz0);
   fChain->SetBranchAddress("exdp0", exdp0, &b_exdp0);
   Notify();

   m_nfilled = 0;
}

void ReadNtpCntTrk::InitEvt(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run",    &m_e_run,    &b_e_run);
   fChain->SetBranchAddress("event",  &m_e_event,  &b_e_event);
   fChain->SetBranchAddress("strig",  &m_e_strig,  &b_e_strig);
   fChain->SetBranchAddress("eseq",   &m_e_eseq,   &b_e_eseq);
   fChain->SetBranchAddress("xctr",   &m_e_xctr,   &b_e_xctr);
   fChain->SetBranchAddress("n0",     &m_e_n0,     &b_e_n0);
   fChain->SetBranchAddress("n1",     &m_e_n1,     &b_e_n1);
   fChain->SetBranchAddress("n2",     &m_e_n2,     &b_e_n2);
   fChain->SetBranchAddress("n3",     &m_e_n3,     &b_e_n3);
   fChain->SetBranchAddress("zvtx",   &m_e_zvtx,   &b_e_zvtx);
   fChain->SetBranchAddress("zbbc",   &m_e_zbbc,   &b_e_zbbc);
   fChain->SetBranchAddress("zzdc",   &m_e_zzdc,   &b_e_zzdc);
   fChain->SetBranchAddress("nseg",   &m_e_nseg,   &b_e_nseg);
   fChain->SetBranchAddress("ntrk",   &m_e_ntrk,   &b_e_trk);
   fChain->SetBranchAddress("bbcq",   &m_e_bbcq,   &b_e_bbcq);
   fChain->SetBranchAddress("zdcq",   &m_e_zdcq,   &b_e_zdcq);
   fChain->SetBranchAddress("xvtxs",  &m_e_xvtxs,  &b_e_xvtxs);
   fChain->SetBranchAddress("yvtxs",  &m_e_yvtxs,  &b_e_yvtxs);
   fChain->SetBranchAddress("zvtxs",  &m_e_zvtxs,  &b_e_zvtxs);
   fChain->SetBranchAddress("zvtxsw", &m_e_zvtxsw, &b_e_zvtxsw);
   fChain->SetBranchAddress("zvtxse", &m_e_zvtxse, &b_e_zvtxse);
   fChain->SetBranchAddress("xvtxp",  &m_e_xvtxp,  &b_e_xvtxp);
   fChain->SetBranchAddress("yvtxp",  &m_e_yvtxp,  &b_e_yvtxp);
   fChain->SetBranchAddress("zvtxp",  &m_e_zvtxp,  &b_e_zvtxp);
   fChain->SetBranchAddress("xvtxpw", &m_e_xvtxpw, &b_e_xvtxpw);
   fChain->SetBranchAddress("yvtxpw", &m_e_yvtxpw, &b_e_yvtxpw);
   fChain->SetBranchAddress("zvtxpw", &m_e_zvtxpw, &b_e_zvtxpw);
   fChain->SetBranchAddress("xvtxpe", &m_e_xvtxpe, &b_e_xvtxpe);
   fChain->SetBranchAddress("yvtxpe", &m_e_yvtxpe, &b_e_yvtxpe);
   fChain->SetBranchAddress("zvtxpe", &m_e_zvtxpe, &b_e_zvtxpe);
   Notify();

   m_nfilled = 0;
}




Bool_t ReadNtpCntTrk::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ReadNtpCntTrk::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ReadNtpCntTrk::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   entry++;
   return 1;
}
#endif // #ifdef ReadNtpCntTrk_cxx
