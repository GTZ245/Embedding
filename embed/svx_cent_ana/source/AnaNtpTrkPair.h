//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 12 16:08:11 2012 by ROOT version 5.30/03
// from TTree ntp_cnttrk/CNTtrack with VTX info tree
// found on file: /phenix/zdata01/phnxreco/VTX/hachiya/source/svx_cent_ana/pack/data_349000_345000/SvxCntQA_file20_0000349679-00.root
//////////////////////////////////////////////////////////

#ifndef AnaNtpTrkPair_h
#define AnaNtpTrkPair_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>

class eTrk;
class TH1F;

class AnaNtpTrkPair {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

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

   Float_t         simdca;
   Float_t         simpt;
   Float_t         simphi0;
   Float_t         simthe0;
   Float_t         simvx;
   Float_t         simvy;
   Float_t         simvz;
   Int_t           simpid;
   Int_t           simpidpa;
   Int_t           simpidpr;
   Float_t         simdcapa;
   Float_t         simptpa;
   Float_t         simphi0pa;
   Float_t         simthe0pa;
   Float_t         simvxpa;
   Float_t         simvypa;
   Float_t         simvzpa;
   Float_t         simdcapr;
   Float_t         simptpr;
   Float_t         simphi0pr;
   Float_t         simthe0pr;
   Float_t         simvxpr;
   Float_t         simvypr;
   Float_t         simvzpr;


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

   TBranch        *b_simdca;   //!
   TBranch        *b_simpt;   //!
   TBranch        *b_simphi0;   //!
   TBranch        *b_simthe0;   //!
   TBranch        *b_simvx;   //!
   TBranch        *b_simvy;   //!
   TBranch        *b_simvz;   //!
   TBranch        *b_simpid;   //!
   TBranch        *b_simpidpa;   //!
   TBranch        *b_simpidpr;   //!
   TBranch        *b_simdcapa;   //!
   TBranch        *b_simptpa;   //!
   TBranch        *b_simphi0pa;   //!
   TBranch        *b_simthe0pa;   //!
   TBranch        *b_simvxpa;   //!
   TBranch        *b_simvypa;   //!
   TBranch        *b_simvzpa;   //!
   TBranch        *b_simdcapr;   //!
   TBranch        *b_simptpr;   //!
   TBranch        *b_simphi0pr;   //!
   TBranch        *b_simthe0pr;   //!
   TBranch        *b_simvxpr;   //!
   TBranch        *b_simvypr;   //!
   TBranch        *b_simvzpr;   //!


   // 
   TFile *m_outfile;
   TTree *m_eepair;
   TTree *m_eepair_bg;
   TTree *m_etree;
   TH1F  *h_pt;

   // variables for pair
   // event
   float m_bbcz, m_bbcq;
   float m_xvtx, m_yvtx, m_zvtx;
   // pair
   float m_mee, m_px, m_py,m_pz, m_pt, m_thev, m_phiv, m_ptep, m_ptem;
   // ep
   float m_momp, m_phi0p, m_the0p, m_n0p, m_ch2npe0p, m_dispp, m_ecorep, m_epp;
   float m_chi2ndfp,  m_d2dcap, m_zdcap;
   int   m_nhitp;
   float m_dphip[4], m_dzp[4];
   int   m_simpaidp;
   float m_simvrp, m_simvzp;
   // em
   float m_momm, m_phi0m, m_the0m, m_n0m, m_ch2npe0m, m_dispm, m_ecorem, m_epm;
   float m_chi2ndfm, m_d2dcam, m_zdcam;
   int m_nhitm;
   float m_dphim[4], m_dzm[4];
   int   m_simpaidm;
   float m_simvrm, m_simvzm;

   // variable for single
   float m_mom, m_phi0, m_the0, m_c;
   int   m_dcqual;
   float m_emcdphi, m_emcdz;
   float m_n0, m_ch2npe0, m_disp, m_ecore, m_ep;
   float m_chi2ndf, m_d2dca, m_zdca;
   int   m_nhit;
   int   m_convtag;
   float m_dphi[4], m_dz[4];
   int   m_simpaid;
   float m_simvr, m_simvz;


   //////////
   std::vector<eTrk*> m_vep;
   std::vector<eTrk*> m_vem;

   int m_simmode;

   AnaNtpTrkPair(TTree *tree=0);
   virtual ~AnaNtpTrkPair();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

////////////
   void init_ana(const char *outname);

   void fill_eTrk(eTrk* etrk);
   void clear_VeTrk(std::vector<eTrk*>& vetrk);
   void analyze_event();
   void fill_pair(eTrk* ep, eTrk* em, TTree *eepair);
   void fill_pair(std::vector<eTrk*>& vep, std::vector<eTrk*>& vem, TTree *eepair); // include loop

   void fill_single(std::vector<eTrk*>& ve); // include loop
   void fill_single(eTrk* etrk);

   void end(){ if(m_outfile!=NULL) { m_outfile->Write();  m_outfile->Close(); } }

   void loop_a_file(const char* listfile="readCntTrk.list");

   void setSimMode(int mode){ m_simmode=mode; }

  private:
   void fill_trkvalue(eTrk *etrk, 
     float& mom, float& phi0, float& the0, 
     float& n0,  float& ch2npe0, float& disp, 
     float& ecore, float& ep,
     int&   nhit,
     float& chi2ndf, float& d2dca, float& zdca,
     float* dphi_ptr, float* dz_ptr,
     int& simpaid, float& simvr, float& simvz);

};

#endif

#ifdef AnaNtpTrkPair_cxx
AnaNtpTrkPair::AnaNtpTrkPair(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/phenix/zdata01/phnxreco/VTX/hachiya/source/svx_cent_ana/pack/data_349000_345000/SvxCntQA_file20_0000349679-00.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/phenix/zdata01/phnxreco/VTX/hachiya/source/svx_cent_ana/pack/data_349000_345000/SvxCntQA_file20_0000349679-00.root");
      }
      f->GetObject("ntp_cnttrk",tree);

   }
   Init(tree);

   m_simmode=0;
}

AnaNtpTrkPair::~AnaNtpTrkPair()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnaNtpTrkPair::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnaNtpTrkPair::LoadTree(Long64_t entry)
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

void AnaNtpTrkPair::Init(TTree *tree)
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

   if(m_simmode==1){
     fChain->SetBranchAddress("simdca", &simdca, &b_simdca);
     fChain->SetBranchAddress("simpt", &simpt, &b_simpt);
     fChain->SetBranchAddress("simphi0", &simphi0, &b_simphi0);
     fChain->SetBranchAddress("simthe0", &simthe0, &b_simthe0);
     fChain->SetBranchAddress("simvx", &simvx, &b_simvx);
     fChain->SetBranchAddress("simvy", &simvy, &b_simvy);
     fChain->SetBranchAddress("simvz", &simvz, &b_simvz);
     fChain->SetBranchAddress("simpid", &simpid, &b_simpid);
     fChain->SetBranchAddress("simpidpa", &simpidpa, &b_simpidpa);
     fChain->SetBranchAddress("simpidpr", &simpidpr, &b_simpidpr);
     fChain->SetBranchAddress("simdcapa", &simdcapa, &b_simdcapa);
     fChain->SetBranchAddress("simptpa", &simptpa, &b_simptpa);
     fChain->SetBranchAddress("simphi0pa", &simphi0pa, &b_simphi0pa);
     fChain->SetBranchAddress("simthe0pa", &simthe0pa, &b_simthe0pa);
     fChain->SetBranchAddress("simvxpa", &simvxpa, &b_simvxpa);
     fChain->SetBranchAddress("simvypa", &simvypa, &b_simvypa);
     fChain->SetBranchAddress("simvzpa", &simvzpa, &b_simvzpa);
     fChain->SetBranchAddress("simdcapr", &simdcapr, &b_simdcapr);
     fChain->SetBranchAddress("simptpr", &simptpr, &b_simptpr);
     fChain->SetBranchAddress("simphi0pr", &simphi0pr, &b_simphi0pr);
     fChain->SetBranchAddress("simthe0pr", &simthe0pr, &b_simthe0pr);
     fChain->SetBranchAddress("simvxpr", &simvxpr, &b_simvxpr);
     fChain->SetBranchAddress("simvypr", &simvypr, &b_simvypr);
     fChain->SetBranchAddress("simvzpr", &simvzpr, &b_simvzpr);
   }
   Notify();
}

Bool_t AnaNtpTrkPair::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnaNtpTrkPair::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnaNtpTrkPair::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AnaNtpTrkPair_cxx
