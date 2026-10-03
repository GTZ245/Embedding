#ifndef __ANACOMPACTCNT_H__
#define __ANACOMPACTCNT_H__

#include "SubsysReco.h"

class PHCompositeNode;
class TFile;
//class TNtuple;
class TTree;
//class TH1F;
//class TH2F;

//class BbcOut;
//class VtxOut;
class PHCentralTrack;
class PHSnglCentralTrack;
class McEvalSingleList;
class SvxCentralTrackRecalList;
class SvxCentralTrackRecal;

class eTrk;

#include <string>
#include <vector>
//#include <map>


////////////

class AnaCompactCNT: public SubsysReco {
public:

  AnaCompactCNT(std::string filename = "anacompactcnt.root");
  virtual ~AnaCompactCNT();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  inline int Reset(PHCompositeNode *topNode) 		{return 0;}
  inline int ResetEvent(PHCompositeNode *topNode) 	{return 0;}
  int InitRun(PHCompositeNode *topNode);
  inline void Print(const char *what) const 		{return;}

protected:

  int CreateNodeTree(PHCompositeNode *topNode) 	{return 0;}

  void initPairTree();

  void analyze_track(PHCentralTrack* cntlist, SvxCentralTrackRecalList *svxcntlist, McEvalSingleList* mctrk);
  void analyze_pair();

  void fill_eTrk(eTrk* etrk, PHSnglCentralTrack *trk, SvxCentralTrackRecal* svxcnt);
  void fill_eTrk_sim(eTrk* etrk, int itrk, McEvalSingleList* mctrklist);
  void fill_pair(eTrk* ep, eTrk* em, TTree *eepair);
  void fill_pair(std::vector<eTrk*>& ep, std::vector<eTrk*>& vem, TTree *eepair);

  void fill_single(std::vector<eTrk*>& ve); // include loop
  void fill_single(eTrk *etrk);
 

  bool isGoodTrack(PHSnglCentralTrack* trk);
  bool isElectron(int n0, float ch2npe0, float disp, float Ep, float prob);

  void clear_VeTrk(std::vector<eTrk*>& vetrk);

  void fill_trkvalue(eTrk *etrk,
    float& mom, float& phi0, float& the0,
    float& n0,  float& ch2npe0, float& disp,
    float& ecore, float& ep,
    int&   nhit,
    float& chi2ndf, float& d2dca, float& zdca,
    float* dphi_ptr, float* dz_ptr,
    int& simpaid, float& simvr, float& simvz);

  float calc_pair(eTrk *ep, eTrk *em,
                float& Mee,  float& px,   float& py,   float& pz, float& pt,
                float& thev, float& ptep, float& ptem, float& phiv);


  void init_idxbuf();
  void calc_iz_icent(float zvtx, float bbcq, int& iz, int& icent);



protected:
  std::vector<eTrk*> m_vep;
  std::vector<eTrk*> m_vem;

  int m_simmode;



protected:

  std::string m_OutputFileName;
  TFile*      m_OutputNtupleFile;

  int init_ana;
  int EventNumber;

  // 
  TTree *m_eepair;
  TTree *m_eepair_bg;
  TTree *m_etree;

  // variables for pair
  int   m_evt;
  float m_bbcq, m_xvtx, m_yvtx, m_zvtx;
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


};

#endif

