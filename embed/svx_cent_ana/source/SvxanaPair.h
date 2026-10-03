#ifndef __SVXANAPAIR_H__
#define __SVXANAPAIR_H__

#include "SubsysReco.h"

class PHCompositeNode;
class TFile;
//class TNtuple;
class TTree;
//class TH1F;
//class TH2F;

//class BbcOut;
//class VtxOut;
class PHPoint;
class PHCentralTrack;
class PHSnglCentralTrack;
class McEvalSingleList;
class SvxCentralTrackRecalList;
class SvxCentralTrackRecal;

class AnaTrk;

#include <string>
#include <vector>
//#include <map>


////////////

class SvxanaPair: public SubsysReco {
public:

  SvxanaPair(std::string filename = "svxanapair.root");
  virtual ~SvxanaPair();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  inline int Reset(PHCompositeNode *topNode) 		{return 0;}
  inline int ResetEvent(PHCompositeNode *topNode) 	{return 0;}
  int InitRun(PHCompositeNode *topNode);
  inline void Print(const char *what) const 		{return;}

  void setSimMode(int flag) { m_simmode = flag;}
  void setPPMode(int flag)  { m_ppmode  = flag;}

protected:

  int CreateNodeTree(PHCompositeNode *topNode) 	{return 0;}

  void initPairTree();

  void analyze_track(PHCentralTrack* cntlist, SvxCentralTrackRecalList *svxcntlist, McEvalSingleList* mctrk);
  void analyze_pair();

  void fill_AnaTrk(AnaTrk* etrk, PHSnglCentralTrack *trk, SvxCentralTrackRecal* svxcnt);
  void fill_AnaTrk_sim(AnaTrk* etrk, int itrk, McEvalSingleList* mctrklist);
  void fill_pair(AnaTrk* ep, AnaTrk* em, float Mep, float Mem, TTree *eepair);
  void fill_pair(std::vector<AnaTrk*>& ep, std::vector<AnaTrk*>& vem, float Mep, float Mem, TTree *eepair);

  void fill_single(std::vector<AnaTrk*>& ve); // include loop
  void fill_single(AnaTrk *etrk);
 

  bool isGoodTrack(PHSnglCentralTrack* trk);
  bool isElectron(PHSnglCentralTrack* trk);
  bool isElectron(int n0, float ch2npe0, float disp, float Ep, float prob);

  void clear_VAnaTrk(std::vector<AnaTrk*>& vetrk);

  void fill_trkvalue(AnaTrk *etrk,
    float& mom, float& phi0, float& the0,
    float& n0,  float& ch2npe0, float& disp,
    float& ecore, float& ep,
    float& emcdp, float& emcdz,
    int&   nhit,
    float& chi2ndf, float& d2dca, float& zdca,
    float* dphi_ptr, float* dz_ptr,
    int& simpaid, float& simvr, float& simvz, float& simptpr);

  float calc_pair(float pxp, float pyp, float pzp, float Mep,
                  float pxm, float pym, float pzm, float Mem,
                float& Mee,  float& px,   float& py,   float& pz, float& pt,
                float& thev, float& ptep, float& ptem, float& phiv);

  float calc_pair(AnaTrk *ep, AnaTrk *em, float Mep, float Mem,
                float& Mee,  float& px,   float& py,   float& pz, float& pt,
                float& thev, float& ptep, float& ptem, float& phiv);


  int  calc_crosspos(AnaTrk* ep, AnaTrk* em, 
                     int& ncross, std::vector<float>& v_crossR,
                     std::vector<float>& v_crossX,
                     std::vector<float>& v_crossY);


  void calc_vector_at_cross(
          float px, float py, float pz, 
          float dcax, float dcay, float dcaz,
          float charge, float R,
          float cross_x, float cross_y,
          float* cpx, float* cpy, float *cpz, float *cdphi
        );

  void calc_dca_with_xvector(
        float xpx, float xpy, float xpz, float xpt,
        float crossx, float crossy,
        float xvtx, float yvtx, float zvtx,
        float* xdca2d, float *xlength, float* lxy
      );



  void init_idxbuf();
  void calc_iz_icent(float zvtx, float bbcq, int& iz, int& icent);


  void getPrimVertex(float* xvtx, float* yvtx, float* zvtx){
    *xvtx = m_xvtx; *yvtx = m_yvtx; *zvtx = m_zvtx;
  }

protected:
  std::vector<AnaTrk*> m_vep;
  std::vector<AnaTrk*> m_vem;

  std::vector<AnaTrk*> m_vhp;
  std::vector<AnaTrk*> m_vhm;

  int m_simmode;
  int m_ppmode;

  float m_fieldScale;



protected:

  std::string m_OutputFileName;
  TFile*      m_OutputNtupleFile;

  int init_ana;
  int EventNumber;

  // 
  TTree *m_eepair;
  TTree *m_eepair_bg;
  TTree *m_etree;

  TTree *m_pipipair;
  TTree *m_pikpair;
  TTree *m_kkpair;

  // variables for pair
  int   m_evt;
  float m_bbcq, m_xvtx, m_yvtx, m_zvtx;
  // pair
  float m_mee, m_px, m_py, m_pz, m_pt, m_thev, m_phiv, m_ptep, m_ptem;
  float m_svxmee, m_svxphiv;
  int   m_ncross;
  float m_crossr[2];
  float m_xmee, m_xpx, m_xpy, m_xpz, m_xpt, m_xdca2d, m_xlength, m_lxy;
  float m_simmee, m_simxpt, m_simxpz,  m_simxdca2d, m_simxlength, m_simlxy;
  // ep
  float m_momp, m_phi0p, m_the0p, m_n0p, m_ch2npe0p, m_dispp, m_ecorep, m_epp, m_emcdpp, m_emcdzp;
  float m_chi2ndfp,  m_d2dcap, m_zdcap;
  int   m_nhitp;
  float m_dphip[4], m_dzp[4];
  int   m_simpaidp;
  float m_simvrp, m_simvzp;
  float m_simptprp;
  // em
  float m_momm, m_phi0m, m_the0m, m_n0m, m_ch2npe0m, m_dispm, m_ecorem, m_epm, m_emcdpm, m_emcdzm;
  float m_chi2ndfm, m_d2dcam, m_zdcam;
  int m_nhitm;
  float m_dphim[4], m_dzm[4];
  int   m_simpaidm;
  float m_simvrm, m_simvzm;
  float m_simptprm;


  // variable for single
  float m_mom, m_phi0, m_the0, m_c;
  int   m_dcqual;
  float m_emcdphi, m_emcdz;
  float m_n0, m_ch2npe0, m_disp, m_ecore, m_ep;
  float m_chi2ndf, m_d2dca, m_zdca;
  int   m_nhit;
  int   m_convtag;
  float m_dphi[4], m_dz[4];
  float m_svxmom, m_svxphi0, m_svxthe0, m_svxc;
  int   m_simpaid;
  float m_simpx, m_simpy, m_simpz;
  float m_simvr, m_simvz;
  float m_simptpr;


};

#endif

