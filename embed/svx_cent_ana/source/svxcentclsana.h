#ifndef __SVXCENTCLSANA_H__
#define __SVXCENTCLSANA_H__

#include "SubsysReco.h"
#include <string>


class TFile;
class TNtuple;
class TTree;
//class SvxCentralClusterReco;
class SvxCentralTrackReco;

class svxcentclsana: public SubsysReco {

public:

  svxcentclsana();
  svxcentclsana(std::string filename);
  virtual ~svxcentclsana();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  inline int Reset(PHCompositeNode *topNode) 		{return 0;}
  inline int ResetEvent(PHCompositeNode *topNode) 	{return 0;}
  int InitRun(PHCompositeNode *topNode);
  inline void Print(const char *what) const 		{return;}

  void setRecoModule(SvxCentralTrackReco *mod) { if(mod!=NULL)  m_mod=mod; }
  void setSimMode(bool mode) { m_simmode=mode; }

private:
  float calcSimD2DCA(float pt, float charge, float phi0, float hx, float hy, float vx, float vy);


protected:

  int CreateNodeTree(PHCompositeNode *topNode) 	{return 0;}

  TFile* m_outFile;
  std::string m_outFileName;
  int m_initana;
  int m_eventNumber;
  bool m_simmode;

  SvxCentralTrackReco *m_mod;

  float* m_datap[36]; // temporary value array for tree filling
  TTree *m_ntp_cntclus;
  TTree *m_ntp_cnttrk;

  float  m_zvtx, m_bbcq, m_t0;
  int    m_vtxid;
  float  m_xvtx, m_yvtx;
  float  m_mom, m_phi0, m_the0, m_zed;
  int    m_dcqual;
  float  m_emcdphi, m_emcdz;
  float  m_tofdphi, m_tofdz;
  float  m_tofwdphi, m_tofwdz;
  float  m_m2emc, m_m2tof, m_m2tofw;
  float  m_temc, m_ttof, m_ttofw;
  float  m_plemc, m_pltof, m_pltofw;
  float  m_n0, m_cch2, m_npe0, m_disp;
  float  m_sn0, m_scch2, m_snpe0, m_sdisp;
  float  m_ecore, m_ecorr;
  int    m_lid; // link id
  float  m_ly3, m_ld3, m_zproj3, m_dproj3, m_bend3, m_zv3, m_phiv3, m_l3, m_size3;
  float  m_ly2, m_ld2, m_zproj2, m_dproj2, m_bend2, m_zv2, m_phiv2, m_l2, m_size2;
  float  m_ly1, m_ld1, m_zproj1, m_dproj1, m_bend1, m_zv1, m_phiv1, m_l1, m_size1;
  float  m_ly0, m_ld0, m_zproj0, m_dproj0, m_bend0, m_zv0, m_phiv0, m_l0, m_size0;
  float  m_chi2;
  int    m_ndf, m_best, m_nhit;
  float  m_d2dca;
  float  m_dL, m_dbend;
  float  m_d2dca0, m_d2dca1, m_d2dca2;
  float  m_dedx1, m_dedx2;
  float  m_simdca, m_simpt, m_simphi0, m_simvx, m_simvy, m_simvz;

};

#endif

