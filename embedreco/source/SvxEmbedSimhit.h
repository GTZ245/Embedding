#ifndef __SVXEMBEDSIMHIT_H__
#define __SVXEMBEDSIMHIT_H__

#include <SubsysReco.h>
#include <phool.h>
#include <TRandom3.h>

class svxDetectorGeo;

class SvxEmbedSimhit : public SubsysReco
{

 public:
  SvxEmbedSimhit(const std::string &name="SVXEMBEDSIMHIT");
  virtual ~SvxEmbedSimhit();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  
  void setRandomVertexShift(bool t, float wshift=0.2) { m_rndm=t; m_wshift=wshift; }
  void set_StripixelNoise(float a);
  void set_mc_zvtx_mean(float a) { mc_zvtx_mean = a; }

  void copyPHCentralTrack(bool flag) { m_copyPHCentralTrack = flag; }
  void copyMcSingle(      bool flag) { m_copyMcSingle       = flag; }

 protected:
  int CreateNodeTree(PHCompositeNode *topNode);

  void copyRunHeader(    PHCompositeNode *topNode);
  void copyEventHeader(  PHCompositeNode *topNode);
  void copyPreviousEvent(PHCompositeNode *topNode);
  void copyTrigLvl1(     PHCompositeNode* topNode);
  void copyTrigRunLvl1(  PHCompositeNode* topNode);
  void copyBbcOut(       PHCompositeNode* topNode);
  void copyPHGlobal(     PHCompositeNode* topNode);

 private:
  int embed_simhit(PHCompositeNode *topNode, float zshift);
  void copy_vtx(PHCompositeNode *topNode);
  
  PHCompositeNode *m_mcnode;
  PHCompositeNode *m_realnode;
  int m_ievt;

  float m_stripixel_sAQ;
  float m_stripixel_sNOISE;
  int m_stripixel_adcthre_zs;
  int m_stripixel_adcthre_rawhit;
  int m_stripixel_adcthre_rawhit_sum;

  float mc_zvtx_mean;
  int m_max_ladder[4];
  int m_max_sensor[4];
  double m_sensorZwhalf[4];

  svxDetectorGeo *m_svxgeo;

  bool     m_rndm;
  float    m_wshift;
  TRandom3 rndmx;
  TRandom3 rndmy;

  bool     m_copyPHCentralTrack; 
  bool     m_copyMcSingle; 

  bool     m_new_runhdr;
  bool     m_new_evthdr;
  bool     m_new_prev;
  bool     m_new_trglvl1;
  bool     m_new_trgrunlvl1;
  bool     m_new_bbcout;
  bool     m_new_phglobal;
};

#endif 
