#ifndef __SVXANAMINI_H__
#define __SVXANAMINI_H__

#include "SubsysReco.h"
#include <string>


class PHCompositeNode;
class TFile;
class TTree;

////////////

class svxanamini: public SubsysReco {
public:

  svxanamini();
  svxanamini(std::string filename);
  virtual ~svxanamini();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  inline int Reset(PHCompositeNode *topNode) 		{return 0;}
  inline int ResetEvent(PHCompositeNode *topNode) 	{return 0;}
  int InitRun(PHCompositeNode *topNode);
  inline void Print(const char *what) const 		{return;}

protected:

  int CreateNodeTree(PHCompositeNode *topNode) 	{return 0;}

  TFile* OutputNtupleFile;
  std::string OutputFileName;
  int init_ana;
  int EventNumber;

  TTree *ntpevt;

private:
  void initEvtTree();

  int m_evt_run;
  int m_evt_event, m_evt_strig;
  int m_evt_eseq;
  float m_evt_zvtx, m_evt_zbbc, m_evt_zzdc;
  float m_evt_bbcq;
  int   m_evt_bbcns, m_evt_bbcnn;
  float m_evt_bbcqs, m_evt_bbcqn;
  float m_evt_xvtxs, m_evt_yvtxs, m_evt_zvtxs;
  float m_evt_xvtxp, m_evt_yvtxp, m_evt_zvtxp;


};

#endif

