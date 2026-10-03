#ifndef __SVXRPQA_H__
#define __SVXRPQA_H__

#include "SubsysReco.h"

class PHCompositeNode;
class TFile;
class TH3F;
class TH2F;

class VtxOut;
class EventHeader;
class PHCentralTrack;
class RpSumXYObject;


#include <string>


////////////

class SvxRpQA: public SubsysReco {
public:
  SvxRpQA(std::string filename="svxrpqa.root");
  virtual ~SvxRpQA();

  virtual int  Init(PHCompositeNode *topNode);
  virtual int  process_event(PHCompositeNode *topNode);
  virtual int  End(PHCompositeNode *topNode);

  virtual int  Reset(PHCompositeNode *topNode) 		{return 0;}
  virtual int  ResetEvent(PHCompositeNode *topNode) 	{return 0;}
  virtual int  InitRun(PHCompositeNode *topNode);
  virtual void Print(const char *what) const 		{return;}

  virtual void enableTickCut(bool flag){ m_istickcut=flag;}

protected:
  virtual bool TickCut(PHCompositeNode *topNode);
  virtual bool EventSelection(PHCompositeNode *topNode, float bbcz);

protected:

  int CreateNodeTree(PHCompositeNode *topNode) 	{return 0;}

private:
  std::string OutputFileName;
  TFile*      OutputNtupleFile;

  int init_ana;
  int EventNumber;

  float m_zvcut;
  bool  m_istickcut;
  int   m_pticks[3];

private:
  void fillRpSum(int eseq, float cent, VtxOut *vtxout, RpSumXYObject* rpsum);
  void initRpSumTree();

 private:

  static const int NRPKIND=63;

  // 
  TH3F *h_qxqy_cent[NRPKIND]; // B0:0-11, B1:12-21, B2:22-29, B3:30-35, 
                         // Ball:36-47, B0-B3(no gap):48-51, Ball(nogap):52
                         // B0-B3,all(eta gap1): 53-57, B0-B3,all(eta gap2): 58-62, 
  TH3F *h_qxqy_zvtx[NRPKIND]; // B0:0-11, B1:12-21, B2:22-29, B3:30-35, 
                         // Ball:36-47, B0-B3(no gap):48-51, Ball(nogap):52
                         // B0-B3,all(eta gap1): 53-57, B0-B3,all(eta gap2): 58-62, 
  TH3F *h_qx_eseq[NRPKIND]; // B0:0-11, B1:12-21, B2:22-29, B3:30-35, 
  TH3F *h_qy_eseq[NRPKIND]; // B0:0-11, B1:12-21, B2:22-29, B3:30-35, 
                         // Ball:36-47, B0-B3(no gap):48-51, Ball(nogap):52
                         // B0-B3,all(eta gap1): 53-57, B0-B3,all(eta gap2): 58-62, 
                         // Ball:36-47, B0-B3(no gap):48-51, Ball(nogap):52
                         // B0-B3,all(eta gap1): 53-57, B0-B3,all(eta gap2): 58-62, 

  TH3F *hbbc_qxqy_cent[3];
  TH3F *hbbc_qxqy_zvtx[3];
  TH3F *hbbc_qx_eseq[3];
  TH3F *hbbc_qy_eseq[3];

};

#endif

