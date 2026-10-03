#ifndef __SVXEMBEDEVAL_H__
#define __SVXEMBEDEVAL_H__

#include <SubsysReco.h>
#include <phool.h>
#include <TRandom3.h>

class PHCompositeNode;
class TFile;
class TNtuple;
class SvxCentralTrack;
class SvxClusterList;

class SvxEmbedEval : public SubsysReco
{

 public:
  SvxEmbedEval();
  SvxEmbedEval(std::string filename);
  virtual ~SvxEmbedEval() {}

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  
 protected:

  void ConversionVeto(float centrality,
                      float pt, 
                      float charge, 
                      bool electron, 
                      float mom,
                      SvxCentralTrack *svxcnt,
                      SvxClusterList* d_svxhit);

  PHCompositeNode *m_mcnode;
  PHCompositeNode *m_realnode;

  TFile* OutputNtupleFile;
  std::string OutputFileName;
  int init_ana;
  int EventNumber;

  TNtuple* ntpvtx;
  TNtuple* ntpeval;

    bool d_conversionVeto;
    bool d_conversionVeto1;

  int n0;
  int n1;
  int n2;
};

#endif 
