#ifndef __EMCEMBEDRECO_H__
#define __EMCEMBEDRECO_H__

#include "SubsysReco.h"

class PHCompositeNode;
class mEmcClusterizerv0;
class EmcMixer;

class EmcEmbedreco: public SubsysReco
{
 public:
  EmcEmbedreco(const char *name = "EmbedEmc");
  virtual ~EmcEmbedreco();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int Reset(PHCompositeNode *topNode);
  int ResetEvent(PHCompositeNode *topNode);
  void Print(const std::string&) const {}

 protected:
  mEmcClusterizerv0* clusterizer;
  EmcMixer*  emcmixer;
};

#endif /* __EMCEMBEDRECO_H__ */
