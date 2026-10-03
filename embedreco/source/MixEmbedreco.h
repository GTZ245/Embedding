#ifndef __MIXEmbedRECO_H__
#define __MIXEmbedRECO_H__

#include "SubsysReco.h"

class PHCompositeNode;

class MixEmbedreco: public SubsysReco
{
 public:
  MixEmbedreco(const char *name = "MixEmbed");
  virtual ~MixEmbedreco();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int Reset(PHCompositeNode *topNode);
  int ResetEvent(PHCompositeNode *topNode);
  void Print(const std::string&) const {}

 protected:
};

#endif 
