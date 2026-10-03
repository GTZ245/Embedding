#ifndef __EVAEMBEDRECO_H__
#define __EVAEMBEDRECO_H__

#include "SubsysReco.h"

class PHCompositeNode;
class PHEmbedMCEvaluator;
class PHEmbededEvent;

class EvaEmbedreco: public SubsysReco
{
 public:
  EvaEmbedreco(const char *name = "Eva");
  virtual ~EvaEmbedreco();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int EndRun(const int runnumber);
  int End(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int Reset(PHCompositeNode *topNode);
  int ResetEvent(PHCompositeNode *topNode);
  void Print(const std::string&) const {}

 protected:
  PHEmbedMCEvaluator* evaluator;
  PHEmbededEvent* event;
};

#endif
