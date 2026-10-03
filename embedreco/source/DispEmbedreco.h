#ifndef __DISPEMBEDRECO_H__
#define __DISPEMBEDRECO_H__

#include "SubsysReco.h"

class PHCompositeNode;

class DispEmbedreco: public SubsysReco
{
 public:
  DispEmbedreco(const char *name = "EmbedDisplay");
  virtual ~DispEmbedreco();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int Reset(PHCompositeNode *topNode);
  int ResetEvent(PHCompositeNode *topNode);
  void Print(const char *what)const;
  int flag;
 protected:
};

#endif 
