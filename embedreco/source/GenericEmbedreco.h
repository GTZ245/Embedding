#ifndef __GENERICSIMRECO_H__
#define __GENERICSIMRECO_H__

#include "SubsysReco.h"

class PHCompositeNode;

class GenericSimreco: public SubsysReco
{
 public:
  GenericSimreco(const char *name = "GENERIC");
  virtual ~GenericSimreco() {}

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int Reset(PHCompositeNode *topNode) {return 0;}
  int ResetEvent(PHCompositeNode *topNode);
  void Print(const char *what) const {return;}

 protected:
  int CreateNodeTree(PHCompositeNode *topNode);

};

#endif /* __GENERICSIMRECO_H__ */
