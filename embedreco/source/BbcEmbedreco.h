#ifndef __BBCEMBEDRECO_H__
#define __BBCEMBEDRECO_H__

#include "SubsysReco.h"

#include <string>

class PHCompositeNode;
class VtxMixer;

class BbcEmbedreco: public SubsysReco
{
 public:
  BbcEmbedreco(const std::string &name = "BBC");
  virtual ~BbcEmbedreco();

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int ResetEvent(PHCompositeNode *topNode);

 protected:
  VtxMixer*     vtxmixer;
};

#endif /* __BBCEMBEDRECO_H__ */
