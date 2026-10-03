#ifndef __CRKEMBEDRECO_H__
#define __CRKEMBEDRECO_H__

#include "SubsysReco.h"

class PHCompositeNode ;
class CrkMixer        ;


class CrkEmbedreco: public SubsysReco
{
 public:
  CrkEmbedreco(const char *name = "CRK");
  virtual ~CrkEmbedreco();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int Reset(PHCompositeNode *topNode) {return 0;}
  int ResetEvent(PHCompositeNode *topNode);
  void Print(const std::string&) const {}

 protected:
  int copyWrapper(PHCompositeNode *);

  CrkMixer        *crkmixer         ;
};

#endif /* __CRKEMBEDRECO_H__ */
