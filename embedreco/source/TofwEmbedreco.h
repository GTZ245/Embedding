#ifndef __TOFWEMBEDRECO_H__
#define __TOFWEMBEDRECO_H__

#include "SubsysReco.h"

class PHCompositeNode;
class TofwMixer;
class TofwGeometry;
class TofwHit;

class TofwEmbedreco: public SubsysReco
{
 public:
  TofwEmbedreco(const char *name = "TOFW");
  virtual ~TofwEmbedreco();
  
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int Reset(PHCompositeNode *topNode) {return 0;}
  int ResetEvent(PHCompositeNode *topNode);

 protected:
  TofwMixer           *tofwmixer;
  TofwGeometry *d_geom;
  int runNumber;

  

};

#endif /* __TOFWEMBEDRECO_H__ */
