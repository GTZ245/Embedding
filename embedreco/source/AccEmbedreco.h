#ifndef __ACCEMBEDRECO_H_
#define __ACCEMBEDRECO_H_


#include "SubsysReco.h"

class AccEvent;
class AccGeometry;
class AccHit;
class AccRaw;

class AccEmbedreco : public SubsysReco
{

 public:
  AccEmbedreco(const char * name = "AccEmbedreco", const int ver=1);
  virtual ~AccEmbedreco();

  int Init(PHCompositeNode* topNode);
  int InitRun(PHCompositeNode* topNode);
  int process_event(PHCompositeNode* topNode);
  int ResetEvent(PHCompositeNode* topNode);

 protected:
  int version;

  AccRaw* d_raw;
  AccHit* d_hit;
  AccEvent* d_acc;
  AccGeometry* d_geo;

  bool CreateNodeTree(PHCompositeNode* topNode);

};

#endif
