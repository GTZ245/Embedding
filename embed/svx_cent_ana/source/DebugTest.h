#ifndef __DEBUGTEST_H__
#define __DEBUGTEST_H__

#include "SubsysReco.h"

class SvxClusterList;
class SvxCentralTrackList;

class DebugTest: public SubsysReco {

public:
  DebugTest();
  virtual ~DebugTest();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  inline int Reset(PHCompositeNode *topNode) 		{return 0;}
  inline int ResetEvent(PHCompositeNode *topNode) 	{return 0;}
  int InitRun(PHCompositeNode *topNode);
  inline void Print(const char *what) const 		{return;}


  bool testSelectedCluster(SvxClusterList      *svxsellist, 
                           SvxCentralTrackList *svxcntlist, 
                           SvxCentralTrackList *svxcntbglist);

};

#endif

