#ifndef __SVXCENTRALTRACKRECODEBUG__
#define __SVXCENTRALTRACKRECODEBUG__

#include <SvxCentralTrackReco.h>

class SvxCentralTrackRecoDebug : public SvxCentralTrackReco {
  public:
    SvxCentralTrackRecoDebug();
    virtual ~SvxCentralTrackRecoDebug(){}

    virtual void testLinkCluster(int nclus, float **ary);

  private:

};

#endif
