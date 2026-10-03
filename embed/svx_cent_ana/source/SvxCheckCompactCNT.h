#ifndef __SVXCHECKCOMPACTCNT_H__
#define __SVXCHECKCOMPACTCNT_H__

#include "SubsysReco.h"

class PHCompositeNode;
class SvxCluster;
class SvxSegment;
class SvxCentralTrack;

#include <string>
#include <vector>

class SvxCheckCompactCNT: public SubsysReco {
  public:
    SvxCheckCompactCNT(const std::string& name="SVXCHECKCOMPACTCNT");
    virtual ~SvxCheckCompactCNT();
  
    int Init(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);
  
    inline int Reset(PHCompositeNode *topNode) 		{return 0;}
    inline int ResetEvent(PHCompositeNode *topNode) 	{return 0;}
    int InitRun(PHCompositeNode *topNode);
    inline void Print(const char *what) const 		{return;}

    void setPrint(bool print) { m_print = print; }

  protected:
  
    int CreateNodeTree(PHCompositeNode *topNode);

    bool compareClusterList(PHCompositeNode *topNode);
    bool compareSegmentList(PHCompositeNode *topNode);
    bool compareCentralTrackList(PHCompositeNode *topNode);

    bool compareCluster(SvxCluster* cls_org, SvxCluster* cls_new);
    bool compareSegment(SvxSegment* seg_org, SvxSegment* seg_new);
    bool compareCentralTrack(SvxCentralTrack* svxcnt_org, SvxCentralTrack* svxcnt_new);

    bool compareInt(int orgval, int newval, const char *err);
    bool compareFloat(float orgval, float newval, float judge, const char *err);
  
    int init_ana;
    int EventNumber;

    bool m_print;
  
  private:
  
};

#endif

