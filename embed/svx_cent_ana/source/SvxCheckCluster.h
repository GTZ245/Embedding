#ifndef __SVXCHECKCLUSTER_H__
#define __SVXCHECKCLUSTER_H__

#include "SubsysReco.h"

class PHCompositeNode;

class SvxClusterList;
class SvxHitMap;

class TFile;
class TNtuple;


////////////

class SvxCheckCluster: public SubsysReco 
{
  public:
    SvxCheckCluster();
    virtual ~SvxCheckCluster();
  
    virtual int Init(PHCompositeNode *topNode);
    virtual int process_event(PHCompositeNode *topNode);
    virtual int End(PHCompositeNode *topNode);
  
    virtual int  Reset(PHCompositeNode *topNode)	{return 0;}
    virtual int  ResetEvent(PHCompositeNode *topNode) 	{return 0;}
    virtual int  InitRun(PHCompositeNode *topNode);
    virtual void Print(const char *what) const 		{return;}
  
  protected:
    int CreateNodeTree(PHCompositeNode *topNode) 	{return 0;}
  
  private:
    TFile   *m_outputFile;
    TNtuple *m_ntpclus;

    int m_eventNumber; 

  
};

#endif

