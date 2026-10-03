#ifndef __DSTMERGER_H__
#define __DSTMERGER_H__

#include "SubsysReco.h"

#include "Riostream.h"

class cglDetectorGeo;
class cglHitAssociate;
class cglTransformDST;
class mPHDchTrackModel;
class mPHLineTrack;
class TriggerHelper;

class PHCompositeNode;
class PHGlobal;

class DstMerger: public SubsysReco {
  
 public:
  
  DstMerger();
  virtual ~DstMerger() {}
  
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  
  int Reset(PHCompositeNode *topNode) 		{return 0;}
  int ResetEvent(PHCompositeNode *topNode) 	{return 0;}
  void Print(const char *what) const 		{return;}
  void setVtxAsciiFileName(const char *filename);

 protected:
  std::string vtxfn;
  ofstream vtxo;
  int runnumber;
  int EventNumber;

};

#endif

