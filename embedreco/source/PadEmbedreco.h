#ifndef __PADEMBEDRECO_H__
#define __PADEMBEDRECO_H__

#include "SubsysReco.h"
#include "TFile.h"
#include "TNtuple.h"

class PHCompositeNode;

class padInclBad            ;
class padDetectorGeo      ;
class padEvtToRaw           ;
class PadRecModule          ;
class PadTrackEvaluate      ;
class pcghitWrapper         ;
class dPadGhitClusWrapper   ;
class PadMixer              ;
class PadEmbedreco: public SubsysReco
{
 public:
  PadEmbedreco(const char *name = "PAD");
  virtual ~PadEmbedreco();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int Reset(PHCompositeNode *topNode) {return 0;}
  int ResetEvent(PHCompositeNode *topNode);
  void Print(const std::string&) const {}

 protected:

  // Utility functions
  int CreateNodeTree(PHCompositeNode *topNode);
  int copyWrapper(PHCompositeNode *);
  void setPadSplitMode(short padMode);
  short getPadSplitMode();


  //  Pointers to data
  padInclBad             *PadInclBad        ;
  padDetectorGeo       *mPadDetGeo        ;
  padEvtToRaw            *PadEvtToRaw       ;
  PadRecModule           *Pc1Rec            ;
  PadRecModule           *Pc2Rec            ;
  PadRecModule           *Pc3Rec            ;
  PadMixer               *padmixer          ;
  short padSplitMode;  // default 0 means no splitting of clusters (Run2 v03 DSTs choice)

};

#endif /* __PADEMBEDRECO_H__ */
