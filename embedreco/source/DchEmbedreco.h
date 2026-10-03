#ifndef __DCHEMBEDRECO_H__
#define __DCHEMBEDRECO_H__

#include "SubsysReco.h"

class PHCompositeNode;
class PHNodeIOManager;
class dDchTracksWrapper;
class dDchTracksExtWrapper;
class mNewDchInitializer;
class mNewDchCandidatory;
class DchMixer;
//class mCentralTrackEvaluator_v1;
class DchEmbedreco: public SubsysReco
{
 public:
  DchEmbedreco(const std::string &name = "DCH");
  virtual ~DchEmbedreco();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int Reset(PHCompositeNode *topNode) {return 0;}
  int ResetEvent(PHCompositeNode *topNode);
  void Print(const std::string&) const {}

 protected:  
  // Pointers to modules...
  dDchTracksWrapper* dDchTrack;
  dDchTracksWrapper* dDchTrackPerf;
  dDchTracksExtWrapper* dDchTrackExt;
  dDchTracksExtWrapper* dDchTrackExtPerf;
  mNewDchInitializer* mDchInitializer;
  mNewDchCandidatory* mDchCandidatory;
  //embedding modules
  DchMixer* dchmixer;
};

#endif /* __DCHEMBEDRECO_H__ */
