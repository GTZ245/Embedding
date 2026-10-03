#ifndef __TECEMBEDRECO_H__
#define __TECEMBEDRECO_H__

#include "SubsysReco.h"

class PHCompositeNode;

class mTecHoughTrackModule  ;
class TecAddressObject      ;
class TecGeometryObject     ;
class TecCalibrationObject  ;
class TecMixer              ;


class TecEmbedreco: public SubsysReco
{
 public:
  TecEmbedreco(const char *name = "TEC");
  virtual ~TecEmbedreco();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int Reset(PHCompositeNode *topNode) {return 0;}
  int ResetEvent(PHCompositeNode *topNode);
  void Print(const std::string&) const {}

 protected:
  mTecHoughTrackModule  * mTecHoughTrack ;
  
  TecAddressObject      * TecAddress     ;
  TecGeometryObject     * TecGeometry    ;
  TecCalibrationObject  * TecCalibration ;
  TecMixer              * tecmixer       ;
};

#endif /* __TECEMBEDRECO_H__ */
