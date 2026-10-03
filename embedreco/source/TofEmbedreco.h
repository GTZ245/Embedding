#ifndef __TOFEMBEDRECO_H__
#define __TOFEMBEDRECO_H__

#include <SubsysReco.h>

class PHCompositeNode;

class TofAddressObject   ;
class TofCalibObject     ;
class TofGeometryObject  ;
class TofMixer           ;

class TofEmbedreco: public SubsysReco
{
 public:
  TofEmbedreco(const std::string &name = "TOF");
  virtual ~TofEmbedreco();

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int ResetEvent(PHCompositeNode *topNode);
  void Print(const std::string&) const {}

 protected:
  int copyWrapper(PHCompositeNode *);

  TofAddressObject   * TofAddress    ;
  TofCalibObject     * TofCalib      ;
  TofGeometryObject  * TofGeometry   ;
  TofMixer           * tofmixer      ;
};

#endif /* __TOFEMBEDRECO_H__ */

