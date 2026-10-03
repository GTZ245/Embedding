#ifndef __SVXSTRIPDISP_H__
#define __SVXSTRIPDISP_H__

#include "SubsysReco.h"

class PHCompositeNode;
class TH1F;
class TH2F;
class TGraph;
class TCanvas;

class SvxRawhitClusterList;
class SvxRawhitList;
class SvxClusterList;

class svxAddress;
class svxDetectorGeo;
class SvxPixelHotDeadMap;

#include <string>

class svxstripdisp: public SubsysReco {

public:

  svxstripdisp();
  virtual ~svxstripdisp();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  inline int Reset(PHCompositeNode *topNode) 		{return 0;}
  inline int ResetEvent(PHCompositeNode *topNode) 	{return 0;}
  int InitRun(PHCompositeNode *topNode);
  inline void Print(const char *what) const 		{return;}

  virtual bool fill();
  virtual void draw(int layer, int ladder, bool isclus=true);
  virtual void drawone(int layer, int ladder, int sensor=0, bool isclus=true);

  virtual void setSizeCut(const int layer, const int val);

protected:

  int CreateNodeTree(PHCompositeNode *topNode) 	{return 0;}

  int init_ana;
  int EventNumber;

  int m_clssizecut[4];

private:
  TCanvas *m_c1;
  TCanvas *m_c2;
  TH2F    *m_h2_sensor[4][24][6]; // layer, ladder, sensor
  TGraph  *m_gcls_sensor[4][24][6]; // layer, ladder, sensor
  TGraph  *m_gcls_sensor_good[4][24][6]; // layer, ladder, sensor

  SvxRawhitList        *m_svxrawhit;
  SvxClusterList       *m_svx;
  SvxRawhitClusterList *m_svxrawhitclus;


private:
  svxAddress*          m_svxAdrObj;
  svxDetectorGeo*      m_svxGeo;
  SvxPixelHotDeadMap*  m_pixelhotdead;

  // copy from SvxPixel1v1
  // Return the local X position from the sensor Ix value
  double get_sensorXpos_pixel(const int layer, const int ladder, const int sens,
                              const int section, const int readout, int ix) const;

  // Return the local Z position from the sensor Iz value
  double get_sensorZpos_pixel(const int layer, const int ladder, const int sens,
                              const int section, const int readout, int ix) const;

};

#endif

