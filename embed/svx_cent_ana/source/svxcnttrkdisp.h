#ifndef __SVXCNTTRKDISP_H__
#define __SVXCNTTRKDISP_H__

#include "SubsysReco.h"

class PHCompositeNode;

class PHCentralTrack;

class SvxClusterList;
class SvxCentralTrackList;
class SvxClusterContainer;

class svxAddress;
class svxDetectorGeo;
class SvxCentralTrackReco;

class TCanvas;
class TH2F;
class TPolyLine;
class TText;
class TMarker;
class TArc;
class TLine;

#include <string>
#include <vector>

class svxcnttrkdisp : public SubsysReco {
public:
  enum {MAXTRK=1000};

public:

  svxcnttrkdisp();
  svxcnttrkdisp(std::string filename);
  virtual ~svxcnttrkdisp();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  inline int Reset(PHCompositeNode *topNode) 		{return 0;}
  inline int ResetEvent(PHCompositeNode *topNode) 	{return 0;}
  int InitRun(PHCompositeNode *topNode);
  inline void Print(const char *what) const 		{return;}

  void setModule(SvxCentralTrackReco *module){ m_module = module;}

  void initCanvas();
  void drawCanvas();
  void drawEvent();
  void drawNextEvent() { drawEvent(); m_itrk++;}

  // 
/*
  void cctestGenerateCluster();
  void cctestDrawPoint(float phi0, float dphi);
  void cctestDrawPointBlock(float phi0, float dphi, int block=-1);
*/


protected:

  int CreateNodeTree(PHCompositeNode *topNode) 	{return 0;}

  int init_ana;
  int EventNumber;
  bool m_isInitCanvas;

private:
  svxAddress*     m_svxAdrObj;
  svxDetectorGeo* m_svxGeo;

  SvxCentralTrackReco *m_module;

  void fillCluster(SvxClusterList *clslist);
//  void fillCntTrack(PHCentralTrack *trk);


//  void fillCentralTrack(float bbcq, float zvtx, SvxCentralTrackList *svxcnttrklist, PHCentralTrack* cntlist);

  void fillModule();

  // use for clustercontainer test

/*
  SvxClusterList      *m_cctestList;
  SvxClusterContainer *m_cctestCont;
*/

private:
  int m_itrk;
  TCanvas   *m_c1;
  TH2F      *m_frame[4];
/*
  TPolyLine *m_pl[4][24];
  TText     *m_txt[4][24];
*/

  std::vector<TMarker*>  m_vClusPos[MAXTRK][4]; // ntrk/nlayer
  std::vector<TMarker*>  m_vClusPosBest[MAXTRK][4]; // ntrk/nlayer
/*
  std::vector<TArc*>     m_vTrkArc;
  std::vector<TLine*>    m_vTrkLine;
  std::vector<TLine*>    m_vRngLine;

  std::vector<TMarker*>  m_vccPos;
*/
};

#endif

