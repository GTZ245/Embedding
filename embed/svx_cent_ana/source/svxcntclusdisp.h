#ifndef __SVXCNTCLUSDISP_H__
#define __SVXCNTCLUSDISP_H__

#include "SubsysReco.h"

class PHCompositeNode;

class PHCentralTrack;

class SvxClusterList;
class SvxCentralTrackList;
class SvxClusterContainer;

class svxAddress;
class svxDetectorGeo;

class TCanvas;
class TH2F;
class TPolyLine;
class TText;
class TMarker;
class TArc;
class TLine;

#include <string>
#include <vector>

class svxcntclusdisp : public SubsysReco {

public:

  svxcntclusdisp();
  svxcntclusdisp(std::string filename);
  virtual ~svxcntclusdisp();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  inline int Reset(PHCompositeNode *topNode) 		{return 0;}
  inline int ResetEvent(PHCompositeNode *topNode) 	{return 0;}
  int InitRun(PHCompositeNode *topNode);
  inline void Print(const char *what) const 		{return;}

  void initCanvas();
  void drawCanvas();
  void drawEvent();

  // 
  void cctestGenerateCluster();
  void cctestDrawPoint(float phi0, float dphi);
  void cctestDrawPointBlock(float phi0, float dphi, int block=-1);


protected:

  int CreateNodeTree(PHCompositeNode *topNode) 	{return 0;}

  int init_ana;
  int EventNumber;
  bool m_isInitCanvas;

private:
  svxAddress*     m_svxAdrObj;
  svxDetectorGeo* m_svxGeo;

  void fillCluster(SvxClusterList *clslist);
  void fillCntTrack(PHCentralTrack *trk);


  void fillCentralTrack(float bbcq, float zvtx, SvxCentralTrackList *svxcnttrklist, PHCentralTrack* cntlist);

  // use for clustercontainer test

  SvxClusterList      *m_cctestList;
  SvxClusterContainer *m_cctestCont;

private:
  TCanvas   *m_c1;
  TH2F      *m_frame;
  TPolyLine *m_pl[4][24];
  TText     *m_txt[4][24];

  std::vector<TMarker*>  m_vClusPos;
  std::vector<TArc*>     m_vTrkArc;
  std::vector<TLine*>    m_vTrkLine;
  std::vector<TLine*>    m_vRngLine;

  std::vector<TMarker*>  m_vccPos;
};

#endif

