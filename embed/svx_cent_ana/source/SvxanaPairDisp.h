#ifndef __SVXANAPAIRDISP_H__
#define __SVXANAPAIRDISP_H__

#include "SubsysReco.h"
#include <vector>

class PHCompositeNode;

class PHCentralTrack;
class SvxCentralTrackList;
class McEvalSingleList;

class TCanvas;
class TH2F;
class TArc;
class TArrow;
class TMarker;

class TrkCont;
class PairCont;

class SvxanaPairDisp : public SubsysReco {

public:

  SvxanaPairDisp();
  virtual ~SvxanaPairDisp();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  inline int Reset(PHCompositeNode *topNode) 		{return 0;}
  inline int ResetEvent(PHCompositeNode *topNode) 	{return 0;}
  int InitRun(PHCompositeNode *topNode);
  inline void Print(const char *what) const 		{return;}

  void initCanvas();
  void drawCanvas();
  void drawEvent(bool drawcross=true, bool drawsim=true);

protected:

  int CreateNodeTree(PHCompositeNode *topNode) 	{return 0;}

  int  m_init_ana;
  int  m_EventNumber;
  bool m_isInitCanvas;

private:
  void fillCntTrack(PHCentralTrack *trk);

  void fillCentralTrack(float bbcq, float zvtx, SvxCentralTrackList *svxcnttrklist, PHCentralTrack* cntlist, McEvalSingleList *mctrklist=NULL);

  void fillpair();

  TArc* makeLine( float cx, float cy, float R); 

  int calc_crosspos(TrkCont *ep, TrkCont *em,
                     int& ncross, 
                     std::vector<float>& v_crossR,
                     std::vector<float>& v_crossX, 
                     std::vector<float>& v_crossY);

  void calc_vector_at_cross(
    float px, float py, float pz,
    float dcax, float dcay, float dcaz,
    float charge, float R,
    float cross_x, float cross_y,
    float* cpx, float* cpy, float *cpz, float *cdphi
   );

  float calc_pair(float pxp, float pyp, float pzp, float Mep,
                  float pxm, float pym, float pzm, float Mem,
                float& Mee,  float& px,   float& py,   float& pz, float& pt,
                float& thev, float& ptep, float& ptem, float& phiv);



private:
  float m_vtx[3];
  std::vector<TrkCont*>  m_vTrkCont;
  std::vector<PairCont*>  m_vPairCont;

private: // canval
  TCanvas   *m_c1;
  TH2F      *m_frame;

  TMarker *m_mvtx;
  std::vector<TArc*>     m_vArc;
  std::vector<TArrow*>   m_vArrow;
  std::vector<TArrow*>   m_vArrowCross;
  std::vector<TArrow*>   m_vArrowMc;
  std::vector<TMarker*>  m_vMarker;
};

#endif

