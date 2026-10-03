#ifndef __EMBEDDISP_H__
#define __EMBEDDISP_H__

#include <SubsysReco.h>
#include <PHPoint.h>

#include <string>
#include <vector>
#include <map>

class PHCompositeNode;
class TFile;
class TNtuple;
class TH2F;
class TCanvas;
class TPad;
class TMarker;
//class TTree;
//class TH1F;
//class TH3F;
//class TH3D;

//class SvxRawhitClusterList;
//class SvxRawhitList;
//class SvxClusterList;
//class SvxCluster;

class PHCentralTrack;
class SvxCentralTrackList;
class SvxCentralTrack;

//class BbcOut;
//class VtxOut;
//class EventHeader;
//class SvxClusterContainer;
//class SvxClusterInfo;
//class RpSumXYObject;
//class McEvalSingleList;
//class SvxGhitList;
//class SvxGhit;
//class SvxGhitRawhitList;
//class SvxHitMap;
//class fkinWrapper;
//class svxAddress;
//class svxDetectorGeo;


////////////

class embeddisp: public SubsysReco {
public:
  embeddisp(std::string filename="embeddisp.root");
  virtual ~embeddisp();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  inline int Reset(PHCompositeNode *topNode) 		{return 0;}
  inline int ResetEvent(PHCompositeNode *topNode) 	{return 0;}
  int InitRun(PHCompositeNode *topNode);
  inline void Print(const char *what) const 		{return;}

  void  drawCanvas();
  void  draw();

private:
  void calcDCA_BCbyCircleProjection
            (
              float pt, float phi, float the,  // pt in xy-plane, phi, theta at inner most layer
              int charge,                      // charge of the track
              float hx, float hy, float hz,    // hit position at inner most layer
              float vx, float vy,              // beam center
              float fieldPolarity,             //
              float* dx, float* dy, float* dz, // dca position
              float* d2dca_bc                  // return
            );

  void remapCntSvx(PHCentralTrack* cnt, SvxCentralTrackList* svxcnt);

private:
  int EventNumber;

  std::string m_outFileName;
  TFile*      m_outfile;

  TNtuple*    m_ntp_embed;

  std::vector<SvxCentralTrack*> m_vsvxcnt;

  // canvas
  TCanvas *m_c1;
  TPad    *m_pad[3];
  TH2F    *m_frame1, *m_frame2, *m_frame3;

  std::vector<PHPoint>  m_vdchhit_real, m_vdchhit_mc, m_vdchhit;
  std::vector<PHPoint>  m_vpc1hit_real, m_vpc1hit_mc, m_vpc1hit;
  std::vector<PHPoint>  m_vpc2hit_real, m_vpc2hit_mc, m_vpc2hit;
  std::vector<PHPoint>  m_vpc3hit_real, m_vpc3hit_mc, m_vpc3hit;
  std::vector<PHPoint>  m_vemchit_real, m_vemchit_mc, m_vemchit;
  std::vector<PHPoint>  m_vcrkhit_real, m_vcrkhit_mc, m_vcrkhit;

  std::vector<TMarker*> m_vdrawmarker;

  PHPoint m_crkpmt[5120];
};

#endif

