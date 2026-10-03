#ifndef __SVXANALYSIS_H__
#define __SVXANALYSIS_H__

#include "SubsysReco.h"

class PHCompositeNode;
class TFile;
class TNtuple;
class TTree;
class TH1F;
class TH2F;
class TH3F;
class TH3D;
class TRandom;

class SvxRawhitClusterList;
class SvxRawhitList;
class SvxClusterList;
class SvxCluster;

class BbcOut;
class VtxOut;
class EventHeader;
class SvxCentralTrackList;
class SvxClusterContainer;
class SvxClusterInfo;
class PHCentralTrack;
class RpSumXYObject;
class McEvalSingleList;
class SvxGhitList;
class SvxGhit;
class SvxGhitRawhitList;

class SvxHitMap;


class fkinWrapper;



class svxAddress;
class svxDetectorGeo;

class residual_data;
class svxana_gTrack;


#include <string>
#include <vector>
#include <map>


////////////
class svxana_gTrack {
  public:
    svxana_gTrack(int mctrack) : m_mctrack(mctrack){
      m_nghit       = 0;      
      m_nrecohit    = 0;   
      m_nrecohitcnt = 0;
      for(int i=0; i<4; i++) {
        m_nghit_layer[i] = 0;
        m_nrecohit_layer[i] = 0;
        m_nrecohitcnt_layer[i] = 0;
      }
    }

    void print();
    void addGhit(SvxGhit *ghit);

    void calcNhit();

  public:
    int m_mctrack;
    std::vector<int> m_vrecoid;

    std::vector<SvxGhit*> m_vghit;
    int m_nghit_layer[4];    // n_ghit in each layer
    int m_nrecohit_layer[4]; // n_cluster associated ghit in each layer
    int m_nrecohitcnt_layer[4]; // n_cluster associated both ghit and SvxCentralTrack in each layer

    int m_nghit;       // n_ghit 0-4. if nhit in each layer must be 1 or 0. then nhit:0-4
    int m_nrecohit;    // n_cluster 0-4
    int m_nrecohitcnt; // n_cluster associated both ghit and SvxCentralTrack in each layer
};

////////////

class svxana: public SubsysReco {
public:
//  enum {MAX=10};
  enum {MAX=3};

public:

  svxana();
  svxana(std::string filename);
  virtual ~svxana();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  inline int Reset(PHCompositeNode *topNode) 		{return 0;}
  inline int ResetEvent(PHCompositeNode *topNode) 	{return 0;}
  int InitRun(PHCompositeNode *topNode);
  inline void Print(const char *what) const 		{return;}

  void setSimMode(const int sim){ m_simmode=sim; }
  void setGoodOnlyFlag(bool flag) { m_goodOnlyFlag = flag; }

public: // relation between gtrack and cnt
  std::map<int, svxana_gTrack> m_vgTrack;      // <mctrack, svxana_gTrack object>
  std::vector<svxana_gTrack*>  m_vreco_gTrack; // reverse map, index is recoid. <svxana_gTrack pointer>
  std::vector<int>             m_vClusterGhit; // m_vClusterGhit[clusterID] = ghitID
  std::vector<int>             m_vGhitCluster; // m_vGhitCluster[ghitID] = clusterID

private:
  static const float GOODZVTX;
  static const float GOODEMCDPHI;
  static const float GOODEMCDZ;
  static const int   GOODDCQUALITY;

protected:

  int CreateNodeTree(PHCompositeNode *topNode) 	{return 0;}

  TRandom* random;
  TFile* OutputNtupleFile;
  std::string OutputFileName;
  int init_ana;
  int EventNumber;

  bool m_goodOnlyFlag; // true: store only good event/track; false: store all events/tracks

  TNtuple* ntpsvxkalfit;
  TNtuple* ntpsvxrawhit;
  TNtuple* ntpsvxcluster;
  TNtuple* ntpsvxsegment;
  TNtuple* ntprawclus;
  TNtuple *ntp_trk;
  TNtuple *ntp_cnt;
  TNtuple *ntp_cnt_clus;
  TNtuple *ntp_cnt_clus1;
  TNtuple *ntp_evt_pxlchip;

  TNtuple *m_ntp_dphidz;

//  TNtuple *ntp_cnttrk;
  TH1F* hdphi;
  TH1F* hdphi2;
  TH1F* hdphi3;
  TH1F* hdphi3a;
  TH1F* hdphi3b;
  TH1F* hdthe;
  TH1F* hdthe2;
  TH2F* hdphidthe;
  TH1F* hsvxsegmentmult;
  TH1F* hsvxclustermult;
  TH1F* hsvxrawhitmult;
  TH1F* hadc2X;
  TH1F* hpix_chip_mult;
  TH1F* hzvtx;
  TH2F* h2pixel;
  TH2F* h2pixel2; //nhit<30
  TH2F* h2pixel3; //nchip<5

  TTree *ntpevt;
  TTree *m_ntp_cnttrk;
  TTree *m_ntp_cnttrk_bg;
  TTree *m_ntpcls;
  TTree *m_ntpclsr;
  TTree *m_ntprp;
  TTree *m_ntprpsum;

  TH3F* hsizebbcqevt[4]; // B0,1,2,3

  TH3D* h3_clsmap;

  int m_fieldScale;

private:
  svxAddress*     m_svxAdrObj;
  svxDetectorGeo* m_svxGeo;

  void fillRawCluster(SvxRawhitClusterList *rawclus, SvxRawhitList* rawhit, SvxClusterList* clus);

  void fillPxlChip(float, float, SvxRawhitClusterList*, SvxRawhitList*, SvxClusterList*);
  void fillCentralTrack(TTree *ntp_trk, BbcOut*, VtxOut*, SvxCentralTrackList*, PHCentralTrack*, 
                        SvxClusterContainer *container, McEvalSingleList *mctrk, int &lowPtCtr);

  void fillCntEval(McEvalSingleList *mctrk, SvxGhitList *svxghit, fkinWrapper* fkinW, PHCentralTrack* cnt, SvxCentralTrackList* svxcntlist);

  void fillGhitCluster(SvxGhitRawhitList *ghitrawlist, SvxRawhitClusterList *rawcluslist, SvxClusterList *svxcluslist, SvxGhitList *ghitlist);

  void fillRpSum(int eseq, float bbcq, float cent, VtxOut *vtxout, RpSumXYObject*, SvxClusterList *clslist);

  void fillQvector(int evt, float bbcq, float cent, VtxOut *vtxout, SvxClusterList *svxlist);
  void fillCompactCluster(int evt, float bbcq, VtxOut *vtxout, SvxHitMap *svxlist);
  void fillRp(int evt, int eseq, float bbcq, float cent, VtxOut *vtxout, SvxClusterList *svxlist);

  void fillTest(PHCompositeNode *topNode);

  void filldphidz(PHCentralTrack *cntlist, SvxCentralTrackList* svxcntlist, SvxClusterList *svxlist, McEvalSingleList* mctrk);

  void fillBadPacket(EventHeader *evthdr);

  void initEvtTree();
  void initCntTree();
  void initClsTree();
  void initRpTree();
  void initRpSumTree();

  // copy from SvxPixel1v1
  // Return the local X position from the sensor Ix value
  double get_sensorXpos_pixel(const int layer, const int ladder, const int sens,
                              const int section, const int readout, int ix) const;

  // Return the local Z position from the sensor Iz value
  double get_sensorZpos_pixel(const int layer, const int ladder, const int sens,
                              const int section, const int readout, int ix) const;


  double get_ROC_pixel(const int layer, const int ladder, const int sens, const float localz) const;

  void calc_dphidz(float xvtx, float yvtx, float zvtx, 
                   float xsvx, float ysvx, float zsvx, 
                   float pt, float charge, float phi0, float the0, 
                   float* dproj, float* magbend, float* zproj
                  );

  void calc_sdphidz(float pt, int layer, float dp, float dz, float* sdp, float* sdz);

  int m_simmode;

  void calc_dca(float mom, float phi0, float the0, float c, float bdir, 
                float x, float y, float xvtx, float yvtx,
                float* dca, float* R, float *L, int v=0);
  void calc_dcafrom2hit(float mom, float phi0, float the0, float c, float bdir, 
                float x0, float y0, float x1, float y1, 
                float xvtx, float yvtx,
                float* dca, float* R, float *L, int v=0);

  void calc_chi2(float mom, int* vlayer, float* vdphi, float *vdz, float *r_dpchi2, float *r_dzchi2, int* r_dpndf, int* r_dzndf);
  void calc_schi2(float mom, int* vlayer, float* vdphi, float *vdz, float *r_dpchi2, float *r_dzchi2, int* r_dpndf, int* r_dzndf);
  void calc_sschi2(float mom, int* vlayer, float* vsdp, float *vsdz, float *r_dpchi2, float *r_dzchi2, int* r_dpndf, int* r_dzndf, float *vdfsdp);
  void calc_chi2p(float mom, float the0, int* vlayer, float* vdphi, float *vdz, float *r_dpchi2, float *r_dzchi2, int* r_dpndf, int* r_dzndf, float *vdfdpp, float *vdfdzz);
  int  getCorrelationIdx(int prevlay, int lay); // prevlay : outer, lay : inner
  void getCorrelationSlope(int prevlay, int lay, float *s_dp, float *s_dz); // prevlay : outer, lay : inner
  void getCorrelationSigma(float mom, int prevlay, int lay, float *sig_dp, float *sig_dz); // prevlay : outer, lay : inner

  int  getCorrelationIdxSub(int prevlay, int lay); // prevlay : outer, lay : inner, sublayer
  void getCorrelationSlopeSub(int prevlay, int lay, float *s_dp, float *s_dz); // prevlay : outer, lay : inner sublayer
  void getCorrelationSigmaSub(float mom, int prevlay, int lay, float *sig_dp, float *sig_dz); // prevlay : outer, lay : inner sublayer

  void getCorrelationSSigma(float mom, int prevlay, int lay, float *sig_dp, float *sig_dz); // prevlay : outer, lay : inner sublayer

  void getCorrelationSlopeP(int prevlay, int lay, float *s_dp, float *s_dz); // prevlay : outer lay : inner layer
  void getCorrelationSigmaP(float mom, float the0, int prevlay, int lay, float *sig_dp, float *sig_dz); // prevlay : outer, lay : inner

  void searchSecondHit(int layer, SvxClusterInfo *hit, SvxClusterContainer *container, 
                       float xvtx, float yvtx, float zvtx,
                       float mom, float phi0, float the0, float c, 
                       std::vector<residual_data>& vlist);

  int get_sublayer(int layer, int ladder);

  float getSigmaDphi(int layer, float mom);

  // simulation
  float calcSimD2DCA(float pt, float charge, float phi0, float hx, float hy, float vx, float vy);

  int findCloseSvxHits(int layer, float phihit, float zhit, int cid, 
                       SvxClusterList* svxlist, std::vector<SvxCluster*>& vhit);


 private:
   float m_beam_x, m_beam_y;
   int   m_lowpt_ctr[2];

 private:
  bool m_isBadPacket[2][60];


  int m_evt_run;
  int m_evt_event, m_evt_strig;
  int m_evt_eseq, m_evt_xctr;
  int m_evt_ptick[3];
  int m_evt_n0, m_evt_n1, m_evt_n2, m_evt_n3;
  float m_evt_zvtx, m_evt_zbbc, m_evt_zzdc;
  int m_evt_ntrk, m_evt_nseg, m_evt_nsvxtrk;
  int m_evt_nsegr, m_evt_nsvxtrkr;
  int m_evt_nraw, m_evt_nclsr;
  float m_evt_bbcq;
  float m_evt_xvtxs, m_evt_yvtxs, m_evt_zvtxs;
  float m_evt_xvtxp, m_evt_yvtxp, m_evt_zvtxp;
  float m_evt_zvtxsw, m_evt_zvtxse;
  float m_evt_xvtxpw, m_evt_yvtxpw, m_evt_zvtxpw;
  float m_evt_xvtxpe, m_evt_yvtxpe, m_evt_zvtxpe;

  ///////////////////
  int   m_cls_event;
  float m_cls_zvtx;
  float m_cls_bbcq;
  float m_cls_cent;
  float m_cls_xyz[3]; // not write
  float m_cls_lxz[2];
  float m_cls_qv[2]; // not write
  float m_cls_adc[2]; // not write
  int   m_cls_lay;
  int   m_cls_lad;
  int   m_cls_sens;
  int   m_cls_size; // not write
  int   m_cls_xzsize[2];

  ///////////////////
  int   m_rp_event[2];
  float m_rp_zvtx;
  float m_rp_bbcq;
  float m_rp_cent;
  float m_rp_qv[2][5][3]; // 2:nosize/size, 5:lay0,1,2,3,all, 2:1=y,0=x, nhit

  ///////////////////
  int   m_rpsum_evtseq;
  float m_rpsum_zvtx;
  float m_rpsum_bbcq;
  float m_rpsum_cent;
  float m_rpsum_qv[6][3]; //B0-3 Ball(all eta), B0 eta(gap2);
  float m_rpsum_bbcqv[3][3]; //Bbc North, South, All;
  int   m_rpsum_nhitlad[10];  //nhit @ b0 

  ///////////////////
  float* m_c_datap[56]; // temporary value array for tree filling
  float  m_c_bbcz, m_c_bbcq, m_c_t0;
  float  m_c_xvtx, m_c_yvtx, m_c_zvtx;
  float  m_c_zvtxp, m_c_zvtxs;
  float  m_c_mom, m_c_phi0, m_c_the0, m_c_zed, m_c_c;
  float  m_c_ephi0[4];
  int    m_c_dcqual;
  float  m_c_emcdphi, m_c_emcdz;
  float  m_c_pc3dphi, m_c_pc3dz;
  float  m_c_pc2dphi, m_c_pc2dz;
  float  m_c_n0, m_c_cch2, m_c_npe0, m_c_disp;
  float  m_c_sn0, m_c_scch2, m_c_snpe0, m_c_sdisp;
  float  m_c_ecore, m_c_ecorr;
  float  m_c_ly3, m_c_ld3, m_c_zproj3, m_c_dproj3, m_c_bend3, m_c_ph3, m_c_zv3, m_c_r3, m_c_sdp3, m_c_sdz3, m_c_bdp3, m_c_bdz3, m_c_fitdp3, m_c_fitdz3;
  float  m_c_ly2, m_c_ld2, m_c_zproj2, m_c_dproj2, m_c_bend2, m_c_ph2, m_c_zv2, m_c_r2, m_c_sdp2, m_c_sdz2, m_c_bdp2, m_c_bdz2, m_c_fitdp2, m_c_fitdz2;
  float  m_c_ly1, m_c_ld1, m_c_zproj1, m_c_dproj1, m_c_bend1, m_c_ph1, m_c_zv1, m_c_r1, m_c_sdp1, m_c_sdz1, m_c_bdp1, m_c_bdz1, m_c_fitdp1, m_c_fitdz1;
  float  m_c_ly0, m_c_ld0, m_c_zproj0, m_c_dproj0, m_c_bend0, m_c_ph0, m_c_zv0, m_c_r0, m_c_sdp0, m_c_sdz0, m_c_bdp0, m_c_bdz0, m_c_fitdp0, m_c_fitdz0;
  float  m_c_chi2;
  int    m_c_ndf, m_c_nhit;
  float  m_c_chi22;
  int    m_c_unique;
  int    m_c_hitptn;
  float  m_c_dpchi2, m_c_dzchi2;
  float  m_c_dpchi2p, m_c_dzchi2p;
  float  m_c_sdpchi2, m_c_sdzchi2;
  float  m_c_ssdpchi2, m_c_ssdzchi2;
  float  m_c_dfsdp[6]; // 3-2,3-1,3-0,2-1,2-0,1-0
  float  m_c_dfdpp[6]; // 3-2,3-1,3-0,2-1,2-0,1-0
  float  m_c_dfdzz[6]; // 3-2,3-1,3-0,2-1,2-0,1-0
  int    m_c_dpndf, m_c_dzndf;
  float  m_c_d2dca, m_c_d2dca0, m_c_d2dca1, m_c_R, m_c_L;
  float  m_c_d2dcab, m_c_d2dca2;
  float  m_c_bd2dca;
  float  m_c_zdca;
  float  m_c_dedx1, m_c_dedx2;
  int    m_c_exn3, m_c_exn2, m_c_exn1, m_c_exn0;
  float  m_c_exdz3[MAX], m_c_exdp3[MAX];
  float  m_c_exdz2[MAX], m_c_exdp2[MAX];
  float  m_c_exdz1[MAX], m_c_exdp1[MAX];
  float  m_c_exdz0[MAX], m_c_exdp0[MAX];
  float  m_c_simdca, m_c_simpt, m_c_simphi0, m_c_simthe0, m_c_simvx, m_c_simvy, m_c_simvz;
  int    m_c_simpid, m_c_simpidpa, m_c_simpidpr;
  float  m_c_simdcapa, m_c_simptpa, m_c_simphi0pa, m_c_simthe0pa, m_c_simvxpa, m_c_simvypa, m_c_simvzpa;
  float  m_c_simdcapr, m_c_simptpr, m_c_simphi0pr, m_c_simthe0pr, m_c_simvxpr, m_c_simvypr, m_c_simvzpr;
  int    m_c_simnghit, m_c_simnreco, m_c_simnrecocnt;

};

#endif

