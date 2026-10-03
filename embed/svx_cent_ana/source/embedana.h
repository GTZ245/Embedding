#ifndef __EMBEDANA_H__
#define __EMBEDANA_H__

#include <SubsysReco.h>
#include <string>
#include <vector>
#include <map>

#include "Run14AuAuLeptonCombyConstants.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <stdlib.h>
#include <cassert>
#include <stdio.h>

#include "TMath.h"
#include "phool.h"
#include "PHTypedNodeIterator.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "Fun4AllReturnCodes.h"
#include "Fun4AllServer.h"
#include "recoConsts.h"

#include "PHCentralTrack.h"
#include "PHSnglCentralTrack.h"
#include "EventHeader.h"
#include "VtxOut.h"
#include "PHGlobal.h"
#include "RunHeader.h"
#include "TrigLvl1.h"

#include <DchHitLineTable.hh>
#include <PadCluster.h>
#include <CrkHit.h>
#include <PHTrackOut.h>
#include <PHDchTrackOut.h>
#include <CglTrack.h>
#include "SvxClusterList.h"
#include "SvxGhitClusterList.h"
#include "SvxCluster.h"
#include "SvxGhitCluster.h"
#include <PHPointerList.h>
#include <PHEmbedMcRecoTrack.hh>

#include "TFile.h"
#include "TNtuple.h"

#include "getClass.h"

#include "EmcCluster.h"
#include "emcClusterContainer.h"
#include "emcClusterContent.h"
#include "emcCalibrationDataHelper.h"
#include "EmcIndexer.h"
#include "emcTowerContainer.h"
#include "emcTowerContent.h"
#include "emcCalibrationData.h"
#include "emcCalibrationDataHelper.h"
#include "emcCalFEM.h"

using namespace std;
using namespace findNode;

#include "MyEvent.h"
#include "TOAD.h"
#include "TMath.h"

class PHCompositeNode;
class TFile;
class TNtuple;

class PHCentralTrack;

struct InData {
    double px, py, pz, vx, vy, vz;
    int id, nn;
};

////////////

class embedana: public SubsysReco {
public:
  embedana(std::string filename="embedana.root", 
           std::string filepath="/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/Dileptons/real/work/output/vertexes.txt",
           std::string oscarpath = "/phenix/plhf/mitran/Simul/Dileptons/output_single/single/00002.oscar.particles.dat",
           const int InID = 1);
  virtual ~embedana();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  int ResetEvent(PHCompositeNode *topNode);

  inline int Reset(PHCompositeNode *topNode) 		{return 0;}
  int InitRun(PHCompositeNode *topNode);
  inline void Print(const char *what) const 		{return;}

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

  //void remapCntSvx(PHCentralTrack* cnt, SvxCentralTrackList* svxcnt);
  void InitParams();
  float get_board(float phi_dc);
  void init_rungroup();
  int get_rungroup(int run_num);
  void init_dead_area();
  bool left_bound(float x, float y, float xx1, float yy1, float xx2, float yy2);
  bool right_bound(float x, float y, float xx1, float yy1, float xx2, float yy2);
  bool dead_region(float x, float y, float xx1, float yy1, float xx2, float yy2, float xx3, float yy3, float xx4, float yy4);
  bool IsCentralSupportCut(const float theta0, const float bbcVertex);
  int applySingleTrackCut(const PHCentralTrack *d_trk, const int itrk, const float vertex, const float centrality, const int run_number);
    

private:
  int EventNumber, InPartNumber;
  int InIDread;

  std::string m_outFileName;
  int fill_TTree, remove_cluster_hits;
  std::string local_filepath;
  std::string local_oscarpath;
  std::vector<double> vertexes;

  int ReadOrigPartMoms();
  InData InData_read[100000];

  int  dcmap_runs[N_RUN_GRP][MAX_RUN];
  float dcmap_xx1[N_RUN_GRP][N_SIDE][N_ARM][MAX_DEAD_AREA];
  float dcmap_yy1[N_RUN_GRP][N_SIDE][N_ARM][MAX_DEAD_AREA];
  float dcmap_xx2[N_RUN_GRP][N_SIDE][N_ARM][MAX_DEAD_AREA];
  float dcmap_yy2[N_RUN_GRP][N_SIDE][N_ARM][MAX_DEAD_AREA];
  float dcmap_xx3[N_RUN_GRP][N_SIDE][N_ARM][MAX_DEAD_AREA];
  float dcmap_yy3[N_RUN_GRP][N_SIDE][N_ARM][MAX_DEAD_AREA];
  float dcmap_xx4[N_RUN_GRP][N_SIDE][N_ARM][MAX_DEAD_AREA];
  float dcmap_yy4[N_RUN_GRP][N_SIDE][N_ARM][MAX_DEAD_AREA];

  int TRIGGERBIT;
  float BBC_VERTEX_CUT;

  // single cut
  int QUALITY[3];
  float Z_GLOBAL; // <|cut value|
  int DC_DEADMAP; // 0--w/o dc dead map, 1--w/ dc dead map

  // eid cut
  float E_PT;
  float MAX_PT;
  float N0;        // >cut value
  float DISP;      // <cut value
  float CHI2_NPE0; //< cut value

  float EOVERP;
  float DEP[2];
  float PROB; // >cut value
  float TEMC;
  float EMCDPHI;
  float EMCDZ;
  float EMCSDPHI;
  float EMCSDZ;
  float RICH_GHOST; // RICH ghost cut



  DileptonAnalysis::MyEventContainer *event_container;

  //std::vector<SvxCentralTrack*> m_vsvxcnt;
};

#endif

