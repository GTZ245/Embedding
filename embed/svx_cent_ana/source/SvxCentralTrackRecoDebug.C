#include <iostream>

#include "SvxCentralTrackRecoDebug.h"

#include "SvxClusterContainer.h"
#include "SvxClusterListv4.h"
#include "PHSnglCentralTrackv24.h"

using namespace std;

SvxCentralTrackRecoDebug::SvxCentralTrackRecoDebug() 
  : SvxCentralTrackReco("SvxCentralTrackRecoDebug") {
}

void SvxCentralTrackRecoDebug::testLinkCluster(int nclus, float **ary) {
  ////////////////////////
  // make dummy data
  cout<<"GenerageClusterList "<<endl;
  SvxClusterList* clslist = new SvxClusterListv4();

  const int NCLS = 6;
  float clspos[NCLS][4] = {
   {3, 15.3658, 0.477552, -1.13551},
   {2, 10.2177, 0.473301, -1.88446},
   {1, 5.18563, 0.353456, -2.54975},
   {0, 2.58232, 0.205352, -2.94392}

/*
    {3, 15.0813, 2.80149, -0.937701},
    {2, 10.013,  1.98752, -2.15},
    {1, 4.93342, 1.05181, -3.31225},
    {0, 2.50278, 0.547906, -3.86475}
*/
/*
    {3, -17.4984, -1.42511,   0.45664},
    {2, -11.4051, -1.09733,   1.75   },
    {1, -5.2,     -0.363796,  2.89397},
    {1, -5.10455, -0.263796,  2.89397},
    {1, -5.0,     -0.163796,  2.89397},
    {0, -2.65076, -0.0490907, 3.34776},
*/
  };

  for(int i=0; i<NCLS; i++){
    SvxCluster *cluster = clslist->addCluster();
    cluster->set_svxSection(0);
    cluster->set_layer((int)clspos[i][0]);
//--    cluster->set_ladder(ladder);
//--    cluster->set_sensor(sensor);
//--    cluster->set_sensorType(sensorType);
//--    cluster->set_adc(0, 1);
//--    cluster->set_adc(1, 1);
//--    cluster->set_size(cluster_list[i].nhits);
//--    cluster->set_xz_size(0, xsize);
//--    cluster->set_xz_size(1, zsize);
//--    cluster->set_circumference(circumference);
//--    cluster->set_edgeflag(edge_flag);
    for ( int j=0; j<3; j++ )
      {
//--        cluster->set_xyz_local (j, clspos[i][j+1]);
        cluster->set_xyz_global(j, clspos[i][j+1]);
      }
  }

  cout<<"fill ClusterContainer "<<endl;

  //float beam_xy[2] = {0.23, 0.08};
  float beam_xy[2] = {0.0, 0.0};

  SvxClusterContainer* container = new SvxClusterContainer();
  container->set_beam_center(beam_xy[0], beam_xy[1]);


  /// load clusters
  container->load_clusters(clslist);


  cout<<"generate CNT track "<<endl;
  ////////////////////////
  // exec LinkCluster
  //Particle : 0.825195, -1, 3.0625, 1.76074, Zvtx= 0.684868,
  //Particle : 0.488525, 1, 0.239502, 1.35449, Zvtx= -4.424, 
  //Particle : 0.402832, 1, 0.059082, 1.43164, Zvtx= -3.3406, 

  float mom = 0.402832;
  float charge =   1;
  float phi0 = 0.059082;
  float the0 = 1.43164;
  //float vtx[3] = {beam_xy[0], beam_xy[1], 0.684868};
  float vtx[3] = {beam_xy[0], beam_xy[1], -3.3406};

  ////////////////////////
  // central track
  PHSnglCentralTrack *sngltrk = new PHSnglCentralTrackv24();
  sngltrk->set_mom   (mom);
  sngltrk->set_charge(charge);
  sngltrk->set_phi0  (phi0);
  sngltrk->set_the0  (the0);


  SvxCentralClusterLink* trkseed = new SvxCentralClusterLink(0, sngltrk);
  m_vtrklist.push_back(trkseed);

  cout<<"linkCluster "<<endl;
  ////////////////////////
  // linkcluster
  SvxClsLink link;

  cout<<"Test listsize  "<<m_vtrklist.size()<<endl;

  SvxTrackPart part(mom, charge, phi0, the0);

  LinkClusters(part, // track info
               vtx[0], vtx[1], vtx[2],  // primary vtx
               link,  container,        // link data
               7, 0.2, 3.0);// search condition


  if(verbosity>1) {
    trkseed->print();
  }

}
