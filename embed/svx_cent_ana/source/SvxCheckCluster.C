#include "SvxCheckCluster.h"

#include <iostream>

#include <TFile.h>
#include <TNtuple.h>

#include <phool.h>
#include <PHTypedNodeIterator.h>
#include <PHCompositeNode.h>
#include <PHIODataNode.h>
#include <Fun4AllReturnCodes.h>
#include <PHPoint.h>
#include <SvxClusterList.h>
#include <SvxCluster.h>
#include <compactCNT/SvxHitMap.h>
#include <VtxOut.h>
#include <getClass.h>

using namespace std;
using namespace findNode;

//==============================================================

SvxCheckCluster::SvxCheckCluster() : 
      SubsysReco("SvxCheckCluster"),
      m_outputFile(NULL),
      m_ntpclus(NULL),
      m_eventNumber(0)
{
}

//==============================================================

SvxCheckCluster::~SvxCheckCluster() {
}

//==============================================================

int SvxCheckCluster::Init(PHCompositeNode *topNode) 
{
  cout << "SvxCheckCluster::Init started..." << endl;

  m_outputFile = new TFile("./test.root", "RECREATE");

  m_ntpclus = new TNtuple("ntpclus","","event:xvtx:yvtx:zvtx:x:y:z:layer:ladder");


  cout << "SvxCheckCluster::Init ended." << endl;
  return 0;
}

//==============================================================
  
int SvxCheckCluster::InitRun(PHCompositeNode *topNode) {
  cout << "SvxCheckCluster::InitRun started..." << endl;

  m_eventNumber = 0;
  cout << "SvxCheckCluster::InitRun ended." << endl;
  return 0;
}

//==============================================================

int SvxCheckCluster::End(PHCompositeNode *topNode) {
  cout << "SvxCheckCluster::End:  Writing out..." << endl;

  m_outputFile->Write();
  m_outputFile->Close();
  delete m_outputFile;
  m_outputFile=0;

  return 0;
}

//==============================================================

int SvxCheckCluster::process_event(PHCompositeNode *topNode) {

  SvxClusterList* svx_list     = getClass<SvxClusterList>(topNode,"SvxClusterList");
  SvxClusterList* svxsel_list  = getClass<SvxClusterList>(topNode,"SvxSelectedClusterList");
  SvxHitMap*      svxrec_list  = getClass<SvxHitMap>(topNode, "SvxHit_comp");

  VtxOut*         vtxout       = getClass<VtxOut>(topNode, "VtxOut");
  PHPoint vtx(-999, -999, -999);
  if(vtxout!=NULL){
    vtx = vtxout->get_Vertex();
  }

  int nsvx    = (svx_list   !=NULL) ? svx_list->get_nClusters() : -9999;
  int nsvxsel = (svxsel_list!=NULL) ? svxsel_list->get_nClusters() : -9999;
  int nsvxrec = (svxrec_list!=NULL) ? svxrec_list->GetNentries() : -9999;

  cout<<"nhit : "<<nsvx<<" "<<nsvxsel<<" "<<nsvxrec<<endl;

  if(svxsel_list!=NULL)
    {
      for(int ic=0; ic<nsvxsel; ic++)
        {
          SvxCluster *cls = svxsel_list->get_Cluster(ic);

          float ntpary[] = {
              (float)m_eventNumber,
              (float)vtx.getX(),
              (float)vtx.getY(),
              (float)vtx.getZ(),
              cls->get_xyz_global(0),
              cls->get_xyz_global(1),
              cls->get_xyz_global(2),
              (float)cls->get_layer(),
              (float)cls->get_ladder()
            };

          m_ntpclus->Fill(ntpary);
        }
    }

  m_eventNumber++;

  return 0;
}



