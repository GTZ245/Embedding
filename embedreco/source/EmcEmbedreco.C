#include "EmcEmbedreco.h"
#include "Fun4AllReturnCodes.h"
#include "Fun4AllServer.h"

#include "mEmcClusterizerv0.h"
#include "mEmcGeometryModule.h"

#include "recoConsts.h"

#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "PHTypedNodeIterator.h"
#include "PHNodeIterator.h"

#include "TObject.h"

#include <string>
#include <vector>
#include "EmcMixer.hh"
#include "Fun4AllServer.h"

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

using namespace findNode;

using namespace std;

EmcEmbedreco::EmcEmbedreco(const char *name){
  ThisName    = name;
  clusterizer = NULL;
  emcmixer    = NULL;
}
EmcEmbedreco::~EmcEmbedreco(){
  delete clusterizer;
  delete emcmixer;
}
int EmcEmbedreco::Init(PHCompositeNode *topNode)
{
  if (verbosity > 0)
    {
      cout << "Calling Init" << endl;
    }
  // Create Histograms here - later you will have to do file magic
  // to make sure they are not deleted when the input file is closed
  return EVENT_OK;
}
int EmcEmbedreco::InitRun(PHCompositeNode *topNode){
  recoConsts *rc = recoConsts::instance();
  if (verbosity > 0){
    // this rc flag is set by the framework
    cout << "Calling InitRun for Run" << rc->get_IntFlag("RUNNUMBER") << endl;
  }
  PHNodeIterator iter(topNode);

  PHCompositeNode *parNode;
  parNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "PAR"));

  mEmcGeometryModule* geometryModule = 
    new mEmcGeometryModule(mEmcGeometryModule::kPISA);
  PHIODataNode<TObject>* geometryNode 
    = new PHIODataNode<TObject>(geometryModule, "mEmcGeometry");
  parNode->addNode(geometryNode);

  float fgTowerThresholdPbSc = 0.010; //GeV
  float fgTowerThresholdPbGl = 0.014;
  float fgMinClusterEnergyPbSc = 0.015;
  float fgMinClusterEnergyPbGl = 0.060;
  clusterizer = new mEmcClusterizerv0(geometryModule);  
  clusterizer->SetTowerThresholdPbSc(fgTowerThresholdPbSc);
  clusterizer->SetTowerThresholdPbGl(fgTowerThresholdPbGl);
  clusterizer->SetMinClusterEnergyPbSc(fgMinClusterEnergyPbSc);
  clusterizer->SetMinClusterEnergyPbGl(fgMinClusterEnergyPbGl);

  //initialize the stuff for embedding
  emcmixer = new EmcMixer;
  //mixer ->setVerbose(rc->get_IntFlag("VERBOSITY"));
  Fun4AllServer* se = Fun4AllServer::instance();
  
  PHCompositeNode* mcnode = se->topNode(rc->get_CharFlag("EMBED_MC_TOPNODE"));
  PHCompositeNode* realnode = se->topNode(rc->get_CharFlag("EMBED_REAL_TOPNODE"));
  PHCompositeNode* mergednode = se->topNode("TOP");
  
  emcmixer->InitRun(mcnode,realnode,mergednode);
  return EVENT_OK;
}

int EmcEmbedreco::process_event(PHCompositeNode *topNode){
  emcmixer->merge();
  clusterizer->process_event(topNode);
  
  if(false)
  {
    Fun4AllServer* se = Fun4AllServer::instance();
    recoConsts *rc = recoConsts::instance();
    PHCompositeNode* realnode = se->topNode(rc->get_CharFlag("EMBED_REAL_TOPNODE"));
    const emcClusterContainer* realemccont =
      findNode::getClass<emcClusterContainer>(realnode, "emcClusterContainer");
    const emcTowerContainer* realtwrcont =
      findNode::getClass<emcTowerContainer>(realnode, "emcTowerContainer");

    if(!realemccont||!realtwrcont) 
    {
      std::cout<<"kek0 "<<std::endl;
      return 0;
    }
    for (unsigned int iemc = 0; iemc < realemccont->size(); iemc++)
    {
      emcClusterContent* emc = realemccont->getCluster(iemc);
      if(!emc||emc->ecore()<0.4)  {
        continue;
      }
      std::vector<double> towerid, twrenergy;
      for (int itower=0; itower < emc->multiplicity(); itower++)
      {
        towerid.push_back(emc->towerid(itower));
	      emcTowerContent* twr = realtwrcont->findTower(emc->towerid(itower));
        if(!twr) continue;
        twrenergy.push_back(twr->Energy());
      }
      double twre = 0, twrecore = 0;
      //std::cout<<towerid.size()<<" "<<twrenergy.size()<<" "<<realtwrcont->size()<<std::endl;
      for (unsigned int itower=0; itower < twrenergy.size(); itower++)
      {
          twre += twrenergy[itower];
      }
      for (unsigned int itower=0; itower < twrenergy.size(); itower++)
      {
          if(twrenergy[itower]>0.02*twre) twrecore+=twrenergy[itower];
      }
      std::cout<<"embed ecore, e, twr ecore, e: "<<emc->ecore()<<" "<<emc->e()<<" "<<twrecore/0.932<<" "<<twre<<" "<<emc->multiplicity()<<std::endl;
    } 
  }

  return EVENT_OK;
}
int EmcEmbedreco::Reset(PHCompositeNode *topNode){
  return EVENT_OK;
}
int EmcEmbedreco::ResetEvent(PHCompositeNode *topNode)
{
  return EVENT_OK;
}
