#include <iostream>
#include <iomanip>
#include <cmath>

/*
#include "McEvalSingleTrack.h"
#include "RecoEvalSingleTrack.h"
#include "TrigLvl1.h"
#include "ErtOut.h"
*/

//#include "cglDetectorGeo.hh"
//#include "cglHitAssociate.hh"
//#include "cglTransformDST.hh"
//#include "mPHDchTrackModel.hh"
//#include "mPHLineTrack.hh"
//#include "PHDetectorGeometry.h"

#include "Fun4AllServer.h"
#include "Fun4AllReturnCodes.h"
#include "PHTypedNodeIterator.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"

#include "phool.h"
#include "PHGlobal.h"
#include "recoConsts.h"
#include "utiCentrality.h"
#include "getClass.h"
#include "BbcOut.h"


//#include "DetectorGeometry.h"
//#include "DetectorGeometryv1.h"

#include "DstMerger.h"

using namespace PhUtilities;
using namespace findNode;

//==============================================================
DstMerger::DstMerger() 
{
  ThisName = "DstMerger";
  std::cout << ThisName << " constructor called " << std::endl;
  EventNumber=0;
  runnumber = 0;
}

//==============================================================
int DstMerger::Init(PHCompositeNode *topNode) 
{
  std::cout << ThisName << "::Init called " << std::endl;
  vtxo.open( vtxfn.c_str());
  if ( ! vtxo.good() ) 
    {
      std::cout << ThisName << "::Init: couldn't open vtx out file " << vtxfn << std::endl;
      return -1;
    }
  return 0;
}

//==============================================================
int DstMerger::InitRun(PHCompositeNode *topNode) 
{
  
  std::cout << ThisName << "::InitRun called " << std::endl;
  recoConsts *rc = recoConsts::instance();
  runnumber = rc->get_IntFlag("RUNNUMBER"); 
  std::cout << "Initializing run " << runnumber << std::endl;
  return 0;
}

//==============================================================

int DstMerger::process_event(PHCompositeNode *topNode) 
{
  
  //std::cout << ThisName << "::process_event called " << std::endl;
  int iret=0;
  

  PHGlobal* evt = getClass<PHGlobal>(topNode,"PHGlobal");
  if (!evt) std::cout << PHWHERE << ThisName << ":: " << "PHGlobal" << " not in Node Tree" << std::endl;
  
  /*
  PHGlobal* evt = 0;
  PHTypedNodeIterator<PHGlobal> evtiter(topNode);
  PHIODataNode<PHGlobal> *EvtNode = evtiter.find("PHGlobal");
  if(EvtNode) { evt = (PHGlobal*)EvtNode->getData(); }
  //    else {std::cerr << "Can't find PHGlobal." << std::endl; return ABORTEVENT;}
  */

  /*
  BbcOut* bbc = 0;
  PHTypedNodeIterator<BbcOut> bbciter(topNode);
  PHIODataNode<BbcOut> *BbcNode = bbciter.find("BbcOut");
  if(BbcNode) { bbc = (BbcOut*)BbcNode->getData(); }
    else {std::cerr << "Can't find BbcOut." << std::endl; return ABORTEVENT;}
  */
  
  int evtnumber = -1;
  float zvtx    = -9999.;
  float cent      = -1;
  if(evt) {
    /*
    runnumber = evt->getRunNumber();
    evtnumber = evt->getEventNumber();
    */
    zvtx = evt->getBbcZVertex();
    //float bbc1 = evt->getBbcChargeN();
    //float bbc2 = evt->getBbcChargeS();
    //float zdc1 = evt->getZdcEnergyN();
    //float zdc2 = evt->getZdcEnergyS();
    cent = evt -> getCentrality();
  } 
  /*
  if(bbc) {
    zvtx = bbc->get_VertexPoint();
  }
  */
  
  if(EventNumber%1000==0) {
    std::cout << "Event: " << EventNumber << " " << runnumber << " " << evtnumber << " "
         << zvtx << " " << cent << std::endl;
  }
  
  if(fabs(zvtx)>30.0) { return ABORTEVENT; } 
  //if(zvtx<5 || zvtx>15.0) { return ABORTEVENT; } 	
  //if(zvtx<-15 || zvtx>-5.0) { return ABORTEVENT; } 	
  //if(zvtx<-25 || zvtx>-15.0) { return ABORTEVENT; } 	
  //if(zvtx<15 || zvtx>25.0) { return ABORTEVENT; } 

  vtxo << EventNumber << "   " << zvtx << "   " << cent << endl;
  EventNumber++;

  return iret;
}

//==============================================================

int DstMerger::End(PHCompositeNode *topNode) {
  vtxo.close();
  return 0;
}

void 
DstMerger::setVtxAsciiFileName(const char *filename)
{
  cout << "DstMerger::setVtxAsciiFileName : Output ascii file: " << filename << endl;
  vtxfn = filename;

}

