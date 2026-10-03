#include "TecEmbedreco.h"
#include "mTecHoughTrackModule.h"
#include "TecAddressObject.hh"
#include "TecCalibrationObject.hh"
#include "TecGeometryObject.hh"
#include "TecOutV3.hh"

#include "dTecFemDataWrapper.h"
#include "dTecGhitRawWrapper.h"
#include "tecghitWrapper.h"

#include "TecGetGEA.h"
#include "TecMixer.hh"

#include "recoConsts.h"

#include "Fun4AllServer.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "PHTypedNodeIterator.h"
#include "PHNodeIterator.h"
#include "PHNodeReset.h"


#include "TObject.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

typedef PHIODataNode<PHObject> PHObjectNode_t;
typedef PHIODataNode<PHTable> PHTableNode_t;
typedef PHIODataNode<TObject> TObjectNode_t;

long TecPutDCM(PHCompositeNode* topNode);

TecEmbedreco::TecEmbedreco(const char *name)
{
  ThisName = name;

  mTecHoughTrack = NULL;
  TecAddress     = NULL;
  TecGeometry    = NULL;
  TecCalibration = NULL;
  tecmixer       = NULL;
}

TecEmbedreco::~TecEmbedreco()
{
  if ( tecmixer)
    {
      delete tecmixer;
    }
  if ( mTecHoughTrack)
    {
      delete mTecHoughTrack;
    }
  if ( TecAddress)
    {
      delete TecAddress;
    }
  if ( TecGeometry)
    {
      delete TecGeometry;
    }
  if ( TecCalibration)
    {
      delete TecCalibration;
    }
}

int TecEmbedreco::Init(PHCompositeNode *topNode)
{

  //  mTecDraw = new mTecDrawModule;
  mTecHoughTrack = new mTecHoughTrackModule;

  TecAddress = new TecAddressObject();
  TecGeometry = new TecGeometryObject();
  TecCalibration = new TecCalibrationObject();
  return 0;
}

int TecEmbedreco::InitRun(PHCompositeNode *topNode)
{
  recoConsts *rc = recoConsts::instance();

  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode;
  dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode) { cout << PHWHERE << "DST Node missing doing nothing" << endl; return -1; }
  
  PHCompositeNode *parNode;
  parNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "PAR"));
  if (!parNode) { cout << PHWHERE << "PAR Node missing doing nothing" << endl; return -1; }

  PHCompositeNode *dcmNode;
  dcmNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DCM"));
  if (!dcmNode) { cout << PHWHERE << "DCM Node missing doing nothing" << endl; return -1; }

  PHCompositeNode *evaNode;
  evaNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVA"));
  if (!evaNode) { cout << PHWHERE << "EVA Node missing doing nothing" << endl; return -1; }

  PHCompositeNode *geaNode;
  geaNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "GEA"));
  if (!geaNode) { cout << PHWHERE << "GEA Node missing doing nothing" << endl; return -1; }

  /* here comes the TEC */
  PHCompositeNode* tecNode = new PHCompositeNode("TEC");
  topNode->addNode(tecNode);

  PHObjectNode_t *TecOutNode =  static_cast<PHObjectNode_t*>(iter.findFirst("PHIODataNode","TecOut"));
  if(!TecOutNode){
    TecOut *tecoutv3 = new TecOutV3();
    TecOutNode = new PHObjectNode_t(tecoutv3, "TecOut", "PHObject");
    dstNode->addNode(TecOutNode);
  }
  PHObjectNode_t *TecHitOutNode =  static_cast<PHObjectNode_t*>(iter.findFirst("PHIODataNode","TecHitOut"));
  if(!TecHitOutNode){
    TecOut *techitout = new TecOutV3();
    TecHitOutNode = new PHObjectNode_t(techitout, "TecHitOut", "PHObject");
    dstNode->addNode(TecHitOutNode);
  }
  int mr;

  mr = 4500;
  tecghitWrapper* tecghit = new tecghitWrapper("tecghit", mr);
  PHIODataNode<PHTable>* tecghitNode = new PHIODataNode<PHTable>(tecghit, "tecghit");
  geaNode->addNode(tecghitNode);
  
  mr = 150000;
  dTecGhitRawWrapper* dTecGhitRaw = new dTecGhitRawWrapper("dTecGhitRaw", mr);
  PHIODataNode<PHTable>* dTecGhitRawNode = new PHIODataNode<PHTable>(dTecGhitRaw, "dTecGhitRaw");
  evaNode->addNode(dTecGhitRawNode);
  
  mr = 8000;
  dTecFemDataWrapper* dTecFemData = new dTecFemDataWrapper("dTecFemData", mr);
  PHIODataNode<PHTable>* dTecFemDataNode = new PHIODataNode<PHTable>(dTecFemData, "dTecFemData");
  tecNode->addNode(dTecFemDataNode);
  
  PHIODataNode<TObject>* TecDetGeoNode = new PHIODataNode<TObject> (TecGeometry, "TecGeometry");
  parNode->addNode(TecDetGeoNode);
  
  PHDataNode<TecCalibrationObject>* TecCalibNode = new PHDataNode<TecCalibrationObject> (TecCalibration, "TecCalibration");
  parNode->addNode(TecCalibNode);
  
//  // In case of simulation the Time Stamp will correspond to system
//  // time. OK for standard simulation, but could cause problems
//  // for non-standard geometry (e.g. retracted arms).
//  PHTimeStamp TimeStp = rc->get_TimeStamp();
// For now set TimeStamp to January 15, 2003 (run3 conditions) SL 12/11/30
// This must be changed later so that simulations for older runs are possible.
  PHTimeStamp TimeStp = PHTimeStamp(2003, 1, 15, 0, 0, 0);
  TecAddress->setTimeStamp(TimeStp);
  TecGeometry->setTimeStamp(TimeStp);
  TecCalibration->setTimeStamp(TimeStp);
  
  TecAddress->UseSimulationDatabase();
  TecGeometry->UseSimulationDatabase();
  TecCalibration->UseSimulationDatabase();

  // If database is not accessible, stop execution.
  PHBoolean status1 = TecAddress->Fetch();
  if (!status1)
    {
      cout << PHWHERE << __FILE__ << ":" << __LINE__ << ":" << endl;
      cout << PHWHERE << "ERROR: Can not access the database (TAO). Execution terminated." << endl;
      exit(1);
    }
  PHBoolean status2 = TecGeometry->Fetch();
  if (!status2)
    {
      cout << PHWHERE << __FILE__ << ":" << __LINE__ << ":" << endl;
      cout << PHWHERE << "ERROR: Can not access the database (TGO). Execution terminated." << endl;
      exit(1);
    }
  PHBoolean status3 = TecCalibration->Fetch();
  if (!status3)
    {
      cout << PHWHERE << __FILE__ << ":" << __LINE__ << ":" << endl;
      cout << PHWHERE << "ERROR: Can not access the database (TCO). Execution terminated." << endl;
      exit(1);
    }
  
  //initialize the stuff for embedding
  tecmixer = new TecMixer;
  //mixer ->setVerbose(rc->get_IntFlag("VERBOSITY"));
  Fun4AllServer* se = Fun4AllServer::instance();
  
  PHCompositeNode* mcnode = se->topNode(rc->get_CharFlag("EMBED_MC_TOPNODE"));
  PHCompositeNode* realnode = se->topNode(rc->get_CharFlag("EMBED_REAL_TOPNODE"));
  PHCompositeNode* mergednode = se->topNode("TOP");
  
  tecmixer->InitRun(mcnode,realnode,mergednode);
  return 0;
}

int TecEmbedreco::process_event(PHCompositeNode *topNode)
{
  tecmixer->merge();
  mTecHoughTrack->set_fillRawTrack(0);
  mTecHoughTrack->event(topNode,TecGeometry, TecCalibration);
  return 0;
}

int TecEmbedreco::ResetEvent(PHCompositeNode *topNode)
{  
  PHNodeIterator mainIter(topNode);
  PHNodeReset reset;
  if (mainIter.cd("TEC"))
    {
      mainIter.forEach(reset);
    }
  return 0;
}
