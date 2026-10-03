#include "TObject.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "PHTypedNodeIterator.h"
#include "PHNodeIterator.h"
#include "PHNodeReset.h"
#include "recoConsts.h"

#include <cmath>

#include "TofEmbedreco.h"
#include "TofAddressObject.hh"
#include "TofCalibObject.hh"
#include "TofGeometryObject.hh"

#include "dTofDCMParWrapper.h"
#include "dTofDCMWrapper.h"
#include "dTofFEMhitGhitWrapper.h"
#include "dTofFEMmapWrapper.h"
#include "dTofFEMWrapper.h"
#include "dTofGdigiRecWrapper.h"
#include "dTofGdigiWrapper.h"
#include "dTofGhitGdigiWrapper.h"
#include "dTofGhitRawWrapper.h"
#include "dTofRawRecWrapper.h"
#include "dTofRawWrapper.h"
#include "dTofReconstructedWrapper.h"
#include "tofghitWrapper.h"

#include "TofPISAHit.h"
#include "TofOutv2.h"
#include "TofGetGEA.h"
#include "TofMixer.hh"
#include "Fun4AllServer.h"

#include <cstdlib>

using namespace std;

typedef PHIODataNode<PHObject> PHObjectNode_t;
typedef PHIODataNode<PHTable> PHTableNode_t;
typedef PHIODataNode<TObject> TObjectNode_t;
typedef PHIODataNode <TofOut> TofOutNode_t;

long TofPutDCM(PHCompositeNode* topNode);

TofEmbedreco::TofEmbedreco(const string &name):
  SubsysReco(name),
  TofAddress(NULL),
  TofCalib(NULL),
  TofGeometry(NULL),
  tofmixer(NULL)
{}

TofEmbedreco::~TofEmbedreco()
{
  delete TofAddress;
  delete TofCalib;
//  delete TofGeometry;
  delete tofmixer;
}


int 
TofEmbedreco::InitRun(PHCompositeNode *topNode)
{
  recoConsts *rc = recoConsts::instance();
  
  size_t mr;
  
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode;
  dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      cout << PHWHERE << "DST Node missing doing nothing" << endl;
      return -1;
    }
  PHCompositeNode *parNode;
  parNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "PAR"));

  PHCompositeNode *evaNode;
  evaNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVA"));

  // here comes the TOF
  PHCompositeNode* tofNode = new PHCompositeNode("TOF");
  topNode->addNode(tofNode);

  TofAddress = new TofAddressObject();
  TofGeometry = new TofGeometryObject();
  TofCalib = new TofCalibObject();

  PHDataNode<TofGeometryObject>* TofDetGeoNode =
    new PHDataNode<TofGeometryObject>(TofGeometry, "TofGeometry");
  parNode->addNode(TofDetGeoNode);
  
  mr = 960;
  dTofRawWrapper* dTofRaw = new dTofRawWrapper("dTofRaw", mr);
  PHIODataNode<PHTable>* dTofRawNode =
    new PHIODataNode<PHTable>(dTofRaw, "dTofRaw");
  tofNode->addNode(dTofRawNode);
  
  mr = 960;
  PHIODataNode<PHTable>* dTofReconstructedNode = 
    (PHIODataNode<PHTable>*)iter.findFirst("PHIODataNode","dTofReconstructed");
  if(!dTofReconstructedNode){
    dTofReconstructedWrapper *dTofReconstructed = 
      new dTofReconstructedWrapper("dTofReconstructed", mr);
    dTofReconstructedNode =
      new PHIODataNode<PHTable>(dTofReconstructed, "dTofReconstructed");
    dstNode->addNode(dTofReconstructedNode);
  }
  mr = 960;
  dTofRawRecWrapper* dTofRawRec = new dTofRawRecWrapper("dTofRawRec", mr);
  PHIODataNode<PHTable>* dTofRawRecNode =
    new PHIODataNode<PHTable>(dTofRawRec, "dTofRawRec");
  evaNode->addNode(dTofRawRecNode);
  
  //===============================================================
  // simulation parameters stored in DB  
  // last update : Oct. 18, 2007, by T.Chujo
  //
  //   Run1 : 2001.1.1.0.0.0 - 2001.1.1.23.59.59
  //   Run2 : 2002.2.1.0.0.0 - 2002.2.1.23.59.59
  //   Run3 : (used Run2 time stamp)
  //   Run4 : 2004.4.4.0.0.0 - 2004.4.4.23.59.59
  //   Run5 : 2005.1.1.0.0.0 - 2005.1.1.23.59.59
  //   Run6 : (will be used Run5 time stamp)
  //   Run7 : (used Run5 time stamp)
  //============================================================
  PHTimeStamp TimeStp; 

  // Run1 simulation (default, but there is little or no Run1 analysis left)
  bool time_stamp_set( false );
  
  if ( rc->FlagExist("RUN2AUAU") && rc->get_IntFlag("RUN2AUAU") )
    {
      // Run2 simulation
      TimeStp = PHTimeStamp(2002, 2, 1, 12, 0, 0);
      time_stamp_set = true;
    }
  else if ( rc->FlagExist("RUN3DAU") && rc->get_IntFlag("RUN3DAU") )
    {
      // Run3 simulation (same as Run2)
      TimeStp = PHTimeStamp(2002, 2, 1, 12, 0, 0);
      time_stamp_set = true;
    }
  else if(
	  ( rc->FlagExist("RUN4AUAU63GEV")  && rc->get_IntFlag("RUN4AUAU63GEV") ) || 
	  ( rc->FlagExist("RUN4AUAU200GEV") && rc->get_IntFlag("RUN4AUAU200GEV") ) )
    { 
      // Run4 simulation
      TimeStp = PHTimeStamp(2004, 4, 4, 12, 0, 0); 
      time_stamp_set = true; 
    }
  else if( rc->FlagExist("RUN5PP200GEV") && rc->get_IntFlag("RUN5PP200GEV") ) 
    {
      // Run5 simulation
      TimeStp = PHTimeStamp(2005, 1, 1, 12, 0, 0); 
      time_stamp_set = true; 
    } 
  else if( rc->FlagExist("RUN7AUAU200GEV") && rc->get_IntFlag("RUN7AUAU200GEV") ) 
    {
      // Run7 simulation (same as Run5) 
      TimeStp = PHTimeStamp(2005, 1, 1, 12, 0, 0); 
      time_stamp_set = true; 
    }

  if( time_stamp_set ) cout << "TofEmbedreco::InitRun - time stamp for MC parameters : " << TimeStp << "\n" << endl;
  else {
    cout << "TofEmbedreco::InitRun - missing Rhic run reconConst (e.g. RUN7AUAU200GEV) to define time stamp." << endl;
    cout << "TofEmbedreco::InitRun - exiting" << endl;
    exit(1);
  }

  TofAddress->setTimeStamp(TimeStp);
  TofAddress->fetch();
  TofGeometry->setTimeStamp(TimeStp);
  TofGeometry->fetch();
  TofCalib->setTimeStamp(TimeStp);
  TofCalib->fetch();

  PHObjectNode_t *TofOutNode =  static_cast<PHObjectNode_t*>(iter.findFirst("PHIODataNode","TofOut"));
  if(!TofOutNode){
    TofOut *tofoutv2 = new TofOutv2();
    TofOutNode =
      new PHIODataNode<PHObject>(tofoutv2, "TofOut", "PHObject");
    dstNode->addNode(TofOutNode);
  }
  //initialize the stuff for embedding
  tofmixer = new TofMixer;
  //mixer ->setVerbose(rc->get_IntFlag("VERBOSITY"));
  Fun4AllServer* se = Fun4AllServer::instance();
  
  PHCompositeNode* mcnode = se->topNode(rc->get_CharFlag("EMBED_MC_TOPNODE"));
  PHCompositeNode* realnode = se->topNode(rc->get_CharFlag("EMBED_REAL_TOPNODE"));
  PHCompositeNode* mergednode = se->topNode("TOP");
  
  tofmixer->InitRun(mcnode,realnode,mergednode);
  return 0;
}

int TofEmbedreco::process_event(PHCompositeNode *topNode)
{

  tofmixer->merge();
  //mTofRawRec->event(topNode, TofAddress, TofGeometry, TofCalib); // comment out by AK 9/17/2002
  //copyWrapper(topNode);
  return 0;
}

int TofEmbedreco::ResetEvent(PHCompositeNode *topNode)
{
  PHNodeIterator mainIter(topNode);
  PHNodeReset reset;
  if (mainIter.cd("TOF"))
    {
      mainIter.forEach(reset);
    }

  return 0;
}

int TofEmbedreco::copyWrapper(PHCompositeNode *topNode)
{
  dTofReconstructedWrapper *dtofout = NULL;
  PHTypedNodeIterator<dTofReconstructedWrapper> dtofout_iter(topNode);
  PHIODataNode <dTofReconstructedWrapper> *dTofOutNode = dtofout_iter.find("dTofReconstructed");
  if (dTofOutNode)
    {
      dtofout = dTofOutNode->getData();
    }
  else
    {
      cout << PHWHERE << "dTofReconstructed Node not found" << endl;
      return -1;
    }

  TofOut *tofout = NULL;
  PHTypedNodeIterator<TofOut> tofoutiter(topNode);
  TofOutNode_t *TofOutNode = tofoutiter.find("TofOut");
  if (TofOutNode)
    {
      tofout = TofOutNode->getData();
    }
  else
    {
      cout << PHWHERE << "TofOut Node not found" << endl;
      return -1;
    }

  if (dtofout && tofout)
    {
      if (!tofout->isValid() && dtofout->RowCount())
        {
          tofout->FillFromWrapper(dtofout);
        }
      if (verbosity > 0)
	{
          tofout->identify();
	}
    }
  return 0;
}
