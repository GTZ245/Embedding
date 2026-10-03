
#include "PHCompositeNode.h"
#include "PHNodeReset.h"
#include "PHIODataNode.h"
#include "phool.h"
#include "RunHeader.h"

#include "AccEmbedreco.h"

#include "AccRawv2.h"
#include "AccHitv1.h"
#include "AccEvent.h"

#include "AccGeometry.h"

typedef PHIODataNode<RunHeader> RunHeaderNode_t;
typedef PHIODataNode <PHObject> PHObjectNode_t;
typedef PHIODataNode <AccRaw>   AccRawNode_t;
typedef PHIODataNode <AccHit>   AccHitNode_t;
typedef PHDataNode <AccGeometry> AccGeometryNode_t;

using namespace std;

AccEmbedreco::AccEmbedreco(const char *name, const int ver)
 : SubsysReco(name)
{
  version = ver; 
  d_acc = 0;
  d_geo = 0;
  d_hit = 0;
  d_raw = 0;
}

AccEmbedreco::~AccEmbedreco()
{ 
  if ( d_acc)
    {
      delete d_acc;
    }
  if ( d_geo)
    {
      delete d_geo;
    }
  if ( d_hit)
    {
      delete d_hit;
    }
  if ( d_raw)
    {
      delete d_raw;
    }
}

int AccEmbedreco::Init(PHCompositeNode* topNode)
{

  d_acc = new AccEvent();
  return 0;
}


int AccEmbedreco::InitRun(PHCompositeNode* topNode)
{
  //In order to run embedding, the nodes have to be put in the InitRu.
  //otherwise phool will not be able to properly read in the events from DSTs 
  CreateNodeTree(topNode);


  RunHeader* d_runhdr = 0;

  PHTypedNodeIterator<RunHeader> iRUN(topNode);
  RunHeaderNode_t *RUN = iRUN.find("RunHeader");
  if (RUN) d_runhdr = RUN->getData();
  if (!d_runhdr) cout << PHWHERE << "AccEmbedreco:: runhdr not in Node Tree" << endl;

  if (d_runhdr){
    int runNumber = d_runhdr->get_RunNumber();
    d_geo->fetch(runNumber);
  }

  return 0;
}

int AccEmbedreco::ResetEvent(PHCompositeNode* topNode)
{

  d_acc->Reset(topNode);

  PHNodeIterator mainIter(topNode);
  PHNodeReset reset;
  if(mainIter.cd("ACC")) {
    mainIter.forEach(reset);
  }

  return 0;
}


int AccEmbedreco::process_event(PHCompositeNode* topNode)
{

  return d_acc->process_event(topNode);
}

bool AccEmbedreco::CreateNodeTree(PHCompositeNode* topNode)
{

  // Find the DST node so we can put objects there
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode;
  dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","DST"));
  PHCompositeNode* parNode;
  parNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","PAR"));

  AccRawNode_t *raw =  static_cast<AccRawNode_t*>(iter.findFirst("PHIODataNode","AccRaw"));
  AccHitNode_t *hit =  static_cast<AccHitNode_t*>(iter.findFirst("PHIODataNode","AccHit"));

  // make the data objects
  switch ( version ){
    case 1:
      if(!raw){
	d_raw = new AccRawv2();
	raw = new PHIODataNode<AccRaw>(d_raw,"AccRaw","PHObject");
	dstNode->addNode(raw);
      }else{
	d_raw = raw->getData();
      }
      if(!hit){
	d_hit = new AccHitv1();
	hit = new PHIODataNode<AccHit>(d_hit,"AccHit","PHObject");
	dstNode->addNode(hit);
      }else{
	d_hit = hit->getData();
      }
      
      break;

    default:
      cout << PHWHERE << Name() << "... Unknown version requested: " << version << endl;
      cout << PHWHERE << Name() << "... You will receive version 1 instead. " << endl;

      if(!raw){
	d_raw = new AccRawv2();
	raw = new PHIODataNode<AccRaw>(d_raw,"AccRaw","PHObject");
	dstNode->addNode(raw);
      }else{
	d_raw = raw->getData();
      }
      if(!hit){
	d_hit = new AccHitv1();
	hit = new PHIODataNode<AccHit>(d_hit,"AccHit","PHObject");
	dstNode->addNode(hit);
      }else{
	d_hit = hit->getData();
      }
      break;
  }

  // Make nodes containing the objects
  AccGeometryNode_t *geo =  static_cast<AccGeometryNode_t*>(iter.findFirst("PHDataNode","AccGeometry"));

  if(!geo){
    // make the geometry object
    d_geo = new AccGeometry();
    geo = new PHDataNode<AccGeometry>(d_geo,"AccGeometry","PHObject");
    // put the nodes in the tree below DST
    parNode->addNode(geo);
  }else{
    d_geo = geo->getData();
  }

  return True;
}
