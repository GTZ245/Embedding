#include "TObject.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "PHTypedNodeIterator.h"
#include "PHNodeIterator.h"
#include "PHNodeReset.h"
#include "recoConsts.h"

#include <cmath>

#include "TofwEmbedreco.h"
#include <TofwGeometry.h>
#include <RunHeader.h>
#include <getClass.h>

#include "TofwHit.h"
#include "TofwHitv1.h"
#include "TofwMixer.hh"
#include "Fun4AllServer.h"

using namespace std;

typedef PHIODataNode<PHObject> PHObjectNode_t;
typedef PHIODataNode <TofwHit> TofwHitNode_t;

TofwEmbedreco::TofwEmbedreco(const char *name)
{
  ThisName = name;

  tofwmixer      = NULL;

  d_geom   = 0;
}

TofwEmbedreco::~TofwEmbedreco()
{
  if (d_geom)
    {
      delete d_geom;
    }
  if ( tofwmixer)
    {
      delete tofwmixer;
    }
}

int TofwEmbedreco::Init(PHCompositeNode *topNode)
{
  return 0;
}

int TofwEmbedreco::InitRun(PHCompositeNode *topNode)
{
  recoConsts *rc = recoConsts::instance();
 
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode;
  dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      cout << PHWHERE << "DST Node missing doing nothing" << endl;
      return -1;
    }

  d_geom = new TofwGeometry();

  PHCompositeNode *parNode;
  parNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "PAR"));
  if (!parNode)
    {
	   cout << PHWHERE << "PAR Node missing doing nothing" << endl;
	   return -1;

    }

  PHDataNode<TofwGeometry> *NewNode = new PHDataNode<TofwGeometry>(d_geom,"TofwGeometry");
  parNode->addNode(NewNode);

  PHObjectNode_t *TofwHitNode =  static_cast<PHObjectNode_t*>(iter.findFirst("PHIODataNode","TofwHit"));
  if(!TofwHitNode){
    TofwHit *tofwhitv1 = new TofwHitv1();
    TofwHitNode = new PHIODataNode<PHObject>(tofwhitv1, "TofwHit", "PHObject");
    dstNode->addNode(TofwHitNode);
  }

  //initialize the stuff for embedding
  tofwmixer = new TofwMixer();
  //mixer ->setVerbose(rc->get_IntFlag("VERBOSITY"));
  Fun4AllServer* se = Fun4AllServer::instance();
  
  PHCompositeNode* mcnode = se->topNode(rc->get_CharFlag("EMBED_MC_TOPNODE"));
  PHCompositeNode* realnode = se->topNode(rc->get_CharFlag("EMBED_REAL_TOPNODE"));
  PHCompositeNode* mergednode = se->topNode("TOP");
  
  tofwmixer->InitRun(mcnode,realnode,mergednode);
  return 0;
}

int TofwEmbedreco::process_event(PHCompositeNode *topNode)
{
  tofwmixer->merge();
  return 0;
}


int TofwEmbedreco::ResetEvent(PHCompositeNode *topNode)
{
  /*
  PHNodeIterator mainIter(topNode);
  PHNodeReset reset;
  if (mainIter.cd("TOFW"))
    {
      mainIter.forEach(reset);
    }
  */
  return 0;
}

