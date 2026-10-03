#include "TObject.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "PHTypedNodeIterator.h"
#include "PHNodeIterator.h"
#include "PHNodeReset.h"
#include "recoConsts.h"

#include <cmath>
#include "CrkEmbedreco.h"



#include "CrkHitv1.h"

#include "dCrkHitWrapper.h"
#include "CrkMixer.hh"
#include "Fun4AllServer.h"

using namespace std;

typedef PHIODataNode<PHObject> PHObjectNode_t;
typedef PHIODataNode<PHTable> PHTableNode_t;
typedef PHIODataNode<TObject> TObjectNode_t;
typedef PHIODataNode<CrkHit> CrkHitNode_t;

CrkEmbedreco::CrkEmbedreco(const char *name)
{
  ThisName = name;

  crkmixer         = NULL;
}

CrkEmbedreco::~CrkEmbedreco()
{
  delete crkmixer;
}

int CrkEmbedreco::Init(PHCompositeNode *topNode)
{
  return 0;
}

int CrkEmbedreco::InitRun(PHCompositeNode *topNode)
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
  
  PHCompositeNode* crkNode = new PHCompositeNode("CRK");
  topNode->addNode(crkNode);
  
  PHIODataNode<PHTable>* dCrkHitNode = (PHIODataNode<PHTable>*)iter.findFirst("PHIODataNode","dCrkHit");
  if(!dCrkHitNode){
    mr = 5120;
    dCrkHitWrapper* dCrkHit = new dCrkHitWrapper("dCrkHit", mr);
    dCrkHitNode = new PHIODataNode<PHTable>(dCrkHit, "dCrkHit");
    dstNode->addNode(dCrkHitNode);
  }

  PHObjectNode_t *CrkHitNode =  static_cast<PHObjectNode_t*>(iter.findFirst("PHIODataNode","CrkHit"));
  if(!CrkHitNode){
    CrkHit *crkhitv1 = new CrkHitv1();
    CrkHitNode =
      new PHIODataNode<PHObject>(crkhitv1, "CrkHit", "PHObject");
    dstNode->addNode(CrkHitNode);
  }
  //initialize the stuff for embedding
  crkmixer = new CrkMixer;
  crkmixer->setVerbose( Verbosity() ); // Verbosity() inherited from Fun4AllBase
  Fun4AllServer* se = Fun4AllServer::instance();
  
  PHCompositeNode* mcnode = se->topNode(rc->get_CharFlag("EMBED_MC_TOPNODE"));
  PHCompositeNode* realnode = se->topNode(rc->get_CharFlag("EMBED_REAL_TOPNODE"));
  PHCompositeNode* mergednode = se->topNode("TOP");
  
  crkmixer->InitRun(mcnode,realnode,mergednode);
  return 0;
}

int CrkEmbedreco::process_event(PHCompositeNode *topNode)
{
  crkmixer->merge();
  //copyWrapper(topNode);
  return 0;
}

int CrkEmbedreco::ResetEvent(PHCompositeNode *topNode)
{
  PHNodeIterator mainIter(topNode);
  PHNodeReset reset;
  if (mainIter.cd("CRK"))
    {
      mainIter.forEach(reset);
    }
  return 0;
}

int CrkEmbedreco::copyWrapper(PHCompositeNode *topNode)
{
  dCrkHitWrapper *dcrkhit = NULL;
  PHTypedNodeIterator<dCrkHitWrapper> dcrkhit_iter(topNode);
  PHIODataNode <dCrkHitWrapper> *dCrkHitNode = dcrkhit_iter.find("dCrkHit");
  if (dCrkHitNode)
    {
      dcrkhit = dCrkHitNode->getData();
    }
  else
    {
      return -1;
    }

  CrkHit *crkhit = NULL;
  PHTypedNodeIterator<CrkHit> crkhititer(topNode);
  CrkHitNode_t *CrkHitNode = crkhititer.find("CrkHit");
  if (CrkHitNode)
    {
      crkhit = CrkHitNode->getData();
    }
  if (dcrkhit && crkhit)
    {
      if (!crkhit->isValid() && dcrkhit->RowCount())
        {
          crkhit->FillFromWrapper(dcrkhit);
        }
    }
  return 0;
}
