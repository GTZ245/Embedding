#include "BbcEmbedreco.h"

#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "PHNodeIterator.h"
#include "PHNodeReset.h"
#include "recoConsts.h"

#include "VtxMixer.hh"
#include "Fun4AllServer.h"

#include <iostream>

using namespace std;

BbcEmbedreco::BbcEmbedreco(const string &name): 
  SubsysReco(name),
  vtxmixer(NULL)
{}

BbcEmbedreco::~BbcEmbedreco()
{
  delete vtxmixer;
}


int BbcEmbedreco::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode;
  dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode) { 
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  
  PHCompositeNode *parNode;
  parNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "PAR"));
  if (!parNode) {
    parNode = new PHCompositeNode("PAR");
    topNode->addNode(parNode);
  }

  PHCompositeNode *dcmNode;
  dcmNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DCM"));
  if (!dcmNode) { 
    dcmNode = new PHCompositeNode("DCM");
    topNode->addNode(dcmNode);
  }

  PHCompositeNode *evaNode;
  evaNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVA"));
  if (!evaNode) { 
    evaNode = new PHCompositeNode("EVA");
    topNode->addNode(evaNode);
  }

  PHCompositeNode *geaNode;
  geaNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "GEA"));
  if (!geaNode) { 
    geaNode = new PHCompositeNode("GEA");
    topNode->addNode(geaNode);
  }
 
  recoConsts *rc = recoConsts::instance();
  vtxmixer = new VtxMixer;
  //mixer ->setVerbose(rc->get_IntFlag("VERBOSITY"));
  Fun4AllServer* se = Fun4AllServer::instance();
  
  PHCompositeNode* mcnode = se->topNode(rc->get_CharFlag("EMBED_MC_TOPNODE"));
  PHCompositeNode* realnode = se->topNode(rc->get_CharFlag("EMBED_REAL_TOPNODE"));
  PHCompositeNode* mergednode = se->topNode("TOP");
  
  vtxmixer->InitRun(mcnode,realnode,mergednode);
  return 0;
}

int BbcEmbedreco::process_event(PHCompositeNode *topNode)
{
  
  vtxmixer->merge();
  return 0;
}

int BbcEmbedreco::ResetEvent(PHCompositeNode *topNode)
{
  PHNodeIterator mainIter(topNode);
  PHNodeReset reset;
  if (mainIter.cd("BBC"))
    {
      mainIter.forEach(reset);
    }

  return 0;
}
