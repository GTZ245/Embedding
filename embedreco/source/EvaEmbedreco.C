
#include "EvaEmbedreco.h"
#include "PHEmbedMCEvaluator.hh"
#include "PHEmbededEvent.hh"
#include "PHEmbedHistogrammer.hh"
#include "Fun4AllReturnCodes.h"
#include "Fun4AllServer.h"

#include "recoConsts.h"

#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "PHTypedNodeIterator.h"
#include "PHNodeIterator.h"

#include "TObject.h"

#include <iostream>

using namespace std;

EvaEmbedreco::EvaEmbedreco(const char *name): SubsysReco(name){
  event     = new PHEmbededEvent;
  evaluator = new PHEmbedMCEvaluator;
}
EvaEmbedreco::~EvaEmbedreco(){
  //delete event; // deleted in dtor of Fun4AllServer
  delete evaluator;
}


int EvaEmbedreco::Init(PHCompositeNode *topNode)
{
  if (verbosity > 0)
    {
      cout << "Calling EvaEmbedreco::Init()" << endl;
    }
  // Create Histograms here - later you will have to do file magic
  // to make sure they are not deleted when the input file is closed
  return EVENT_OK;
}
int EvaEmbedreco::InitRun(PHCompositeNode *topNode){

  recoConsts *rc = recoConsts::instance();
  //set the output evaluation ntuple name
  PHEmbedHistogrammer::instance()->setFileName(rc->get_CharFlag("EMBED_CHARGED_EVAOUT")); 

  if (verbosity > 0){
    // this rc flag is set by the framework
    cout << "Calling EvaEmbedReco::InitRun() for Run " 
	 << rc->get_IntFlag("RUNNUMBER") << endl;
  }
  evaluator->setVerbose(rc->get_IntFlag("VERBOSITY"));

  Fun4AllServer* se = Fun4AllServer::instance();
  PHCompositeNode* mcnode = se->topNode(rc->get_CharFlag("EMBED_MC_TOPNODE"));
  PHCompositeNode* realnode = se->topNode(rc->get_CharFlag("EMBED_REAL_TOPNODE"));
  PHCompositeNode* mergednode = se->topNode("TOP");
  PHNodeIterator iter3(mergednode);  

  PHCompositeNode*embedNode = static_cast<PHCompositeNode*>(iter3.findFirst("PHCompositeNode","EMBED"));
  if(!embedNode){
    embedNode = new PHCompositeNode("EMBED");
    mergednode->addNode(embedNode);
  }
  PHDataNode<PHEmbededEvent>* embedEventNode = new PHDataNode<PHEmbededEvent>(event,"PHEmbededEvent");
  embedNode->addNode(embedEventNode);

  cout<<"MC node structure"<<endl;
  mcnode->print();
  cout<<"Real node structure"<<endl;
  realnode->print();
  cout<<"Merged node structure"<<endl;  
  mergednode->print();
  return EVENT_OK;
}

int EvaEmbedreco::process_event(PHCompositeNode *topNode){
  Fun4AllServer* se = Fun4AllServer::instance();

  recoConsts *rc = recoConsts::instance();
  PHCompositeNode* mcnode = se->topNode(rc->get_CharFlag("EMBED_MC_TOPNODE"));
  PHCompositeNode* realnode = se->topNode(rc->get_CharFlag("EMBED_REAL_TOPNODE"));
  PHCompositeNode* mergednode = topNode;
  event->event(mcnode, realnode, mergednode);

  evaluator->setVerbose(rc->get_IntFlag("VERBOSITY"));
  evaluator->event(mcnode, realnode, mergednode);
  //evaluator->setEventNum(currentEvtSingle);

  evaluator->associateDC(); // calls main contrib. fns. and fill fns.

  evaluator->associatePC(1);
  evaluator->associatePC(2);
  evaluator->associatePC(3);
  evaluator->associateTOF();

  evaluator->associateTOFW();
  //evaluator->associateTEC();

  evaluator->associateCRK();
  evaluator->associateEMC();
  // Can't find AccRaw anyway...
  //  evaluator->associateACC();
  static int counter = 0;
  counter ++;
  evaluator->associateMatching();
  evaluator->fillPHEmbedMcRecoTrack();
  PHEmbedHistogrammer* histo = PHEmbedHistogrammer::instance();
  if (counter % 50 == 1) {
    histo->flush();
  }

  return EVENT_OK;
}

int EvaEmbedreco::EndRun(const int runnumber){
  return EVENT_OK;
}
int EvaEmbedreco::End(PHCompositeNode *topNode){
  PHEmbedHistogrammer* histo = PHEmbedHistogrammer::instance();  
  histo->flush();
  histo->CloseOutputFile();
  return EVENT_OK;
}

int EvaEmbedreco::Reset(PHCompositeNode *topNode){
  return EVENT_OK;
}
int EvaEmbedreco::ResetEvent(PHCompositeNode *topNode)
{
  return EVENT_OK;
}
