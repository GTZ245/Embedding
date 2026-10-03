#include <string>
#include <vector>
#include "TObject.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "PHTypedNodeIterator.h"
#include "PHNodeIterator.h"
#include "recoConsts.h"
#include "DispEmbedreco.h"
#include "Fun4AllReturnCodes.h"
#include "Fun4AllServer.h"
#include "PhenixRun.hh"
#include "PhenixDisplay.hh"
#include "PhDchDisplay.hh"
#include "PhEventDisplay.hh"
#include "PhPadDisplay.hh"
DispEmbedreco::DispEmbedreco(const char *name){
  ThisName = name;
  flag =1;
}
DispEmbedreco::~DispEmbedreco(){
}


int DispEmbedreco::Init(PHCompositeNode *topNode)
{
  if (verbosity > 0)
    {
      cout << "Calling Init" << endl;
    }

  Fun4AllServer* se = Fun4AllServer::instance();
  //removes the tables that we are going to read in from DST file.
  //DetectorGeometry, VtxOut,CglTrack,CglTrackBack,PHTrackOut,PHTrackOutBack,PHDchTrackOut,AccRaw;
  PHCompositeNode* dstNode = se->getNode("DST","TOP");
  PHNodeIterator iter(dstNode);
  PHCompositeNode *geaNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "GEA"));
  if (!geaNode){
    geaNode = new PHCompositeNode("GEA");
    dstNode->addNode(geaNode);       
  }
  PHCompositeNode *evaNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVA"));
  if (!evaNode){
    evaNode = new PHCompositeNode("EVA");
    dstNode->addNode(evaNode);         
  }

  return EVENT_OK;
}
int DispEmbedreco::InitRun(PHCompositeNode *topNode){
  recoConsts *rc = recoConsts::instance();
  if (verbosity > 0){
    // this rc flag is set by the framework
    cout << "Calling InitRun for Run" << rc->get_IntFlag("RUNNUMBER") << endl;
  }
  Fun4AllServer* se = Fun4AllServer::instance();

  PHCompositeNode* mcnode = se->topNode(rc->get_CharFlag("EMBED_MC_TOPNODE"));
  PHCompositeNode* realnode = se->topNode(rc->get_CharFlag("EMBED_REAL_TOPNODE"));
  PHCompositeNode* mergednode = se->topNode("TOP");
  PhenixRun* phenix = new PhenixRun("Phenix-Run","Phenix-Run",mergednode);
  PhDchDisplay  *dchdisplay = new PhDchDisplay("Dch-Display","Dch-Display");
  gPhenix->GetDisplay()->AddSubsystem(dchdisplay);
  gPhenix->GetDisplay()->AddSubsystem(new PhPadDisplay("Pad-Display","Pad-Display"));
  //gPhenix->GetDisplay()->AddSubsystem(dchdisplay);
  //gPhenix->GetDisplay()->AddSubsystem(new PhPadDisplay("Pad-Display","Pad-Display"));
  Int_t TrackDraw_Mask   = 0x7;//BIT(0)|BIT(1)|BIT(2)
  Int_t HitBaseDraw_Mask = 0x3<<5;//BIT(5)|BIT(6)
  Int_t HitTypeDraw_Mask = 0x3<<11;//BIT(11)|BIT(12)|BIT(13)|BIT(14)|BIT(15)|BIT(16);
  dchdisplay->SetControlMode(11);
  dchdisplay->GetCurEventAttr()->SetTrackDrawMode(TrackDraw_Mask|HitBaseDraw_Mask|HitTypeDraw_Mask);
  gPhenix->StartRun();
  PhEventDisplay* dch = new PhEventDisplay("Main",gClient->GetRoot(),1000,1000,gPhenix->GetDisplay()->GetList());
  dch->Popup();
  //gPhenix->event(mergednode);
  return EVENT_OK;
}

int DispEmbedreco::process_event(PHCompositeNode *topNode){
  Fun4AllServer* se = Fun4AllServer::instance();

  recoConsts *rc = recoConsts::instance();
  PHCompositeNode* mcnode = se->topNode(rc->get_CharFlag("EMBED_MC_TOPNODE"));
  PHCompositeNode* realnode = se->topNode(rc->get_CharFlag("EMBED_REAL_TOPNODE"));
  PHCompositeNode* mergednode = topNode;
  cout<<mergednode<<endl;


  if(flag==1){
    gPhenix->event(mergednode);
  }else if(flag==2){ 
    gPhenix->event(mcnode);
  }else{
    gPhenix->event(realnode);
  }
  gPhenix->Draw();  
  return EVENT_OK;
}
int DispEmbedreco::Reset(PHCompositeNode *topNode){
  return EVENT_OK;
}
int DispEmbedreco::ResetEvent(PHCompositeNode *topNode)
{
  return EVENT_OK;
}
void DispEmbedreco::Print(const char *what) const{
}
