#include "MixEmbedreco.h"
#include "getClass.h"
#include "Fun4AllReturnCodes.h"
#include "Fun4AllServer.h"

#include "recoConsts.h"

#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "PHTypedNodeIterator.h"
#include "PHNodeIterator.h"
#include "PHGlobal.h"
#include "PHGlobalv11.h"

#include "TObject.h"
#include "BbcOutv1.h"
#include "dBbcOutWrapper.h"

#include "T0Out.h"
#include "T0Outv1.h"

#include "EventHeader.h"
#include "EventHeaderv2.h"

#include "VtxOutv7.h"

#include "RunHeader.h"
#include "RunHeaderv3.h"

#include "TrigLvl1v3.h"
#include "TrigRunLvl1v3.h"
#include "TrigRunLvl2v3.h"


#include <string>
#include <vector>

using namespace std;
typedef PHIODataNode<PHObject> PHObjectNode_t;
typedef PHIODataNode<PHTable> PHTableNode_t;

MixEmbedreco::MixEmbedreco(const char *name): SubsysReco(name){
}
MixEmbedreco::~MixEmbedreco(){
}

int MixEmbedreco::Init(PHCompositeNode *topNode)
{
  if (verbosity > 0)
    {
      cout << "MixEmbedReco: Calling Init" << endl;
    }

  //removes the tables that we are going to read in from DST file.
  //DetectorGeometry, VtxOut,CglTrack,CglTrackBack,PHTrackOut,PHTrackOutBack,PHDchTrackOut,AccRaw;
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode){
      dstNode = new PHCompositeNode("DST");
      topNode->addNode(dstNode);
  }



  PHCompositeNode* dstNode_h = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!dstNode_h){
      dstNode_h = new PHCompositeNode("RUN");
      topNode->addNode(dstNode_h);
  }

  PHGlobal *globalout = findNode::getClass<PHGlobal>(topNode, "PHGlobal");
  if(!globalout){
    globalout = new PHGlobalv11();
    PHObjectNode_t* GlobalOutNode = new PHObjectNode_t(globalout, "PHGlobal", "PHObject");
    dstNode->addNode(GlobalOutNode);
  }
  
  TrigLvl1 *triglvl1out = findNode::getClass<TrigLvl1>(topNode, "TrigLvl1");
  if(!triglvl1out){
    triglvl1out = new TrigLvl1v3();
    PHObjectNode_t* TrigLvl1OutNode = new PHObjectNode_t(triglvl1out, "TrigLvl1", "PHObject");
    dstNode->addNode(TrigLvl1OutNode);
  }

  TrigRunLvl1 *trigrunlvl1out = findNode::getClass<TrigRunLvl1>(topNode, "TrigRunLvl1");
  if(!trigrunlvl1out){
    trigrunlvl1out = new TrigRunLvl1v3();
    PHObjectNode_t* TrigRunLvl1OutNode = new PHObjectNode_t(trigrunlvl1out, "TrigRunLvl1", "PHObject");
    dstNode_h->addNode(TrigRunLvl1OutNode);
  }

  TrigRunLvl2 *trigrunlvl2out = findNode::getClass<TrigRunLvl2>(topNode, "TrigRunLvl2");
  if(!trigrunlvl2out){
    trigrunlvl2out = new TrigRunLvl2v3();
    PHObjectNode_t* TrigRunLvl2OutNode = new PHObjectNode_t(trigrunlvl2out, "TrigRunLvl2", "PHObject");
    dstNode_h->addNode(TrigRunLvl2OutNode);
  }
  

  //  ReactionPlaneObject* rpOut = findNode::getClass<ReactionPlaneObject>(topNode,"ReactionPlaneObject");
  //if(!rpOut){
  //rpOut = new ReactionPlaneObject("ReactionPlaneObject"
				  


  RunHeader *runout = findNode::getClass<RunHeader>(topNode,"RunHeader");
  if(!runout){
    runout = new RunHeaderv3();
    PHObjectNode_t* RunOutNode = new PHObjectNode_t(runout,"RunHeader","PHObject");
    dstNode_h->addNode(RunOutNode);
  }

  
  EventHeader *eventheader = findNode::getClass<EventHeader>(topNode,"EventHeader");
  if(!eventheader){
    eventheader = new EventHeaderv2();
    PHObjectNode_t* EventHeaderNode = new PHObjectNode_t(eventheader,"EventHeader","PHObject");
    dstNode->addNode(EventHeaderNode);
  }


  BbcOut* bbcout = findNode::getClass<BbcOut>(topNode, "BbcOut");
  if(!bbcout){
    bbcout = new BbcOutv1();
    PHObjectNode_t* BbcOutNode = new PHObjectNode_t(bbcout, "BbcOut", "PHObject");
    dstNode->addNode(BbcOutNode);
  }
  int  mr = 1;
  dBbcOutWrapper* dBbcOut = findNode::getClass<dBbcOutWrapper>(topNode, "dBbcOut");
  if(!dBbcOut){
    dBbcOut= new dBbcOutWrapper("dBbcOut", mr);
    PHTableNode_t* dBbcOutNode =new PHTableNode_t(dBbcOut, "dBbcOut");
    dstNode->addNode(dBbcOutNode);
  }

  VtxOut *vtxout = findNode::getClass<VtxOut>(topNode, "VtxOut");
  if(!vtxout){
    vtxout = new VtxOutv7();
    PHObjectNode_t* VtxOutNode = new PHObjectNode_t(vtxout, "VtxOut","PHObject");
    dstNode->addNode(VtxOutNode);
  }
 
  
  T0Out* t0out = findNode::getClass<T0Out>(topNode,"T0Out");
  if(!t0out){
    t0out = new T0Outv1();
    PHObjectNode_t *T0OutNode = new PHObjectNode_t(t0out,"T0Out","PHObject");
    dstNode->addNode(T0OutNode);
  }
  
  if (verbosity > 0)
    {
      cout << "MixEmbedReco: printing topNode" << endl;
    }
  topNode->print();
  return EVENT_OK;
}
int MixEmbedreco::InitRun(PHCompositeNode *topNode){
  return EVENT_OK;
}

int MixEmbedreco::process_event(PHCompositeNode *topNode){

  
  PHCompositeNode *m_realnode;

  Fun4AllServer *se = Fun4AllServer::instance();

  m_realnode = se->topNode("REAL");


  PHGlobal *global_real = findNode::getClass<PHGlobal>(m_realnode,"PHGlobal");
  PHGlobal *global = findNode::getClass<PHGlobal>(topNode,"PHGlobal");

  if(global_real)
    *global = *global_real;
  else
    cout << "hold up" << endl;

  
  //cout << endl;
  //cout << global_real->getBbcZVertex() << endl;

  global->setBbcZVertex(global_real->getBbcZVertex());
  global->setBbcZVertexError(global_real->getBbcZVertexError());

  //cout << global->getBbcZVertex() << endl;



  // copy over runheader from real data node
  RunHeader *head_real = findNode::getClass<RunHeader>(m_realnode,"RunHeader");
  RunHeader *head = findNode::getClass<RunHeader>(topNode,"RunHeader");

  if(head_real)
    *head = *head_real;
  else
    cout << PHWHERE << "hold up" << endl;

  head->set_RunNumber(head_real->get_RunNumber());
  head->set_TimeStart(head_real->get_TimeStart());
  head->set_TimeStop(head_real->get_TimeStop());
  head->set_currentNorth(head_real->get_currentNorth());
  head->set_currentSouth(head_real->get_currentSouth());
  head->set_currentCentral(head_real->get_currentCentral());
  head->set_currentInner(head_real->get_currentInner());

  //cout << PHWHERE << " runumber = " << head->get_RunNumber() << endl;
  //cout << PHWHERE << " runumber_real = " << head_real->get_RunNumber() << endl;
  //cout << PHWHERE << " timestart = " << head->get_TimeStart() << endl;
  //cout << PHWHERE << " timestart_real = " << head_real->get_TimeStart() << endl;


  // copy over EventHeader from real data node
  EventHeader *ehead_real = findNode::getClass<EventHeader>(m_realnode,"EventHeader");
  EventHeader *ehead = findNode::getClass<EventHeader>(topNode,"EventHeader");

  if(ehead_real)
    *ehead = *ehead_real;
  else
    cout << "hold up...no existing EventHeader in the topNode" << endl;
  
  ehead->set_EvtSequence(ehead_real->get_EvtSequence());
  ehead->set_EvtType(ehead_real->get_EvtType());

  ehead->set_TimeStamp(ehead_real->get_TimeStamp());

  TrigLvl1 *trig_real = findNode::getClass<TrigLvl1>(m_realnode,"TrigLvl1");
  TrigLvl1 *trig = findNode::getClass<TrigLvl1>(topNode,"TrigLvl1");

  //  *trig = *trig_real;
  
  //trig = trig_real->clone();

  trig->set_lvl1_trigraw(trig_real->get_lvl1_trigraw());  
  trig->set_lvl1_triglive(trig_real->get_lvl1_triglive());
  trig->set_lvl1_trigscaled(trig_real->get_lvl1_trigscaled());
  trig->set_lvl1_clock_cross(trig_real->get_lvl1_clock_cross());
  // source code says there are 130 bits
  //  but in DumpTrigLvl1, its seems to only have 5 values
  for(int i=0; i<5; i++)
    trig->set_lvl1_rbits(i, trig_real->get_lvl1_rbits(i));
  
  //trig->set_lvl1_beam_clk(trig_real->get_lvl1_beam_clk());
  

  TrigRunLvl1 *trigrun_real = findNode::getClass<TrigRunLvl1>(m_realnode,"TrigRunLvl1");
  TrigRunLvl1 *trigrun = findNode::getClass<TrigRunLvl1>(topNode,"TrigRunLvl1");

  *trigrun = *trigrun_real;

  //trigrun = (TrigRunLvl1v3*)trigrun_real->clone();
  
  trigrun->set_lvl1_trigger_version(trigrun_real->get_lvl1_trigger_version());
  trigrun->set_lvl1_trigger_description(trigrun_real->get_lvl1_trigger_description());
  trigrun->set_lvl1_bbcll1_description(trigrun_real->get_lvl1_bbcll1_description());
  trigrun->set_lvl1_bbcll1_version(trigrun_real->get_lvl1_bbcll1_version());
  trigrun->set_lvl1_partition_name(trigrun_real->get_lvl1_partition_name());
  
  for(int ilvl1=0; ilvl1<32; ilvl1++)
    {
      trigrun->set_lvl1_trig_name(trigrun_real->get_lvl1_trig_name(ilvl1),ilvl1);//
      trigrun->set_lvl1_rbit_name(trigrun_real->get_lvl1_rbit_name(ilvl1),ilvl1);//
      trigrun->set_lvl1_trigger_enable(trigrun_real->get_lvl1_trigger_enable(ilvl1),ilvl1);//
      trigrun->set_lvl1_trig_bit(trigrun_real->get_lvl1_trig_bit(ilvl1),ilvl1);//
      ////      trigrun->set_lvl1_trig_index(trigrun_real->get_lvl1_trig_index(ilvl1),ilvl1);
      trigrun->set_lvl1_trig_bitmask(trigrun_real->get_lvl1_trig_bitmask(ilvl1),ilvl1);
      trigrun->set_lvl1_trig_scale_down(trigrun_real->get_lvl1_trig_scale_down(ilvl1),ilvl1);//
      trigrun->set_lvl1_lvl2_reject_enable(trigrun_real->get_lvl1_lvl2_reject_enable(ilvl1),ilvl1);//
      trigrun->set_lvl1_trig_rate_begin(trigrun_real->get_lvl1_trig_rate_begin(ilvl1),ilvl1);
      
    }
  
  trigrun->set_start_time(trigrun_real->get_start_time());
  trigrun->set_run_number(trigrun_real->get_run_number());

  


  TrigRunLvl2 *trigrun2_real = findNode::getClass<TrigRunLvl2>(m_realnode,"TrigRunLvl2");
  TrigRunLvl2 *trigrun2 = findNode::getClass<TrigRunLvl2>(topNode,"TrigRunLvl2");

  *trigrun2 = *trigrun2_real;

  //  trigrun2 = (TrigRunLvl2v3*)trigrun2_real->clone();
  trigrun2->set_lvl2_version(trigrun2_real->get_lvl2_version());
  trigrun2->set_lvl2_description(trigrun2_real->get_lvl2_description());
  trigrun2->set_lvl2_run_enable(trigrun2_real->get_lvl2_run_enable());
  trigrun2->set_lvl2_reject_enable(trigrun2_real->get_lvl2_reject_enable());
  
  for(int ilvl2=0; ilvl2<32; ilvl2++)
    {
      trigrun2->set_lvl1_lvl2_reject_enable(trigrun2_real->get_lvl1_lvl2_reject_enable(ilvl2),ilvl2);
      trigrun2->set_lvl1_lvl2_force_accept(trigrun2_real->get_lvl1_lvl2_force_accept(ilvl2),ilvl2);
      trigrun2->set_lvl2_trig_name(trigrun2_real->get_lvl2_trig_name(ilvl2),ilvl2);
      trigrun2->set_lvl2_trig_bit(trigrun2_real->get_lvl2_trig_bit(ilvl2),ilvl2);
      trigrun2->set_lvl2_trig_version(trigrun2_real->get_lvl2_trig_version(ilvl2),ilvl2);
      for(int ilvl1=0; ilvl1<32; ilvl1++)
	{
	  trigrun2->set_lvl2_lvl1_assoc(trigrun2_real->get_lvl2_lvl1_assoc(ilvl2,ilvl1),ilvl2,ilvl1);
	  trigrun2->set_lvl2_lvl1_prescale(trigrun2_real->get_lvl2_lvl1_prescale(ilvl2,ilvl1),ilvl2,ilvl1);
	}
    }
  
  BbcOut *bbcout = findNode::getClass<BbcOut>(topNode,"BbcOut");
  BbcOut *bbcout_real = findNode::getClass<BbcOut>(m_realnode,"BbcOut");

  bbcout->set_TimeZero(bbcout_real->get_TimeZero());
  bbcout->set_Vertex(bbcout_real->get_VertexPoint(), bbcout_real->get_dVertexPoint());



  VtxOut *vtxout = findNode::getClass<VtxOut>(topNode,"VtxOut");

  float vtx[3], vtxerr[3];
  vtx[2] = bbcout->get_VertexPoint();
  vtx[0] = 0.;
  vtx[1] = 0.;
  
  vtxerr[2] = bbcout->get_dVertexPoint();
  vtxerr[0] = 0.;
  vtxerr[1] = 0.;

  vtxout->AddVtx("BBC",vtx,vtxerr,VTX::BBCORDER);


  
  //cout << PHWHERE << "bbc time zero = " << bbcout->get_TimeZero() << endl;
  //cout << PHWHERE << "bbc_real time zero = " << bbcout_real->get_TimeZero() << endl;

  
  T0Out *t0Out = findNode::getClass<T0Out>(topNode,"T0Out");
  t0Out->set_BbcT0(bbcout_real->get_TimeZero(), bbcout_real->get_dTimeZero());
  //  t0Out->set_NtcT0(t0Out_real->get_NtcT0());
  // t0Out->set_TzrT0(t0Out_real->get_TzrT0());
  //t0Out->set_ZdcT0(t0Out_real->get_ZdcT0());

  return EVENT_OK;
}
int MixEmbedreco::Reset(PHCompositeNode *topNode){
  return EVENT_OK;
}
int MixEmbedreco::ResetEvent(PHCompositeNode *topNode)
{
  return EVENT_OK;
}
