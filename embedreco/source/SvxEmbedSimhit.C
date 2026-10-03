#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

#include <SvxEmbedSimhit.h>

#include <SvxGhitListv1.h>
#include <SvxGhitv1.h>
#include <SvxRawhitListv4.h>
#include <SvxRawhitv4.h>
#include <SvxGhitRawhitListv1.h>
#include <SvxGhitRawhitv1.h>
#include <RunHeader.h>
#include <RunHeaderv3.h>
#include <EventHeader.h>
#include <EventHeaderv2.h>
#include <PreviousEvent.h>
#include <PreviousEventv1.h>
#include <TrigLvl1.h>
#include <TrigLvl1v3.h>
#include <TrigRunLvl1.h>
#include <TrigRunLvl1v3.h>
#include <BbcOut.h>
#include <BbcOutv1.h>
#include <PHGlobal.h>
#include <PHGlobalv11.h>
#include <PHCentralTrack.h>
#include <McEvalSingleList.h>
#include <VtxOutv7.h>
#include <PHPoint.h>
#include <svxAddress.hh>
#include <svxDetectorGeo.hh>
#include <SvxSensor.h>
#include <SvxParameters.h>

#include <Fun4AllReturnCodes.h>

#include <Fun4AllServer.h>
#include <recoConsts.h>
#include <PHNodeIterator.h>
#include <PHCompositeNode.h>
#include <PHIODataNode.h>
#include <getClass.h>

using namespace std;

SvxEmbedSimhit::SvxEmbedSimhit(const string &name) : SubsysReco(name)
{
  m_ievt = 0;

  /// following values are the same as those in SvxSimulator
  m_stripixel_sAQ                = 0.109; // based on beam test
  m_stripixel_sNOISE             = 10.2; // based on measurement (2008 June)
  m_stripixel_adcthre_zs         = 31; // 3 sigma of sNOISE
  m_stripixel_adcthre_rawhit     = 31; // 3 sigma of sNOISE
  m_stripixel_adcthre_rawhit_sum = 50; // a speculated value

  m_max_ladder[0] = SVXLADDERSLAYER0*2;
  m_max_ladder[1] = SVXLADDERSLAYER1*2;
  m_max_ladder[2] = SVXLADDERSLAYER2*2;
  m_max_ladder[3] = SVXLADDERSLAYER3*2;
  m_max_sensor[0] = SVXSENSORSLAYER0;
  m_max_sensor[1] = SVXSENSORSLAYER1;
  m_max_sensor[2] = SVXSENSORSLAYER2;
  m_max_sensor[3] = SVXSENSORSLAYER3;

  m_sensorZwhalf[0] = 2.836;
  m_sensorZwhalf[1] = 2.836;
  m_sensorZwhalf[2] = 3.1877;
  m_sensorZwhalf[3] = 3.1877;

  m_wshift = 0.;
  m_rndm = false;
  rndmx.SetSeed(0);
  rndmy.SetSeed(0);

  m_mcnode = NULL;
  m_realnode = NULL;
  m_svxgeo = NULL;

  mc_zvtx_mean = 0.;

  m_copyPHCentralTrack = true; 
  m_copyMcSingle       = true; 

  m_new_runhdr    = false;
  m_new_evthdr    = false;
  m_new_prev      = false;
  m_new_trglvl1   = false;
  m_new_trgrunlvl1= false;
  m_new_bbcout    = false;
  m_new_phglobal  = false;
}

SvxEmbedSimhit::~SvxEmbedSimhit() {}

//----------------------------------------------------------------------------------

void SvxEmbedSimhit::set_StripixelNoise(float a) 
{
  m_stripixel_sNOISE = a;
}

//----------------------------------------------------------------------------------

int SvxEmbedSimhit::Init(PHCompositeNode *topNode)
{
  if ( verbosity>0 ) {
    cout << PHWHERE << "SvxEmbedSimhit::Init started..." << endl;
  }

  m_svxgeo = findNode::getClass<svxDetectorGeo>(topNode, "svxDetectorGeo");
  if ( m_svxgeo == NULL) {
    if(verbosity>0) { cout << PHWHERE<< "Can't find svxDetectorGeo. " << endl; }
    return ABORTRUN;
  }
  SvxSensor *svxsen;
  for ( int ilr=2; ilr<4; ilr++)  {
    for ( int ild=0; ild<m_max_ladder[ilr]; ild++) {
      for ( int isn=0; isn<m_max_sensor[ilr]; isn++) {
        svxsen = m_svxgeo->GetSensorPtr(ilr, ild, isn);
	svxsen->set_sAQ(m_stripixel_sAQ);
	svxsen->set_sNOISE(m_stripixel_sNOISE);
	svxsen->set_adcthre_zs(m_stripixel_adcthre_zs);
	svxsen->set_adcthre_rawhit(m_stripixel_adcthre_rawhit);
	svxsen->set_adcthre_rawhit_sum(m_stripixel_adcthre_rawhit_sum);
      }
    }
  }

  return EVENT_OK;
}

//----------------------------------------------------------------------------------

int SvxEmbedSimhit::InitRun(PHCompositeNode *topNode)
{
  if ( verbosity>0 ) {
    cout << PHWHERE << "SvxEmbedSimhit::InitRun started..." << endl;
  }

  ///
  /// get nodes for real data and simulation which you want to embed
  ///
  Fun4AllServer* se = Fun4AllServer::instance();
  recoConsts* rc = recoConsts::instance();
  /// simulation
  if ( rc->FlagExist("EMBED_MC_TOPNODE") ) {
    m_mcnode = se->topNode(rc->get_CharFlag("EMBED_MC_TOPNODE"));
  } else {
    cout << PHWHERE << "EMBED_MC_TOPNODE char flag should be set" << endl;
    return EVENT_OK;
  }
  /// real data
  if ( rc->FlagExist("EMBED_REAL_TOPNODE") ) {
    m_realnode = se->topNode(rc->get_CharFlag("EMBED_REAL_TOPNODE"));
  } else {
    cout << PHWHERE << "EMBED_REAL_TOPNODE char flag should be set" << endl;
    return EVENT_OK;
  }

  if ( verbosity>0 && m_ievt==0 ) {
    m_mcnode->print();
    m_realnode->print();
    topNode->print();
  }

  int done = CreateNodeTree(topNode);
  if ( done==0 ) return ABORTEVENT;

  cout<<(m_copyPHCentralTrack ? "copy" : "not copy")<<" PHCentralTrack"<<endl;
  cout<<(m_copyMcSingle       ? "copy" : "not copy")<<" McSingle"      <<endl;

  return EVENT_OK;
}

//----------------------------------------------------------------------------------

int SvxEmbedSimhit::process_event(PHCompositeNode *topNode)
{
  if ( verbosity>0 ) {
    cout << PHWHERE << "SvxEmbedSimhit::process_event started..." << endl;
  }

  // copy prev
  if(m_new_runhdr){     copyRunHeader(topNode);     }
  if(m_new_evthdr){     copyEventHeader(topNode);   }
  if(m_new_prev)  {     copyPreviousEvent(topNode); }
  if(m_new_trglvl1){    copyTrigLvl1(topNode);      }
  if(m_new_trgrunlvl1){ copyTrigRunLvl1(topNode);   }
  if(m_new_bbcout){     copyBbcOut(topNode);        }
  if(m_new_phglobal){   copyPHGlobal(topNode);      }

  ///
  /// get primary vertex information
  /// 
  VtxOut *vtxout      = findNode::getClass<VtxOut>(topNode,   "VtxOut");
  VtxOut *mc_vtxout   = findNode::getClass<VtxOut>(m_mcnode,  "VtxOut");
  VtxOut *real_vtxout = findNode::getClass<VtxOut>(m_realnode,"VtxOut");
  if ( !vtxout || !mc_vtxout || !real_vtxout ) {
    cerr << "VtxOut node not found!" << endl;
    return EVENT_OK;
  }
  ///
  /// following vertices are stored in output node.
  /// - event vertex in simulation
  ///   + stored as "SIM"
  /// - reconstructed vertex in real data
  ///   + stored as "SVX_PRECISE"
  /// - reconstructed beam center in real data
  ///   + stored as "SVX"
  /// - combination of simulation & real data
  ///   + (x,y)=simulation,
  ///   + z=real data
  ///   + stored as "FORCED"
  ///   + this vertex is used to calculate DCA
  ///
  //cout << "----------------------------------------------------------------------" << endl;
  PHPoint mc_vtx     = mc_vtxout->get_Vertex("SIM"); 
  PHPoint real_vtx   = real_vtxout->get_Vertex();
  PHPoint real_vtxbc = real_vtxout->get_Vertex("SVX");
  float mc_vtx_xyz[3] = {(float)mc_vtx.getX(), (float)mc_vtx.getY(), (float)mc_vtx.getZ()};  // NO shift
  //float mc_vtx_xyz[3] = {mc_vtx.getX(), mc_vtx.getY(), real_vtx.getZ()+mc_vtx.getZ()}; 
  //float mc_vtx_xyz[3] = {0., 0., (float) real_vtx.getZ() }; 
  float mc_vtx_xyz_err[3] = {0.01, 0.01, 0.01};
  vtxout->AddVtx("SIM", mc_vtx_xyz, mc_vtx_xyz_err, VTX::SIMORDER);
  //cout << "===============================================" << endl;
  //cout << "LEBEDEV SIM VTX: " << mc_vtx_xyz[0] << " " << mc_vtx_xyz[1] << " " << mc_vtx_xyz[2] << endl;
  //float real_vtx_xyz[3] = {real_vtx.getX(), real_vtx.getY(), real_vtx.getZ()};
  //cout << "LEBEDEV REAL VTX: " << real_vtx_xyz[0] << " " << real_vtx_xyz[1] << " " << real_vtx_xyz[2] << endl;
  //float real_vtx_xyz_err[3] = {0.01, 0.01, 0.01};
  //vtxout->AddVtx("SVX_PRECISE", real_vtx_xyz, real_vtx_xyz_err, VTX::SVX_PRECISEORDER);
  //float real_vtxbc_xyz[3] = {real_vtxbc.getX(), real_vtxbc.getY(), real_vtxbc.getZ()};
  //vtxout->AddVtx("SVX", real_vtxbc_xyz, real_vtx_xyz_err, VTX::SVXORDER);
  //float vtx_xyz[3] = {mc_vtx.getX(), mc_vtx.getY(), real_vtx.getZ()};
  //if ( m_rndm ) {
  //  vtx_xyz[0] += rndmx.Uniform(-m_wshift, m_wshift);
  //  vtx_xyz[1] += rndmy.Uniform(-m_wshift, m_wshift);
 //}
  //float vtx_xyz_err[3] = {0.01, 0.01, 0.01};
//  vtxout->AddVtx("FORCED", vtx_xyz, vtx_xyz_err, VTX::FORCEDORDER);


//  if ( verbosity>0 ) {
//    cout << "forced vertex : " << m_ievt <<" | "<< m_wshift <<", "
//	 << vtx_xyz[0] <<" "<< vtx_xyz[1] <<" "<< vtx_xyz[2] << endl;
//  }
  /// vertex difference in z-direction is used to shift rawhit in simulation.
  float diff_vtxz = real_vtx.getZ()*0 - mc_zvtx_mean;
  //if(fabs(diff_vtxz)>0.5) return ABORTEVENT;

//  diff_vtxz = 0.; // don't shift
//  cout << "LEBEDEV WARNING: NO SHIFT!!!" << endl;

  ///
  /// start embedding
  ///
  int nraw_embed = 0*diff_vtxz;//embed_simhit(topNode, diff_vtxz);

  ////////////////////////////////////////////////////////
  /////Takashi and Sasha used this for embed sim//////////
  /////It converts convert ghit to rawhit/////////////////
  /////My code below just merge rawhit in mc and real/////
  /////Choose your fighter////////////////////////////////
  /////Mine produces a lot of replicas////////////////////
  /////Yuri Mitrankov mitrankovy@gmail.com////////////////
  ////////////////////////////////////////////////////////
  int Takashi_Sasha_embed = 1;
  int Yura_embed = 0;
  if(Takashi_Sasha_embed) nraw_embed = embed_simhit(topNode, diff_vtxz);

  ///
  /// save rawhits of real data
  ///
  SvxRawhitList *rawhitlist  = findNode::getClass<SvxRawhitList>(topNode,   "SvxRawhitList");
  SvxRawhitList *real_rawhit = findNode::getClass<SvxRawhitList>(m_realnode,"SvxRawhitList");
  SvxRawhitList *mc_rawhit   = findNode::getClass<SvxRawhitList>(m_mcnode,  "SvxRawhitList");
  SvxGhitRawhitList *g2rlist = findNode::getClass<SvxGhitRawhitList>(topNode,"SvxGhitRawhitList");
  if ( !rawhitlist ) {
    cerr << "SvxRawhitList node (embed) not found : " << m_ievt << endl;
    return EVENT_OK;
  }
  if ( !real_rawhit ) {
    cerr << "SvxRawhitList node (real) not found : " << m_ievt << endl;
    return EVENT_OK;
  }
  if ( !mc_rawhit ) {
    cerr << "SvxRawhitList node (mc) not found : " << m_ievt << endl;
    return EVENT_OK;
  }
  if ( !g2rlist ) {
    cerr << "SvxGhitRawhitList node (geant) not found : " << m_ievt << endl;
    return EVENT_OK;
  }

//  real_rawhit->Reset();  // embed in empty event

  int nrealhit = real_rawhit->get_nRawhits();
  for ( int i=0; i<nrealhit; i++ ) {
    SvxRawhit *realhit = real_rawhit->get_Rawhit(i);
    SvxRawhit *tmphit  = rawhitlist->addRawhit();
    tmphit->set_layer(        realhit->get_layer());
    tmphit->set_ladder(       realhit->get_ladder());
    tmphit->set_sensor(       realhit->get_sensor());
    tmphit->set_sensorSection(realhit->get_sensorSection());
    tmphit->set_channel(      realhit->get_channel());
    tmphit->set_sensorType(   realhit->get_sensorType());
    tmphit->set_sensorReadout(realhit->get_sensorReadout());
    tmphit->set_adc(          realhit->get_adc());
    tmphit->set_pixelROC(     realhit->get_pixelROC());
    tmphit->set_pixelModule(  realhit->get_pixelModule());
  }
  int nsimhit = mc_rawhit->get_nRawhits()*Yura_embed;
  int k = 0;
  for ( int i=0; i<nsimhit; i++ ) {
    SvxRawhit *simhit = mc_rawhit->get_Rawhit(i);
    int adc = simhit->get_adc();//(int) floor(0.0019*static_cast<float>(simhit->get_adc()));
    //std::cout<<simhit->get_adc()<<" "<<adc<<std::endl;
    if(adc<1)///adcthreshold is 1; adcgain = 0.0019 counts/electron
    {
      continue;
    }
    if(adc>255) continue;//adc = 255; ///adctop value is 255
    SvxRawhit *tmphit  = rawhitlist->addRawhit();
    tmphit->set_layer(        simhit->get_layer());
    tmphit->set_ladder(       simhit->get_ladder());
    tmphit->set_sensor(       simhit->get_sensor());
    tmphit->set_sensorSection(simhit->get_sensorSection());
    tmphit->set_channel(      simhit->get_channel());
    tmphit->set_sensorType(   simhit->get_sensorType());
    tmphit->set_sensorReadout(simhit->get_sensorReadout());
    tmphit->set_adc(          adc);
    tmphit->set_pixelROC(     simhit->get_pixelROC());
    tmphit->set_pixelModule(  simhit->get_pixelModule());
    SvxGhitRawhit *ghit  = g2rlist->addGhitRawhit();
    ghit->set_ghitID(k);
    ghit->set_rawhitID(tmphit->get_hitID());
    k++;
  }

  if ( verbosity>0 ) {
    cout << "number of rawhit : " 
	 << rawhitlist->get_nRawhits() <<" = "
	 << nraw_embed <<" + "<< nrealhit << endl;
  }

  m_ievt++;
  return EVENT_OK;
}

//----------------------------------------------------------------------------------

int SvxEmbedSimhit::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = NULL;
  dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","DST"));
  if ( !dstNode ) {
    cerr << PHWHERE << " DST node missing, doing nothing." << endl;
    return 0;
  }
  PHCompositeNode *svxNode = NULL;
  svxNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","SVX"));
  if ( !svxNode ) {
    svxNode = new PHCompositeNode("SVX");
    dstNode->addNode(svxNode);
  }

  /// SvxGhitList
  PHIODataNode<PHObject> *ghitNode = NULL;
  ghitNode = (PHIODataNode<PHObject>*)iter.findFirst("PHIODataNode","SvxGhitList");
  if ( !ghitNode ) {
    SvxGhitList *ghitlist = new SvxGhitListv1();
    ghitNode = new PHIODataNode<PHObject>(ghitlist,"SvxGhitList","PHObject");
    svxNode->addNode(ghitNode);
  }

  /// SvxRawhitList
  PHIODataNode<PHObject> *rawhitNode = NULL;
  rawhitNode = (PHIODataNode<PHObject>*)iter.findFirst("PHIODataNode","SvxRawhitList");
  if ( !rawhitNode ) {
    SvxRawhitList *rawhitlist = new SvxRawhitListv4();
    rawhitNode = new PHIODataNode<PHObject>(rawhitlist,"SvxRawhitList","PHObject");
    svxNode->addNode(rawhitNode);
  }

  /// SvxGhitRawhitList
  PHIODataNode<PHObject> *g2rNode = NULL;
  g2rNode = (PHIODataNode<PHObject>*)iter.findFirst("PHIODataNode","SvxGhitRawhitList");
  if ( !g2rNode ) {
    SvxGhitRawhitList *g2rlist = new SvxGhitRawhitListv1();
    g2rNode = new PHIODataNode<PHObject>(g2rlist,"SvxGhitRawhitList","PHObject");
    svxNode->addNode(g2rNode);
  }

  /// VtxOut
  PHIODataNode<PHObject> *vtxoutNode = NULL;
  vtxoutNode = (PHIODataNode<PHObject>*)iter.findFirst("PHIODataNode","VtxOut");
  if ( !vtxoutNode ) {
    VtxOut *vtxout = new VtxOutv7();
    vtxoutNode = new PHIODataNode<PHObject>(vtxout,"VtxOut","PHObject");
    dstNode->addNode(vtxoutNode);
    cout<<"SvxEmbedSimhit : no VtxOut node. newly added"<<endl;
  }

  /// RunHeader, EventHeader, PreviouseEvent, TrigLvl1, BbcOut, and PHGlobal nodes
  /// of real data are saved in topNode
  /// RunHeader
  RunHeader *run = findNode::getClass<RunHeader>(m_realnode,"RunHeader");
  if ( !run ) {
    cerr << "RunHeader can not be found in node tree of real data" << endl;
    return 0;
  }
  RunHeader *run_current = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if(run_current==NULL){
    run = new RunHeaderv3();
    PHIODataNode<PHObject> *runNode = new PHIODataNode<PHObject>(run,"RunHeader","PHObject");
    dstNode->addNode(runNode);
    m_new_runhdr = true;
  } else {
    cout<<"SvxEmbedSimHit : RunHeader exists"<<endl;
  }

  /// EventHeader
  EventHeader *event = findNode::getClass<EventHeader>(m_realnode,"EventHeader");
  if ( !event ) {
    cerr << "EventHeader can not be found in node tree of real data" << endl;
    return 0;
  }
  EventHeader *event_current = findNode::getClass<EventHeader>(topNode,"EventHeader");
  if(event_current==NULL){
    event = new EventHeaderv2();
    PHIODataNode<PHObject> *eventNode = new PHIODataNode<PHObject>(event,"EventHeader","PHObject");
    dstNode->addNode(eventNode);
    m_new_evthdr = true;
  } else {
    cout<<"SvxEmbedSimHit : EventHeader exists"<<endl;
  }

  /// PreviousEvent
  PreviousEvent *prev = findNode::getClass<PreviousEvent>(m_realnode,"PreviousEvent");
  if ( !prev ) {
    cerr << "PreviousEvent can not be found in node tree of real data" << endl;
    return 0;
  }
  PreviousEvent *prev_current = findNode::getClass<PreviousEvent>(topNode,"PreviousEvent");
  if(prev_current==NULL){
    prev = new PreviousEventv1();
    PHIODataNode<PHObject> *prevNode = new PHIODataNode<PHObject>(prev,"PreviousEvent","PHObject");
    dstNode->addNode(prevNode);
    m_new_prev = true;
    cout<<"SvxEmbedSimHit : PreviousEvent newly added"<<endl;
  } else {
    cout<<"SvxEmbedSimHit : PreviousEvent exists"<<endl;
  }

  /// TrigLvl1
  TrigLvl1 *trig = findNode::getClass<TrigLvl1>(m_realnode,"TrigLvl1");
  if ( !trig ) {
    cerr << "TrigLvl1 can not be found in node tree of real data" << endl;
    return 0;
  }
  TrigLvl1 *trig_current = findNode::getClass<TrigLvl1>(topNode,"TrigLvl1");
  if(trig_current==NULL){
    trig = new TrigLvl1v3();
    PHIODataNode<PHObject> *trigNode = new PHIODataNode<PHObject>(trig,"TrigLvl1","PHObject");
    dstNode->addNode(trigNode);
    m_new_trglvl1 = true;
  } else {
    cout<<"SvxEmbedSimHit : TrigLvl1 exists"<<endl;
  }

  /// TrigRunLvl1
  TrigRunLvl1 *trigrun = findNode::getClass<TrigRunLvl1>(m_realnode,"TrigRunLvl1");
  if ( !trigrun ) {
    cerr << "TrigRunLvl1 can not be found in node tree of real data" << endl;
    return 0;
  }
  TrigRunLvl1 *trigrun_current = findNode::getClass<TrigRunLvl1>(topNode,"TrigRunLvl1");
  if(trigrun_current==NULL){
    trigrun = new TrigRunLvl1v3();
    PHIODataNode<PHObject> *trigrunNode = new PHIODataNode<PHObject>(trigrun,"TrigRunLvl1","PHObject");
    dstNode->addNode(trigrunNode);
    m_new_trgrunlvl1 = true;
  } else {
    cout<<"SvxEmbedSimHit : TrigRunLvl1 exists"<<endl;
  }

  /// BbcOut
  BbcOut *bbcout = findNode::getClass<BbcOut>(m_realnode,"BbcOut");
  if ( !bbcout ) {
    cerr << "BbcOut can not be found in node tree of real data" << endl;
    return 0;
  }
  BbcOut *bbcout_current = findNode::getClass<BbcOut>(topNode,"BbcOut");
  if(bbcout_current==NULL){
    bbcout = new BbcOutv1();
    PHIODataNode<PHObject> *bbcoutNode = new PHIODataNode<PHObject>(bbcout,"BbcOut","PHObject");
    dstNode->addNode(bbcoutNode);
    m_new_bbcout = true;
  } else {
    cout<<"SvxEmbedSimHit : BbcOut exists"<<endl;
  }

  /// PHGlobal
  PHGlobal *global = findNode::getClass<PHGlobal>(m_realnode,"PHGlobal");
  if ( !global ) {
    cerr << "PHGlobal can not be found in node tree of real data" << endl;
    return 0;
  }
  PHGlobal *global_current = findNode::getClass<PHGlobal>(topNode,"PHGlobal");
  if(global_current==NULL){
    global = new PHGlobalv11();
    PHIODataNode<PHObject> *globalNode = new PHIODataNode<PHObject>(global,"PHGlobal","PHObject");
    dstNode->addNode(globalNode);
    m_new_phglobal = true;
  } else {
    cout<<"SvxEmbedSimHit : PHGlobal exists"<<endl;
  }
  
  /// McSingle & PHCentralTrack nodes of simulation are saved in topNode
  /// McSingle
  // Do not _require_ McSingle and PHCentralTrack for embedding so that
  // svx embedding can be called even if the central arm reconstruction is not.

  if(m_copyMcSingle){
    PHIODataNode<PHObject> *mcNode = NULL;
    McEvalSingleList *mctrk = findNode::getClass<McEvalSingleList>(m_mcnode, "McSingle");
    if ( !mctrk ) {
      cerr << "McSingle can not be found in node tree of simulation" << endl;
    }
    else{
      mcNode = new PHIODataNode<PHObject>(mctrk,"McSingle","PHObject");
      dstNode->addNode(mcNode);
    }
  }
  
  if(m_copyPHCentralTrack){
    /// PHCentralTrack
    PHIODataNode<PHObject> *cntNode = NULL;
    PHCentralTrack *cnt = findNode::getClass<PHCentralTrack>(m_mcnode,"PHCentralTrack");
    if ( !cnt ) {
      cerr << "PHCentralTrack can not be found in node tree of simulation" << endl;
    }
    else{
      cntNode = new PHIODataNode<PHObject>(cnt,"PHCentralTrack","PHObject");
      dstNode->addNode(cntNode);
    }
  }

  return 1;
}

//----------------------------------------------------------------------------------

int SvxEmbedSimhit::embed_simhit(PHCompositeNode *topNode, float zshift)
{
  if ( verbosity>0 ) {
    cout << "SvxEmbedSimhit : start embedding..." <<endl;
  }

  SvxGhitList *ghitlist = findNode::getClass<SvxGhitList>(topNode,"SvxGhitList");
  SvxGhitList *mc_ghit = findNode::getClass<SvxGhitList>(m_mcnode,"SvxGhitList");
  SvxRawhitList *rawhitlist = findNode::getClass<SvxRawhitList>(topNode,"SvxRawhitList");
  SvxGhitRawhitList *g2rlist = findNode::getClass<SvxGhitRawhitList>(topNode,"SvxGhitRawhitList");
  svxAddress* address = findNode::getClass<svxAddress>(topNode, "svxAddress");
  if ( address==NULL ) {
    if ( verbosity>0 ) cerr << PHWHERE<< "Can't find svxAddress." << endl;
    return 0;
  } 
  svxAddress &SvxAddObj = *address;

  if ( !ghitlist || !mc_ghit || !rawhitlist ) {
    cerr << "SvxGhitList or SvxRawhitList node not found!" << endl;
    return 0;
  }
  
  ///
  /// copy ghits from simulation, shift z position, and save in output node
  ///
  int nmchit = mc_ghit->get_nGhits();
  if ( verbosity>0 ) cout << "ghit in MC : " << nmchit << endl;
  SvxSensor *svxsen;
  for ( int imchit=0; imchit<nmchit; imchit++ ) {
    SvxGhit *mchit = mc_ghit->get_Ghit(imchit);
    int layer = mchit->get_layer();
    if(layer<0 || layer>3) continue; // FVTX
    int ladder = mchit->get_ladder();
    int sensor = mchit->get_sensor();
    svxsen = m_svxgeo->GetSensorPtr(layer, ladder, sensor);
    if(!svxsen) continue;
    double pos_localout[3];
    double pos_globalin[3];
    double pos_localin[3];
    double pos_globalout[3];
    pos_localin[0] = mchit->get_xyzlocalin(0);
    pos_localin[1] = mchit->get_xyzlocalin(1);
    pos_localin[2] = mchit->get_xyzlocalin(2);
    pos_localout[0] = mchit->get_xyzlocalout(0);
    pos_localout[1] = mchit->get_xyzlocalout(1);
    pos_localout[2] = mchit->get_xyzlocalout(2);
    svxsen->position_local2global(pos_localin, pos_globalin);
    svxsen->position_local2global(pos_localout, pos_globalout);
    pos_globalin[2] += zshift;
    pos_globalout[2] += zshift;
    bool found = false;
    for ( int isn=0; isn<m_max_sensor[layer]; isn++ ) {
      svxsen = m_svxgeo->GetSensorPtr(layer, ladder, isn);
      svxsen->position_global2local(pos_globalin, pos_localin);
      if ( fabs(pos_localin[2])<m_sensorZwhalf[layer] ) {
	sensor = isn;
	svxsen->position_global2local(pos_globalout, pos_localout);
	found = true;
	break;
      }
    }
    if ( found ) {  /// if shifted ghit is still in VTX acceptance, save it.
      SvxGhit *ghit = ghitlist->addGhit();
      ghit->set_svxSection(mchit->get_svxSection());
      ghit->set_layer(layer);
      ghit->set_ladder(ladder);
      ghit->set_sensor(sensor);
      ghit->set_mctrack(mchit->get_mctrack());
      ghit->set_idPart(mchit->get_idPart());
      ghit->set_track(mchit->get_track());
      ghit->set_dele(mchit->get_dele());
      ghit->set_tof(mchit->get_tof());
      ghit->set_isubevent(mchit->get_isubevent());
      ghit->set_nfile(mchit->get_nfile());
      ghit->set_xyzglobal(0, mchit->get_xyzglobal(0));
      ghit->set_xyzglobal(1, mchit->get_xyzglobal(1));
      ghit->set_xyzglobal(2, mchit->get_xyzglobal(2)+zshift);
      ghit->set_xyzlocalin(0, pos_localin[0]);
      ghit->set_xyzlocalin(1, pos_localin[1]);
      ghit->set_xyzlocalin(2, pos_localin[2]);
      ghit->set_xyzlocalout(0, pos_localout[0]);
      ghit->set_xyzlocalout(1, pos_localout[1]);
      ghit->set_xyzlocalout(2, pos_localout[2]);
      ghit->set_pmomxyz(0, mchit->get_pmomxyz(0));
      ghit->set_pmomxyz(1, mchit->get_pmomxyz(1));
      ghit->set_pmomxyz(2, mchit->get_pmomxyz(2));
    } /// if : found
  } /// for : imchit

  ///
  /// convert ghit to rawhit
  ///
  int nghit = ghitlist->get_nGhits();
  int nrawhit = rawhitlist->get_nRawhits();
  if ( verbosity>0 ) {
    cout << "number of ghit after shifting vertex : " << nghit << " " << nrawhit << endl;
  }
  if ( nghit>0 ) {
    ghitlist->sort_sensorID();
    for ( int ilr=0; ilr<4; ilr++ ) {
      for ( int ild=0; ild<m_max_ladder[ilr]; ild++ ) {
	for ( int isn=0; isn<m_max_sensor[ilr]; isn++ ) {
	  int ng2r = g2rlist->get_nGhitRawhits();
	  svxsen = m_svxgeo->GetSensorPtr(ilr, ild, isn);
	  svxsen->Reset();
	  int nreal = svxsen->makeRawhits(ghitlist,rawhitlist,g2rlist,&SvxAddObj);///defined in SvxStripixel.h
	  int nnoise = svxsen->AddNoise(nrawhit,nrawhit+nreal,rawhitlist,g2rlist,ng2r);
    if(false)std::cout<<svxsen->get_adcthresh()<<std::endl;
    if(false)svxsen->printPar();
	  nrawhit += nreal+nnoise;
	} // isn
      } // ild
    } // ilr
  }
  
  if ( verbosity>0 ) {
    cout << "number of rawhit after shifting vertex : " << nrawhit << endl;
  }
  return nrawhit;
}

//----------------------------------------------------------------------------------

int SvxEmbedSimhit::End(PHCompositeNode *topNode)
{
  return EVENT_OK;
}

void SvxEmbedSimhit::copyRunHeader(PHCompositeNode *topNode)
{
  RunHeader *run      = findNode::getClass<RunHeader>(topNode,   "RunHeader");
  RunHeader *run_real = findNode::getClass<RunHeader>(m_realnode,"RunHeader");

  // copy
  //*run = *run_real;
  run->set_RunNumber(     run_real->get_RunNumber());
  run->set_Bfield(        run_real->get_Bfield());
  run->set_TimeStart(     run_real->get_TimeStart());
  run->set_TimeStop(      run_real->get_TimeStop());
  run->set_currentNorth(  run_real->get_currentNorth());
  run->set_currentSouth(  run_real->get_currentSouth());
  run->set_currentCentral(run_real->get_currentCentral());
  run->set_currentInner(  run_real->get_currentInner());
}

void SvxEmbedSimhit::copyEventHeader(PHCompositeNode *topNode)
{
  EventHeader *evt      = findNode::getClass<EventHeader>(topNode,   "EventHeader");
  EventHeader *evt_real = findNode::getClass<EventHeader>(m_realnode,"EventHeader");

  *evt = *evt_real;

  evt->set_EvtSequence(evt_real->get_EvtSequence());
  evt->set_EvtType(    evt_real->get_EvtType()    );
  evt->set_TimeStamp(  evt_real->get_TimeStamp()  );

  //const set<unsigned int>* badpkt_p = evt_real->GetBadPacketSet();
  //set<unsigned int>::iterator itr;
  //for(itr=badpkt_p->begin(); itr!=badpkt_p->end(); ++itr){
  //  evt->AddBadPacket( *itr );
  //}
}


void SvxEmbedSimhit::copyPreviousEvent(PHCompositeNode *topNode)
{
  PreviousEvent *prev      = findNode::getClass<PreviousEvent>(topNode,   "PreviousEvent");
  PreviousEvent *prev_real = findNode::getClass<PreviousEvent>(m_realnode,"PreviousEvent");

  *prev = *prev_real;

  for(int itick=0; itick<3; itick++){ prev->set_clockticks(prev_real->get_clockticks(itick), itick); }
}

void SvxEmbedSimhit::copyTrigLvl1(     PHCompositeNode* topNode)
{
  TrigLvl1 *prev      = findNode::getClass<TrigLvl1>(topNode,   "TrigLvl1");
  TrigLvl1 *prev_real = findNode::getClass<TrigLvl1>(m_realnode,"TrigLvl1");

  *prev = *prev_real;

  prev->set_lvl1_trigraw(    prev_real->get_lvl1_trigraw()     );
  prev->set_lvl1_triglive(   prev_real->get_lvl1_triglive()    );
  prev->set_lvl1_trigscaled( prev_real->get_lvl1_trigscaled()  );
  prev->set_lvl1_clock_cross(prev_real->get_lvl1_clock_cross() );
  for(int i=0; i<5; i++){
    prev->set_lvl1_rbits(prev_real->get_lvl1_rbits(i), i);
  }
  for(int i=0; i<2; i++){
    prev->set_lvl1_beam_clk(prev_real->get_lvl1_beam_clk(i), i);
  }

}
void SvxEmbedSimhit::copyTrigRunLvl1(  PHCompositeNode* topNode)
{
  TrigRunLvl1 *trgrun      = findNode::getClass<TrigRunLvl1>(topNode,   "TrigRunLvl1");
  TrigRunLvl1 *trgrun_real = findNode::getClass<TrigRunLvl1>(m_realnode,"TrigRunLvl1");

  *trgrun = *trgrun_real;

  trgrun->set_lvl1_trigger_description( trgrun_real->get_lvl1_trigger_description() );
  trgrun->set_lvl1_trigger_version(     trgrun_real->get_lvl1_trigger_version() );

  trgrun->set_lvl1_bbcll1_description( trgrun_real->get_lvl1_bbcll1_description() );
  trgrun->set_lvl1_bbcll1_version(     trgrun_real->get_lvl1_bbcll1_version() );
  trgrun->set_lvl1_partition_name(     trgrun_real->get_lvl1_partition_name() );

  for(int i=0; i<32; i++){
    trgrun->set_lvl1_trig_name(         trgrun_real->get_lvl1_trig_name(i),         i);
    trgrun->set_lvl1_trigger_enable(    trgrun_real->get_lvl1_trigger_enable(i),    i);//
    trgrun->set_lvl1_trig_bitmask(      trgrun_real->get_lvl1_trig_bitmask(i),      i);
    trgrun->set_lvl1_trig_bit(          trgrun_real->get_lvl1_trig_bit(i),          i);//
    trgrun->set_lvl1_trig_scale_down(   trgrun_real->get_lvl1_trig_scale_down(i),   i);//
    trgrun->set_lvl1_lvl2_reject_enable(trgrun_real->get_lvl1_lvl2_reject_enable(i),i);//
    trgrun->set_lvl1_trig_rate_begin(   trgrun_real->get_lvl1_trig_rate_begin(i),   i);
  }

  for(int i=0; i<130; i++){
    trgrun->set_lvl1_rbit_name( trgrun_real->get_lvl1_rbit_name(i), i);
  }

  trgrun->set_start_time(trgrun_real->get_start_time());
  trgrun->set_run_number(trgrun_real->get_run_number());
}
void SvxEmbedSimhit::copyBbcOut(       PHCompositeNode* topNode)
{
  BbcOut *bbcout = findNode::getClass<BbcOut>(topNode,"BbcOut");
  BbcOut *bbcout_real = findNode::getClass<BbcOut>(m_realnode,"BbcOut");

  *bbcout = *bbcout_real;

  bbcout->set_TimeZero(bbcout_real->get_TimeZero());
  bbcout->set_Vertex(  bbcout_real->get_VertexPoint(), bbcout_real->get_dVertexPoint());

  //for(int i=0; i<2; i++){
  //  bbcout->AddBbcNS(bbcout_real->get_nPmt(i), bbcout_real->get_ChargeSum(i), bbcout_real->get_Timing(i), i);
  //}
}
void SvxEmbedSimhit::copyPHGlobal(     PHCompositeNode* topNode)
{
  PHGlobal *global      = findNode::getClass<PHGlobal>(topNode,   "PHGlobal");
  PHGlobal *global_real = findNode::getClass<PHGlobal>(m_realnode,"PHGlobal");

  *global = *global_real;

  global->setBbcZVertex(     global_real->getBbcZVertex());
  global->setBbcZVertexError(global_real->getBbcZVertexError());
}

