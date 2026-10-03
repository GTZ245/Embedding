#include <iostream>
#include <iomanip>
#include "phool.h"
#include "PHTypedNodeIterator.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "Fun4AllReturnCodes.h"
#include "PHPoint.h"

#include "svxanamini.h"

#include "PHGlobal.h"
#include "VtxOut.h"
#include "TrigLvl1.h"
#include "RunHeader.h"
#include "EventHeader.h"
#include "Bbc.hh"
#include "BbcOut.h"
#include <EventHeader.h>

#include "TFile.h"
#include "TTree.h"

#include "getClass.h"

using namespace std;
using namespace findNode;



//==============================================================

svxanamini::svxanamini() {
  ThisName = "svxanaminimini";
  init_ana=0;
  EventNumber=0;
  OutputFileName="svxanalysis.root";
}

//==============================================================

svxanamini::svxanamini(string filename) : OutputFileName(filename)
{
  ThisName = "svxanamini";
  init_ana=0;
  EventNumber=0;

}

//==============================================================

svxanamini::~svxanamini() {
}

//==============================================================

int svxanamini::Init(PHCompositeNode *topNode) {

  cout << "svxanamini::Init started..." << endl;
  OutputNtupleFile = new TFile(OutputFileName.c_str(),"RECREATE");
  cout << "svxanamini::Init: output file " << OutputFileName << " opened." << endl;


  initEvtTree();

  cout << "svxanamini::Init ended." << endl;
  return 0;
}

//==============================================================
  
int svxanamini::InitRun(PHCompositeNode *topNode) {
  cout << "svxanamini::InitRun started..." << endl;

  cout << "svxanamini::InitRun ended." << endl;
  return 0;
}

//==============================================================


int svxanamini::process_event(PHCompositeNode *topNode) {

  //float ntp1[99],ntp2[99],ntp3[99],ntp4[99],ntp5[99],ntp6[99],ntp7[99];

  RunHeader   *run     = getClass<RunHeader>(  topNode, "RunHeader");
  EventHeader *evthdr  = getClass<EventHeader>(topNode, "EventHeader");
  PHGlobal    *global  = getClass<PHGlobal>(   topNode, "PHGlobal");
  BbcOut      *bbc     = getClass<BbcOut>(     topNode, "BbcOut");
  VtxOut      *vtxout  = getClass<VtxOut>(     topNode, "VtxOut");
  TrigLvl1    *trglvl1 = getClass<TrigLvl1>(   topNode, "TrigLvl1");

  int   RunNumber = (run!=NULL)    ? run->get_RunNumber()      : -9999;
  int   EventNum  = (evthdr!=NULL) ? evthdr->get_EvtSequence() : -9999;

  float zvtx = (bbc!=NULL) ? bbc->get_VertexPoint() : -999.0;
  if(zvtx<-900){
    zvtx = (vtxout!=NULL) ? vtxout->get_ZVertex("BBC") : -999.0;
  }

  int   nbbcs = -999,  nbbcn = -999;
  float qbbcs = -999., qbbcn = -999.;
  if(bbc!=NULL)         { 
    nbbcs = bbc->get_nPmt(Bbc::South); 
    nbbcn = bbc->get_nPmt(Bbc::North); 
    qbbcs = bbc->get_ChargeSum(Bbc::South); 
    qbbcn = bbc->get_ChargeSum(Bbc::North); 
  }
  else if(global!=NULL) { 
    nbbcs = global->getBbcMultS();     
    nbbcn = global->getBbcMultN();     
    qbbcs = global->getBbcChargeS();     
    qbbcn = global->getBbcChargeN();     
  }

  float bbcq = (bbc==NULL&&global==NULL) ? -999 : qbbcs+qbbcn;

  // most precise primary vertex 
  float vtx_z = (vtxout!=NULL) ? vtxout->get_ZVertex() : -999.;

  

  if(init_ana==0) {
    cout<<"init ana"<<endl;
    topNode->print(); 
    init_ana=1;
 
    cout<<endl<<endl;
    cout<<"init ana ends"<<endl;
  }

  if((EventNumber%1000)==0){
    cout << "------- Event # " << EventNumber << endl;
  }



  /////////////////////////////
  // event by event info
  float zzdc = -9999.0;
  PHPoint vtxs_pos(-9999.0, -9999.0, -9999.0);
  PHPoint vtxp_pos(-9999.0, -9999.0, -9999.0);
  PHPoint vtxbbc_pos(-9999.0, -9999.0, -9999.0);

  if(vtxout!=NULL) {
    zzdc     = vtxout->get_ZVertex("ZDC");
    vtxs_pos = vtxout->get_Vertex("SVX");
    vtxp_pos = vtxout->get_Vertex("SVX_PRECISE");
    vtxbbc_pos = vtxout->get_Vertex("BBC");

    if(verbosity>2){
      cout<<"SVXW :  "<<(vtxout->isVtx("SVXW")?"OK":"Fail")<<endl;
      cout<<"SVXE :  "<<(vtxout->isVtx("SVXE")?"OK":"Fail")<<endl;
    }
  }

  m_evt_run   = RunNumber;
  m_evt_event = EventNumber;
  m_evt_strig = (trglvl1!=NULL) ? trglvl1->get_lvl1_trigscaled()  : -1;
  m_evt_eseq  = EventNum;
  m_evt_zvtx  = vtx_z;
  m_evt_zbbc  = zvtx;
  m_evt_zzdc  = zzdc;
  m_evt_bbcq  = bbcq;
  m_evt_bbcns = nbbcs;
  m_evt_bbcnn = nbbcn;
  m_evt_bbcqs = qbbcs;
  m_evt_bbcqn = qbbcn;
  m_evt_xvtxs = vtxs_pos.getX();
  m_evt_yvtxs = vtxs_pos.getY();
  m_evt_zvtxs = vtxs_pos.getZ();
  m_evt_xvtxp = vtxp_pos.getX();
  m_evt_yvtxp = vtxp_pos.getY();
  m_evt_zvtxp = vtxp_pos.getZ();

  if(verbosity>2){
    cout<<"x_vtx_s : x_vtx_p : "<< vtxs_pos.getZ()<<" "<<vtxp_pos.getZ()<<endl;
    cout<<"x_vtx_s : y_vtx_s : "<< vtxs_pos.getX()<<" "<<vtxs_pos.getY()<<endl;
    cout<<"zbbc_x : zbbc_y : zbbc_z "
        << vtxbbc_pos.getX()<<" "<<vtxbbc_pos.getY()<<" "<<vtxbbc_pos.getZ()<<endl;
  }

  ntpevt->Fill();

  EventNumber++;
  return 0;
}

//==============================================================

int svxanamini::End(PHCompositeNode *topNode) {
  cout << "svxanamini::End:  Writing out..." << endl;
  OutputNtupleFile->Write();
  cout << "svxanamini::End:  Closing output file..." << endl;
  OutputNtupleFile->Close();
  delete OutputNtupleFile;
  OutputNtupleFile=0;
  return 0;
}

void svxanamini::initEvtTree(){
  ntpevt = new TTree("ntpevt", "Event Info Tree");
  ntpevt->Branch("run",    &m_evt_run,     "run/I");
  ntpevt->Branch("event",  &m_evt_event,   "event/I");
  ntpevt->Branch("strig",  &m_evt_strig,   "strig/I");
  ntpevt->Branch("eseq",   &m_evt_eseq,    "eseq/I");
  ntpevt->Branch("zvtx",   &m_evt_zvtx,    "zvtx/F");
  ntpevt->Branch("zbbc",   &m_evt_zbbc,    "zbbc/F");
  ntpevt->Branch("zzdc",   &m_evt_zzdc,    "zzdc/F");
  ntpevt->Branch("bbcq",   &m_evt_bbcq,    "bbcq/F");
  ntpevt->Branch("bbcns",  &m_evt_bbcns,   "bbcns/I");
  ntpevt->Branch("bbcnn",  &m_evt_bbcnn,   "bbcnn/I");
  ntpevt->Branch("bbcqs",  &m_evt_bbcqs,   "bbcns/F");
  ntpevt->Branch("bbcqn",  &m_evt_bbcqn,   "bbcnn/F");
  ntpevt->Branch("xvtxs",  &m_evt_xvtxs,   "xvtxs/F");
  ntpevt->Branch("yvtxs",  &m_evt_yvtxs,   "yvtxs/F");
  ntpevt->Branch("zvtxs",  &m_evt_zvtxs,   "zvtxs/F");
  ntpevt->Branch("xvtxp",  &m_evt_xvtxp,   "xvtxp/F");
  ntpevt->Branch("yvtxp",  &m_evt_yvtxp,   "yvtxp/F");
  ntpevt->Branch("zvtxp",  &m_evt_zvtxp,   "zvtxp/F");
}

