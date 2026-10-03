#include "SvxRpQA.h"

#include <iostream>

#include <phool.h>
#include <PHTypedNodeIterator.h>
#include <PHCompositeNode.h>
#include <PHIODataNode.h>
#include <Fun4AllReturnCodes.h>


#include <PHGlobal.h>
#include <VtxOut.h>
#include <RunHeader.h>
#include <EventHeader.h>
#include <PreviousEvent.h>
#include <TriggerHelper.h>

#include <RpConst.h>
#include <RpSnglSumXY.h>
#include <RpSumXYObject.h>

#include "TFile.h"
#include "TH3.h"
#include "TH2.h"

#include "getClass.h"

using namespace std;
using namespace findNode;

//==============================================================

SvxRpQA::SvxRpQA(string filename) : 
  SubsysReco("SVXRPQA"),
  OutputFileName(filename),
  OutputNtupleFile(NULL),
  init_ana(0),
  EventNumber(0),
  m_zvcut(10), // 10cm
  m_istickcut(true)
{
  for(int i=0; i<3; i++){ m_pticks[i]=-1; }
}

//==============================================================

SvxRpQA::~SvxRpQA() {
}

//==============================================================

int SvxRpQA::Init(PHCompositeNode *topNode) {

  cout << "SvxRpQA::Init started..." << endl;
  OutputNtupleFile = new TFile(OutputFileName.c_str(),"RECREATE");
  cout << "SvxRpQA::Init: output file " << OutputFileName << " opened." << endl;

  initRpSumTree();


  cout << "SvxRpQA::Init ended." << endl;
  return 0;
}

//==============================================================
  
int SvxRpQA::InitRun(PHCompositeNode *topNode) {
  cout << "SvxRpQA::InitRun started..." << endl;

  // check magnet current 
  RunHeader* runheader = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if (runheader==NULL) {
    cout << PHWHERE<< "Can't find runheader. " << endl;
    return ABORTRUN;
  }

  cout << "SvxRpQA::InitRun ended." << endl;
  return 0;
}

//==============================================================


int SvxRpQA::process_event(PHCompositeNode *topNode) {
  PHGlobal      *global = getClass<PHGlobal>(topNode,"PHGlobal");
  RunHeader     *run    = getClass<RunHeader>(topNode,"RunHeader");
  EventHeader   *evthdr = getClass<EventHeader>(topNode,"EventHeader");
  PreviousEvent *peve   = getClass<PreviousEvent>(topNode,"PreviousEvent");
  VtxOut        *vtxout = getClass<VtxOut>(topNode,"VtxOut");
  RpSumXYObject *rpsum  = getClass<RpSumXYObject>(topNode,"RpSumXYObject");

  if(EventNumber==0)
    if(peve==NULL){ cout<<"no PreviousEvent"<<endl; /* return 0; */ }


  int   RunNumber = (run!=NULL)    ? run->get_RunNumber()      : -9999;
  int   EvtSeq    = (evthdr!=NULL) ? evthdr->get_EvtSequence() : -9999;
  float cent      = (global!=NULL) ? global->getCentrality()   : -9999;

  
  EventNumber++;

  if(init_ana==0) {
    cout<<"init ana"<<endl;
    topNode->print(); 
    init_ana=1;
    cout << "run = "<<RunNumber<<" event = "<<EvtSeq<<endl;
    cout<<endl<<"init ana ends"<<endl;
  }

  if((EventNumber%1000)==0){
    cout << "------- Event # " << EventNumber << endl;
  }

  /////////////////////////////////////
  //float zbbc = vtxout->get_ZVertex("BBC");
  float zvtx = vtxout->get_ZVertex();
  if(!EventSelection(topNode,zvtx)) 
    return EVENT_OK;


  fillRpSum(EvtSeq, cent, vtxout, rpsum);


  return 0;
}

//==============================================================

int SvxRpQA::End(PHCompositeNode *topNode) {
  cout << "svxana::End:  Writing out..." << endl;
  OutputNtupleFile->Write();
  cout << "svxana::End:  Closing output file..." << endl;

  OutputNtupleFile->Close();
  delete OutputNtupleFile;
  OutputNtupleFile=0;
  return 0;
}

void SvxRpQA::initRpSumTree(){
  //
  //     -3.0 -2.5 -2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0 2.5 3.0
  // B0      0    1    2    3    4    5   6   7   8   9   10  11
  // B1          12   13   14   15   16  17  18  19  20   21
  // B2               22   23   24   25  26  27  28  29 
  // B3                    30   31   32  33  34  35 
  // Ball   36   37   38   39   40   41  42  43  44  45   46  47
  // B0-3,all(nogap)   48, 49, 50, 51, 52
  // B0-3,all(egap 1)  53, 54, 55, 56, 57
  // B0-3,all(egap 2)  58, 59, 60, 61, 62
  static const char clayer[5][256] = {"B0", "B1", "B2", "B3", "Ball"};
  static const char ceta[15][256] = {"-3:-2.5", "-2.5:-2", "-2:-1.5", "-1.5:-1",
                                     "-1:-0.5", "-0.5:0", "0:0.5", "0.5:1",
                                     "1:1.5", "1.5:2", "2:2.5", "2.5:3",
                                     "no gap", "gap1", "gap2"
                                    };
                                    //                    1,1
                                    //0,1,2,3,4,5,6,7,8,9,0,1
  static const int ilayer[NRPKIND] = {0,0,0,0,0,0,0,0,0,0,0,0, // B0 12 kinds 
                                      1,1,1,1,1,1,1,1,1,1,     // B1 10 kinds
                                      2,2,2,2,2,2,2,2,         // B2  8 kinds
                                      3,3,3,3,3,3,             // B3  6 kinds
                                      4,4,4,4,4,4,4,4,4,4,4,4, // Ball 12 kinds
                                      0,1,2,3,4,               // B0-3, all, nogap
                                      0,1,2,3,4,               // B0-3, all, gap 1
                                      0,1,2,3,4,               // B0-3, all, gap 2
                                     };
  static const int iegap[NRPKIND] = {0,1,2,3,4,5,6,7,8,9,10,11, // B0 12 kinds 
                                       1,2,3,4,5,6,7,8,9,10,    // B1 10 kinds
                                         2,3,4,5,6,7,8,9,       // B2  8 kinds
                                           3,4,5,6,7,8,         // B3  6 kinds
                                     0,1,2,3,4,5,6,7,8,9,10,11, // Ball 12 kinds
                                      12,12,12,12,12,           // B0-3, all, nogap
                                      13,13,13,13,13,           // B0-3, all, gap 1
                                      14,14,14,14,14,           // B0-3, all, gap 2
                                     };

  for(int ikind=0; ikind<NRPKIND; ikind++){
    // B0:0-11, B1:12-21, B2:22-29, B3:30-35,
    // Ball:36-47, B0-B3(no gap):48-51, Ball(nogap):52
    // B0-B3,all(eta gap1): 53-57, B0-B3,all(eta gap2): 58-62,
    h_qxqy_cent[ikind] = new TH3F(
                         Form("h_qxqy_cent_%d", ikind),
                         Form("qy vs qx vs cent eta=%s at %s (%d)", 
                               ceta[iegap[ikind]], clayer[ilayer[ikind]], ikind),
                          50,0,100, 100,-1,1, 100,-1,1
                       );

    h_qxqy_zvtx[ikind] = new TH3F(
                         Form("h_qxqy_zvtx_%d", ikind),
                         Form("qy vs qx vs zvtx eta=%s at %s (%d)", 
                               ceta[iegap[ikind]], clayer[ilayer[ikind]], ikind),
                          50,-10,10, 100,-1,1, 100,-1,1
                       );

    h_qx_eseq[ikind] = new TH3F(
                         Form("h_qx_eseq_%d", ikind),
                         Form("qx vs eseq eta=%s at %s (%d)", 
                               ceta[iegap[ikind]], clayer[ilayer[ikind]], ikind),
                          500,0,60*1e6, 100,-1,1, 5,0,5
                       );

    h_qy_eseq[ikind] = new TH3F(
                         Form("h_qy_eseq_%d", ikind),
                         Form("qy vs eseq eta=%s at %s (%d)", 
                               ceta[iegap[ikind]], clayer[ilayer[ikind]], ikind),
                          500,0,60*1e6, 100,-1,1, 5,0,5
                       );

  }

  static const char cbbc[3][256] = {"North", "South", "All"};
  for(int ikind=0; ikind<3; ikind++){
    hbbc_qxqy_cent[ikind] = new TH3F(
                            Form("hbbc_qxqy_cent_%d", ikind),
                            Form("qy vs qx vs cent BBC %s", cbbc[ikind]),
                            50,0,100, 100,-1,1, 100,-1,1);
    hbbc_qxqy_zvtx[ikind] = new TH3F(
                            Form("hbbc_qxqy_zvtx_%d", ikind),
                            Form("qy vs qx vs zvtx BBC %s", cbbc[ikind]),
                            50,-10,10, 100,-1,1, 100,-1,1);
    hbbc_qx_eseq[ikind]   = new TH3F(
                            Form("hbbc_qx_eseq_%d", ikind),
                            Form("qx vs eseq vtx BBC %s", cbbc[ikind]),
                            500,0,60*1e6, 100,-1,1, 5,0,5);
    hbbc_qy_eseq[ikind]   = new TH3F(
                            Form("hbbc_qy_eseq_%d", ikind),
                            Form("qy vs eseq vtx BBC %s", cbbc[ikind]),
                            500,0,60*1e6, 100,-1,1, 5,0,5);
  }

}

void SvxRpQA::fillRpSum(int eseq, float cent, VtxOut *vtxout, RpSumXYObject* rpsum)
{

  if(verbosity>0) cout<<"fillRpSum"<<endl;

  if(vtxout==NULL){
    if(EventNumber==0){
      cout<<PHWHERE<<" VtxOut is NULL "<<endl;
    }
    return;
  }

  if(rpsum==NULL){
    if(EventNumber==0){
      cout<<PHWHERE<<" RpSumXYObject is NULL "<<endl;
    }
    return;
  }

  if(verbosity>1){
    rpsum->identify();
  }

  //float zvtx = vtxout->get_ZVertex("BBC");  // should use BBCZ?
  float zvtx = vtxout->get_ZVertex();  // should use BBCZ?

  int icent = (0<=cent&&cent<100) ? (int) cent/20. : -1;


  // fill histogram
  for(int id=0; id<NRPKIND-10; id++){ // 0-52
    int rpID = RP::calcIdCode(RP::ID_SVX, id, 1); // Ball alleta 2nd order
    RpSnglSumXY* obj = rpsum->getRpSumXY(rpID);
    if(obj==NULL) 
      {
        //cout<<"No RpSnglSumXY object event : kind:"<<id<<" event "<<eseq<<endl;
        continue;
      }

    float qx  = obj->QVector(0);
    float qy  = obj->QVector(1);
    float qw  = obj->Weight();
    if(qw>0) 
      {
        h_qxqy_cent[id]->Fill(cent,  qx/qw, qy/qw);
        h_qxqy_zvtx[id]->Fill(zvtx,  qx/qw, qy/qw);

        h_qx_eseq[id]->Fill(eseq,  qx/qw, icent);
        h_qy_eseq[id]->Fill(eseq,  qy/qw, icent);
      }
  }

  // etagap1, etagap2
  // B0      0    1    2    3    4    5   6   7   8   9   10  11
  // B1          12   13   14   15   16  17  18  19  20   21
  // B2               22   23   24   25  26  27  28  29 
  // B3                    30   31   32  33  34  35 
  // Ball   36   37   38   39   40   41  42  43  44  45   46  47
  static const int nid_etagap[10] = {10,8,6,4,10, // egap1
                                      8,6,4,2,8}; // egap2
  static const int id_etagap[10][12] = {
      { 0, 1, 2, 3, 4,     7, 8, 9,10,11}, // egap 1jB0  
      {   12,13,14,15,    18,19,20,21   }, // egap 1jB1  
      {      22,23,24,    27,28,29      }, // egap 1jB2  
      {         30,31,    34,35         }, // egap 1jB3  
      {36,37,38,39,40,    43,44,45,46,47}, // egap 1jBall
      { 0, 1, 2, 3,           8, 9,10,11}, // egap 2jB0     
      {   12,13,14,          19,20,21   }, // egap 2jB1     
      {      22,23,          28,29      }, // egap 2jB2     
      {         30,          35         }, // egap 2jB3     
      {36,37,38,39,          44,45,46,47}, // egap 2jBall   
    };


  for(int ikind=0; ikind<10; ikind++)
    {
      float qv[3]={0,0,0};

      // B0 eta delta +-1
      for(int id=0; id<nid_etagap[ikind]; id++){
        int rpID = RP::calcIdCode(RP::ID_SVX, id_etagap[ikind][id], 1); // Ball alleta 2nd order
        RpSnglSumXY* obj = rpsum->getRpSumXY(rpID);
        if(obj==NULL) {
          continue;
        }
        qv[0] += obj->QVector(0);
        qv[1] += obj->QVector(1);
        qv[2] += obj->Weight();
      }

      if(qv[2]>0) 
        {
          int id_kind = ikind + 53;
          h_qxqy_cent[id_kind]->Fill(cent,  qv[0]/qv[2], qv[1]/qv[2]);
          h_qxqy_zvtx[id_kind]->Fill(zvtx,  qv[0]/qv[2], qv[1]/qv[2]);

          h_qx_eseq[id_kind]->Fill(eseq,  qv[0]/qv[2], icent);
          h_qy_eseq[id_kind]->Fill(eseq,  qv[1]/qv[2], icent);
        }
    }

  /////////////
  // BBC
  for(int id=0; id<3; id++){
    int bbcid = RP::calcIdCode(RP::ID_BBC, id, 1); // Bbc north, south, all

    RpSnglSumXY* obj = rpsum->getRpSumXY(bbcid);
    if(obj==NULL) {
      continue;
    }

    float qx = obj->QVector(0);
    float qy = obj->QVector(1);
    float qw = obj->Weight();

    if(qw>0)
      {
        hbbc_qxqy_cent[id]->Fill(cent, qx/qw, qy/qw);
        hbbc_qxqy_zvtx[id]->Fill(zvtx, qx/qw, qy/qw);
        if(cent>80)
          hbbc_qx_eseq[id]->Fill(eseq, qx/qw,0.);
          hbbc_qy_eseq[id]->Fill(eseq, qy/qw,0.);
      }
  }
}


bool SvxRpQA::TickCut(PHCompositeNode *topNode) 
{ 
  PreviousEvent *peve = getClass<PreviousEvent>(topNode,"PreviousEvent");
  if(peve==NULL)
    {
      cerr<<"No PreviousEvent Object"<<endl;
      return false;
    }

  int pticks[3]={-1,-1,-1};
  for( int i = 0; i < 3; i++ ) pticks[i] = peve->get_clockticks(i);

  return ( ( 50<pticks[0]&&pticks[0]<120)|| 
           (700<pticks[1]&&pticks[1]<780) ); 
}


//event selection
bool SvxRpQA::EventSelection(PHCompositeNode *topNode, float zbbc)
{
  //pixel tick cut
  if(m_istickcut){
    if( TickCut(topNode) != EVENT_OK ) return false;
  }
  //bbcz vertex cut
  if(fabs(zbbc)>m_zvcut) return false;

   
/*
  //trigger selection
  TriggerHelper d_trghelp(topNode);
  bool isL1     = d_trghelp.didLevel1TriggerGetScaled( "BBCLL1(>1 tubes)" ); //BBCL1
  bool isnarrow = d_trghelp.didLevel1TriggerGetScaled( "BBCLL1(>1 tubes) narrowvtx" ); //BBCL1 narrow
  bool isnarrow_copyA = d_trghelp.didLevel1TriggerGetScaled( "BBCLL1(>1 tubes) narrowvtx Copy A" ); // BBCL1 narrow A
  bool isnarrow_copyB = d_trghelp.didLevel1TriggerGetScaled( "BBCLL1(>1 tubes) narrowvtx Copy B" ); // BBCL1 narrow B
  bool isnovtx  = d_trghelp.didLevel1TriggerGetScaled( "BBCLL1(>1 tubes) novertex" ); // BBCL1 novertex

  //store trigger info
  if( isL1 || isnarrow || isnarrow_copyA || isnarrow_copyB || isnovtx){
    //OK. This is MB event
    return true;
  }

  return false;
*/
  return true;
}

