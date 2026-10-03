#include <iostream>
//#include <iomanip>
//#include <vector>
//#include <map>
//#include <algorithm>
//#include <cstdlib>
//#include "gsl/gsl_rng.h"
//#include <math.h>

#include "phool.h"
#include "PHTypedNodeIterator.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "Fun4AllReturnCodes.h"
#include "Fun4AllServer.h"
#include "recoConsts.h"

#include "embeddisp.h"

#include "CrkGeometryObject.hh"
#include "PHCentralTrack.h"
#include "PHSnglCentralTrack.h"
#include "EventHeader.h"
#include "VtxOut.h"
//#include "PHGlobal.h"
//#include "TriggerHelper.h"
//#include "RunHeader.h"
//#include <PreviousEvent.h>
//#include "TrigLvl1.h"
//#include "TrigRunLvl1.h"
//#include "Bbc.hh"
//#include "BbcOut.h"

#include <DchHitLineTable.hh>
#include <PadCluster.h>
#include <CrkHit.h>
#include <PHTrackOut.h>
#include <PHDchTrackOut.h>
#include <emcClusterContainer.h>
#include <emcClusterContent.h>

//#include "fkinWrapper.h"

//#include "McEvalSingleList.h"

//#include "SvxRawhitList.h"
//#include "SvxRawhit.h"
//#include "SvxSegmentList.h"
//#include "SvxSegment.h"
//#include "SvxTracker.h"
//#include "SvxRawhitClusterList.h"
//#include "SvxRawhitCluster.h"
#include "SvxClusterList.h"
#include "SvxCluster.h"
#include "SvxCentralTrackList.h"
#include "SvxCentralTrack.h"
//#include "SvxClusterInfo.h"
//#include "SvxResidualInfo.h"
//#include "SvxBeamCenterPar.h"
//#include "SvxEventInfo.h"

//#include "SvxGhitList.h"
//#include "SvxGhit.h"
//#include "SvxGhitRawhitList.h"
//#include "SvxGhitRawhit.h"

//#include "SvxClusterContainer.h"

//#include "svxAddress.hh"
//#include "svxDetectorGeo.hh"
//#include "SvxSensor.h"
//
//#include "compactCNT/SvxHitMap.h"
//#include "compactCNT/SvxTrackMap.h"
//#include <SvxCentralTrackRecalList.h>

//#include <TofwHit.h>
//#include <HbdBlobList.h>
#include <PHPointerList.h>
#include <PHEmbedMcRecoTrack.hh>

#include "TFile.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TMarker.h"
#include "TMath.h"
//#include "TTree.h"
//#include "TF1.h"
//#include "TH1.h"
//#include "TH3.h"
//#include "TLorentzVector.h"
//#include "TVector3.h"
//#include "TRandom2.h"
//#include "TMath.h"

#include "getClass.h"

using namespace std;
using namespace findNode;

//==============================================================

embeddisp::embeddisp(string filename) : m_outFileName(filename),
  m_c1(NULL),
  m_frame1(NULL), m_frame2(NULL), m_frame3(NULL)
{
  ThisName   = "embeddisp";
  EventNumber=0;

  m_pad[0] = m_pad[1] = m_pad[2] = NULL;
}

//==============================================================

embeddisp::~embeddisp() {
  if(m_c1    !=NULL) delete m_c1;
  for(int ipad=0; ipad<3; ipad++) {if(m_pad[ipad]!=NULL) delete m_pad[ipad];}
  if(m_frame1!=NULL) delete m_frame1;
  if(m_frame2!=NULL) delete m_frame2;
  if(m_frame3!=NULL) delete m_frame3;
}

//==============================================================

int embeddisp::Init(PHCompositeNode *topNode) {

  cout << "embeddisp::Init started..." << endl;
  //OutputNtupleFile = new TFile(OutputFileName.c_str(),"RECREATE");
  //cout << "embeddisp::Init: output file " << OutputFileName << " opened." << endl;

  cout << "embeddisp::Init ended." << endl;
  return 0;
}

//==============================================================
  
int embeddisp::InitRun(PHCompositeNode *topNode) {
  cout << "embeddisp::InitRun started..." << endl;

/*
  // check magnet current 
  RunHeader* runheader = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if (runheader==NULL) {
    cout << PHWHERE<< "Can't find runheader. " << endl;
    return ABORTRUN;
  }
  if(runheader->get_currentCentral()>0){
   m_fieldScale = 1.0;
  } else {
   m_fieldScale = -1.0;
  }
  cout<<"SvxCentralTrackReco::InitRun  fieldScale="<<m_fieldScale<<endl;
*/


  recoConsts *rc = recoConsts::instance();
  int runnumber = rc->get_IntFlag("RUNNUMBER");
  
  if (rc->FlagExist("SIMULATIONFLAG"))
    {
      if (rc->get_IntFlag("SIMULATIONFLAG") >= 1)
        {
          runnumber = -1;
        }
    }

  if( rc->get_IntFlag("LOOKATSIMFROMEMBED",1) == 0)
    {
      runnumber = rc->get_IntFlag("RUNNUMBER");
    }

  
  cout << "embeddisp using Run number " << runnumber << endl;
  {
    CrkGeometryObject cgo;

    int fUseSurvey = 0;
    switch (fUseSurvey)
      {
      case 0:             //Select to use survey or not using run number
        if(runnumber>0){
          cgo.UseSurvey();
        }
        break;
      case 1:             //force to use survey
        cgo.UseSurvey();
        break;
      case -1:            //force not to use survey
        break;
      default:
        if(runnumber>0){
          cgo.UseSurvey();
        }
      }

    cgo.Fetch(runnumber);

    for(int i=0; i<5120; i++){
      m_crkpmt[i] = 
        cgo.GetPmtPosition(cgo.IdToArm(i), cgo.IdToSide(i),
                           cgo.IdToSm(i),  cgo.IdToPmt( i));
    }
  }



  /////////////////////

  //TofwHit        *tofw = findNode::getClass<TofwHit>(       topNode, "TofwHit");
  //SvxClusterList *svx  = findNode::getClass<SvxClusterList>(topNode, "SvxClusterList");
  //HbdBlobList    *hbd  = findNode::getClass<HbdBlobList>(   topNode, "HbdBlobList");
  //cout<<"tofw "<<(tofw==NULL?"doesnot":"")<<" exist in merged topNode"<<endl;
  //cout<<"svx  "<<(svx ==NULL?"doesnot":"")<<" exist in merged topNode"<<endl;
  //cout<<"hbd  "<<(hbd ==NULL?"doesnot":"")<<" exist in merged topNode"<<endl;

  m_outfile = new TFile(m_outFileName.c_str(), "RECREATE");

  m_ntp_embed = new TNtuple("ntp_embed","embedding",
      "zvtx"
      ":pt:charge:mom:phi0:the0:dcqual"
      ":emcdz:emcdphi:pc3dz:pc3dphi"
      ":ecore:prob:n0:npe0:ch2npe0:disp"
      ":dca2d:dcaz:chi2ndf:hitptn"
      ":simpt:simmom:simphi0:simthe0"
      ":simecore:simprob:simn0:simnpe0:simch2npe0:simdisp"
      ":genpt:genmom:genphi0:genthe0:genvx:genvy:genvz:genpid:gendca2d"
    );


  cout << "embeddisp::InitRun ended." << endl;
  return 0;
}

//==============================================================

int embeddisp::process_event(PHCompositeNode *topNode) {

  static int first_event = 0;
  if(first_event==0){
    cout<<"embeddisp::process_event firstevent"<<endl;
    topNode->print();
    first_event++;
  }

  // event info
  EventHeader*  evthdr  = getClass<EventHeader>(topNode,"EventHeader");
  VtxOut*       vtxout  = getClass<VtxOut>(     topNode,"VtxOut");

  // track info
  DchHitLineTable     *dchit   = getClass<DchHitLineTable>(    topNode,"DchHitLineTable");
  PadCluster          *pc1     = getClass<PadCluster>(         topNode,"Pc1Cluster");
  PadCluster          *pc2     = getClass<PadCluster>(         topNode,"Pc2Cluster");
  PadCluster          *pc3     = getClass<PadCluster>(         topNode,"Pc3Cluster");
  CrkHit              *crkhit  = getClass<CrkHit>(             topNode,"CrkHit");
  emcClusterContainer *emc     = getClass<emcClusterContainer>(topNode,"emcClusterContainer");
  PHCentralTrack      *trk     = getClass<PHCentralTrack>(     topNode,"PHCentralTrack");
  PHTrackOut          *phtrk   = getClass<PHTrackOut>(         topNode,"PHTrackOut");
  PHDchTrackOut       *phdchtrk= getClass<PHDchTrackOut>(      topNode,"PHDchTrackOut");
  SvxCentralTrackList *svxcnt  = getClass<SvxCentralTrackList>(topNode,"SvxCentralTrackList");
  SvxClusterList      *svx     = getClass<SvxClusterList>(     topNode,"SvxClusterList");

  cout<<"event : "<<EventNumber<<"  "<<((evthdr!=NULL) ? evthdr->get_EvtSequence() : -1 ) <<endl;
  
  // topnode
  recoConsts*    rc = recoConsts::instance();
  Fun4AllServer* se = Fun4AllServer::instance();
  PHCompositeNode* mcnode   = se->topNode(rc->get_CharFlag("EMBED_MC_TOPNODE"));
  PHCompositeNode* realnode = se->topNode(rc->get_CharFlag("EMBED_REAL_TOPNODE"));

    
  //EventHeader    *evthdr_mc = getClass<EventHeader>(   mcnode,"EventHeader");
  DchHitLineTable     *dchit_mc  = getClass<DchHitLineTable>(    mcnode,"DchHitLineTable");
  PadCluster          *pc1_mc    = getClass<PadCluster>(         mcnode,"Pc1Cluster");
  PadCluster          *pc2_mc    = getClass<PadCluster>(         mcnode,"Pc2Cluster");
  PadCluster          *pc3_mc    = getClass<PadCluster>(         mcnode,"Pc3Cluster");
  CrkHit              *crkhit_mc = getClass<CrkHit>(             mcnode,"CrkHit");
  emcClusterContainer *emc_mc    = getClass<emcClusterContainer>(mcnode,"emcClusterContainer");
  PHCentralTrack      *trk_mc    = getClass<PHCentralTrack>(     mcnode,"PHCentralTrack");
  PHTrackOut          *phtrk_mc  = getClass<PHTrackOut>(         mcnode,"PHTrackOut");
    
  //EventHeader    *evthdr_real = getClass<EventHeader>(   realnode,"EventHeader");
  DchHitLineTable     *dchit_real  = getClass<DchHitLineTable>(    realnode,"DchHitLineTable");
  PadCluster          *pc1_real    = getClass<PadCluster>(         realnode,"Pc1Cluster");
  PadCluster          *pc2_real    = getClass<PadCluster>(         realnode,"Pc2Cluster");
  PadCluster          *pc3_real    = getClass<PadCluster>(         realnode,"Pc3Cluster");
  CrkHit              *crkhit_real = getClass<CrkHit>(             realnode,"CrkHit");
  emcClusterContainer *emc_real    = getClass<emcClusterContainer>(realnode,"emcClusterContainer");
  PHCentralTrack      *trk_real    = getClass<PHCentralTrack>(     realnode,"PHCentralTrack");
  PHTrackOut          *phtrk_real  = getClass<PHTrackOut>(         realnode,"PHTrackOut");
    
    
  if(dchit_mc   !=NULL) cout<<"dchit nhit : "<<dchit_mc->Entries()<<flush;
  if(dchit_real !=NULL) cout<<          " : "<<dchit_real->Entries()<<flush;
  if(dchit      !=NULL) cout<<          " : "<<dchit->Entries()<<flush;
  cout<<endl;
  if(pc1_mc   !=NULL) cout<<"pc1 nhit : "<<pc1_mc->get_PadNCluster()<<flush;
  if(pc1_real !=NULL) cout<<        " : "<<pc1_real->get_PadNCluster()<<flush;
  if(pc1      !=NULL) cout<<        " : "<<pc1->get_PadNCluster()<<flush;
  cout<<endl;
  if(pc2_mc   !=NULL) cout<<"pc2 nhit : "<<pc2_mc->get_PadNCluster()<<flush;
  if(pc2_real !=NULL) cout<<        " : "<<pc2_real->get_PadNCluster()<<flush;
  if(pc2      !=NULL) cout<<        " : "<<pc2->get_PadNCluster()<<flush;
  cout<<endl;
  if(pc3_mc   !=NULL) cout<<"pc3 nhit : "<<pc3_mc->get_PadNCluster()<<flush;
  if(pc3_real !=NULL) cout<<        " : "<<pc3_real->get_PadNCluster()<<flush;
  if(pc3      !=NULL) cout<<        " : "<<pc3->get_PadNCluster()<<flush;
  cout<<endl;
  if(crkhit_mc  !=NULL) cout<<"crk nhit : "<<crkhit_mc->get_CrkNHit()<<flush;
  if(crkhit_real!=NULL) cout<<        " : "<<crkhit_real->get_CrkNHit()<<flush;
  if(crkhit     !=NULL) cout<<        " : "<<crkhit->get_CrkNHit()<<flush;
  cout<<endl;
  cout<<"cnt ntrk : "<<((trk_mc  !=NULL) ? trk_mc->get_npart()   : 999999)<<flush;
  cout<<        " : "<<((trk_real!=NULL) ? trk_real->get_npart() : 999999)<<flush;
  cout<<        " : "<<((trk     !=NULL) ? trk->get_npart()      : 999999)<<flush;
  cout<<endl;
  if(phtrk_mc  !=NULL) cout<<"phtrk ntrk : "<<phtrk_mc->get_PHNTrack()<<flush;
  if(phtrk_real!=NULL) cout<<          " : "<<phtrk_real->get_PHNTrack()<<flush;
  if(phtrk     !=NULL) cout<<          " : "<<phtrk->get_PHNTrack()<<flush;
  cout<<endl;
  cout<<        "svxcnttrk : "<<((svxcnt  !=NULL) ? svxcnt->get_nCentralTracks()      : 999999)<<endl;
  cout<<        "svxcls    : "<<((svx     !=NULL) ? svx->get_nClusters()      : 999999)<<endl;

  //////////////////////////////////////////////////////////////////////
  PHPointerList<PHEmbedMcRecoTrack>  *embedtrk = getClass< PHPointerList<PHEmbedMcRecoTrack> >(topNode, "PHEmbedMcRecoTrack");

  // mapping between PHCentralTrack and SvxCentralTrack
  remapCntSvx(trk, svxcnt);

  if(embedtrk==NULL) { cout<<" no PHEmbedMcRecoTrack"<<endl;}
  else {
    cout<<"PHEmbedMcRecoTrack: length = "<<embedtrk->length()<<endl;
    for(unsigned int i=0; i<embedtrk->length(); i++){
      PHEmbedMcRecoTrack *embed = (*embedtrk)[i];
      cout<<"   id_g, _s, _r = "<<embed->get_dctrkidG()<<" "<<embed->get_dctrkidS()<<" "<<embed->get_dctrkidR()<<" ";
      cout<<" mom="<<embed->get_momG()<<" "<<embed->get_momR()<<endl;

      float momr = -9999, momr1=-9999;
      int dcidr = embed->get_dctrkidR();
      if(0<=dcidr&&dcidr<(int)trk->get_npart()){
         PHSnglCentralTrack* sngl = trk->get_track(dcidr);
         momr = sngl->get_mom();
         momr1= phdchtrk->get_momentum(dcidr);
         cout<<"  reco mom= "<<momr<<" "<<momr1<<endl;
      }
      else {
        cout<<"out of range : "<<dcidr<<endl;
      }
    }

    int Ncnt = trk->get_npart();

    // ntp fill
    unsigned int Nembed = (unsigned int)(0.5*embedtrk->length());

    for(unsigned int itrk=0; itrk<Nembed; itrk++){
      PHEmbedMcRecoTrack *embed = (*embedtrk)[itrk];
      //int idG = embed->get_dctrkidG();
      int idR = embed->get_dctrkidR();
      //int idS = embed->get_dctrkidS();

      if(idR<0||Ncnt<=idR) {
        cout<<" RealtrackID is out of range = "<<idR<<" : "<<Ncnt<<endl;
        continue;
      }

      PHSnglCentralTrack* sngl       = trk->get_track(idR);
      SvxCentralTrack*    svxcntsngl = m_vsvxcnt[idR];

      float ntp[100];
      for(int i=0; i<100; i++) ntp[i] = -9999.;

      float charge = sngl->get_charge();

      ntp[ 0] = vtxout->get_ZVertex(); //zvtx
      ntp[ 1] = sngl->get_mom()*sin(sngl->get_the0());
      ntp[ 2] = charge;
      ntp[ 3] = sngl->get_mom();
      ntp[ 4] = sngl->get_phi0();
      ntp[ 5] = sngl->get_the0();
      ntp[ 6] = sngl->get_quality();
      ntp[ 7] = sngl->get_emcdz();
      ntp[ 8] = sngl->get_emcdphi();
      ntp[ 9] = sngl->get_pc3dz();
      ntp[10] = sngl->get_pc3dphi();
      ntp[11] = sngl->get_ecore();
      ntp[12] = sngl->get_prob();
      ntp[13] = sngl->get_n0();
      ntp[14] = sngl->get_npe0();
      ntp[15] = sngl->get_chi2()/sngl->get_npe0();
      ntp[16] = sngl->get_disp();
      
      if(svxcntsngl!=NULL){
        float ndf        = svxcntsngl->getNDF();
        float chisqr     = svxcntsngl->getChiSquare();
        ntp[17] = svxcntsngl->getDCA2D(); // dca2d
        ntp[18] = svxcntsngl->getDCAZ();  // dcaz
        ntp[19] = (ndf>0) ? chisqr/ndf : -9999.; // chi2ndf
        ntp[20] = svxcntsngl->getHitPattern(); // hitptn
      }

      ntp[21] = embed->get_ptS();
      ntp[22] = embed->get_momS();
      ntp[23] = embed->get_phi0S();
      ntp[24] = embed->get_the0S();
      ntp[25] = embed->get_emcecoreS();
      ntp[26] = embed->get_emcprobphotS();
      ntp[27] = embed->get_crknpmt0S();
      ntp[28] = embed->get_crknpe0S();
      ntp[29] = embed->get_crkchi2S()/embed->get_crknpe0S();
      ntp[30] = embed->get_crkdispS();

      // use gen info
      float genpt  = embed->get_ptG();
      float genphi = embed->get_phi0G();
      float genthe = embed->get_the0G();
      float genvx  = embed->get_xvtxG();
      float genvy  = embed->get_yvtxG();
      float genvz  = embed->get_zvtxG();

      float beamxy[2] = {0., 0.};
      float fieldPol  = 1.0;
      float gendca[3], gendca2d;

      calcDCA_BCbyCircleProjection
            (
              genpt, genphi, genthe,  // pt in xy-plane, phi, theta at inner most layer
              charge,                 // charge of the track
              genvx, genvy, genvz,    // hit position at inner most layer
              beamxy[0], beamxy[1],   // beam center
              fieldPol,               //
              &gendca[0], &gendca[1], &gendca[2], // dca position
              &gendca2d               // return
            );

      ntp[31] = embed->get_ptG();
      ntp[32] = embed->get_momG();
      ntp[33] = embed->get_phi0G();
      ntp[34] = embed->get_the0G();
      ntp[35] = embed->get_xvtxG();
      ntp[36] = embed->get_yvtxG();
      ntp[37] = embed->get_zvtxG();
      ntp[38] = embed->get_partidG();
      ntp[39] = gendca2d;

      m_ntp_embed->Fill(ntp);
    }
  }

  ///////////////////////////////////
  // dchit
  {
    vector<PHPoint>* vhit[3] = {&m_vdchhit, &m_vdchhit_real, &m_vdchhit_mc};
    DchHitLineTable* dch[3]  = { dchit,      dchit_real,      dchit_mc};

    for(int i=0; i<3; i++){
      vhit[i]->clear(); 
      for( int ihit=0; ihit<dch[i]->Entries(); ihit++){
        PHPoint pos = dch[i]->getXYZ(ihit);
        //cout<<"pos : "<<pos.getX()<<" "<<pos.getY()<<" "<<pos.getZ()<<endl;
        vhit[i]->push_back(pos);
      }
    }
  }

  ///////////////////////////////////
  // pc
  {
    for(int ipc=0; ipc<3; ipc++){
      vector<PHPoint>* vhit[3];
      if(ipc==0)      { vhit[0] = &m_vpc1hit; vhit[1] = &m_vpc1hit_real; vhit[2] = &m_vpc1hit_mc;}
      else if(ipc==1) { vhit[0] = &m_vpc2hit; vhit[1] = &m_vpc2hit_real; vhit[2] = &m_vpc2hit_mc;}
      else            { vhit[0] = &m_vpc3hit; vhit[1] = &m_vpc3hit_real; vhit[2] = &m_vpc3hit_mc;}

      PadCluster*      pc[3]  ;
      if(ipc==0)      {pc[0] = pc1; pc[1] = pc1_real; pc[2] = pc1_mc;}
      else if(ipc==1) {pc[0] = pc2; pc[1] = pc2_real; pc[2] = pc2_mc;}
      else            {pc[0] = pc3; pc[1] = pc3_real; pc[2] = pc3_mc;}

      for(int i=0; i<3; i++){
        vhit[i]->clear(); 
        for(unsigned int ihit=0; ihit<pc[i]->get_PadNCluster(); ihit++){
          PHPoint pos(pc[i]->get_xyz(ihit, 0), pc[i]->get_xyz(ihit, 1), pc[i]->get_xyz(ihit, 2));
          //cout<<"pos : "<<pos.getX()<<" "<<pos.getY()<<" "<<pos.getZ()<<endl;
          vhit[i]->push_back(pos);
        }
      }
    }
  }
  //
  ///////////////////////////////////
  // emc
  {
    vector<PHPoint>*     vhit[3]   = {&m_vemchit, &m_vemchit_real, &m_vemchit_mc};
    emcClusterContainer* emchit[3] = { emc,        emc_real,        emc_mc};

    for(int i=0; i<3; i++){
      vhit[i]->clear(); 
      for(unsigned int ihit=0; ihit<emchit[i]->size(); ihit++){
        emcClusterContent* clus = emchit[i]->getCluster(ihit);
        PHPoint pos(clus->x(), clus->y(), clus->z());
        //cout<<"pos : "<<pos.getX()<<" "<<pos.getY()<<" "<<pos.getZ()<<endl;
        vhit[i]->push_back(pos);
      }
    }
  }
  //
  ///////////////////////////////////
  // crk
  {
    vector<PHPoint>* vhit[3] = {&m_vcrkhit, &m_vcrkhit_real, &m_vcrkhit_mc};
    CrkHit*          crk[3]  = { crkhit,     crkhit_real,     crkhit_mc};

    for(int i=0; i<3; i++){
      vhit[i]->clear(); 
      for(unsigned int ihit=0; ihit<crk[i]->get_CrkNHit(); ihit++){
        PHPoint pos = m_crkpmt[crk[i]->get_pmt(ihit)];
        //cout<<"pos : "<<pos.getX()<<" "<<pos.getY()<<" "<<pos.getZ()<<endl;
        vhit[i]->push_back(pos);
      }
    }
  }

/*
  //  cout<<"aa"<<endl;
  //  PHTrackOut         *proj = getClass<PHTrackOut>(topNode,"PHTrackOut");
  //  CrkRing             *crk = getClass<CrkRing>(topNode,"CrkRing");
  //  emcClusterContainer *emc = getClass<emcClusterContainer>(topNode,"emcClusterContainer");
  //  DchTrack            *dch = getClass<DchTrack>(topNode,"DchTrack");
  //  CglTrack            *cgl = getClass<CglTrack>(topNode,"CglTrack");
  //TrigRunLvl1 *lvl1    = getClass<TrigRunLvl1>(topNode,"TrigRunLvl1");
  PHGlobal    *global       = getClass<PHGlobal>(topNode,"PHGlobal");
  RunHeader   *run          = getClass<RunHeader>(topNode,"RunHeader");
  TrigLvl1    *trglvl1      = getClass<TrigLvl1>(topNode,"TrigLvl1");
  BbcOut *bbc               = getClass<BbcOut>(topNode,"BbcOut");

  PreviousEvent *peve = getClass<PreviousEvent>(topNode, "PreviousEvent");
  if(EventNumber==0)
    if(peve==NULL){ cout<<"no PreviousEvent"<<endl; }


  int   RunNumber = (run!=NULL)          ? run->get_RunNumber()            : -9999;
  int   EventNum  = (event_header!=NULL) ? event_header->get_EvtSequence() : -9999;
  float cent      = (global!=NULL)       ? global->getCentrality()         : -9999;

  if(EventNumber==0) {
    cout << "run = "<<RunNumber<<" event = "<<EventNum<<endl;
  }
*/


/*

  SvxRawhitList  *svxrawhit = getClass<SvxRawhitList>(topNode,"SvxRawhitList");
  SvxClusterList *svxsel    = getClass<SvxClusterList>(topNode,"SvxSelectedClusterList");
  SvxSegmentList *svxtracks = getClass<SvxSegmentList>(topNode,"SvxSegmentList");
  SvxRawhitClusterList 
             *svxrawhitclus = getClass<SvxRawhitClusterList>(topNode,"SvxRawhitClusterList");

  SvxCentralTrackList
             *svxcnttrklist = getClass<SvxCentralTrackList>(topNode,"SvxCentralTrackList");


  McEvalSingleList *mctrk   = getClass<McEvalSingleList>(topNode,"McSingle");
  SvxGhitList      *svxghit = getClass<SvxGhitList>(topNode,"SvxGhitList");
  fkinWrapper      *fkinW   = getClass<fkinWrapper>(topNode, "fkin");
  SvxGhitRawhitList *svxghitrawlist = getClass<SvxGhitRawhitList>(topNode,"SvxGhitRawhitList");
  if(EventNumber==0){
    if(m_simmode){
      if(mctrk==NULL){
        cout<<"No McSingle Info in Nodetree"<<endl;
      }
      if(svxghit==NULL){
        cout<<"No SvxGhitList Info in Nodetree"<<endl;
      }
      if(fkinW==NULL){
        cout<<"No fkinWrapper Info in Nodetree"<<endl;
      }
      if(svxghitrawlist==NULL){
        cout<<"No SvxGhitRawhitList Info in Nodetree"<<endl;
      }
    }
  }

  int ntrk       = (trk!=NULL) ? trk->get_npart() : 0;
  int nsvxtracks = (svxtracks!=NULL) ? svxtracks->get_nSegments() : -9999;
  int nsvx       = (svx!=NULL)? svx->get_nClusters() : -9999;
  int nsvxrawhit = (svxrawhit!=NULL) ? svxrawhit->get_nRawhits() : -9999;
  //int nsvxrawhitclus = 0;
  int nkalfit = 0;
  float zvtx = (bbc!=NULL) ? bbc->get_VertexPoint() : -999.0;
  if(zvtx<-900){
    zvtx = (vtxout!=NULL) ? vtxout->get_ZVertex("BBC") : -999.0;
  }
  float bbcq = 50;
  if(bbc!=NULL) {
    bbcq = bbc->get_ChargeSum(Bbc::North) + bbc->get_ChargeSum(Bbc::South);
  }
  else if(global!=NULL) {
    bbcq = global->getBbcChargeN() + global->getBbcChargeS();
  }
  //if(bbcq>200 || bbcq<10)return 1; 

  float vtx_z = (vtxout!=NULL) ? vtxout->get_ZVertex() : -999.;
  if(m_goodOnlyFlag&&(fabs(vtx_z)>GOODZVTX)){ // skip this event
    if(verbosity>0) cout << "embeddisp  skip this event : run = "<<RunNumber<<" event = "<<EventNum<<endl;
    return EVENT_OK;
  }

  

  if(init_ana==0) {
    cout<<"init ana"<<endl;
    topNode->print(); 
    init_ana=1;
 
    if(svx!=NULL)          svx->identify();          else cout<<"no SvxCluster object"<<endl; 
    if(svxtracks!=NULL)    svxtracks->identify();    else cout<<"no SvxSegment object"<<endl;
    if(svxrawhitclus!=NULL)svxrawhitclus->identify();else cout<<"no SvxRawhitClus object"<<endl;
    if(svxrawhit!=NULL)    svxrawhit->identify();    else cout<<"no SvxRawhit object"<<endl;
    if(trglvl1!=NULL)      trglvl1->identify();      else cout<<"no TrigLvl1 object"<<endl;

    cout<<endl<<endl;
    cout<<"init ana ends"<<endl;
  }

  if(verbosity>1){
    if(EventNumber%1==0) {
      cout << "------- Event # " << EventNumber << " nsvxrawhit= " << nsvxrawhit << " nsvx= " << nsvx << " ntrack= " << ntrk << " nkalfit= " << nkalfit << " nsvxtracks=" << nsvxtracks << endl;
    }
  }

  if((EventNumber%1000)==0){
    cout << "------- Event # " << EventNumber << endl;
  }
*/

  EventNumber++;
  return 0;
}

int embeddisp::End(PHCompositeNode *topNode)
{
  cout << "Writing out..." << endl;
  m_outfile->Write();
  cout << "Closing output file..." << endl;
  m_outfile->Close();
  delete m_outfile;
  
  return 0;
}

void embeddisp::calcDCA_BCbyCircleProjection
(
 float pt, float phi, float the,  // pt in xy-plane, phi, theta at inner most layer
 int charge,                      // charge of the track
 float hx, float hy, float hz,    // hit position at inner most layer
 float vx, float vy,              // beam center
 float fieldPolarity,             //
 float* dx, float* dy, float* dz, // dca position
 float* d2dca_bc)                 // return
{
  static const float B = 0.90;
  static const float b = 0.003*B;

  // vector to rotation center at hit1 
  float R  = pt/b;
  float pz = pt/tan(the);

  float b_sign = (fieldPolarity>0) ? -1.0 : 1.0;
  float dir = ( (b_sign)*charge>0. ) ? -1.0 : 1.0;
  float cx = hx - dir*R*sin(phi);
  float cy = hy + dir*R*cos(phi);

  // L is a distance btw the rotation center and primary vtx  
  // psi is a angle of the vector starting from the rotation center to primary vertex
  float L = sqrt((vx-cx)*(vx-cx) + (vy-cy)*(vy-cy));
  float psi = atan2((vy-cy), (vx-cx));

  // DCA point
  *dx = cx + R*cos(psi);
  *dy = cy + R*sin(psi);

  float dphi = phi - dir*0.5*M_PI - psi;
  if ( dphi>M_PI*1.5 ) dphi -= M_PI*2.;
  else if ( dphi<-M_PI*1.5 ) dphi += M_PI*2.;
  float dzdphi = dir*pz*R/pt;
  *dz = hz + dphi*dzdphi;

  // DCA value
  *d2dca_bc = b_sign*(charge*(R - L));
}

void embeddisp::remapCntSvx(PHCentralTrack* cnt, SvxCentralTrackList* svxcnt)
{
  m_vsvxcnt.clear();
  if(cnt==NULL || svxcnt==NULL){
    if(cnt==NULL)    cout<<"No PHCentralTrack available"<<endl;
    if(svxcnt==NULL) cout<<"No SvxCentralTrackList available"<<endl;
    return;
  }

  int npart = cnt->get_npart();
  m_vsvxcnt.assign(npart, (SvxCentralTrack*)NULL);
  
  if( svxcnt!=NULL ) {
    for(int i=0; i<svxcnt->get_nCentralTracks(); i++){
      SvxCentralTrack *svxcnttrk = svxcnt->getCentralTrack(i);
      if(svxcnttrk==NULL) {
        cout<<"no svxcnt"<<endl;
        continue;
      }
      
      int itrk = svxcnttrk->getDchIndex();
      m_vsvxcnt[itrk] = svxcnttrk;
    }
  }
}
  
void embeddisp::drawCanvas()
{
  if(m_c1==NULL){
    m_c1 = new TCanvas("c1", "c1", 800, 800);
    m_pad[0] = new TPad("pad_0", "", 0.21, 0.21, 0.80, 0.79);
    m_pad[0]->SetLeftMargin(0.05);
    m_pad[0]->SetRightMargin(0.05);
    m_pad[1] = new TPad("pad_1", "", 0.01, 0.01, 0.20, 0.99);
    m_pad[2] = new TPad("pad_2", "", 0.81, 0.01, 0.99, 0.99);

    m_frame1 = new TH2F("frame1", "", 100, -600, 600, 100, -600, 600);
    //m_frame2 = new TH2F("frame2", "", 100, -400, 400, 100,  1.7, 4.8);
    m_frame2 = new TH2F("frame2", "", 100, -400, 400, 100, -1.5, 1.5);
    m_frame3 = new TH2F("frame3", "", 100, -400, 400, 100, -1.5, 1.5);
  }

  m_c1->cd();
  m_pad[0]->Draw();
  m_pad[1]->Draw();
  m_pad[2]->Draw();

  m_pad[0]->cd(); m_frame1->Draw();
  m_pad[1]->cd(); m_frame2->Draw();
  m_pad[2]->cd(); m_frame3->Draw();
}

void embeddisp::draw()
{

  if(m_c1==NULL){ drawCanvas(); }

  // reset
  vector<TMarker*>::iterator itr_m;
  for(itr_m=m_vdrawmarker.begin(); itr_m!=m_vdrawmarker.end(); ++itr_m){
    TMarker *m = *itr_m;
    delete m;
  }
  m_vdrawmarker.clear();
  

  // draw
  { // dchit
    m_pad[0]->cd();
    vector<PHPoint>::iterator itr_dch;

    vector<PHPoint>* vdch[3] = {&m_vdchhit_real, &m_vdchhit_mc, &m_vdchhit}; 
    int   color[3] = {3,2,1};
    float size[3]  = {0.5, 0.5, 0.2};
    for(int i=0;i<3; i++){
      for(itr_dch=vdch[i]->begin(); itr_dch!=vdch[i]->end(); ++itr_dch){
        PHPoint pos = *itr_dch;
        TMarker *m = new TMarker(pos.getX(), pos.getY(), 20);
        m->SetMarkerColor(color[i]);
        m->SetMarkerSize(size[i]);
        m->Draw();
        m_vdrawmarker.push_back(m);
      }
    }
  }

  { // pc1hit
    m_pad[0]->cd();
    vector<PHPoint>::iterator itr_hit;

    vector<PHPoint>* vhit[3];
    for(int ipc=0; ipc<3; ipc++){
      if(     ipc==0) { vhit[0] = &m_vpc1hit_real; vhit[1] = &m_vpc1hit_mc; vhit[2] = &m_vpc1hit; }
      else if(ipc==1) { vhit[0] = &m_vpc2hit_real; vhit[1] = &m_vpc2hit_mc; vhit[2] = &m_vpc2hit; }
      else            { vhit[0] = &m_vpc3hit_real; vhit[1] = &m_vpc3hit_mc; vhit[2] = &m_vpc3hit; }

      int   color[3] = {3,2,1};
      float size[3]  = {1.0, 1.0, 0.5};
      for(int i=0;i<3; i++){
        for(itr_hit=vhit[i]->begin(); itr_hit!=vhit[i]->end(); ++itr_hit){
          PHPoint pos = *itr_hit;
          TMarker *m = new TMarker(pos.getX(), pos.getY(), 24);
          m->SetMarkerColor(color[i]);
          m->SetMarkerSize(size[i]);
          m->Draw();
          m_vdrawmarker.push_back(m);
        }
      }
    }
  }

  { // emc
    m_pad[0]->cd();
    vector<PHPoint>::iterator itr_hit;

    vector<PHPoint>* vhit[3] = {&m_vemchit_real, &m_vemchit_mc, &m_vemchit}; 
    int   color[3] = {3,2,1};
    float size[3]  = {1.0, 1.0, 0.7};
    for(int i=0;i<3; i++){
      for(itr_hit=vhit[i]->begin(); itr_hit!=vhit[i]->end(); ++itr_hit){
        PHPoint pos = *itr_hit;
        TMarker *m = new TMarker(pos.getX(), pos.getY(), 28);
        m->SetMarkerColor(color[i]);
        m->SetMarkerSize(size[i]);
        m->Draw();
        m_vdrawmarker.push_back(m);
      }
    }
  }

  { // crk
    vector<PHPoint>::iterator itr_hit;

    vector<PHPoint>* vhit[3] = {&m_vcrkhit_real, &m_vcrkhit_mc, &m_vcrkhit}; 
    int   color[3] = {3,2,1};
    float size[3]  = {1.0, 1.0, 0.7};
    for(int i=0;i<3; i++){
      for(itr_hit=vhit[i]->begin(); itr_hit!=vhit[i]->end(); ++itr_hit){
        PHPoint pos = *itr_hit;
        float phi = atan2(pos.getY(), pos.getX());
        if(phi< - 0.5*TMath::Pi()) phi+= 2*TMath::Pi();

        int ipad = (phi<1.5) ? 2 : 1;
        m_pad[ipad]->cd();

        if(pos.getX()<0) { phi = atan2(pos.getY(), -pos.getX()); }

        TMarker *m = new TMarker(pos.getZ(), phi, 28);
        m->SetMarkerColor(color[i]);
        m->SetMarkerSize(size[i]);
        m->Draw();
        m_vdrawmarker.push_back(m);
      }
    }
  }
}
