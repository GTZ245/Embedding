#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include "gsl/gsl_rng.h"
#include <math.h>

#include "phool.h"
#include "PHTypedNodeIterator.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "Fun4AllReturnCodes.h"

#include "svxcntclusdisp.h"

#include "PHCentralTrack.h"
#include "VtxOut.h"
#include "Bbc.hh"
#include "BbcOut.h"

#include "SvxClusterList.h"
#include "SvxCluster.h"
#include "SvxCentralTrackList.h"
#include "SvxCentralTrack.h"
#include "SvxCentralTrackv1.h" // for ClusterInfo

#include "SvxClusterv4.h"
#include "SvxClusterListv4.h"
#include "SvxClusterContainer.h"

#include "svxAddress.hh"
#include "svxDetectorGeo.hh"
#include "SvxSensor.h"

#include "getClass.h"

#include "TCanvas.h"
#include "TH2.h"
#include "TPolyLine.h"
#include "TText.h"
#include "TMarker.h"
#include "TArc.h"
#include "TLine.h"

using namespace std;
using namespace findNode;


//==============================================================

svxcntclusdisp::svxcntclusdisp() {
  ThisName = "svxcntclusdisp";
  init_ana=0;
  EventNumber=0;

  m_svxAdrObj = NULL;
  m_svxGeo    = NULL;

  m_isInitCanvas = false;

  m_c1 = NULL;
  m_frame = NULL;
  for(int i=0; i<4; i++){
    for(int j=0; j<24; j++) {
      m_pl[i][j] = NULL; 
      m_txt[i][j] = NULL; 
    }
  }

  //
  m_cctestList = NULL;
  m_cctestCont = NULL;

}

//==============================================================

svxcntclusdisp::svxcntclusdisp(string filename) {
  ThisName = "svxcntclusdisp";
  init_ana=0;
  EventNumber=0;

  m_svxAdrObj = NULL;
  m_svxGeo    = NULL;

  m_isInitCanvas = false;

  m_c1 = NULL;
  m_frame = NULL;
  for(int i=0; i<4; i++){
    for(int j=0; j<24; j++) {
      m_pl[i][j] = NULL; 
      m_txt[i][j] = NULL; 
    }
  }

  //
  m_cctestList = NULL;
  m_cctestCont = NULL;
}

//==============================================================

svxcntclusdisp::~svxcntclusdisp() {
  if(m_c1!=NULL) delete m_c1;
  if(m_frame!=NULL) delete m_frame;
  for(int i=0; i<4; i++){
    for(int j=0; j<24; j++) {
      if(m_pl[i][j]!=NULL) delete m_pl[i][j]; 
      if(m_txt[i][j]!=NULL) delete m_txt[i][j]; 
    }
  }

  if(m_cctestList!=NULL) delete m_cctestList;
  if(m_cctestCont!=NULL) delete m_cctestCont;
}

//==============================================================

int svxcntclusdisp::Init(PHCompositeNode *topNode) {

  cout << "svxcntclusdisp::Init started..." << endl;
  cout << "svxcntclusdisp::Init ended." << endl;
  return 0;
}

//==============================================================
  
int svxcntclusdisp::InitRun(PHCompositeNode *topNode) {
  cout << "svxcntclusdisp::InitRun started..." << endl;
  svxAddress* address = findNode::getClass<svxAddress>(topNode, "svxAddress");
  if ( address == NULL) {
    if(verbosity>0) { cout << PHWHERE<< "Can't find svxAddress. " << endl; }
    return ABORTRUN;
  }
  m_svxAdrObj = address;

  svxDetectorGeo* detgeo = findNode::getClass<svxDetectorGeo>(topNode, "svxDetectorGeo");
  if ( detgeo == NULL) {
    if(verbosity>0) { cout << PHWHERE<< "Can't find svxDetectorGeo. " << endl; }
    return ABORTRUN;
  }
  m_svxGeo = detgeo;

  initCanvas();

  cout << "svxcntclusdisp::InitRun ended." << endl;
  return 0;
}

//==============================================================

int svxcntclusdisp::process_event(PHCompositeNode *topNode) {


  //BbcOut *bbc    = getClass<BbcOut>(topNode,"BbcOut");
  //VtxOut *vtxout = getClass<VtxOut>(topNode,"VtxOut");

  //int RunNumber = (run!=NULL) ? run->get_RunNumber() : -9999;
  //int EventNum = (event_header!=NULL) ? event_header->get_EvtSequence() : -9999;
  //cout << "run = "<<RunNumber<<" event = "<<EventNum<<endl;


  PHCentralTrack *trk = getClass<PHCentralTrack>(topNode,"PHCentralTrack");
  int ntrk = (trk!=NULL) ? trk->get_npart() : 0;

  SvxClusterList *svx = getClass<SvxClusterList>(topNode,"SvxClusterList");

  //SvxCentralTrackList *svxcnttrklist = getClass<SvxCentralTrackList>(topNode,"SvxCentralTrackList");

  int nsvx = svx->get_nClusters();
  //int nsvxrawhitclus = 0;
  //int nkalfit = 0;
  //float zvtx = (bbc!=NULL) ? bbc->get_VertexPoint() : -999.0;
  //float bbcq = (bbc!=NULL) ? (bbc->get_ChargeSum(Bbc::North) + bbc->get_ChargeSum(Bbc::South)) : 50.0;
  

  if(init_ana==0) {
    cout<<"init ana"<<endl;
    topNode->print(); 
    init_ana=1;
 
    if(svx!=NULL)          svx->identify();          else cout<<"no SvxCluster object"<<endl; 

    cout<<endl<<endl;
    cout<<"init ana ends"<<endl;
  }

  if(EventNumber%1==0) {
    cout << "------- Event # " << EventNumber << " nsvx= " << nsvx << " ntrack= " << ntrk << endl;
  }


//------------------------ VTX information --------------------------

  fillCluster(svx);
  fillCntTrack(trk);
  //fillCentralTrack(bbcq, zvtx, svxcnttrklist, trk);


  EventNumber++;
  return 0;
}

void svxcntclusdisp::fillCluster(SvxClusterList* clslist){
  // reset vector

  vector<TMarker*>::iterator itr;
  for(itr=m_vClusPos.begin(); itr!=m_vClusPos.end(); ++itr){
    TMarker *tmp = *itr;
    delete tmp;
  }
  m_vClusPos.clear();

  int ncluster = clslist->get_nClusters();
  for(int icls=0; icls<ncluster; icls++){
    SvxCluster *cls = clslist->get_Cluster(icls);
    float x = cls->get_xyz_global(0);
    float y = cls->get_xyz_global(1);
    //cout<<" clus "<<icls<<" "<<x<<" "<<y<<endl;
    TMarker *m = new TMarker(x, y, 20);
    m->SetMarkerSize(1.0);
    m->SetMarkerColor(2);
    m_vClusPos.push_back(m);
  }
}

void svxcntclusdisp::fillCntTrack(PHCentralTrack *trk){
  // reset vector
  vector<TArc*>::iterator itr;
  for(itr=m_vTrkArc.begin(); itr!=m_vTrkArc.end(); ++itr){
    TArc *tmp = *itr;
    delete tmp;
  }
  m_vTrkArc.clear();

  vector<TLine*>::iterator itr1;
  for(itr1=m_vTrkLine.begin(); itr1!=m_vTrkLine.end(); ++itr1){
    TLine *tmp = *itr1;
    delete tmp;
  }
  m_vTrkLine.clear();

  vector<TLine*>::iterator itr2;
  for(itr2=m_vRngLine.begin(); itr2!=m_vRngLine.end(); ++itr2){
    TLine *tmp = *itr2;
    delete tmp;
  }
  m_vRngLine.clear();


  //
  static const float B = 0.9; // T
  static const float L = 0.17; // m Radius of B3

  float pStart[2]={0,0};

  int ntrk = (trk!=NULL) ? trk->get_npart() : 0;
  for(int itrk=0; itrk<ntrk; itrk++){
    float mom  = trk->get_mom(itrk);
    float charge  = trk->get_charge(itrk);
    float phi0 = trk->get_phi0(itrk);
    float the0 = trk->get_the0(itrk);
    float dcarm = trk->get_dcarm(itrk);

    cout<<"itrk : "<<itrk<<" pt="<<mom*sin(the0)<<" phi0="<<phi0<<endl;

    float pt = mom*sin(the0);

    float R = 100.*(pt/0.3*B);  // cm

    float cx = pStart[0] + charge * sin(phi0)*R;
    float cy = pStart[1] - charge * cos(phi0)*R;

    float draw_phi = atan2(-cy, -cx);

    //cout<<"  "<<cx<<" "<<cy<<" "<<R<<" "<<draw_phi<<endl;

    TArc *a = new TArc(cx, cy, R, 0, draw_phi*180./3.1415);
    a->SetLineColor(2);
    a->SetFillStyle(1);
    a->SetNoEdges();
    m_vTrkArc.push_back(a);

    float lx = pStart[0] + 20. * cos(phi0);
    float ly = pStart[0] + 20. * sin(phi0);
    TLine *l = new TLine(pStart[0], pStart[1], lx, ly);
    l->SetLineColor(4);
    m_vTrkLine.push_back(l);

    // range
    static const float Bdir = -1.0;
    static const float factor = 0.3;
    float dphi = factor*Bdir*charge*0.3*B*L/(2*pt); //= L/2R = Lx0.3B/2pT
    cout<<pt<<" "<<dcarm<<" "<<phi0<<" "<<dphi<<" "<<charge<<endl;

    for(int i=0; i<2; i++){
      float rng = (i==0) ? 0.2 : -0.2;
      lx = pStart[0] + 20. * cos(phi0+rng+dphi);
      ly = pStart[0] + 20. * sin(phi0+rng+dphi);
      l = new TLine(pStart[0], pStart[1], lx, ly);
      l->SetLineColor(6);
      m_vRngLine.push_back(l);
    }
  }

}

void svxcntclusdisp::fillCentralTrack(float bbcq, float zvtx, SvxCentralTrackList* svxcnttrklist, PHCentralTrack *cntlist){
  if(svxcnttrklist==NULL || cntlist==NULL){
    return;
  }

  float ary[100];
  for(int i=0; i<100; i++){ ary[i] = -9999.0; }

  int ntrk = svxcnttrklist->get_nCentralTracks();
  for(int itrk=0; itrk<ntrk; itrk++){
    SvxCentralTrack *svxtrk = svxcnttrklist->getCentralTrack(itrk);

    int cntidx=svxtrk->getDchIndex();
    
    ary[ 0] = EventNumber;
    ary[ 1] = bbcq;
    ary[ 2] = zvtx;
    ary[ 3] = cntlist->get_mom(cntidx);
    ary[ 4] = cntlist->get_phi0(cntidx);
    ary[ 5] = cntlist->get_the0(cntidx);
    ary[ 6] = cntlist->get_zed(cntidx);
    ary[ 7] = svxtrk->getQuality();
    ary[ 8] = svxtrk->getChiSquare();
    ary[ 9] = svxtrk->getNDF();
    ary[10] = svxtrk->getDCA2D();
    ary[11] = svxtrk->getClosestApproach(0);
    ary[12] = svxtrk->getClosestApproach(1);
    ary[13] = svxtrk->getClosestApproach(2);
    ary[14] = svxtrk->getNhits();
    cout<<" NassociatedCluster : "<<svxtrk->getNhits()<<ary[0]<<endl;

    bool isN4 = (svxtrk->getNhits()>3);
    if(isN4){
      cout<<"Particle : ";
      cout<<cntlist->get_mom(cntidx)<<", ";
      cout<<cntlist->get_charge(cntidx)<<", ";
      cout<<cntlist->get_phi0(cntidx)<<", ";
      cout<<cntlist->get_the0(cntidx)<<", ";
      cout<<"Zvtx= "<<zvtx<<", ";
      cout<<endl;
      
    }
    for(int ihit=0; ihit<svxtrk->getNhits(); ihit++){
      //float *info = svxtrk->getClusterInfo(ihit);
      //cout<<"ihit : "<<ihit<<" "<<*info<<endl;

      SvxClusterInfo *info = svxtrk->getClusterInfo(ihit);
      int idx=0;
      if     (info->getLayer()==3){ idx=15; }
      else if(info->getLayer()==2){ idx=19; }
      else if(info->getLayer()==1){ idx=23; }
      else                        { idx=27; }
     
      ary[idx+0] = info->getLayer();
      ary[idx+1] = info->getPosition(0);
      ary[idx+2] = info->getPosition(1);
      ary[idx+3] = info->getPosition(2);

      if(isN4){
        cout<<"   layer : x:y:z : "<<(int)info->getLayer()<<" ";
        cout<<info->getPosition(0)<<" ";
        cout<<info->getPosition(1)<<" ";
        cout<<info->getPosition(2)<<" "<<endl;
      }

   
    }
  }
}



//==============================================================

int svxcntclusdisp::End(PHCompositeNode *topNode) {
  cout << "svxcntclusdisp::End" << endl;
  return 0;
}

void svxcntclusdisp::initCanvas(){

  if(m_svxGeo==NULL){
    cout<<"initCanvas svxDetectorGeo is not initialized"<<endl;
    return ;
  }

  m_c1 = new TCanvas("c1","c1", 800, 800);

  m_frame = new TH2F("frame", "", 100, -20, 20, 100, -20, 20);


  int isen=0;
  for(int ily=0; ily<4; ily++){
    int nLadder = m_svxGeo->get_nBarLadder(ily);
    for(int ildr=0; ildr<nLadder; ildr++){
      SvxSensor *sens_p = m_svxGeo->GetSensorPtr(ily, ildr, isen);
      float xyz[3];
      float xwid = sens_p->get_xhalfWidth();
      float ywid = 0.1 ; // 500um
      float zwid = 0;

      for(int i=0; i<3; i++){ xyz[i] = sens_p->get_transVector(i); }

      float rot[3][3];
      for(int i=0; i<3; i++){ 
        for(int j=0; j<3; j++){ 
          rot[i][j] = sens_p->get_rotMatrix(j, i);
        }
      }

      float xyz_lorg[4][3]={{-xwid, -ywid, -zwid},
                            { xwid, -ywid, -zwid},
                            { xwid,  ywid, -zwid},
                            {-xwid,  ywid, -zwid}};

      //  ratation
      //         | (0,0), (0,1), (0,2) |   |cosphi, -sinphicostheta,  sinphisintheta|
      //  Mrot = | (1,0), (1,1), (1,2) | = |sinphi,  cosphicostheta, -cosphisintheta|
      //         | (2,0), (2,1), (2,2) |   |     0,        sintheta,        costheta|
      float xyz_lo[4][3];
      for(int ipos=0; ipos<4; ipos++){
        xyz_lo[ipos][0] =rot[0][0]*xyz_lorg[ipos][0] + rot[0][1]*xyz_lorg[ipos][1] + rot[0][2]*xyz_lorg[ipos][2];
        xyz_lo[ipos][1] =rot[1][0]*xyz_lorg[ipos][0] + rot[1][1]*xyz_lorg[ipos][1] + rot[1][2]*xyz_lorg[ipos][2];
        xyz_lo[ipos][2] =rot[2][0]*xyz_lorg[ipos][0] + rot[2][1]*xyz_lorg[ipos][1] + rot[2][2]*xyz_lorg[ipos][2];
//        cout<<ily<<" "<<ildr<<" ipos ="<<ipos<<", "<<xyz_lo[ipos][0]<<" "<<xyz_lo[ipos][1]<<" "<<xyz_lo[ipos][2]<<" --- ";
//        cout<<rot[0][0]<<" "<<xyz_lorg[ipos][0]<<" "<<rot[1][0]<<" "<<xyz_lorg[ipos][1]<<" "<<rot[2][0]<<" "<<xyz_lorg[ipos][2]<<endl;
      }

      float x_pol[5], y_pol[5], z_pol[5];
      for(int i=0; i<5; i++){
        if(i<4){
          x_pol[i] = xyz[0] + xyz_lo[i][0];
          y_pol[i] = xyz[1] + xyz_lo[i][1];
          z_pol[i] = xyz[2] + xyz_lo[i][2];
        } 
        else {
          x_pol[i] = xyz[0] + xyz_lo[0][0];
          y_pol[i] = xyz[1] + xyz_lo[0][1];
          z_pol[i] = xyz[2] + xyz_lo[0][2];
        }
      }
 
      m_pl[ily][ildr] = new TPolyLine(5, x_pol, y_pol);
      m_pl[ily][ildr]->SetFillStyle(1001);
      m_pl[ily][ildr]->SetFillColor(6);
	if(0==1) std::cout<<z_pol[0]<<std::endl;
/*
      TPolyLine *plse = new TPolyLine(2, &x_pol[2], &y_pol[2]);
      plse->SetLineColor(2);
      plse->SetLineWidth(1.5);
      plse->Draw("same");

      TPolyLine *plsex = new TPolyLine(2, &x_pol[1], &y_pol[1]);
      plsex->SetLineColor(4);
      plsex->SetLineWidth(3.0);
      plsex->Draw("same");
*/
      
      m_txt[ily][ildr] = new TText(xyz[0], xyz[1], TString::Format("%d", ildr));
      m_txt[ily][ildr]->SetTextSize(0.02);
      m_txt[ily][ildr]->Draw("same");
     
    }
  }

  m_isInitCanvas = true;
}

void svxcntclusdisp::drawCanvas(){
  if(!m_isInitCanvas){
    cout<<"Canvas is Not initialized"<<endl;
    return;
  }

  m_c1->cd();
  m_frame->Draw();
  for(int ily=0; ily<4; ily++){
    int nLadder = m_svxGeo->get_nBarLadder(ily);
    for(int ildr=0; ildr<nLadder; ildr++){
      m_pl[ily][ildr]->Draw("same");
      m_txt[ily][ildr]->Draw("same");
    }
  }
}

void svxcntclusdisp::drawEvent(){
  if(!m_isInitCanvas){
    cout<<"Canvas is Not initialized"<<endl;
    return;
  }

  m_c1->cd();

  vector<TMarker*>::iterator itr;
  for(itr=m_vClusPos.begin(); itr!=m_vClusPos.end(); ++itr){
    TMarker *tmp = *itr;
    tmp->Draw("same");
  }


  vector<TArc*>::iterator itr1;
  for(itr1=m_vTrkArc.begin(); itr1!=m_vTrkArc.end(); ++itr1){
    TArc *tmp = *itr1;
    tmp->Draw("same");
  }

  vector<TLine*>::iterator itr2;
  for(itr2=m_vTrkLine.begin(); itr2!=m_vTrkLine.end(); ++itr2){
    TLine *tmp = *itr2;
    tmp->Draw(">same");
  }

  vector<TLine*>::iterator itr3;
  for(itr3=m_vRngLine.begin(); itr3!=m_vRngLine.end(); ++itr3){
    TLine *tmp = *itr3;
    tmp->Draw(">same");
  }
}

void svxcntclusdisp::cctestGenerateCluster(){
  if(m_cctestList!=NULL){
    return;
  }

  if(m_svxGeo==NULL) m_svxGeo = new svxDetectorGeo();

  if(!m_svxGeo->Read_svxPISApar()){
    cout<<"Failed to open the svxPISA.par"<<endl;
    return;
  }

  /////////////////////
  m_cctestList = new SvxClusterListv4();

  // genrate data
  static const int NBIN = 10;

  // 10 cluster per sensor
  int nlayer = 4;
  for(int ilay=0; ilay<nlayer; ilay++){
    int nladder = m_svxGeo->get_nBarLadder(ilay);
    int nsensor = m_svxGeo->get_nBarSensor(ilay);
    //cout<<nladder<<" "<<nsensor<<endl;

    for(int ilad=0; ilad<nladder; ilad++){
      for(int isen=0; isen<nsensor; isen++){
        SvxSensor *sens = m_svxGeo->GetSensorPtr(ilay, ilad, isen);

        int nsec = sens->get_nSection();
        float xhwid_total = sens->get_xhalfWidth(0, 0);
        float zhwid_total = 0;
        for(int isec=0; isec<nsec; isec++){
          zhwid_total += sens->get_zhalfWidth(isec, 0);
        }
        //cout<<ilay<<" "<<ilad<<" "<<isen<<" "<<xhwid_total<<" "<<zhwid_total<<endl;

        float binwid = 2.*xhwid_total/NBIN;
        for(int i=0; i<NBIN+1; i++){
          double xpos = i*binwid - xhwid_total;
          double pos[3] = {xpos, 0.0, 0.0};
          double gpos[3];
          sens->position_local2global(pos, gpos);
          //cout<<ilay<<" "<<ilad<<" "<<isen<<" "<<i<<" "<<pos[0]<<" "<<gpos[0]<<endl;

          SvxCluster *cluster = m_cctestList->addCluster();
          cluster->set_svxSection(0);
          cluster->set_layer(ilay);
          cluster->set_ladder(ilad);
          cluster->set_sensor(isen);
        //--    cluster->set_sensorType(sensorType);
        //--    cluster->set_adc(0, 1);
        //--    cluster->set_adc(1, 1);
        //--    cluster->set_size(cluster_list[i].nhits);
        //--    cluster->set_xz_size(0, xsize);
        //--    cluster->set_xz_size(1, zsize);
        //--    cluster->set_circumference(circumference);
        //--    cluster->set_edgeflag(edge_flag);
          for ( int j=0; j<3; j++ )
            {
               cluster->set_xyz_local (j, pos[j]);
               cluster->set_xyz_global(j, gpos[j]);
            }
        }
      }
    }
  }

  m_cctestCont = new SvxClusterContainer();
  m_cctestCont->load_clusters(m_cctestList);
}

void svxcntclusdisp::cctestDrawPoint(float phi0, float dphi){
  if(m_cctestList==NULL){
    cout<<"No ClusterList object. generate first"<<endl;
    return ;
  }
  // clear current marker
  vector<TMarker*>::iterator itr;
  for(itr=m_vccPos.begin(); itr!=m_vccPos.end(); ++itr){
    TMarker *m = *itr;
    delete m;
  }
  m_vccPos.clear();

  // fill & draw
  int nclus = m_cctestList->get_nClusters();
  for(int i=0; i<nclus; i++){
    SvxCluster *clus = m_cctestList->get_Cluster(i);
    float x = clus->get_xyz_global(0);
    float y = clus->get_xyz_global(1);

    TMarker *m = new TMarker(x, y, 20);
    m->SetMarkerColor(2);
    m->SetMarkerSize(1.0);
    m->Draw("same");
    m_vccPos.push_back(m);
  }

  int nlayer=4;
  int nsub[4] = {1,1,3,3};
  int sub[4][3] = {{0,0,0},{1,0,0},{2,3,4},{5,6,7}};

  vector<SvxCluster*> vtmp;
  vector<int> vsublayer;
  for(int isub=0; isub<nlayer; isub++){

    vtmp.clear();
    int nclus = m_cctestCont->find_clusters_multi(vtmp, vsublayer, nsub[isub], sub[isub], phi0, dphi, 0, 5);
    //int nclus = m_cctestCont->find_clusters_layer(vtmp, vsublayer, isub, phi0, dphi, 0, 5);
    cout<<"isub : "<<isub<<" "<<nclus<<endl;
    for(int i=0; i<nclus; i++){
      SvxCluster *clus = vtmp[i];
      float x = clus->get_xyz_global(0);
      float y = clus->get_xyz_global(1);
      //if(i<10){
        cout<<isub<<" : "<<x<<" "<<y<<" "<<vsublayer[i]<<endl;
      //}

      TMarker *m = new TMarker(x, y, 20);
      m->SetMarkerColor(4+isub);
      m->SetMarkerSize(1.0);
      m->Draw("same");
      m_vccPos.push_back(m);
    }
  }

}

void svxcntclusdisp::cctestDrawPointBlock(float phi0, float dphi, int block){
  if(m_cctestList==NULL){
    cout<<"No ClusterList object. generate first"<<endl;
    return ;
  }
  // clear current marker
  vector<TMarker*>::iterator itr;
  for(itr=m_vccPos.begin(); itr!=m_vccPos.end(); ++itr){
    TMarker *m = *itr;
    delete m;
  }
  m_vccPos.clear();

  // fill & draw
  int nclus = m_cctestList->get_nClusters();
  for(int i=0; i<nclus; i++){
    SvxCluster *clus = m_cctestList->get_Cluster(i);
    float x = clus->get_xyz_global(0);
    float y = clus->get_xyz_global(1);

    TMarker *m = new TMarker(x, y, 20);
    m->SetMarkerColor(2);
    m->SetMarkerSize(1.0);
    m->Draw("same");
    m_vccPos.push_back(m);
  }

  vector<SvxCluster*> vtmp;
  vector<int> vsublayer;
  nclus = m_cctestCont->find_clusters_block(vtmp, vsublayer, block, phi0, dphi, 0, 5);

  for(int i=0; i<nclus; i++){
    SvxCluster *clus = vtmp[i];
    float x = clus->get_xyz_global(0);
    float y = clus->get_xyz_global(1);
    if(i<10){
      cout<<i<<" : "<<x<<" "<<y<<endl;
    }

    TMarker *m = new TMarker(x, y, 20);
    m->SetMarkerColor(4);
    m->SetMarkerSize(1.0);
    m->Draw("same");
    m_vccPos.push_back(m);
  }
  

}
