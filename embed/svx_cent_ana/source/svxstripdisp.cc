#include <iostream>
#include <iomanip>
#include <string>

#include "phool.h"
#include "PHTypedNodeIterator.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "Fun4AllReturnCodes.h"

#include "svxstripdisp.h"

#include "VtxOut.h"
#include "RunHeader.h"
#include "EventHeader.h"
#include "Bbc.hh"
#include "BbcOut.h"

#include "SvxRawhitList.h"
#include "SvxRawhit.h"
#include "SvxClusterList.h"
#include "SvxRawhitClusterList.h"
#include "SvxRawhitCluster.h"
#include "SvxCluster.h"
#include "SvxPixelHotDeadMap.h"

#include "svxAddress.hh"
#include "svxDetectorGeo.hh"
#include "SvxSensor.h"

#include "getClass.h"

#include "TCanvas.h"
#include "TH2.h"
#include "TGraph.h"

using namespace std;
using namespace findNode;
//==============================================================

svxstripdisp::svxstripdisp() {
  ThisName = "svxstripdisp";
  init_ana=0;
  EventNumber=0;

  m_c1        = NULL;
  m_c2        = NULL;
  for(int ilay=0; ilay<4; ilay++){
    for(int ilad=0; ilad<24; ilad++){
      for(int isens=0; isens<6; isens++){
        m_h2_sensor[ilay][ilad][isens] = NULL;
        m_gcls_sensor[ilay][ilad][isens] = NULL; // layer, ladder, sensor
        m_gcls_sensor_good[ilay][ilad][isens] = NULL; // layer, ladder, sensor
      }
    }
  }

  m_clssizecut[0] = 0;
  m_clssizecut[1] = 0;
  m_clssizecut[2] = 0;
  m_clssizecut[3] = 0;

  m_svxrawhit     = NULL;
  m_svx           = NULL;
  m_svxrawhitclus = NULL; 

  m_svxAdrObj = NULL;
  m_svxGeo    = NULL;
  m_pixelhotdead = NULL;
}

//==============================================================

svxstripdisp::~svxstripdisp() {
  if(m_c1!=NULL)        delete m_c1;
  if(m_c2!=NULL)        delete m_c2;
  for(int ilay=0; ilay<4; ilay++){
    for(int ilad=0; ilad<24; ilad++){
      for(int isens=0; isens<6; isens++){
        if(m_h2_sensor[ilay][ilad][isens]!=NULL)   delete m_h2_sensor[ilay][ilad][isens];
        if(m_gcls_sensor[ilay][ilad][isens]!=NULL) delete m_gcls_sensor[ilay][ilad][isens];
        if(m_gcls_sensor_good[ilay][ilad][isens]!=NULL) delete m_gcls_sensor_good[ilay][ilad][isens];
      }
    }
  }

}

//==============================================================
TH2F *inithist_pixelmap(){
  float zbin[131];
  int idx=1;

  int   zsecid[9] = {3,0,2,1,2,1,2,0,3};
  float zwidth[4] = {0.0425, 0.0425, 0.0625, 0.2};
  int   nzsec[4]  = {31, 30, 2, 1};
  float zoffset = -2.78-0.2;

  zbin[0] = zoffset;

  for(int zid=0; zid<9; zid++){
    for(int i=0; i<nzsec[zsecid[zid]]; i++){
      zoffset += zwidth[zsecid[zid]];
      zbin[idx] = zoffset;
//      cout<<zid<<" "<<i<<" "<<idx<<" "<<zbin[idx]<<" "<<zwidth[zsecid[zid]]<<endl;
      idx++;
    }
  }

  idx=2;
  float xbin[259];
  float xoffset = -0.640;
  xbin[0] = xoffset-0.2;
  xbin[1] = xoffset;
  for(int i=0; i<256; i++){
    xoffset += 0.005;
    xbin[idx] = xoffset;
//    cout<<i<<" "<<idx<<" "<<xbin[idx]<<endl;
    idx++;
  }
  xbin[258] = xoffset+0.2;

  return new TH2F("h", "h", 130, zbin, 258, xbin);
}


int svxstripdisp::Init(PHCompositeNode *topNode) {

  cout << "svxstripdisp::Init started..." << endl;

  m_c1 = new TCanvas("c1", "c1", 800, 900);
  m_c1->Divide(2,3);

  
  for(int ilay=0; ilay<4; ilay++){
    for(int ilad=0; ilad<24; ilad++){
      for(int isens=0; isens<6; isens++){
        TString s_name = TString::Format("h2_sensor_%d_%d_%d", ilay, ilad, isens);
        TString s_tit  = TString::Format("Layer:%d, Ladder:%d, Sensor:%d", ilay, ilad, isens);
        if(ilay<2){
          m_h2_sensor[ilay][ilad][isens] = inithist_pixelmap();
          m_h2_sensor[ilay][ilad][isens]->SetName(s_name);
          m_h2_sensor[ilay][ilad][isens]->SetTitle(s_tit);
        } else {
          m_h2_sensor[ilay][ilad][isens] = new TH2F(s_name, s_tit, 60, -3, 3, 768, -1.5360, 1.5360);
        }
        m_h2_sensor[ilay][ilad][isens]->SetMaximum(3);

        m_gcls_sensor[ilay][ilad][isens] = new TGraph();
        m_gcls_sensor[ilay][ilad][isens]->SetMarkerStyle(20);
        m_gcls_sensor[ilay][ilad][isens]->SetMarkerSize(0.5);

        m_gcls_sensor_good[ilay][ilad][isens] = new TGraph();
        m_gcls_sensor_good[ilay][ilad][isens]->SetMarkerStyle(20);
        m_gcls_sensor_good[ilay][ilad][isens]->SetMarkerColor(2);
        m_gcls_sensor_good[ilay][ilad][isens]->SetMarkerSize(0.5);
      }
    }
  }

  cout << "svxstripdisp::Init ended." << endl;
  return 0;
}

//==============================================================
  
int svxstripdisp::InitRun(PHCompositeNode *topNode) {
  cout << "svxstripdisp::InitRun started..." << endl;
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

  SvxPixelHotDeadMap* pixhotdead = findNode::getClass<SvxPixelHotDeadMap>(topNode, "SvxPixelHotDeadMap");
  if ( pixhotdead == NULL) {
    if(verbosity>0) { cout << PHWHERE<< "Can't find SvxPixelHotDeadMap. " << endl; }
    return ABORTRUN;
  }
  m_pixelhotdead = pixhotdead;

  cout << "svxstripdisp::InitRun ended." << endl;
  return 0;
}

//==============================================================


int svxstripdisp::process_event(PHCompositeNode *topNode) {

  //float ntp1[99],ntp2[99],ntp3[99],ntp4[99],ntp5[99],ntp6[99],ntp7[99];

  //  cout<<"aa"<<endl;
  //RunHeader *run = getClass<RunHeader>(topNode,"RunHeader");
  EventHeader *event_header = getClass<EventHeader>(topNode,"EventHeader");
  //BbcOut *bbc = getClass<BbcOut>(topNode,"BbcOut");

  //int RunNumber = (run!=NULL) ? run->get_RunNumber() : -9999;
  int EventNum = (event_header!=NULL) ? event_header->get_EvtSequence() : -9999;
  //cout << "run = "<<RunNumber<<" event = "<<EventNum<<endl;


  m_svxrawhit     = getClass<SvxRawhitList>(topNode,"SvxRawhitList");
  m_svx           = getClass<SvxClusterList>(topNode,"SvxClusterList");
  m_svxrawhitclus = getClass<SvxRawhitClusterList>(topNode,"SvxRawhitClusterList");

  if( m_svxrawhit    ==NULL ) {cout<<"No SvxRawhitList"<<endl; return EVENT_OK; }
  if( m_svx          ==NULL ) {cout<<"No SvxClusterList"<<endl; return EVENT_OK; }
  if( m_svxrawhitclus==NULL ) {cout<<"No SvxRawhitClusterList"<<endl; return EVENT_OK; }

  int nsvx       = m_svx->get_nClusters();
  int nsvxrawhit = m_svxrawhit->get_nRawhits();

  if(init_ana==0) {
    cout<<"init ana"<<endl;
    topNode->print(); 
    init_ana=1;
 
    m_svx->identify();          
    m_svxrawhitclus->identify();
    m_svxrawhit->identify();    

    cout<<endl<<endl;
    cout<<"init ana ends"<<endl;
  }

  if(EventNumber%1==0) {
    cout << "------- Event # " << EventNumber << " (EvtSeq= "<<EventNum<<" )  ";
    cout<< "nsvxrawhit= " << nsvxrawhit << " nsvx= " << nsvx << endl;
  }

  EventNumber++;

  if( fill() ){
    // found . terminate
    cout<<"large size cluster found. terminate"<<endl;
    return ABORTRUN;
  }

  return 0;
}

bool svxstripdisp::fill(){
  if(m_svxrawhit==NULL) {
    cout<<"No rawhit object"<<endl;
    return false;
  }
 

  for(int ilay=0; ilay<4; ilay++){
    for(int ilad=0; ilad<24; ilad++){
      for(int isens=0; isens<6; isens++){
        m_h2_sensor[ilay][ilad][isens]->Reset();
        m_gcls_sensor[ilay][ilad][isens]->Set(0);
        m_gcls_sensor_good[ilay][ilad][isens]->Set(0);
      }
    }
  }

  int nsvxrawhit = m_svxrawhit->get_nRawhits();
  //cout<<"nrawhit : "<<nsvxrawhit<<endl;
  for(int ihit=0; ihit<nsvxrawhit; ihit++) {
    SvxRawhit* tmp = m_svxrawhit->get_Rawhit(ihit);

    //int raw_hitID         = tmp->get_hitID();
    //int raw_section       = tmp->get_svxSection();
    int raw_layer         = tmp->get_layer();
    int raw_ladder        = tmp->get_ladder();
    int raw_sensor        = tmp->get_sensor();
    int raw_SS      = tmp->get_sensorSection();
    int raw_readout = tmp->get_sensorReadout();
    int raw_channel = tmp->get_channel();
    int hotdead     = tmp->get_HotDeadFlag();

    // check 
    
    if(raw_layer<2){
      //int raw_pixelModule = tmp->get_pixelModule();
      int raw_pixelROC    = tmp->get_pixelROC();
      //int  raw_ix  = m_svxAdrObj->getPixelRocIX0(raw_channel);
      //int  raw_iz  = m_svxAdrObj->getPixelRocIZ0(raw_channel);
      int  raw_ixs = m_svxAdrObj->getPixelSensorIX0(raw_pixelROC, raw_channel);
      int  raw_izs = m_svxAdrObj->getPixelSensorIZ0(raw_pixelROC, raw_channel);
      float lx = get_sensorXpos_pixel(raw_layer, raw_ladder, raw_sensor, raw_SS, raw_readout, raw_ixs);
      float lz = get_sensorZpos_pixel(raw_layer, raw_ladder, raw_sensor, raw_SS, raw_readout, raw_izs);
      //if(raw_layer==0&&raw_ladder==0&&raw_sensor==0){
      //  cout<<raw_SS<<" "<<raw_readout<<" "<<raw_channel<<" "<<lx<<" "<<lz<<endl;
      //}
      //if(raw_pixelModule==0&&raw_pixelROC==0){
      //  cout<<ihit<<" "<<raw_layer<<" "<<raw_ladder<<" "<<raw_sensor<<" ";
      //  cout<<raw_channel<<" "<<raw_ixs<<" "<<raw_izs<<" "<<lx<<" "<<lz<<endl;
      //}

      //int status = m_pixelhotdead->getStatus(raw_pixelModule, raw_pixelROC, raw_iz, raw_ix);
      //int status = m_pixelhotdead->getPixelStatus(raw_pixelModule, raw_pixelROC, raw_iz, raw_ix);
      //if(status!=0 || hotdead!=0){
      //  int mapflag = (status==0) ? 1 : 0;
      //  cout<<"HotChannel!! : "<<raw_layer<<" "<<raw_ladder<<" "<<raw_sensor<<" ";
      //  cout<<raw_pixelModule<<" "<<raw_pixelROC<<" "<<raw_iz<<" "<<raw_ix<<" : "<<status<<" "<<hotdead<<endl;
      //}
     

      if(hotdead==0){
      //if(status==0){}
        m_h2_sensor[raw_layer][raw_ladder][raw_sensor]->Fill(lz, lx, 2);
      } else {
        m_h2_sensor[raw_layer][raw_ladder][raw_sensor]->Fill(lz, lx, 3);
      }
    }
    else
    {
      int raw_adc     = tmp->get_adc();
    
      if(raw_adc<30) continue;

      float raw_pos = raw_channel*0.008 - 1.536 + 0.004; //get_sensorXpos_pixel(raw_layer, raw_ladder, raw_sensor, raw_SS, raw_readout, raw_channel);

      //if(raw_layer==2&&raw_ladder==0&&raw_sensor==0) cout<<raw_channel<<" "<<raw_pos<<endl;

      float zoffset = (raw_SS==0) ? -3.0 : 0.0; // z-strip
      for(int i=0; i<30; i++){
        float xstrip  = (i+0.5)*0.1 + zoffset; // x-strip
        //cout<<xstrip<<" "<<raw_lx<<endl;

        if(raw_readout==0) {
          m_h2_sensor[raw_layer][raw_ladder][raw_sensor]->Fill(xstrip,  raw_pos-0.002, raw_adc); // x-strip
        }
        else {
          float xoffset = (raw_SS==0) ? 0.0 : -29*0.008;
          float upos = raw_pos + i*0.008 + 0.002 + xoffset;
          if(upos>=1.536)  upos -= 3.072;
          if(upos<=-1.536) upos += 3.072;

          //if(raw_layer==2&&raw_ladder==0&&raw_sensor==0) cout<<xstrip<<" "<<upos<<" "<<raw_channel<<endl;
          m_h2_sensor[raw_layer][raw_ladder][raw_sensor]->Fill(xstrip,  upos, raw_adc);         // u-strip
        }
      }
    }
  }

  // Graph
  bool islarge=false;
  int ncls = m_svx->get_nClusters();
  for(int icls=0; icls<ncls; icls++){
    SvxCluster* cluster = m_svx->get_Cluster(icls);

    int cls_layer  = cluster->get_layer();
    int cls_ladder = cluster->get_ladder();
    int cls_sensor = cluster->get_sensor();
    int cls_size   = cluster->get_size();
    int cls_xsize  = cluster->get_xz_size(0);
    int cls_zsize  = cluster->get_xz_size(1);

    //if(cls_layer>=2){
      float cls_x = cluster->get_xyz_local(0);
      float cls_z = cluster->get_xyz_local(2);

      //cout<<"cls : "<<cls_layer<<" "<<cls_ladder<<" "<<cls_sensor<<endl;

      int n = m_gcls_sensor[cls_layer][cls_ladder][cls_sensor]->GetN();
      m_gcls_sensor[cls_layer][cls_ladder][cls_sensor]->SetPoint(n, cls_z, cls_x);
    //}

    bool isOK=false;
    if( 
        (cls_layer==0 && m_clssizecut[0]>0 && cls_size>m_clssizecut[0]) ||
        (cls_layer==1 && m_clssizecut[1]>0 && cls_size>m_clssizecut[1]) ||
        (cls_layer==2 && m_clssizecut[2]>0 && cls_size>m_clssizecut[2]) ||
        (cls_layer==3 && m_clssizecut[3]>0 && cls_size>m_clssizecut[3]) 
      )
    {
      islarge = true;
      isOK = true;
    }

    if((verbosity>0&&cls_size>40)|| isOK) { 
      cout<<"layer : ladder : sensor : size = ";
      cout<<cls_layer<<" : ";
      cout<<cls_ladder<<" : ";
      cout<<cls_sensor<<" : ";
      cout<<cls_size<<" : ";
      cout<<cls_xsize<<" : ";
      cout<<cls_zsize<<" : ";
      cout<<endl; 
    }
  }

//  cout<<endl<<endl;;
/*
  bool isHot = false;
  int prev_clsid=-1;
  int nclsraw = m_svxrawhitclus->get_nRawhitClusters();
  for(int icls=0; icls<nclsraw; icls++){
    SvxRawhitCluster *rcls = m_svxrawhitclus->get_RawhitCluster(icls);
    int rawid = rcls->get_rawhitID();
    int clsid = rcls->get_clusterID();

    ////////////////
    {
      SvxRawhit*  raw     = m_svxrawhit->get_Rawhit(rawid);
      SvxCluster* cluster = m_svx->get_Cluster(clsid);

      int raw_layer       = raw->get_layer();
      int raw_pixelmodule = raw->get_pixelModule();
      int raw_pixelROC    = raw->get_pixelROC();
      int c_size   = cluster->get_size();

      bool isSelected = (raw_layer<2&&raw_pixelmodule==0&&raw_pixelROC==0);

      if(clsid!=prev_clsid){
        if(isSelected) {
          cout<<endl<<"cluster "<<clsid<<"("<<c_size<<") : "<<endl; //("<<nrhit<<") : "<<c_layer<<" "<<raw_pixelmodule<<" "<<raw_pixelROC<<endl;
        }

        prev_clsid = clsid;
      }

      if(isSelected){
        cout<<rawid<<" "<<flush;
      }
    }
    ////////////////


    // reset routine
    if(clsid!=prev_clsid){
      //cout<< ((isHot) ?  "Hot": "N/A")<<" "<<endl;

      if(prev_clsid!=-1&&!isHot){
        SvxCluster* cluster = m_svx->get_Cluster(prev_clsid);

        int cls_layer  = cluster->get_layer();
        int cls_ladder = cluster->get_ladder();
        int cls_sensor = cluster->get_sensor();

        float cls_x = cluster->get_xyz_local(0);
        float cls_z = cluster->get_xyz_local(2);

        int n = m_gcls_sensor_good[cls_layer][cls_ladder][cls_sensor]->GetN();
        m_gcls_sensor_good[cls_layer][cls_ladder][cls_sensor]->SetPoint(n, cls_z, cls_x);
      }

      // reset
      prev_clsid = clsid;
      isHot=false;
    }

    // this cluster
    SvxRawhit *rawhit = m_svxrawhit->get_Rawhit(rawid);
    bool isHotChan=false;
    if(rawhit!=NULL){
      isHotChan = (rawhit->get_HotDeadFlag()!=0);
      isHot |= isHotChan;
    }

    //cout<<"rawclus : "<<icls<<" "<<clsid<<" "<<rawid<<" : ";
    //cout<< ((isHotChan) ?  "Hot": "N/A")<<" : ";

    //cout<<endl;
  }
*/

  return islarge;
}

void svxstripdisp::drawone(int layer, int ladder, int sensor, bool isclus){
  if(m_c2==NULL) {
    m_c2 = new TCanvas("c2", "c2", 500, 500);
  }
  if( 0<=layer&&layer<4 &&
      0<=ladder&&ladder<24&&
      0<=sensor&&sensor<6)
  {
      m_c2->cd();

      m_h2_sensor[layer][ladder][sensor]->Draw("colz");

      if(isclus) {
        m_gcls_sensor[layer][ladder][sensor]->Draw("p");
        m_gcls_sensor_good[layer][ladder][sensor]->Draw("p");
      }
  }
  else {
    cout<<"out of range "<<layer<<" "<<ladder<<endl;
  }
}

void svxstripdisp::draw(int layer, int ladder, bool isclus){
  if(m_c1==NULL) {
    cout<<"No canvas is initialized "<<endl;
    return;
  }

  if( 0<=layer&&layer<4 &&
      0<=ladder&&ladder<24) 
      //0<=sensor&&sensor<6)
  {
    for(int i=0; i<6; i++){
      m_c1->cd(i+1);
      m_h2_sensor[layer][ladder][i]->Draw("colz");

      if(isclus) {
        m_gcls_sensor[layer][ladder][i]->Draw("p");
        m_gcls_sensor_good[layer][ladder][i]->Draw("p");
      }
      //cout<<"i : "<<i<<endl;
      //m_gcls_sensor[layer][ladder][i]->Print();
    }
  }
  else {
    cout<<"out of range "<<layer<<" "<<ladder<<endl;
  }
}

void svxstripdisp::setSizeCut(const int layer, const int val)
{
  if(layer<0||3<layer){
    cerr<<__FUNCTION__<<" layer out of range "<< layer <<endl;
    return;
  }

  m_clssizecut[layer] = val;
}




//==============================================================

int svxstripdisp::End(PHCompositeNode *topNode) {
  cout << "svxstripdisp::End:  Writing out..." << endl;
  cout << "svxstripdisp::End:  Closing output file..." << endl;
  return 0;
}


// copy from SvxPixel1v1
// Return the local X position from the sensor Ix value
double
svxstripdisp::get_sensorXpos_pixel(const int layer, const int ladder, const int sens, 
                       const int section, const int readout, int ix) const
{
  SvxSensor *sensor = m_svxGeo->GetSensorPtr(layer, ladder, sens);

  double xPitch = sensor->get_xPitch(section, readout);
  double xhalfWidth = sensor->get_xhalfWidth(section, readout);
  double xpos = (ix+0.5)*xPitch - xhalfWidth;
  xpos += sensor->get_secXpos(section, readout);
  if ( sensor->get_xcntRvrs(section, readout) ) xpos *= -1.0;
  return xpos;
}

// copy from SvxPixel1v1
// Return the local Z position from the sensor Iz value
double
svxstripdisp::get_sensorZpos_pixel(const int layer, const int ladder, const int sens, 
                       const int section, const int readout, int iz) const
{
  SvxSensor *sensor = m_svxGeo->GetSensorPtr(layer, ladder, sens);

  // First adjust the Iz to its value within the section
  for (int isec=0; isec<section; isec++)
    {
      int npitch = sensor->get_nZpitch(isec, readout);
      iz -= npitch;
    }
  double pitch = sensor->get_zPitch(section, readout);
  double halfWidth = sensor->get_zhalfWidth(section, readout);
  double zpos = (iz+0.5)*pitch - halfWidth;
  zpos += sensor->get_secZpos(section, readout);
  if ( sensor->get_zcntRvrs(section, readout) ) zpos *= -1.0;
  return zpos;
}

