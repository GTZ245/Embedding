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

#include "svxana.h"

#include "PHCentralTrack.h"
#include "PHTrackOut.h"
#include "PHGlobal.h"
#include "emcClusterContainer.h"
#include "emcClusterContent.h"
#include "TriggerHelper.h"
#include "VtxOut.h"
#include "ErtOut.h"
#include "PadCluster.h"
#include "CglTrack.h"
#include "DchTrack.h"
#include "CrkRing.h"
#include "RunHeader.h"
#include "EventHeader.h"
#include "TrigRunLvl1.h"
#include "Bbc.hh"
#include "BbcOut.h"
//#include "KalFitOut.h"

#include "SvxRawhitList.h"
#include "SvxRawhit.h"
#include "SvxClusterList.h"
#include "SvxSegmentList.h"
#include "SvxRawhitClusterList.h"
#include "SvxRawhitCluster.h"
#include "SvxClusterList.h"
#include "SvxCluster.h"
#include "SvxSegmentv1.h"
#include "SvxTracker.h"

#include "svxAddress.hh"

#include "TFile.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom2.h"


#include "getClass.h"

using namespace std;
using namespace findNode;
//==============================================================

svxana::svxana() {
  ThisName = "svxana";
  init_ana=0;
  EventNumber=0;
  OutputFileName="svxanalysis.root";

  m_svxAdrObj = NULL;
}

//==============================================================

svxana::svxana(string filename) {
  ThisName = "svxana";
  init_ana=0;
  EventNumber=0;
  OutputFileName=filename;

//   m_svxAdrObj = NULL;
//   // hotchannel 
//   for(int imod=0; imod<NMODULE; imod++){
//     for(int ichip=0; ichip<NCHIP; ichip++){
//       for(int icol=0; icol<NCOLUMN; icol++){
//         for(int irow=0; irow<NROW; irow++){
//           m_hotchannel[imod][ichip][icol][irow]=false;
//         }
//       }
//     }
//   }
  //cluscont = new SvxClusterContainer();

}

//==============================================================

svxana::~svxana() {
//  if(m_svxAdrObj!=NULL) delete m_svxAdrObj;
}

//==============================================================

int svxana::Init(PHCompositeNode *topNode) {

  cout << "svxana::Init started..." << endl;
  OutputNtupleFile = new TFile(OutputFileName.c_str(),"RECREATE");
  cout << "svxana::Init: output file " << OutputFileName << " opened." << endl;

  random = new TRandom2();


  ntpsvxrawhit = new  TNtuple("ntpsvxrawhit","","event:hitID:section:layer:ladder:sensor:SS:readout:type:channel:adc:pixModule:pixROC");
  ntpsvxcluster = new  TNtuple("ntpsvxcluster","","event:x:y:z:adc1:adc2:layer:hitID:section:ladder:sensor:type:edgeflag:lx:ly:lz:phi");
  ntpsvxsegment = new  TNtuple("ntpsvxsegment","","event:mom:px:py:pz:L0phi:L1phi:L2phi:L3phi:charge:nhit0:nhit1:nhit2:nhit3:quality:dedx1:dedx2:dca:dca2d:nseg:bbcz:vtxz:vtxx:vtxy:phi0:phi0b:vtxxc:vtxyc:dca0:bend:L0ladder:L1ladder:L0sensor:L1sensor:dproj2:bend2:dproj3:bend3:dproj32:bend32:L2ladder:L3ladder:x0:y0:z0:x1:y1:z1:x2:y2:z2:x3:y3:z3");
  ntprawclus = new TNtuple("ntprawclus","","event:rlayer:rladder:rsensor:rsection:rreadout:rchannel:pixROC:pixMod:clayer:cladder:csensor:lx:lz:x:y:z:ix:iz:ixs:izs");
  ntpevt = new TNtuple("ntpevt","","event:n0:n1:n2:n3:zvtx:nseg:ntrk:bbcq");
  ntp_trk = new TNtuple("ntp_trk","","mom_cnt:phi0_cnt:mom:phi0:phi00:the0_cnt:the0:dphi0:n0:emcdphi:emcdz:ecore:dca0");
  ntp_cnt = new TNtuple("ntp_cnt","","mom:emcdphi:emcdz:the0:phi0:zed:zvtx");
  ntp_cnt_clus = new TNtuple("ntp_cnt_clus","","mom:emcdphi:emcdz:the0:phi0:zed:layer:zproj:dproj:bend:zv:phiv:zvtx");
  ntp_cnt_clus1 = new TNtuple("ntp_cnt_clus1","","mom:emcdphi:emcdz:the0:phi0:zed:zproj3:dproj3:zv3:zvtx:layer1:dproj1:bend1:zproj1:zv1");

  hdphi = new TH1F("hdphi","dphi(svx-central)",400,-4,4);
  hdphi2 = new TH1F("hdphi2","dphi(svx-central)",400,-4,4);
  hdphi3 = new TH1F("hdphi3","dphi(svx-central)",400,-4,4);
  hdphi3a = new TH1F("hdphi3a","dphi(svx-central)",400,-4,4);
  hdphi3b = new TH1F("hdphi3b","dphi(svx-central)",400,-4,4);
  hdthe = new TH1F("hdthe","dtheta(svx-central)",100,-0.3,0.3);
  hdthe2 = new TH1F("hdthe2","dtheta(svx-central)",100,-0.3,0.3);

  hsvxsegmentmult = new TH1F("hsvxsegmentmult","",1000,0.,1000.);
  hsvxclustermult = new TH1F("hsvxclustermult","",1000,0.,1000.);
  hsvxrawhitmult = new TH1F("hsvxrawhitmult","",1000,0.,1000.);
  hadc2X = new TH1F("hadc2X","",200,0.,200);

  hzvtx = new TH1F("hzvtx","bbc vertex",100,-50,50);

  hpix_chip_mult = new TH1F("hpix_chip_mult","pix chip multiplicity",1000,0.,1000.);

  h2pixel  = new TH2F("h2pixel","",18,-1.,17.0,32,-1.0,31.0);
  h2pixel2 = new TH2F("h2pixel2","nhit<30",18,-1.,17.0,32,-1.0,31.0);
  h2pixel3 = new TH2F("h2pixel3","nchip<5",18,-1.,17.0,32,-1.0,31.0);

  m_svxAdrObj = new svxAddress();
  m_svxAdrObj->set_usedatabase(0);
  m_svxAdrObj->Initialize();

  cout << "svxana::Init ended." << endl;
  return 0;
}

//==============================================================
  
int svxana::InitRun(PHCompositeNode *topNode) {
  cout << "svxana::InitRun started..." << endl;
  cout << "svxana::InitRun ended." << endl;
  return 0;
}

//==============================================================


bool bad_sensor(int layer, int ladder, int sensor) {
  return false;
}


bool bad_chip(int chip, int pix_ladder) {
  return false;
}


int svxana::process_event(PHCompositeNode *topNode) {

  float ntp1[99],ntp2[99],ntp3[99],ntp4[99],ntp5[99],ntp6[99],ntp7[99];

  //   cout<<"00a"<<endl;

  //  cout<<"aa"<<endl;
  //  PHTrackOut         *proj = getClass<PHTrackOut>(topNode,"PHTrackOut");
  //  CrkRing             *crk = getClass<CrkRing>(topNode,"CrkRing");
  //  PadCluster          *pc3 = getClass<PadCluster>(topNode,"Pc3Cluster");
  //  emcClusterContainer *emc = getClass<emcClusterContainer>(topNode,"emcClusterContainer");
  //  DchTrack            *dch = getClass<DchTrack>(topNode,"DchTrack");
  //  CglTrack            *cgl = getClass<CglTrack>(topNode,"CglTrack");
  //  PHGlobal *global = getClass<PHGlobal>(topNode,"PHGlobal");
  RunHeader *run = getClass<RunHeader>(topNode,"RunHeader");
  EventHeader *event_header = getClass<EventHeader>(topNode,"EventHeader");
  TrigRunLvl1 *lvl1 = getClass<TrigRunLvl1>(topNode,"TrigRunLvl1");
  BbcOut *bbc = getClass<BbcOut>(topNode,"BbcOut");

  int RunNumber = run->get_RunNumber();
  int EventNum = event_header->get_EvtSequence();
  cout << "run = "<<RunNumber<<" event = "<<EventNum<<endl;


  PHCentralTrack      *trk = getClass<PHCentralTrack>(topNode,"PHCentralTrack");
  int ntrk = (trk!=NULL) ? trk->get_npart() : 0;

  SvxRawhitList *svxrawhit = getClass<SvxRawhitList>(topNode,"SvxRawhitList");
  SvxClusterList
      *svx = getClass<SvxClusterList>(topNode,"SvxClusterList");
  SvxSegmentList *svxtracks= getClass<SvxSegmentList>(topNode,"SvxSegmentList");
  SvxRawhitClusterList 
            *svxrawhitclus = getClass<SvxRawhitClusterList>(topNode,"SvxRawhitClusterList");

  int nsvxtracks = svxtracks->get_nSegments();
  int nsvx = svx->get_nClusters();
  int nsvxrawhit = svxrawhit->get_nRawhits();
  int nsvxrawhitclus = 0;
  int nkalfit = 0;
  float zvtx = bbc->get_VertexPoint();
  float bbcq = bbc->get_ChargeSum(Bbc::North) + bbc->get_ChargeSum(Bbc::South);
  if(bbcq>200 || bbcq<10)return 1; 
  

  if(init_ana==0) {
    cout<<"init ana"<<endl;
    topNode->print(); 
    init_ana=1;
 
    if(svx!=NULL)          svx->identify();          else cout<<"no SvxCluster object"<<endl; 
    if(svxtracks!=NULL)    svxtracks->identify();    else cout<<"no SvxSegment object"<<endl;
    if(svxrawhitclus!=NULL)svxrawhitclus->identify();else cout<<"no SvxRawhitClus object"<<endl;
    if(svxrawhit!=NULL)    svxrawhit->identify();    else cout<<"no SvxRawhit object"<<endl;

    cout<<endl<<endl;
    cout<<"init ana ends"<<endl;
  }

  if(EventNumber%1==0) {
    cout << "------- Event # " << EventNumber << " nsvxrawhit= " << nsvxrawhit << " nsvx= " << nsvx << " ntrack= " << ntrk << " nkalfit= " << nkalfit << " nsvxtracks=" << nsvxtracks << endl;
  }


//------------------------ VTX information --------------------------
/*
  hsvxrawhitmult->Fill(float(nsvxrawhit));

  
  int L2X[16][5][2][384];
  int L2U[16][5][2][384];
  int L3X[24][6][2][384];
  int L3U[24][6][2][384];

  fill_n(&L2X[0][0][0][0],16*5*2*384,0);
  fill_n(&L2U[0][0][0][0],16*5*2*384,0);
  fill_n(&L3X[0][0][0][0],24*6*2*384,0);
  fill_n(&L3U[0][0][0][0],24*6*2*384,0);
  
  for(int i=0; i<nsvxrawhit; i++) {
    SvxRawhit* tmp = svxrawhit->get_Rawhit(i);

    int hitID         = tmp->get_hitID();
    int section       = tmp->get_svxSection();
    int layer         = tmp->get_layer();
    int ladder        = tmp->get_ladder();
    int sensor        = tmp->get_sensor();
    int sensorSection = tmp->get_sensorSection();
    int readout       = tmp->get_sensorReadout();
    int type          = tmp->get_sensorType();
    int channel       = tmp->get_channel();
    int adc           = tmp->get_adc();
    int pixelModule   = tmp->get_pixelModule();
    int pixelROC      = tmp->get_pixelROC();

    int strip_ladder = -1;

    ntp4[0]=(float)EventNumber;
    //if(readout==0) { ntp4[1]=adc; ntp4[2]=0.; }
    //if(readout==1) { ntp4[1]=0.; ntp4[2]=adc; }
    ntp4[1]=(float)hitID;
    ntp4[2]=(float)section;
    ntp4[3]=(float)layer;
    ntp4[4]=(float)ladder;
    ntp4[5]=(float)sensor;
    ntp4[6]=(float)sensorSection;
    ntp4[7]=(float)readout;
    ntp4[8]=(float)type;
    ntp4[9]=(float)channel;
    ntp4[10]=(float)adc;
    ntp4[11]=(float)pixelModule;
    ntp4[12]=(float)pixelROC;
    
    if(EventNumber<100) ntpsvxrawhit->Fill(ntp4);

    if(ladder == 2) {
      if(readout==0) L2X[ladder][sensor][section][channel] = adc;
      if(readout==1) L2U[ladder][sensor][section][channel] = adc;
    }
    if(ladder == 3) {
      if(readout==0) L3X[ladder][sensor][section][channel] = adc;
      if(readout==1) L3U[ladder][sensor][section][channel] = adc;
    }

  }//for(nrawhit)

  if(fabs(zvtx)<10) {
    for(int iladder=0;iladder<16;iladder++)
      for(int isensor=0;isensor<5;isensor++)
	for(int isection=0;isection<2;isection++) {
	  int ichannel=1;
	  while(ichannel<384) {
	    int adc00 = L2X[iladder][isensor][isection][ichannel-1];
	    int adc0 = L2X[iladder][isensor][isection][ichannel];
	    if(adc0>15 && adc00<10) {
	      ichannel++;
	      int adc1 = L2X[iladder][isensor][isection][ichannel];
	      if(adc1 > 15) {
		ichannel++;
		if(L2X[iladder][isensor][isection][ichannel]<10) {
		  hadc2X->Fill(adc1+adc0);
		}
	      }
	    }
	    ichannel++;
	  }
	}
  }//if(zvtx)



  hsvxclustermult->Fill(float(nsvx));

  //chick if this is a noisy event or not
  int nhit_layer[4];
  for(int i=0;i<4;i++) nhit_layer[i]=0;

  int pix_chip_mult[16][30];
  for(int i=0;i<16;i++)
    for(int j=0;j<30;j++)
      pix_chip_mult[i][j]=0;

  for(int i=0;i<nsvx;i++) {
     SvxCluster* clus0 = svx->get_Cluster(i);
     int layer  = clus0->get_layer();
     int ladder = clus0->get_ladder();
     int sensor = clus0->get_sensor();
     float lz = clus0->get_xyz_local(2);

     if(layer<2) {
       int pix_ladder;
       int chip = sensor*4 + (lz/1.3)+2;
       if(!bad_chip(chip,pix_ladder)) {
	 nhit_layer[layer]++;
       }
     } else if(!bad_sensor(layer,ladder,sensor)) {
       nhit_layer[layer]++;
     }

     if(0<=layer&&layer<2) {
       int pix_ladder;
       int chip = sensor*4 + (lz/1.3)+2;

       if(layer==0) pix_ladder = ladder;
       else pix_ladder = ladder + 10;

       if(0<= chip && chip<16 &&
	  0<= pix_ladder && pix_ladder < 30) {
	 pix_chip_mult[chip][pix_ladder]++;
       }
     }
  }

  for(int i=0;i<16;i++)
    for(int j=0;j<30;j++)
      hpix_chip_mult->Fill(pix_chip_mult[i][j]);

  for(int i=0; i<nsvx; i++) {
     SvxCluster* clus0 = svx->get_Cluster(i);
     float svxx = clus0->get_xyz_global(0);
     float svxy = clus0->get_xyz_global(1);
     float svxz = clus0->get_xyz_global(2);
     int   layer = clus0->get_layer();

     int chip_mult = -1;

     if(layer < 2) {//pixel
       int sensor = clus0->get_sensor();
       int ladder = clus0->get_ladder();
       float lz = clus0->get_xyz_local(2);
       int pix_ladder;
       int chip = sensor*4 + (lz/1.3)+2;

       if(layer==0) pix_ladder = ladder;
       else pix_ladder = ladder + 10;

       if(0<=chip&&chip<16&&0<=pix_ladder&&pix_ladder<30) {
	 h2pixel->Fill(chip,pix_ladder);
	 if(!bad_chip(chip,pix_ladder)) {
	   h2pixel2->Fill(chip,pix_ladder);
	   chip_mult = pix_chip_mult[chip][pix_ladder];
	   if(pix_chip_mult[chip][pix_ladder]<5) {
	     h2pixel3->Fill(chip,pix_ladder);
	   }
	 }
       }
     }//if(layer)
     
     if(!bad_sensor(layer, clus0->get_ladder(), clus0->get_sensor())) {
       float phi = atan2(svxy,svxx);
       ntp5[0]=(float)EventNumber;
       ntp5[1]=svxx;
       ntp5[2]=svxy;
       ntp5[3]=svxz;
       ntp5[4]=clus0->get_adc(0);
       ntp5[5]=clus0->get_adc(1);
       ntp5[6]=(float)layer;
       ntp5[7]=(float) clus0->get_hitID();
       ntp5[8]=(float) clus0->get_svxSection();
       ntp5[9]=(float) clus0->get_ladder();
       ntp5[10]=(float) clus0->get_sensor();
       ntp5[11]=(float) clus0->get_sensorType();
       ntp5[12]=(float) clus0->get_edgeflag();
       ntp5[13]= clus0->get_xyz_local(0);
       ntp5[14]= clus0->get_xyz_local(1);
       ntp5[15]= clus0->get_xyz_local(2);
       ntp5[16]= phi;
       ntpsvxcluster->Fill(ntp5);
     }
  }

  ntpevt->Fill(EventNumber,nhit_layer[0],nhit_layer[1],nhit_layer[2],nhit_layer[3],zvtx,nsvxtracks,ntrk,bbcq);

  hzvtx->Fill(zvtx);
*/

  /*
   * offset of the beam center
   */
  float xoffset1 = 0.23;
  float yoffset1 = 0.08;

  /*
   * Central tracks
   */
  vector<float> v_mom;
  vector<float> v_the0;
  vector<float> v_phi0;
  vector<float> v_n0;
  vector<float> v_emcdphi;
  vector<float> v_emcdz;
  vector<float> v_ecore;
  vector<short> v_quality;
  if(ntrk>0) {
    for(int itrk=0;itrk<ntrk;itrk++) {
      v_mom.push_back(trk->get_mom(itrk)*trk->get_charge(itrk));
      v_the0.push_back(trk->get_the0(itrk));
      float phi0_trk = trk->get_phi0(itrk);
      if(phi0_trk > 3.141592) phi0_trk = phi0_trk - 6.283184;
      v_phi0.push_back(phi0_trk);
      v_n0.push_back(trk->get_n0(itrk));
      v_emcdphi.push_back(trk->get_emcdphi(itrk));
      v_emcdz.push_back(trk->get_emcdz(itrk));
      v_quality.push_back(trk->get_quality(itrk));
      v_ecore.push_back(trk->get_ecore(itrk));

      if(trk->get_mom(itrk)>0.5&&trk->get_mom(itrk)<4.0) {
	float mom_cnt    = trk->get_mom(itrk);
	int   charge_cnt = trk->get_charge(itrk);
	float the0_cnt   = trk->get_the0(itrk);
	float phi0_cnt   = trk->get_phi0(itrk);
	if(phi0_cnt > 3.141592) phi0_cnt = phi0_cnt - 6.283184;
	float emcdphi    = trk->get_emcdphi(itrk);
	float emcdz      = trk->get_emcdz(itrk);
	float pt_cnt     = mom_cnt*sin(the0_cnt);
	float zed        = trk->get_zed(itrk);
	//	float zproj0     = zvtx+R0/tan(the0_cnt);
		    
	float ux = cos(phi0_cnt);
	float uy = sin(phi0_cnt);

	if(fabs(emcdz-1.4)<15 && fabs(emcdphi+0.00015)<0.02) {
	  ntp_cnt->Fill(mom_cnt,emcdphi,emcdz,the0_cnt,phi0_cnt,zed,zvtx);
	  for(int i=0;i<nsvx;i++) {
	    SvxCluster *clus0 = svx->get_Cluster(i);
	    int   layer = clus0->get_layer();
	    float svxx = clus0->get_xyz_global(0);
	    float svxy = clus0->get_xyz_global(1);
	    float svxz = clus0->get_xyz_global(2);
	    float phiv = atan2(svxy,svxx);
	    
	    float dx = svxx - xoffset1;
	    float dy = svxy - yoffset1;
	    float rhit_2 = dx*dx+dy*dy;
	    float rhit = sqrt(rhit_2);
	    
	    // vector product (ux,uy) x (dx, dy)
	    // this is signed difference between straight line projection
	    // from beam spot (xoffset1,yoffset1) and a svx hit position.
	    float u_cross_dx = ux*dy - uy*dx;
	    float mag_bend   = 0.0013*charge_cnt*rhit_2/pt_cnt;
	    
	    float dproj = u_cross_dx;
	    float zproj = zvtx + rhit/tan(the0_cnt);
	    
	    if(layer == 0) {
	      if(fabs(phi0_cnt)<1.5) dproj = dproj + 0.01;
	      else dproj = dproj - 0.04;
	    }
	    if(layer == 1) {
	      if(fabs(phi0_cnt)<1.5) dproj = dproj + 0.04;
	      else dproj = dproj - 0.065;
	    }
	    if(layer == 2) {
	      if(fabs(phi0_cnt)<1.5) dproj = dproj + 0.06;
	      else dproj = dproj - 0.14;
	    }
	    if(layer == 3) {
	      if(fabs(phi0_cnt)<1.5) dproj = dproj + 0.1;
	      else dproj = dproj - 0.2;
	    }

	    if(fabs(zproj - svxz)<3) {
	      if(fabs(dproj+mag_bend)<0.1+0.1*layer){
		ntp_cnt_clus->Fill(mom_cnt,emcdphi,emcdz,the0_cnt,phi0_cnt,zed,layer,zproj,dproj,mag_bend,svxz,phiv,zvtx);
	      }
	    }
	    if(layer==3 && fabs(zproj-svxz)<2&& fabs(dproj+mag_bend)<0.4) {
	      for(int j=0;j<nsvx;j++) {
		SvxCluster *clus1 = svx->get_Cluster(j);
		int   layer1 = clus1->get_layer();
		if(layer1 < 3) {
		  float svxx1 = clus1->get_xyz_global(0);
		  float svxy1 = clus1->get_xyz_global(1);
		  float svxz1 = clus1->get_xyz_global(2);
		  float phiv1 = atan2(svxy1,svxx1);
		  
		  float dx1 = svxx1 - xoffset1;
		  float dy1 = svxy1 - yoffset1;
		  float rhit1_2 = dx1*dx1+dy1*dy1;
		  float rhit1 = sqrt(rhit1_2);
		  
		  // vector product (ux,uy) x (dx, dy)
		  // this is signed difference between straight line projection
		  // from beam spot (xoffset1,yoffset1) and a svx hit position.
		  float u_cross_dx1 = ux*dy1 - uy*dx1;
		  float mag_bend1   = 0.0013*charge_cnt*rhit1_2/pt_cnt;
		  
		  float dproj1 = u_cross_dx1;
		  float zproj1 = zvtx + rhit1/tan(the0_cnt);
		  
		  if(layer1 == 0) {
		    if(fabs(phi0_cnt)<1.5) dproj1 = dproj1 + 0.01;
		    else dproj1 = dproj1 - 0.04;
		  }
		  if(layer1 == 1) {
		    if(fabs(phi0_cnt)<1.5) dproj1 = dproj1 + 0.04;
		    else dproj1 = dproj1 - 0.065;
		  }
		  if(layer1 == 2) {
		    if(fabs(phi0_cnt)<1.5) dproj1 = dproj1 + 0.06;
		    else dproj1 = dproj1 - 0.14;
		  }

		  if(fabs(zproj1-svxz1)<2.5&& fabs(dproj1+mag_bend1)<0.1+0.1*layer1) {
		    ntp_cnt_clus1->Fill(mom_cnt,emcdphi,emcdz,the0_cnt,phi0_cnt,zed,zproj,dproj+mag_bend,svxz,zvtx,layer1,dproj1,mag_bend1,zproj1,svxz1);
		  }
		}//if(layer1<3)
	      }//for(nsvx)j
	    }//if(layer==3)
	  }//for(svx)
	}//if(|emcdz|<15 && |emcdphi|<0.02)
      }//if(0.5<mom_cnt<4)
    }
  }

  /*
   * Segment (StandAloneTracking)
   */
  /*
  if(fabs(zvtx)<10) {
    hsvxsegmentmult->Fill(float(nsvxtracks));
    for(int i=0; i<nsvxtracks; i++) {
      SvxSegment* seg0 = svxtracks->get_segment(i);
      float dca3d = seg0->getDCA();  // 3d dca
      float dca2d = seg0->getDCA2D();  // 2d dca (x-y)
      float charge, primary;
      if(seg0->IsPositive()) charge = 1; else charge = -1;
      if(seg0->getPrimary()) primary= 1; else primary = 0;
      
      // print out reconstructed track
      cout << "track "<<i<<endl;
      cout << "mom ="<<seg0->getMomentum()<<endl;
      for(int jj=0;jj<4;jj++) {
	cout << jj<<":(";
	cout <<seg0->getProjectedPosition(jj,0)<<",";
	cout <<seg0->getProjectedPosition(jj,1)<<",";
	cout <<seg0->getProjectedPosition(jj,2)<<")  ";
	
	int clusterID = seg0->getClusterID(jj,0);
	cout << "clusterID="<<clusterID<<endl;
	if(clusterID>=0) {
	  SvxCluster *cluster = svx->get_Cluster(clusterID);
	  if(cluster>0) {
	    cout << "  (";
	    cout << cluster->get_xyz_global(0)<<",";
	    cout << cluster->get_xyz_global(1)<<",";
	    cout << cluster->get_xyz_global(2)<<")"<<endl;
	  } else {
	    cout <<"error..cluster is empty. skip this event"<<endl;
	    return 0;
	  }
	} else {
	  cout<< " ( blanck )" <<endl;
	}
      }
      // calculate projected vertex from L0 and L1 hit
      float x0 = seg0->getProjectedPosition(0,0);
      float y0 = seg0->getProjectedPosition(0,1);
      float z0 = seg0->getProjectedPosition(0,2);
      float r0 = sqrt(x0*x0+y0*y0);
      float x1 = seg0->getProjectedPosition(1,0);
      float y1 = seg0->getProjectedPosition(1,1);
      float z1 = seg0->getProjectedPosition(1,2);
      float r1 = sqrt(x1*x1+y1*y1);

      int pixel0id = seg0->getClusterID(0,0);
      int pixel0ladder = svx->get_Cluster(pixel0id)->get_ladder();
      int pixel0sensor = svx->get_Cluster(pixel0id)->get_sensor();
      float L0x = svx->get_Cluster(pixel0id)->get_xyz_global(0);
      float L0y = svx->get_Cluster(pixel0id)->get_xyz_global(1);
      float L0z = svx->get_Cluster(pixel0id)->get_xyz_global(2);
      float L0phi = atan2(L0y,L0x);

      int pixel1id = seg0->getClusterID(1,0);
      int pixel1ladder = svx->get_Cluster(pixel1id)->get_ladder();
      int pixel1sensor = svx->get_Cluster(pixel1id)->get_sensor();
      float L1x = svx->get_Cluster(pixel1id)->get_xyz_global(0);
      float L1y = svx->get_Cluster(pixel1id)->get_xyz_global(1);
      float L1z = svx->get_Cluster(pixel1id)->get_xyz_global(2);
      float L1phi = atan2(L1y,L1x);

      int strip0id = seg0->getClusterID(2,0);
      int strip0ladder = -1;
      if(strip0id>0) strip0ladder = svx->get_Cluster(strip0id)->get_ladder();
      float L2x=0;
      float L2y=0;
      float L2z=0;
      float r2 = 11.0; //nominal radius of L2
      float L2phi= -10;
      if(strip0id>0) {
	L2x = svx->get_Cluster(strip0id)->get_xyz_global(0);
	L2y = svx->get_Cluster(strip0id)->get_xyz_global(1);
	L2z = svx->get_Cluster(strip0id)->get_xyz_global(2);
	L2phi = atan2(L2y,L2x);
	r2 = sqrt(L2x*L2x + L2y*L2y);
      }
      
      int strip1id = seg0->getClusterID(3,0);
      int strip1ladder = -1;
      if(strip1id>0) strip1ladder = svx->get_Cluster(strip1id)->get_ladder();
      float L3x=0;
      float L3y=0;
      float L3z=0;
      float r3 = 17.0;
      float L3phi= -10;
      if(strip1id>0) {
	L3x = svx->get_Cluster(strip1id)->get_xyz_global(0);
	L3y = svx->get_Cluster(strip1id)->get_xyz_global(1);
	L3z = svx->get_Cluster(strip1id)->get_xyz_global(2);
	L3phi = atan2(L3y,L3x);
	r3 = sqrt(L3x*L3x + L3y*L3y);
      }
      
      float x0c = x0 - xoffset1;
      float x1c = x1 - xoffset1;
      float y0c = y0 - yoffset1;
      float y1c = y1 - yoffset1;
      
      float zvtx_proj = z0 - r0*(z1-z0)/(r1-r0);
      float phi0 = atan2(y0,x0);
      float phi0b = atan2(y0,-x0);

      float yvtx_proj = y0 - x0*(y1-y0)/(x1-x0);
      float xvtx_proj = -1.0;
      if(fabs(y1-y0)>0.0001) xvtx_proj = x0 - y0*(x1-x0)/(y1-y0);

      float yvtx_projc = y0c - x0c*(y1c-y0c)/(x1c-x0c);
      float xvtx_projc = -1.0;
      if(fabs(y1-y0)>0.01) xvtx_projc = x0c - y0c*(x1c-x0c)/(y1c-y0c);

      // DCA in (x,y) plane
      float ux = x1 - x0;
      float uy = y1 - y0;
      float ul = sqrt(ux*ux + uy*uy);
      ux = ux/ul;
      uy = uy/ul;
      // now (ux, uy) is a unit vector.

      float u_dot_dx = ux*x0c + uy*y0c;
      float u_cross_dx = ux*y0c - uy*x0c;

      float dca0 = sqrt(x0c*x0c + y0c*y0c - u_dot_dx*u_dot_dx);
      if(u_cross_dx < 0) dca0 = -dca0;
      float mom = seg0->getMomentum();
      float px  = seg0->get3Momentum(0);
      float py  = seg0->get3Momentum(1);
      float pz  = seg0->get3Momentum(2);
      float pt  = sqrt(px*px+py*py); 
      float magnetic_bend = 0.017*charge/pt;
      float theta0 = atan2(pt,pz);
      float phi00  = atan2(py,px);


      // projection to L2
      //
      float dproj2 = -10;
      float mag_bend2 = -20;
      if(strip0id>0) {
	float dx21 = L2x - x1;
	float dy21 = L2y - y1;
	float u_cross_dx2 = ux*dy21 - uy*dx21;
	// u_cross_dx2 is a vector product of (ux,uy) and (dx21,dy21).
	// its absolute value is the distance between the (L2x,L2y) and the
	// straight line projecton point of the track from (x0,y0) to (x1,y1).
	dproj2 = u_cross_dx2;
	mag_bend2 = 1.2E-3*r2*(r2-r1)*charge/pt;
      }
      
      // projection to L3
      float dproj3 = -10;
      float mag_bend3 = -20;
      if(strip1id>0) {
	float dx31 = L3x - x1;
	float dy31 = L3y - y1;
	float u_cross_dx3 = ux*dy31 - uy*dx31;
	// u_cross_dx3 is a vector product of (ux,uy) and (dx31,dy31).
	// its absolute value is the distance between the (L3x,L3y) and the
	// straight line projecton point of the track from (x0,y0) to (x1,y1).
	dproj3 = u_cross_dx3;
	mag_bend3 = 1.2E-3*r3*(r3-r1)*charge/pt;
      }

      // projection from L2 and L3 (if both exists)
      float dproj32 = -10;
      float mag_bend32= -20;
      if( strip0id>0 && strip1id>0) {
	float dx32 = L3x - L2x;
	float dy32 = L3y - L2y;
	float dx21 = L2x - x1;
	float dy21 = L2y - y1;
	float l21 = sqrt(dx21*dx21+dy21*dy21);
	float ux21 = dx21/l21;
	float uy21 = dy21/l21;
	dproj32 = ux21*dy32 - uy21*dx32;
	mag_bend32 = 1.2E-3*r3*(r3-r2)*charge/pt;
      }

      // search for the closest central
      if(ntrk>0 && mom>0.2 &&abs(theta0-1.57)<0.35 ) {
	for(int itrk=0;itrk<ntrk;itrk++) {
	  hdphi->Fill(phi00-v_phi0[itrk]);
	  hdthe->Fill(theta0-v_the0[itrk]);

	  if(fabs(phi00 - v_phi0[itrk])<0.05) hdthe2->Fill(theta0-v_the0[itrk]);

	  if(fabs(theta0 - v_the0[itrk])<0.03) hdphi2->Fill(phi00-v_phi0[itrk]);
	  if(fabs(theta0 - v_the0[itrk])<0.03) {
	    hdphi3->Fill(phi00-v_phi0[itrk]);
	    float dphi0 = phi00 - v_phi0[itrk] + 0.006/v_mom[itrk];
	    if(fabs(phi00)<1.5) dphi0 = dphi0 + 7.3E-3;
	    else dphi0 = dphi0 - 1.06E-2;

	    if(fabs(dphi0)<0.05) {
	      if(v_mom[itrk]>0) {
		hdphi3a->Fill(phi00-v_phi0[itrk]);
		ntp_trk->Fill(v_mom[itrk],v_phi0[itrk],mom,phi0,phi00,theta0,v_the0[itrk],dphi0,
			      v_n0[itrk],v_emcdphi[itrk],v_emcdz[itrk],v_ecore[itrk],dca0+magnetic_bend);
	      } else {
		hdphi3b->Fill(phi00-v_phi0[itrk]);
		ntp_trk->Fill(v_mom[itrk],v_phi0[itrk],-mom,phi0,phi00,theta0,v_the0[itrk],dphi0,
			      v_n0[itrk],v_emcdphi[itrk],v_emcdz[itrk],v_ecore[itrk],dca0+magnetic_bend);
	      }
	    }
	  }
	}
      }


      ntp6[0]=(float)EventNumber;
      ntp6[1]=mom;
      ntp6[2]=px;
      ntp6[3]=py;
      ntp6[4]=pz;
      ntp6[5]=L0phi;
      ntp6[6]=L1phi;
      ntp6[7]=L2phi;
      ntp6[8]=L3phi;
      ntp6[9]=charge;
      ntp6[10]=(float) seg0->getNhits(0);
      ntp6[11]=(float) seg0->getNhits(1);
      ntp6[12]=(float) seg0->getNhits(2);
      ntp6[13]=(float) seg0->getNhits(3);
      ntp6[14]= seg0->getQuality();
      ntp6[15]= seg0->get_dEdX1();
      ntp6[16]= seg0->get_dEdX2();
      ntp6[17]= dca3d;
      ntp6[18]= dca2d;
      ntp6[19]= nsvxtracks;
      ntp6[20]= zvtx;
      ntp6[21]= zvtx_proj;
      ntp6[22]= xvtx_proj;
      ntp6[23]= yvtx_proj;
      ntp6[24]= phi0;
      ntp6[25]= phi0b;
      ntp6[26]= xvtx_projc;
      ntp6[27]= yvtx_projc;
      ntp6[28]= dca0;
      ntp6[29]=magnetic_bend;
      ntp6[30]= pixel0ladder;
      ntp6[31]= pixel1ladder;
      ntp6[32]= pixel0sensor;
      ntp6[33]= pixel1sensor;
      ntp6[34]= dproj2;
      ntp6[35]= mag_bend2;
      ntp6[36]= dproj3;
      ntp6[37]= mag_bend3;
      ntp6[38]= dproj32;
      ntp6[39]= mag_bend32;
      ntp6[40] = strip0ladder;
      ntp6[41] = strip1ladder;
      ntp6[42] = L0x;
      ntp6[43] = L0y;
      ntp6[44] = L0z;
      ntp6[45] = L1x;
      ntp6[46] = L1y;
      ntp6[47] = L1z;
      ntp6[48] = L2x;
      ntp6[49] = L2y;
      ntp6[50] = L2z;
      ntp6[51] = L3x;
      ntp6[52] = L3y;
      ntp6[53] = L3z;

      ntpsvxsegment->Fill(ntp6);
      }
    
    }//if(zvtx)
  */
  EventNumber++;
  return 0;
}

//==============================================================

int svxana::End(PHCompositeNode *topNode) {
  cout << "svxana::End:  Writing out..." << endl;
  OutputNtupleFile->Write();
  cout << "svxana::End:  Closing output file..." << endl;
  OutputNtupleFile->Close();
  delete OutputNtupleFile;
  OutputNtupleFile=0;
  return 0;
}

