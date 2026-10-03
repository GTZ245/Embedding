
#include <iostream>
#include <iomanip>
#include "gsl/gsl_rng.h"

#include "phool.h"
#include <getClass.h>
#include "PHTypedNodeIterator.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "Fun4AllReturnCodes.h"
#include "Fun4AllServer.h"
#include <recoConsts.h>

#include "SvxEmbedEval.h"

#include "PHCentralTrack.h"
#include "PHTrackOut.h"
#include "PHGlobal.h"
#include "TriggerHelper.h"
#include "VtxOut.h"
#include "ErtOut.h"

#include "SvxGhitList.h"
#include "SvxGhit.h"
#include "SvxRawhitList.h"
#include "SvxRawhit.h"
#include "SvxClusterList.h"
#include "SvxClusterInfo.h"
#include "SvxSegmentList.h"
#include "SvxGhitClusterList.h"
#include "SvxGhitRawhitList.h"
#include "SvxRawhitClusterList.h"
#include "SvxClusterv1.h"
#include "SvxSegmentv1.h"
#include "SvxCentralTrackList.h"
#include "SvxCentralTrack.h"

#include "TFile.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom2.h"

using namespace std;

static SvxCluster* searchSvxHit(SvxClusterList *d_svxhit, int clusterID);
static int findCloseSvxHits(SvxClusterList*,int,float,float,float[],float[]);
static int layerNumber(float);
static const int MAXHIT=50;

//==============================================================

SvxEmbedEval::SvxEmbedEval() {
  ThisName = "SvxEmbedEval";
  init_ana=0;
  EventNumber=0;
  OutputFileName="svxembedeval.root";
  OutputNtupleFile = NULL;
  d_conversionVeto = false;
  d_conversionVeto1 = false;
  m_mcnode = NULL;
  m_realnode = NULL;
  ntpvtx = NULL;
  ntpeval = NULL;
  n0=0;
  n1=0;
  n2=0;
}

//==============================================================

SvxEmbedEval::SvxEmbedEval(std::string filename) {
  ThisName = "SvxEmbedEval";
  init_ana=0;
  EventNumber=0;
  OutputFileName=filename;
  OutputNtupleFile = NULL;
  d_conversionVeto = false;
  d_conversionVeto1 = false;
  m_mcnode = NULL;
  m_realnode = NULL;
  ntpvtx = NULL;
  ntpeval = NULL;
}

//==============================================================

int SvxEmbedEval::Init(PHCompositeNode *topNode) {

  if ( verbosity>0 ) std::cout << "SvxEmbedEval::Init started..." << std::endl;
  OutputNtupleFile = new TFile(OutputFileName.c_str(),"RECREATE");
  if ( verbosity>0 ) std::cout << "SvxEmbedEval::Init: output file " << OutputFileName << " opened." << std::endl;

  ntpvtx = new TNtuple("ntpvtx","","bbccharge:vtxx:vtxy:xtvz:vtxtype:cent");
  ntpeval = new TNtuple("ntpeval","","cent:bbccharge:charge:pt:mom:qual:n0:emce:snglpt:snglnhit:sngldca2d:snglchi2phi:snglchi2z:snglchi2:snglndfdphi:snglndfdz:snglndf:embedpt:embednhit:embeddca2d:embedchi2phi:embedchi2z:embedchi2:embedndfdphi:embedndfdz:embedndf:snglveto:embedveto:snglhitpattern:embedhitpattern:zed:emcdphi:emcdz:emcsdphi:emcsdz:dep:sngldcaz:embeddcaz:sngllinkscore:embedlinkscore");

  if ( verbosity>0 ) std::cout << "SvxEmbedEval::Init ended." << std::endl;
  return 0;
}

//==============================================================
  
int SvxEmbedEval::InitRun(PHCompositeNode *topNode) {
  if ( verbosity>0 ) std::cout << "SvxEmbedEval::InitRun started..." << std::endl;

  Fun4AllServer* se = Fun4AllServer::instance();
  //recoConsts* rc = recoConsts::instance();

  m_mcnode = se->topNode("SINGLE");
  m_realnode = se->topNode("REAL");

  if ( verbosity>0 ) std::cout << "SvxEmbedEval::InitRun ended." << std::endl;
  return 0;
}

//==============================================================

int SvxEmbedEval::process_event(PHCompositeNode *topNode) {
EventNumber++;
float ntp[99]; for(int i=0; i<99; i++) {ntp[i]=-9999.;} 
if(verbosity>1 && EventNumber%1==0) cout << "=================================== Event # " << EventNumber << endl;

  PHGlobal* evt = findNode::getClass<PHGlobal>(topNode,"PHGlobal");
  //PHGlobal* evt = findNode::getClass<PHGlobal>(m_realnode,"PHGlobal");
  if (!evt) {std::cout << PHWHERE << " ERROR: PHGlobal not in Node Tree" << std::endl; return DISCARDEVENT;}
  float cent = (float)evt->getCentrality();
  if ( verbosity>1 ) cout << "centrality = " << cent << " " << evt->getBbcChargeN() << " " << evt->getBbcChargeS() << endl;
  //float bbczvtx = evt->getBbcZVertex();
  //cout << "bbcZvertex = " << bbczvtx << endl;

/*
  VtxOut *d_vtx = findNode::getClass<VtxOut>(topNode,"VtxOut");
  if(!d_vtx) { cerr << PHWHERE << "VtxOut node not found." << endl; return DISCARDEVENT;}
  float xvtx0 =  (d_vtx->get_Vertex()).getX();
  float yvtx0 =  (d_vtx->get_Vertex()).getY();
  float zvtx0 =  (d_vtx->get_Vertex()).getZ();
  float xvtx  =  (d_vtx->get_Vertex("SIM")).getX();
  float yvtx  =  (d_vtx->get_Vertex("SIM")).getY();
  float zvtx  =  (d_vtx->get_Vertex("SIM")).getZ();
  float xvtx1  =  (d_vtx->get_Vertex("BBC")).getX();
  float yvtx1  =  (d_vtx->get_Vertex("BBC")).getY();
  float zvtx1  =  (d_vtx->get_Vertex("BBC")).getZ();
  float xvtx2  =  (d_vtx->get_Vertex("SVX")).getX();
  float yvtx2  =  (d_vtx->get_Vertex("SVX")).getY();
  float zvtx2  =  (d_vtx->get_Vertex("SVX")).getZ();
  float xvtx3  =  (d_vtx->get_Vertex("SVX_PRECISE")).getX();
  float yvtx3  =  (d_vtx->get_Vertex("SVX_PRECISE")).getY();
  float zvtx3  =  (d_vtx->get_Vertex("SVX_PRECISE")).getZ();
  string s_vtx = d_vtx->which_Vtx();
  if ( verbosity>1 ) {
    cout << "SIM Event vertex = " << xvtx << " " << yvtx << " " << zvtx << endl;
    cout << "BBC Event vertex = " << xvtx1 << " " << yvtx1 << " " << zvtx1 << endl;
    cout << "SVX Event vertex = " << xvtx2 << " " << yvtx2 << " " << zvtx2 << endl;
    cout << "SVX_PRECISE Event vertex = " << xvtx3 << " " << yvtx3 << " " << zvtx3 << endl;
    cout << "Default Event vertex = " << xvtx0 << " " << yvtx0 << " " << zvtx0 << " " << s_vtx << endl;
  }
*/

  PHCentralTrack* trkembed = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack");
  //PHCentralTrack* trkembed = findNode::getClass<PHCentralTrack>(m_realnode, "PHCentralTrack");
  if(!trkembed) {cerr << PHWHERE << "ERROR: No REAL PHCentralTrack object" << endl; return DISCARDEVENT;}
  int ntrkembed = trkembed->get_npart();
  if ( verbosity>1 ) cout << "Number of REAL CNT tracks = " << ntrkembed << endl;

  SvxCentralTrackList* svxcentraltracksembed = findNode::getClass<SvxCentralTrackList>(topNode, "SvxCentralTrackList");
  //SvxCentralTrackList* svxcentraltracksembed = findNode::getClass<SvxCentralTrackList>(m_realnode, "SvxCentralTrackList");
  if (!svxcentraltracksembed) { cerr << PHWHERE << " No REAL SvxCentralTrackList object found." << endl; return DISCARDEVENT; }
  int nsvxcenttrkembed = svxcentraltracksembed->get_nCentralTracks();
  if ( verbosity>1 ) cout << "# of REAL SvxCentralTracks = " << nsvxcenttrkembed << endl;

  SvxClusterList* d_svxclustembed = findNode::getClass<SvxClusterList>(topNode, "SvxClusterList");
  if (!d_svxclustembed) { cerr << PHWHERE << " No REAL SvxCluserList object found." << endl; return DISCARDEVENT; }
  int nsvxclustembed = d_svxclustembed->get_nClusters();
  if ( verbosity>1 ) cout << "# of REAL SvxClusters = " << nsvxclustembed << endl;


// single simulated particle

  PHCentralTrack* trksngl = findNode::getClass<PHCentralTrack>(m_mcnode, "PHCentralTrack");
  if(!trksngl) {cerr << PHWHERE << "ERROR: No SINGLE PHCentralTrack object" << endl; return DISCARDEVENT;}
  int ntrksngl = trksngl->get_npart();
  if ( verbosity>1 ) cout << "Number of SINGLE CNT tracks = " << ntrksngl << endl;

  SvxCentralTrackList* svxcentraltrackssngl = findNode::getClass<SvxCentralTrackList>(m_mcnode, "SvxCentralTrackList");
  if (!svxcentraltrackssngl) { cerr << PHWHERE << " No SINGLE SvxCentralTrackList object found." << endl; return DISCARDEVENT; }
  int nsvxcenttrksngl = svxcentraltrackssngl->get_nCentralTracks();
  if ( verbosity>1 ) cout << "# of SINGLE SvxCentralTracks = " << nsvxcenttrksngl << endl;

  SvxClusterList* d_svxclustsngl = findNode::getClass<SvxClusterList>(m_mcnode, "SvxClusterList");
  if (!d_svxclustsngl) { cerr << PHWHERE << " No SINGLE SvxCluserList object found." << endl; return DISCARDEVENT; }
  int nsvxclustsngl = d_svxclustsngl->get_nClusters();
  if ( verbosity>1 ) cout << "# of SINGLE SvxClusters = " << nsvxclustsngl << endl;


  if(ntrksngl!=1) return 0;

  ntp[0] = cent;
  ntp[1] = evt->getBbcChargeN()+evt->getBbcChargeS();
 

// BEFORE embedding

  n0 += trksngl->get_npart();

  int charge = 0;

  for(unsigned int i=0; i<trksngl->get_npart(); i++) {
    double mom = trksngl->get_mom(i);
    double px = trksngl->get_px(i);
    double py = trksngl->get_py(i);
    //double pz = trksngl->get_pz(i);
    double pt = sqrt(px*px+py*py);
    float n0 = trksngl->get_n0(i);
    float emce = trksngl->get_emce(i);
    float dep = trksngl->get_dep(i);
    float zed = trksngl->get_zed(i);
    float emcdphi = trksngl->get_emcdphi(i);
    float emcdz = trksngl->get_emcdz(i);
    float emcsdphi = trksngl->get_emcsdphi(i);
    float emcsdz = trksngl->get_emcsdz(i);
    int qual = trksngl->get_quality(i);
    //int emcptr = trksngl->get_emcid(i);
    //int pc1ptr = trksngl->get_pc1id(i);
    charge = trksngl->get_charge(i);
    //cout << "SINGLE PHcentralTrack: " << charge << " " << pt << " " << mom << " " << qual << " " << n0 << " " << emce << endl;
    ntp[2] = charge;
    ntp[3] = pt;
    ntp[4] = mom;
    ntp[5] = qual;
    ntp[6] = n0;
    ntp[7] = emce;
    ntp[30] = zed;
    ntp[31] = emcdphi;
    ntp[32] = emcdz;
    ntp[33] = emcsdphi;
    ntp[34] = emcsdz;
    ntp[35] = dep;
  }

  n1 += nsvxcenttrksngl;

  for(int i=0; i<nsvxcenttrksngl; i++) {
    SvxCentralTrack* tmp = svxcentraltrackssngl->getCentralTrack(i);
    double px = tmp->get3MomentumAtPrimaryVertex(0);
    double py = tmp->get3MomentumAtPrimaryVertex(1);
    double pz = tmp->get3MomentumAtPrimaryVertex(2);
    double pt = sqrt(px*px+py*py);
    double mom = sqrt(px*px+py*py+pz*pz);;
    //int dchid = tmp->getDchIndex();
    //int charge = trksngl->get_charge(dchid);
    double dca2d = tmp->getDCA2D();
    double dcaz = tmp->getDCAZ();
    //double dcax = tmp->getClosestApproach(0);
    //double dcay = tmp->getClosestApproach(1);
    //double dcaz = tmp->getClosestApproach(2);
    double chi2 = tmp->getChiSquare();
    double chi2phi = tmp->getChiSquareDPHI();
    double chi2z = tmp->getChiSquareDZ();
    double ndf = tmp->getNDF();
    double ndfdphi = tmp->getNDFDZ();
    double ndfdz = tmp->getNDFDPHI();
    int svxnhits = tmp->getNhits();
    int hitpattern = tmp->getHitPattern();
    float linkscore = tmp->getLinkScore();
    //cout << "SINGLE: " << pt << " " << mom << " " << svxnhits << " " << dca2d << " " << chi2phi << " " << chi2z << " " << chi2 << " " << ndf << endl;

    bool electron = true;
    d_conversionVeto = false;
    d_conversionVeto1 = false;
    ConversionVeto(cent, pt, charge, electron, mom,
                   tmp,
                   d_svxclustsngl);
    if( verbosity>1 ) cout << "SINGLE VETO: " << d_conversionVeto << " " << d_conversionVeto1 << endl;

    ntp[8] = pt; 
    ntp[9] = svxnhits;
    ntp[10] = dca2d;
    ntp[11] = chi2phi;
    ntp[12] = chi2z;
    ntp[13] = chi2;
    ntp[14] = ndfdphi;
    ntp[15] = ndfdz;
    ntp[16] = ndf;
    ntp[26] = 0.; if(d_conversionVeto) ntp[26] = 1.;
    ntp[28] = (float)hitpattern;
    ntp[36] = dcaz;
    ntp[38] = linkscore;
  }

// AFTER embedding

  n2 += nsvxcenttrkembed;

  for(int i=0; i<nsvxcenttrkembed; i++) {
    SvxCentralTrack* tmp = svxcentraltracksembed->getCentralTrack(i);
    double px = tmp->get3MomentumAtPrimaryVertex(0);
    double py = tmp->get3MomentumAtPrimaryVertex(1);
    double pz = tmp->get3MomentumAtPrimaryVertex(2);
    double pt = sqrt(px*px+py*py);
    double mom = sqrt(px*px+py*py+pz*pz);;
    //int dchid = tmp->getDchIndex();
    //int charge = trksngl->get_charge(dchid);
    double dca2d = tmp->getDCA2D();
    double dcaz = tmp->getDCAZ();
    //double dcax = tmp->getClosestApproach(0);
    //double dcay = tmp->getClosestApproach(1);
    //double dcaz = tmp->getClosestApproach(2);
    double chi2 = tmp->getChiSquare();
    double chi2phi = tmp->getChiSquareDPHI();
    double chi2z = tmp->getChiSquareDZ();
    double ndf = tmp->getNDF();
    double ndfdphi = tmp->getNDFDPHI();
    double ndfdz = tmp->getNDFDZ();
    int svxnhits = tmp->getNhits();
    int hitpattern = (int)tmp->getHitPattern();
    float linkscore = tmp->getLinkScore();
    ntp[17] = pt; 
    ntp[18] = svxnhits;

    bool electron = true;
    d_conversionVeto = false;
    d_conversionVeto1 = false;
    ConversionVeto(cent, pt, charge, electron, mom,
                   tmp,
                   d_svxclustembed);
    if( verbosity>1 ) cout << "EMBED VETO: " << d_conversionVeto << " " << d_conversionVeto1 << endl;

    ntp[19] = dca2d;
    ntp[20] = chi2phi;
    ntp[21] = chi2z;
    ntp[22] = chi2;
    ntp[23] = ndfdphi;
    ntp[24] = ndfdz;
    ntp[25] = ndf;
    ntp[27] = 0.; if(d_conversionVeto) ntp[27] = 1.;
    ntp[29] = (float)hitpattern;
    ntp[37] = dcaz;
    ntp[39] = linkscore;
  }

  ntpeval->Fill(ntp);

cout << "n: " << n0 << " " << n1 << " " << n2 << endl;
return 0;
}

//==============================================================

int SvxEmbedEval::End(PHCompositeNode *topNode) {
  if ( verbosity>1 ) std::cout << "SvxEmbedEval::End:  Writing out..." << std::endl;
  OutputNtupleFile->Write();
  if ( verbosity>1 ) std::cout << "SvxEmbedEval::End:  Closing output file..." << std::endl;
  OutputNtupleFile->Close();
  delete OutputNtupleFile;
  OutputNtupleFile=0;
  cout << "final n: " << n0 << " " << n1 << " " << n2 << endl;
  return 0;
}

//==============================================================

void SvxEmbedEval::ConversionVeto(float centrality,
                 float pt, float charge, bool electron, float mom,
                 SvxCentralTrack *svxcnt,
                 SvxClusterList* d_svxhit)
{

  int d_fieldPolarity=1; // ???

  float vphi[MAXHIT],vz[MAXHIT];

  int d_B0p_near_hit=0;
  int d_B0m_near_hit=0;
  int d_B1p_near_hit=0;
  int d_B1m_near_hit=0;
  int d_B2_near_hit=0;
  int d_B3p_near_hit=0;
  int d_B3m_near_hit=0;

  int d_B0m1hit=0;
  int d_B0m2hit=0;
  int d_B0p1hit=0;
  int d_B0p2hit=0;
  int d_B0p3hit=0;

  int d_B1m1hit=0;
  int d_B1m2hit=0;
  int d_B1p1hit=0;
  int d_B1p2hit=0;
  int d_B1p3hit=0;

  int d_B2m1hit=0;
  int d_B2p1hit=0;
  int d_B2p2hit=0;

  int d_B3m1hit=0;
  int d_B3p1hit=0;
  int d_B3p2hit=0;

    int nclust = svxcnt->getNhits();

    for(int ic=0;ic<nclust;ic++) {
      SvxClusterInfo *info = svxcnt->getClusterInfo(ic);
      if(info==NULL) continue;

      int layer          = info->getLayer();
      int clusterID      = info->getClusterId();
      SvxCluster *svxhit = searchSvxHit(d_svxhit,clusterID);

      if(svxhit == NULL) {
        cout << "svxhit is not found"<<endl;
      } else {
        float xhit = svxhit->get_xyz_global(0);
        float yhit = svxhit->get_xyz_global(1);
        float zhit = svxhit->get_xyz_global(2);

        float phihit  = atan2(yhit,xhit);
        if(phihit< - 1.5708) phihit += 6.283184;
        int nfound = findCloseSvxHits(d_svxhit,layer,phihit,zhit,vphi,vz);

        for(int ifound=0;ifound<nfound;ifound++) {

          float dphi=vphi[ifound]-phihit;
          float dz = vz[ifound]-zhit;
          float cdphi = charge*dphi*(-d_fieldPolarity);

          if(layer==0&&fabs(dz)<0.05&&0.001<cdphi&&cdphi<0.04) {
            d_B0p_near_hit++;
          }
          if(layer==0&&fabs(dz)<0.05&&-0.02<cdphi&&cdphi<-0.001) {
            d_B0m_near_hit++;
          }

          if(layer==0&&fabs(dz)<0.05) {
            if(-0.04<=cdphi&&cdphi<-0.02) d_B0m1hit++;
            if(-0.02<=cdphi&&cdphi<-0.001) d_B0m2hit++;
            if(0.001<=cdphi&&cdphi<0.02)  d_B0p1hit++;
            if(0.02<=cdphi&&cdphi<0.04)   d_B0p2hit++;
            if(0.04<=cdphi&&cdphi<0.08)   d_B0p3hit++;
          }

          if(layer==1&&fabs(dz)<0.05&&0.001<cdphi&&cdphi<0.06) {
            d_B1p_near_hit++;
          }
          if(layer==1&&fabs(dz)<0.05&&-0.02<cdphi&&cdphi<-0.001) {
            d_B1m_near_hit++;
          }

          if(layer==1&&fabs(dz)<0.05) {
            if(-0.04<=cdphi&&cdphi<-0.02) d_B1m1hit++;
            if(-0.02<=cdphi&&cdphi<-0.001) d_B1m2hit++;
            if(0.001<=cdphi&&cdphi<0.02)   d_B1p1hit++;
            if(0.02 <=cdphi&&cdphi<0.05)   d_B1p2hit++;
            if(0.05 <=cdphi&&cdphi<0.08)   d_B1p3hit++;
          }
          if(layer==2&&fabs(dz)<0.1&&
             ((0.002<cdphi&&cdphi<0.08)||(-0.04<cdphi&&cdphi<-0.002))) {
            d_B2_near_hit++;
          }

          if(layer==2&&fabs(dz)<0.1) {
            if(-0.04<=cdphi&&cdphi<-0.001) d_B2m1hit++;
            if(0.001<=cdphi&&cdphi<0.04)   d_B2p1hit++;
            if(0.04 <=cdphi&&cdphi<0.08)   d_B2p2hit++;
          }

          if(layer==3&&fabs(dz)<0.1&&0.001<cdphi&&cdphi<0.08) {
          d_B3p_near_hit++;
          }
          if(layer==3&&fabs(dz)<0.1&&-0.02<cdphi&&cdphi<-0.001) {
            d_B3m_near_hit++;
          }

          if(layer==3&&fabs(dz)<0.1) {
            if(-0.03<=cdphi&&cdphi<-0.001) d_B3m1hit++;
            if(0.001<=cdphi&&cdphi<0.04)   d_B3p1hit++;
            if(0.04 <=cdphi&&cdphi<0.08)   d_B3p2hit++;
          }

        }//for(ifound)
      }//if(svxhit)
    }//for(ic)

  d_conversionVeto  =
              (d_B0m_near_hit==0 && d_B0p_near_hit==0 &&
               d_B1m_near_hit==0 && d_B1p_near_hit==0 &&
               d_B2_near_hit==0  &&
               d_B3m_near_hit==0 && d_B3p_near_hit==0);

  d_conversionVeto1 =
              (d_B2_near_hit==0  &&
               d_B3m_near_hit==0 && d_B3p_near_hit==0);

}

//================================================================

static SvxCluster* searchSvxHit(SvxClusterList *d_svxhit, int clusterID) {
  if(d_svxhit == NULL) {
    cout << "searchSvxHit: d_svxhit is NULL!"<<endl;
    return NULL;
  } else {
    int nhit = d_svxhit->get_nClusters();
    for(int i=0; i<nhit; i++) {
      SvxCluster *svxhit = d_svxhit->get_Cluster(i);
      if(svxhit==NULL) {
         cout<<"cluster NULL : "<<i<<endl;
         return NULL;
      }

      if(svxhit->get_hitID() == clusterID) return svxhit;
    }
    cout<<"no cluster found in SvxHitMapEntry"<<endl;
    return NULL;
  }
}

//======================================================================

static int findCloseSvxHits(SvxClusterList *d_svxhit,
                            int layer0, float phi0,float zhit0,
                            float vphi[],float vz[])
{
  const static float philimit[4]={0.12,0.12,0.12,0.12};
  if(d_svxhit==NULL) return 0;
  int nfound = 0;

  int nhit = d_svxhit->get_nClusters();
  for(int ihit=0; ihit<nhit; ihit++) {
    SvxCluster *svxhit = d_svxhit->get_Cluster(ihit);
    float xhit = svxhit->get_xyz_global(0);
    float yhit = svxhit->get_xyz_global(1);
    float zhit = svxhit->get_xyz_global(2);
    float rhit = sqrt(xhit*xhit+yhit*yhit);
    float phi  = atan2(yhit,xhit);
    if(phi < -1.5708) phi += 6.283184;
    int layer  = layerNumber(rhit);
    if(layer == layer0) {
      if(fabs(phi-phi0)<philimit[layer] && fabs(zhit-zhit0)<0.4) {
        vphi[nfound]=phi;
        vz[nfound]=zhit;
        nfound++;
      }
   }
    if(nfound>=MAXHIT) return nfound;
  }
  return nfound;
}

//================================================================

static int layerNumber(float rhit)
{
  if(rhit<3.0) return 0;
  else if(rhit<6.0) return 1;
  else if(rhit<14.0) return 2;
  else return 3;
}







