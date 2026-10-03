#include <iostream>
#include <iomanip>

#include "Fun4AllReturnCodes.h"

#include "svxcentclsana.h"
#include "SvxCentralTrackReco.h"

#include "PHSnglCentralTrack.h"
#include "VtxOut.h"
#include "Bbc.hh"
#include "BbcOut.h"
#include "PHPoint.h"

#include "SvxClusterList.h"
#include "SvxCluster.h"

#include "McEvalSingleList.h"

#include "getClass.h"

#include "TFile.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TTree.h"

using namespace std;
using namespace findNode;
//==============================================================

svxcentclsana::svxcentclsana() {
  ThisName = "svxcentclsana";

  m_outFile    =NULL;
  m_outFileName="svxanalysis.root";
  m_initana    =0;
  m_eventNumber=0;
  m_mod=NULL;
  m_simmode=false;
}

//==============================================================

svxcentclsana::svxcentclsana(string filename) {
  ThisName = "svxcentclsana";

  m_outFile    =NULL;
  m_outFileName=filename;
  m_initana    =0;
  m_eventNumber=0;
  m_mod=NULL;
  m_simmode=false;
}

//==============================================================

svxcentclsana::~svxcentclsana() {
}

//==============================================================

int svxcentclsana::Init(PHCompositeNode *topNode) {

  cout << "svxcentclsana::Init started..." << endl;
  m_outFile = new TFile(m_outFileName.c_str(),"RECREATE");
  cout << "svxcentclsana::Init: output file " << m_outFileName << " opened." << endl;

/*
  m_ntp_cntclus = new TNtuple("ntp_cntclus", "","zvtx:mom:the0:phi0:zed:emcdphi:emcdz:" 
                                                "ly3:ld3:zproj3:dproj3:bend3:zv3:phiv3:"
                                                "ly2:ld2:zproj2:dproj2:bend2:zv2:phiv2:"
                                                "ly1:ld1:zproj1:dproj1:bend1:zv1:phiv1:"
                                                "ly0:ld0:zproj0:dproj0:bend0:zv0:phiv0"
*/
                                                                                                     
  m_datap[0] = &m_ly0;    m_datap[ 9] = &m_ly1;    m_datap[18] = &m_ly2;    m_datap[27] = &m_ly3;
  m_datap[1] = &m_ld0;    m_datap[10] = &m_ld1;    m_datap[19] = &m_ld2;    m_datap[28] = &m_ld3;
  m_datap[2] = &m_zproj0; m_datap[11] = &m_zproj1; m_datap[20] = &m_zproj2; m_datap[29] = &m_zproj3;
  m_datap[3] = &m_dproj0; m_datap[12] = &m_dproj1; m_datap[21] = &m_dproj2; m_datap[30] = &m_dproj3;
  m_datap[4] = &m_bend0;  m_datap[13] = &m_bend1;  m_datap[22] = &m_bend2;  m_datap[31] = &m_bend3;
  m_datap[5] = &m_zv0;    m_datap[14] = &m_zv1;    m_datap[23] = &m_zv2;    m_datap[32] = &m_zv3;
  m_datap[6] = &m_phiv0;  m_datap[15] = &m_phiv1;  m_datap[24] = &m_phiv2;  m_datap[33] = &m_phiv3;
  m_datap[7] = &m_l0;     m_datap[16] = &m_l1;     m_datap[25] = &m_l2;     m_datap[34] = &m_l3;
  m_datap[8] = &m_size0;  m_datap[17] = &m_size1;  m_datap[26] = &m_size2;  m_datap[35] = &m_size3;
                                  
  m_ntp_cntclus = new TTree("ntp_cntclus", "association tree");
  m_ntp_cnttrk  = new TTree("ntp_cnttrk",  "CNTtrack with VTX info tree");

  TTree *t[2] = {m_ntp_cntclus, m_ntp_cnttrk};
  for(int i=0; i<2; i++){
    t[i]->Branch("zvtx",    &m_zvtx,     "zvtx/F");
    t[i]->Branch("vtxid",   &m_vtxid,    "vtxid/I");
    t[i]->Branch("xvtx",    &m_xvtx,     "xvtx/F");
    t[i]->Branch("yvtx",    &m_yvtx,     "yvtx/F");
    t[i]->Branch("bbcq",    &m_bbcq,     "bbcq/F");
    t[i]->Branch("t0",      &m_t0,       "t0/F");
    t[i]->Branch("mom",     &m_mom,      "mom/F");
    t[i]->Branch("the0",    &m_the0,     "the0/F");
    t[i]->Branch("phi0",    &m_phi0,     "phi0/F");
    t[i]->Branch("zed",     &m_zed,      "zed/F");
    t[i]->Branch("dcqual",  &m_dcqual,   "dcqual/I");
    t[i]->Branch("emcdphi", &m_emcdphi,  "emcdphi/F");
    t[i]->Branch("emcdz",   &m_emcdz,    "emcdz/F");
    t[i]->Branch("tofdphi", &m_tofdphi,  "tofdphi/F");
    t[i]->Branch("tofdz",   &m_tofdz,    "tofdz/F");
    t[i]->Branch("tofwdphi",&m_tofwdphi, "tofwdphi/F");
    t[i]->Branch("tofwdz",  &m_tofwdz,   "tofwdz/F");
    t[i]->Branch("m2emc",   &m_m2emc,    "m2emc/F");
    t[i]->Branch("m2tof",   &m_m2tof,    "m2tof/F");
    t[i]->Branch("m2tofw",  &m_m2tofw,   "m2tofw/F");
    t[i]->Branch("temc",   &m_temc,    "temc/F");
    t[i]->Branch("ttof",   &m_ttof,    "ttof/F");
    t[i]->Branch("ttofw",  &m_ttofw,   "ttofw/F");
    t[i]->Branch("plemc",   &m_plemc,    "plemc/F");
    t[i]->Branch("pltof",   &m_pltof,    "pltof/F");
    t[i]->Branch("pltofw",   &m_pltofw,  "pltofw/F");
    t[i]->Branch("n0",     &m_n0,    "n0/F");
    t[i]->Branch("cch2",   &m_cch2,  "cch2/F");
    t[i]->Branch("npe0",   &m_npe0,  "npe0/F");
    t[i]->Branch("disp",   &m_disp,  "disp/F");
    t[i]->Branch("sn0",    &m_sn0,   "sn0/F");
    t[i]->Branch("scch2",  &m_scch2, "scch2/F");
    t[i]->Branch("snpe0",  &m_snpe0, "snpe0/F");
    t[i]->Branch("sdisp",  &m_sdisp, "sdisp/F");
    t[i]->Branch("ecore",  &m_ecore,  "ecore/F");
    t[i]->Branch("ecorr",  &m_ecorr,  "ecorr/F");

    t[i]->Branch("lid",     &m_lid,      "lid/I");
    t[i]->Branch("ly3",     &m_ly3,      "ly3/F");
    t[i]->Branch("ld3",     &m_ld3,      "ld3/F");
    t[i]->Branch("zproj3",  &m_zproj3,   "zproj3/F");
    t[i]->Branch("dproj3",  &m_dproj3,   "dproj3/F");
    t[i]->Branch("bend3",   &m_bend3,    "bend3/F");
    t[i]->Branch("zv3",     &m_zv3,      "zv3/F");
    t[i]->Branch("phiv3",   &m_phiv3,    "phiv3/F");
    t[i]->Branch("l3",      &m_l3,       "l3/F");
    t[i]->Branch("size3",   &m_size3,    "size3/F");
    t[i]->Branch("ly2",     &m_ly2,      "ly2/F");
    t[i]->Branch("ld2",     &m_ld2,      "ld2/F");
    t[i]->Branch("zproj2",  &m_zproj2,   "zproj2/F");
    t[i]->Branch("dproj2",  &m_dproj2,   "dproj2/F");
    t[i]->Branch("bend2",   &m_bend2,    "bend2/F");
    t[i]->Branch("zv2",     &m_zv2,      "zv2/F");
    t[i]->Branch("phiv2",   &m_phiv2,    "phiv2/F");
    t[i]->Branch("l2",      &m_l2,       "l2/F");
    t[i]->Branch("size2",   &m_size2,    "size2/F");
    t[i]->Branch("ly1",     &m_ly1,      "ly1/F");
    t[i]->Branch("ld1",     &m_ld1,      "ld1/F");
    t[i]->Branch("zproj1",  &m_zproj1,   "zproj1/F");
    t[i]->Branch("dproj1",  &m_dproj1,   "dproj1/F");
    t[i]->Branch("bend1",   &m_bend1,    "bend1/F");
    t[i]->Branch("zv1",     &m_zv1,      "zv1/F");
    t[i]->Branch("phiv1",   &m_phiv1,    "phiv1/F");
    t[i]->Branch("l1",      &m_l1,       "l1/F");
    t[i]->Branch("size1",   &m_size1,    "size1/F");
    t[i]->Branch("ly0",     &m_ly0,      "ly0/F");
    t[i]->Branch("ld0",     &m_ld0,      "ld0/F");
    t[i]->Branch("zproj0",  &m_zproj0,   "zproj0/F");
    t[i]->Branch("dproj0",  &m_dproj0,   "dproj0/F");
    t[i]->Branch("bend0",   &m_bend0,    "bend0/F");
    t[i]->Branch("zv0",     &m_zv0,      "zv0/F");
    t[i]->Branch("phiv0",   &m_phiv0,    "phiv0/F");
    t[i]->Branch("l0",      &m_l0,       "l0/F");
    t[i]->Branch("size0",   &m_size0,    "size0/F");
    t[i]->Branch("chi2",    &m_chi2,     "chi2/F");
    t[i]->Branch("ndf",     &m_ndf,      "ndf/I");
    t[i]->Branch("best",    &m_best,     "best/I");
    t[i]->Branch("nhit",    &m_nhit,     "nhit/I");
    t[i]->Branch("d2dca",   &m_d2dca,    "d2dca/F");
    t[i]->Branch("dL",      &m_dL,       "dL/F");
    t[i]->Branch("dbend",   &m_dbend,    "dbend/F");
    t[i]->Branch("d2dca0",  &m_d2dca0,   "d2dca0/F");
    t[i]->Branch("d2dca1",  &m_d2dca1,   "d2dca1/F");
    t[i]->Branch("d2dca2",  &m_d2dca2,   "d2dca2/F");
    t[i]->Branch("dedx1",   &m_dedx1,    "dedx1/F");
    t[i]->Branch("dedx2",   &m_dedx2,    "dedx2/F");

    // forsim
    if(m_simmode){
      t[i]->Branch("simdca",  &m_simdca,  "simdca/F");
      t[i]->Branch("simpt",   &m_simpt,   "simpt/F");
      t[i]->Branch("simphi0", &m_simphi0, "simphi0/F");
      t[i]->Branch("simvx",   &m_simvx,   "simvx/F");
      t[i]->Branch("simvy",   &m_simvy,   "simvy/F");
      t[i]->Branch("simvz",   &m_simvz,   "simvz/F");
    }
  }

  return 0;
}

//==============================================================
  
int svxcentclsana::InitRun(PHCompositeNode *topNode) {
  cout << "svxcentclsana::InitRun started..." << endl;
  cout << "svxcentclsana::InitRun ended." << endl;
  return 0;
}

//==============================================================


int svxcentclsana::process_event(PHCompositeNode *topNode) {
  cout<<"svxcentclsana process_event "<<endl;
  cout<<"svxcentclsana   modname : "<<( m_mod!=NULL ? m_mod->Name() : "NULL") <<endl;

  VtxOut *vtxout = getClass<VtxOut>(topNode,"VtxOut");
  if(vtxout==NULL){
    cout<<"no Vtxout, Skip"<<endl;
    return ABORTEVENT;
  }

  BbcOut *bbcout = getClass<BbcOut>(topNode,"BbcOut");
  if(bbcout==NULL){
    cout<<"no Bbcout, Skip"<<endl;
//    return ABORTEVENT;
  }
  m_bbcq = (bbcout!=NULL) ? (bbcout->get_ChargeSum(Bbc::North) + bbcout->get_ChargeSum(Bbc::South)) : 50.0;
  m_t0   = (bbcout!=NULL) ? bbcout->get_TimeZero() : -9999.0;

  
  McEvalSingleList *mctrk = getClass<McEvalSingleList>(topNode,"McSingle");
  if(m_simmode){
    if(mctrk==NULL){
      cout<<"No McSingle Info in Nodetree"<<endl;
    }
  }


  m_outFile->cd();

  std::vector<SvxCentralClusterLink*>& trklist = m_mod->getTrackList();
  cout<<"svxcentclsana   ntrk : "<<trklist.size() <<endl;

  int ntrk = trklist.size();

  if(ntrk>0) {
    PHPoint vtx = vtxout->get_Vertex();



    for(int itrk=0;itrk<ntrk;itrk++) {
      SvxCentralClusterLink *trklink = trklist[itrk];
      PHSnglCentralTrack    *trk     = trklink->m_central;
      vector<SvxClsLink>    &vlink   = trklink->m_vlink;

//      if(trk->get_mom()>0.5&&trk->get_mom()<4.0) {
	float mom_cnt    = trk->get_mom();
	int   charge_cnt = trk->get_charge();
	float the0_cnt   = trk->get_the0();
	float phi0_cnt   = trk->get_phi0();
	//if(phi0_cnt > TMath::Pi()) phi0_cnt -= (2.*TMath::Pi());
	float emcdphi    = trk->get_emcdphi();
	float emcdz      = trk->get_emcdz();
	//float pt_cnt     = mom_cnt*sin(the0_cnt);
	float zed        = trk->get_zed();

//	if(fabs(emcdz-1.4)<15 && fabs(emcdphi+0.00015)<0.02) {
          
          m_zvtx    = vtxout->get_ZVertex();
          m_xvtx    = vtx.getX();
          m_yvtx    = vtx.getY();
          string s_vtx = vtxout->which_Vtx();
          if     (s_vtx=="SVX_PRECISE") m_vtxid=2;
          else if(s_vtx=="SVX")         m_vtxid=3;
          else if(s_vtx=="BBC")         m_vtxid=4;
          else                          m_vtxid=0;
          

          m_mom     = mom_cnt*charge_cnt;
          m_the0    = the0_cnt;
          m_phi0    = phi0_cnt;
          m_zed     = zed;
          m_dcqual  = trk->get_quality();
          m_emcdphi = emcdphi;
          m_emcdz   = emcdz;

          m_tofdphi = trk->get_tofdphi();
          m_tofdz   = trk->get_tofdz();
//          m_tofwdphi= trk->get_tofwdphi();
//          m_tofwdz  = trk->get_tofwdz();

          m_temc = trk->get_temc();
          m_ttof = trk->get_ttof();
          //m_ttofw = trk->get_ttofw();
          m_plemc = trk->get_plemc();
          m_pltof = trk->get_pltof();
          //m_pltofw = trk->get_pltofw();
           
          //  m2 =(::pow(t*c/L,2.0)-1)*p*p;
          static const float C =  29.9792458; // cm/n
          m_m2emc = (m_plemc>0) ? ((TMath::Power(m_temc*C/m_plemc, 2.0) -1.0)*m_mom*m_mom) : -9999.;
          m_m2tof = (m_plemc>0) ? ((TMath::Power(m_ttof*C/m_pltof, 2.0) -1.0)*m_mom*m_mom) : -9999.;
          //m_m2tofw = trk->get_m2tofw();

          m_n0    = trk->get_n0();
          m_npe0  = trk->get_npe0();
          m_cch2  = trk->get_chi2();
          m_disp  = trk->get_disp();
          m_sn0   = trk->get_sn0();
          m_snpe0 = trk->get_snpe0();
          m_scch2 = trk->get_schi2();
          m_sdisp = trk->get_sdisp();
          m_ecore = trk->get_ecore();
          m_ecorr = 0; //trk->get_ecorr();

          int nlink = vlink.size();    

          if(mctrk!=NULL){
            float simpx = mctrk->get_momentumx(itrk);
            float simpy = mctrk->get_momentumy(itrk);
            
            m_simpt   =  sqrt(simpx*simpx+simpy*simpy);
            m_simphi0 =  atan2(simpy, simpx);
            m_simvx =  mctrk->get_vertexx(itrk);
            m_simvy =  mctrk->get_vertexy(itrk);
            m_simvz =  mctrk->get_vertexz(itrk);
            m_simdca = calcSimD2DCA(m_simpt, charge_cnt, m_simphi0, 
                                    m_simvx, m_simvy, vtx.getX(),  vtx.getY());
          } else {
            m_simdca=-999.0;
          }
        
          // initial value for svx value
          m_lid    = 0;
          m_chi2   = -9999.0;
          m_ndf    = 0;
          m_nhit   = 0;
          m_d2dca  = -9999.0;
          m_dL     = -9999.0;
          m_dbend  = -9999.0;
          m_d2dca0 = -9999.0;
          m_d2dca1 = -9999.0;
          m_d2dca2 = -9999.0;
          m_dedx1 = 0.0;
          m_dedx2 = 0.0;
          for(int ient=0; ient<36; ient++){ *(m_datap[ient]) = -999.0; }; // initial value for each association
          
          //check link
          for(int ilink=0;ilink<nlink;ilink++) {
            // init
            for(int ient=0; ient<36; ient++){ *(m_datap[ient]) = -999.0; }; // initial value for each association
            m_dedx1 = 0.0;
            m_dedx2 = 0.0;

            // fill
            SvxClsLink&             link = vlink[ilink];
            vector<SvxClsLinkNode>& nodelink = link.m_nodelink;

            m_best   = link.m_bestlink ? 1 : 0;

            //if(m_best==1){
              m_lid    = ilink;
              m_chi2   = link.m_chi2;
              m_ndf    = link.m_ndf;
              m_nhit   = link.getNHit();
              m_d2dca  = link.m_d2dca;
              m_dL     = link.m_dL;
              m_dbend  = link.m_dbend;
              m_d2dca0 = link.m_d2dca0;
              m_d2dca1 = link.m_d2dca1;
              m_d2dca2 = link.m_d2dca2;

              //cout<<itrk<<" "<<ilink<<" "<<nodelink.size()<<endl;
              int nnode = nodelink.size();

              //float xyz[8] = {-999,-999,-999,-999,-999,-999,-999,-999};
              for(int inode=0;inode<nnode;inode++) {
                SvxClsLinkNode& node = nodelink[inode];
                bool found = node.found;
                int layer=-1;
                if     (node.sublayer==0){ layer=0; }
                else if(node.sublayer==1){ layer=1; }
                else if(2<=node.sublayer&&node.sublayer<5){ layer=2; }
                else                     { layer=3; }

                if(found){
                  SvxCluster *cls = node.cluster;
                  float xsvx = cls->get_xyz_global(0);
                  float ysvx = cls->get_xyz_global(1);
                  float zsvx = cls->get_xyz_global(2);
                  float phiv = atan2(ysvx, xsvx);
                  float lt   = sqrt((xsvx-vtx.getX())*(xsvx-vtx.getX())+(ysvx-vtx.getY())*(ysvx-vtx.getY()));

                  *(m_datap[9*layer  ]) = (float) node.sublayer;
                  *(m_datap[9*layer+1]) = (float) cls->get_ladder();
                  *(m_datap[9*layer+2]) = node.zproj;
                  *(m_datap[9*layer+3]) = node.dproj;
                  *(m_datap[9*layer+4]) = node.mag_bend;
                  *(m_datap[9*layer+5]) = zsvx;
                  *(m_datap[9*layer+6]) = phiv;
                  *(m_datap[9*layer+7]) = lt;
                  *(m_datap[9*layer+8]) = (float) cls->get_size();

                   //xyz[2*layer  ] = xsvx;
                   //xyz[2*layer+1] = ysvx;

                  // cout<<"sublayer : "<<(*(m_datap[7*layer]))<<" "<<m_ly3<<" "<<zsvx<<endl;

                  if(layer==2) m_dedx1 += (cls->get_adc(0)+cls->get_adc(1));
                  if(layer==3) m_dedx2 += (cls->get_adc(0)+cls->get_adc(1));
                }
              } // inode loop


/*
              if(m_best==1){
                cout<<"trk : ";
                cout<<m_mom*sin(m_the0)<<" ";
                cout<<m_phi0<<", ";
                cout<<vtx.getX()<<", "<<vtx.getY()<<", ";
                for(int ixyz=0; ixyz<8; ixyz++){
                  cout<<xyz[ixyz]<<", ";
                }
                cout<<" "<<m_d2dca0<<", "<<m_d2dca1<<", "<<m_d2dca2<<endl;
              }
*/


              m_ntp_cntclus->Fill();
              //break;
            //} // m_best==1
          }//for(ilink)

          m_ntp_cnttrk->Fill();

//	}//if(|emcdz|<15 && |emcdphi|<0.02)
//      }//if(0.5<mom_cnt<4)
    }
  }

  m_eventNumber++;
  return 0;
}

//==============================================================

int svxcentclsana::End(PHCompositeNode *topNode) {
  cout << "svxcentclsana::End:  Writing out..." << endl;
  m_outFile->cd();
  m_outFile->Write();
  cout << "svxcentclsana::End:  Closing output file..." << endl;
  m_outFile->Close();
  delete m_outFile;
  m_outFile=NULL;
  return 0;
}


float svxcentclsana::calcSimD2DCA(float pt, float charge, float phi0, float hx, float hy, float vx, float vy){
  static const float B = 0.92;
  static const float b = 0.003*B;

  float R = pt/b;
  float cx = hx + charge*R*sin(phi0);
  float cy = hy - charge*R*cos(phi0);

  float L = sqrt((vx-cx)*(vx-cx) + (vy-cy)*(vy-cy));
  //float psi = atan2((vy-cy), (vx-cx));

  // DCA point
  //float xx = cx + R*cos(psi);
  //float xy = cy + R*sin(psi);

  // DCA value
  float d2dca = R - L;
  return d2dca;
}


