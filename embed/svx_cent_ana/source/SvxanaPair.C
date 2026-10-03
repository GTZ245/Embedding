#include <iostream>
#include <string>

#include "SvxanaPair.h"

#include "phool.h"
#include "PHPoint.h"
#include "PHTypedNodeIterator.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "Fun4AllReturnCodes.h"


#include "RunHeader.h"
#include "EventHeader.h"
#include "PHGlobal.h"
#include "VtxOut.h"
#include "PHCentralTrack.h"
#include "PHSnglCentralTrack.h"

#include "SvxCentralTrackRecalList.h"
#include "SvxCentralTrackRecal.h"

#include "McEvalSingleList.h"


/*
#include "RpConst.h"
#include "RpSnglSumXY.h"
#include "RpSumXYObject.h"

*/

#include "TFile.h"
#include "TTree.h"


#include "getClass.h"

using namespace std;
using namespace findNode;

//==============================================================
static const float M_e  = 0.000511; // GeV/c^2
static const float M_pi = 0.1395;   // GeV/c^2
static const float M_K  = 0.4936;   // GeV/c^2

static const int NBUF  = 20;
static const int NZ    = 16;
static const int NCENT = 18;

static int           idxbuf[NZ][NCENT];
static vector<AnaTrk*> vepbuf[NBUF][NZ][NCENT];
static vector<AnaTrk*> vembuf[NBUF][NZ][NCENT];


class AnaTrk { // only dcqual==31 or 63 + n0>2 + E/p>0.7
  public:
    AnaTrk() : mom(-9999.), phi0(-9999.), the0(-9999.), c(-9999.),
             dcqual(-999),
             emcdphi(-9999.), emcdz(-9999.), 
             n0(-9999), cch2npe0(-9999.), disp(-9999.),
             ecore(-9999.), ep(-9999.),
             svxchi2ndf(-9999.), 
             nhit(-9999), 
             d2dca(-9999.), zdca(-9999.)
      {
        svxmom=0.0; svxphi0=0.0; svxthe0=0.0; svxc=0.0;
        dcapos[0]=0.0; dcapos[1]=0.0; dcapos[2]=0.0;
        for(int i=0; i<4; i++){
          svxdphi[i] = -9999.;
          svxdz[i]   = -9999.;
        }
        convtag=0;
        cnt    = NULL;
        svxcnt = NULL;
        evalid = -1;

        R   = 0.0;
        phi = 0.0;
        cx  = 0.0;
        cy  = 0.0;

        simpaid = -1;
        simpx=0; simpy=0; simpz=0;
        simvx=0; simvy=0; simvz=0;
        simvr=0;
        simptpr=0;
      }
    virtual ~AnaTrk(){}

    void fillTrackInfo(PHSnglCentralTrack *trk, SvxCentralTrackRecal* svxtrk, float m_fieldScale=1.0);

    void calc_circleinfo(float m_fieldScale);

    float pt() { return mom*sin(the0); }
    float px() { return pt()*cos(phi0); }
    float py() { return pt()*sin(phi0); }
    float pz() { return mom*cos(the0); }

    float svxpt() { return svxmom*sin(svxthe0); }
    float svxpx() { return svxpt()*cos(svxphi0); }
    float svxpy() { return svxpt()*sin(svxphi0); }
    float svxpz() { return svxmom*cos(svxthe0); }

    void print(){
      cout<<"AnaTrk : "<<mom<<" "<<phi0<<" "<<the0<<" "<<c<<" : ";
      cout<<nhit<<" ";
      for(int i=0; i<4; i++){
        cout<<"("<<svxdphi[i]<<" "<<svxdz[i]<<") ";
      }
      cout<<d2dca<<endl;
    }

  public:
    float mom, phi0, the0, c;
    int   dcqual;
    float emcdphi, emcdz;
    int   n0;
    float cch2npe0, disp;
    float ecore, ep;
 
    float dcapos[3];
    float svxmom, svxphi0, svxthe0, svxc;
    float svxdphi[4], svxdz[4];
    float svxchi2ndf;
    int   nhit;
    float d2dca, zdca;

    int convtag;

    PHSnglCentralTrack   *cnt;
    SvxCentralTrackRecal *svxcnt;
    int evalid;


    int   simpaid;
    float simpx, simpy, simpz;
    float simvx, simvy, simvz;
    float simvr;
    float simptpr;

    // circle info
    float R;
    float phi;
    float cx;
    float cy;
};


void AnaTrk::fillTrackInfo(PHSnglCentralTrack *trk, SvxCentralTrackRecal *svxtrk, float m_fieldScale){

  if(trk!=NULL){
    mom      = trk->get_mom();
    phi0     = trk->get_phi0()     ;
    the0     = trk->get_the0()     ;
    c        = trk->get_charge()   ;
    dcqual   = trk->get_quality()   ;
    emcdphi  = trk->get_emcdphi()  ;
    emcdz    = trk->get_emcdz()    ;
    n0       = trk->get_n0()       ;
    cch2npe0 = trk->get_chi2()/trk->get_npe0();
    disp     = trk->get_disp()     ;
    ecore    = trk->get_ecore()    ;
    ep       = trk->get_ecore()/trk->get_mom();
    cnt      = trk;
  }

  if(svxtrk!=NULL){
    float dphi[4]={-9999., -9999., -9999., -9999.}, dz[4]={-9999., -9999., -9999., -9999.};

    for(int ilay=0; ilay<4; ilay++){
      int nhit = svxtrk->get_nhit(ilay);
      if(nhit>0){
        dphi[ilay] = svxtrk->get_ClusterDPhi(ilay, 0);
        dz[ilay]   = svxtrk->get_ClusterDZ(ilay, 0);
      }
    }

    for(int i=0; i<4; i++){
      svxdphi[i] = dphi[i];
      svxdz[i]   = dz[i];
    }
    int nhitall = svxtrk->get_nhit();
    int ndf = (nhitall==0) ? 0 : (nhitall-1)*2;
    svxchi2ndf = svxtrk->get_Chisquare()/ndf;
    nhit       = nhitall;
    d2dca      = svxtrk->get_DCA2D();
    zdca       = svxtrk->get_DCAZ();

    svxmom     = svxtrk->get_mom();
    svxphi0    = svxtrk->get_phi0()     ;
    svxthe0    = svxtrk->get_the0()     ;
    svxc       = svxtrk->get_charge()   ;

    dcapos[0]  = svxtrk->get_ClosestApproach(0);
    dcapos[1]  = svxtrk->get_ClosestApproach(1);
    dcapos[2]  = svxtrk->get_ClosestApproach(2);

    calc_circleinfo(m_fieldScale);

    svxcnt     = svxtrk;
  }
}

//------------------------------------------

void SvxanaPair::init_idxbuf(){
 for(int iz=0; iz<NZ; iz++){
   for(int icent=0; icent<NCENT; icent++){
     idxbuf[iz][icent] = 0;
   }
 }
}

void SvxanaPair::calc_iz_icent(float zvtx, float bbcq, int& iz, int& icent){
  iz = (int)((zvtx+8.0)/(1.0)); // +-8cm 1bin=1.cm 16bin max
  if(iz<=0)  iz=0;
  if(16<=iz) iz=15;

  static const float bbcqrng[19] = {0, 13.5, 22, 34, 50, 
                                    72, 102, 141, 192, 254, 
                                   330, 421, 527, 646, 793, 
                                   960, 160, 1400, 2000};
  for(int i=0; i<NCENT; i++){
    if(bbcqrng[i]<=bbcq&&bbcq<bbcqrng[i+1]) { icent = i; break; }
  }

  if(m_simmode>0|| m_ppmode>0){
    //iz = 0; // this is a bug
    icent=0;
  }

//  cout<<"iz : icent = "<<zvtx<<" "<<iz<<" : "<<bbcq<<" : "<<icent<<endl;
}

float SvxanaPair::calc_pair(AnaTrk *ep, AnaTrk *em, float Mep, float Mem,
                float& Mee,  float& px,   float& py,   float& pz, float& pt,
                float& thev, float& ptep, float& ptem, float& phiv)
{
  return calc_pair(ep->px(), ep->py(), ep->pz(), Mep,
                   em->px(), em->py(), em->pz(), Mem,
                   Mee,  px,  py,  pz, pt,
                   thev, ptep, ptem, phiv);
}

float SvxanaPair::calc_pair(float pxp, float pyp, float pzp, float Mep,
                            float pxm, float pym, float pzm, float Mem,
                float& Mee,  float& px,   float& py,   float& pz, float& pt,
                float& thev, float& ptep, float& ptem, float& phiv)
{
  float pxep = pxp;
  float pyep = pyp;
  float pzep = pzp;
  float momp = sqrt(pxep*pxep+pyep*pyep+pzep*pzep);
  float Eep  = sqrt(momp*momp+Mep*Mep);

  float pxem = pxm;
  float pyem = pym;
  float pzem = pzm;
  float momm = sqrt(pxem*pxem+pyem*pyem+pzem*pzem);
  float Eem  = sqrt(momm*momm+Mem*Mem);

  px = pxep + pxem;
  py = pyep + pyem;
  pz = pzep + pzem;
  float E  = Eep + Eem;

  pt = sqrt(px*px + py*py);
  Mee = sqrt(E*E - (px*px+py*py+pz*pz));

  //////////////
  //unit vector of (pep+pem)
  float pl = sqrt(px*px+py*py+pz*pz);
  float ux = px/pl;
  float uy = py/pl;
  float uz = pz/pl;

  //axis defined by (ux,uy,ux)X(0,0,1).
  // this is the axis that is perpendicular to the direction of
  // pair, and also perpendicular to the Z axis (field direction).
  // If the pair is conversion at R!=0, it must have (apparent)
  // momentum component in this axis (caused by field intergral from the
  // vertex point to the conversion point).
  // The sign of the component is opposite for e+ and e-.
  //
  // (ux,uy,ux)X(0,0,1)=(uy,-ux,0)
  //
  //cout<<sqrt(ux*ux+uy*uy)<<endl;
  float ax =  uy/sqrt(ux*ux+uy*uy);
  float ay = -ux/sqrt(ux*ux+uy*uy);

  //momentum of e+ and e- in (ax,ay,az) axis. Note that az=0 by
  //definition.
  ptep = pxep*ax + pyep*ay;
  ptem = pxem*ax + pyem*ay;

  //vector product of pep X pem
  float vpx = pyep*pzem - pzep*pyem;
  float vpy = pzep*pxem - pxep*pzem;
  float vpz = pxep*pyem - pyep*pxem;
  float vp = sqrt(vpx*vpx+vpy*vpy+vpz*vpz);
  thev = acos(vpz/vp);

  //unit vector of pep X pem
  float vx = vpx/vp;
  float vy = vpy/vp;
  float vz = vpz/vp;

  //The third axis defined by vector product (ux,uy,uz)X(vx,vy,vz)
  float wx = uy*vz - uz*vy;
  float wy = uz*vx - ux*vz;
  float wz = ux*vy - uy*vx;
  float wl = sqrt(wx*wx+wy*wy+wz*wz);
  // by construction, (wx,wy,wz) must be a unit vector.
  if(fabs(wl - 1.0) > 0.00001) cout << "Calculation error in W vector"<<endl;

  // measure angle between (wx,wy,wz) and (ax,ay,0). The angle between them
  // should be small if the pair is conversion
  //
  float cosPhiV = wx*ax + wy*ay;
  phiv = acos(cosPhiV);

  //


  return Mee;
}


/*
struct circle_info {
  float pt;
  float R;
  float phi;
  float cx;
  float cy;
};
*/

void AnaTrk::calc_circleinfo(float m_fieldScale){
  static const float B = 0.9;
  static const float b = 0.003*B;

  float b_sign = (m_fieldScale>0) ? 1.0 : -1.0;

  R   = svxpt()/b;
  phi = atan2(svxpy(), svxpx());
  cx  = dcapos[0] + (b_sign)*c*R*(-sin(phi));
  cy  = dcapos[1] + (b_sign)*c*R*( cos(phi));

  //cout<<"AnaTrk : "<<R<<" "<<phi<<" "<<cx<<" "<<cy<<endl;
}

int SvxanaPair::calc_crosspos(AnaTrk *ep, AnaTrk *em,
                               int& ncross, vector<float>& v_crossR, 
                   vector<float>& v_crossX, vector<float>& v_crossY)
{

//  circle_info cp, cm;
//  calc_circleinfo(pxp, pyp, chargep, vxp, vyp, cp, m_fieldScale);
//  calc_circleinfo(pxm, pym, chargem, vxm, vym, cm, m_fieldScale);

  float dx = em->cx - ep->cx;
  float dy = em->cy - ep->cy;
  float L  = sqrt(dx*dx + dy*dy);

  float Rdiff = fabs(ep->R - em->R);
  float Radd  = ep->R+em->R;
  if(L<Rdiff) {
    //cout<<" No crossing point"<<endl;
    ncross = 0;
    return ncross;
  }


  // if 2 circle are crossing  (or touching)
  if(Rdiff<=L && L<=Radd) {
    // law of cosine
    //         R2
    //     /----------
    //    /      _--
    // R1/    _--      
    //  /a _--          
    // /--    L
    //     
    
    float phi      = atan2(dy, dx);
    float cosAlpha = (L*L + ep->R*ep->R - em->R*em->R)/(2*L*ep->R);
    float alpha    = acos(cosAlpha);
    
    //cout<<"calc_cross"<<endl;
    //cout<<ep->cx<<" "<<ep->cy<<" "<<ep->R<<endl;
    //cout<<em->cx<<" "<<em->cy<<" "<<em->R<<endl;
    //cout<<phi<<" "<<cosAlpha<<" "<<alpha<<endl;
    
    float cross_x[2], cross_y[2], cross_r[2];
    for(int i=0; i<2; i++){
      float sign = (i==0) ? 1:-1;
      float x_x = ep->cx + ep->R*cos(phi+sign*alpha);
      float x_y = ep->cy + ep->R*sin(phi+sign*alpha);
    
      cross_x[i] = x_x;
      cross_y[i] = x_y;
      cross_r[i] = sqrt(x_x*x_x + x_y*x_y);
      //cout<<" cross : "<<i<<" "<<x_x<<" "<<x_y<<" "<<x_r<<endl;
    }
    
    // sort if x_r1<x_r0
    if(cross_r[1]<cross_r[0]){
      float tmp_x = cross_x[0];
      float tmp_y = cross_y[0];
      float tmp_r = cross_r[0];
      cross_x[0] = cross_x[1];
      cross_y[0] = cross_y[1];
      cross_r[0] = cross_r[1];
      cross_x[1] = tmp_x;
      cross_y[1] = tmp_y;
      cross_r[1] = tmp_r;
    }
    
    for(int i=0; i<2; i++){
      v_crossX.push_back(cross_x[i]);
      v_crossY.push_back(cross_y[i]);
      v_crossR.push_back(cross_r[i]);
    }
    
    ncross = (int)v_crossR.size();
  }
  else if(Radd<L) { // 2 circle are separated
    float vx = em->cx - ep->cx;
    float vy = em->cy - ep->cy;
    float phi = atan2(vy, vx);
    float d   = 0.5*(L -Radd);

    float x_x = ep->cx + d*cos(phi);
    float x_y = ep->cy + d*sin(phi);

    v_crossX.push_back(x_x);
    v_crossY.push_back(x_y);
    v_crossR.push_back(sqrt(x_x*x_x + x_y*x_y));

    ncross = 1;
  }
  else {
    ncross = 0;
  }

  return ncross;
  
}

void SvxanaPair::calc_vector_at_cross(
  float px, float py, float pz, 
  float dcax, float dcay, float dcaz,
  float charge, float R,
  float cross_x, float cross_y,
  float* cpx, float* cpy, float *cpz, float *cdphi
)
{

  float lx = cross_x - dcax;
  float ly = cross_y - dcay;
  float L  = sqrt(lx*lx + ly*ly);

  float sign = (charge>0) ? -1 : 1;
  float dphi = sign*(L/R); // assume angle is small = phi = L(length of dca to cross)/R(bending radius)

  float npx = px*cos(dphi) - py*sin(dphi);
  float npy = px*sin(dphi) + py*cos(dphi);
  float npz = pz;

  // fill
  *cpx = npx;
  *cpy = npy;
  *cpz = npz;
  *cdphi = dphi;
}

void SvxanaPair::calc_dca_with_xvector(
  float xpx, float xpy, float xpz, float xpt,
  float crossx, float crossy,
  float xvtx, float yvtx, float zvtx,
  float* xdca2d, float* xlength, float* lxy
){
  // at this moment, consider 2d  

  float u_xpx = xpx/xpt;
  float u_xpy = xpy/xpt;

  // vector Cross to PrimVertex
  float v_x = xvtx - crossx;
  float v_y = yvtx - crossy;

  // (u x v)
  float v_cross_u = u_xpx*v_y - u_xpy*v_x;
 
  *xdca2d = v_cross_u;

  // vector L which is from prim_vtx to crossing point
  float lx  = crossx - xvtx;
  float ly  = crossy - yvtx;

  float dlen    = sqrt(lx*lx  + ly*ly);
  float l_dot_u = lx*u_xpx + ly*u_xpy;
  float dsign   = (l_dot_u>0) ? 1 : -1;

  //float cosTheta = (lx*xpx + ly*xpy)/(dlen*xpt);
  //float dsign = (cosTheta>=0) ? 1 : -1; // if vector L and its daugher x goes to same dir, sign is 1
  *xlength = dsign * dlen;

  // scalar product of L.u u(u_xpx, u_xpy) is unit vector of =(xpx, xpy)
  *lxy = l_dot_u;
}




//==============================================================

SvxanaPair::SvxanaPair(string filename) : m_OutputFileName(filename)
{
  ThisName = "SvxanaPair";
  init_ana=0;
  EventNumber=0;

  m_fieldScale =  1.0; 

  m_simmode = 0;
  m_ppmode  = 0;

  init_idxbuf();
}

//==============================================================

SvxanaPair::~SvxanaPair() {
}

//==============================================================

int SvxanaPair::Init(PHCompositeNode *topNode) {

  cout << "SvxanaPair::Init started..." << endl;
  m_OutputNtupleFile = new TFile(m_OutputFileName.c_str(),"RECREATE");
  cout << "SvxanaPair::Init: output file " << m_OutputFileName << " opened." << endl;


  initPairTree();

  cout << "SvxanaPair::Init ended." << endl;
  return 0;
}

//==============================================================
  
int SvxanaPair::InitRun(PHCompositeNode *topNode) {
  cout << "SvxanaPair::InitRun started..." << endl;

  // check magnet current 
  RunHeader* runheader = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if (runheader==NULL) {
    cout << PHWHERE<< "Can't find runheader. " << endl;
    return ABORTRUN;
  }
  if(runheader->get_currentCentral()>0){ m_fieldScale =  1.0; } 
  else                                 { m_fieldScale = -1.0; }
  cout<<"SvxanaPair::InitRun  fieldScale="<<m_fieldScale<<endl;


  cout << "SvxanaPair::InitRun ended." << endl;
  return 0;
}

//==============================================================


int SvxanaPair::process_event(PHCompositeNode *topNode) {

  RunHeader   *run    = getClass<RunHeader>(topNode,"RunHeader");
  EventHeader *evthdr = getClass<EventHeader>(topNode,"EventHeader");
  VtxOut      *vtxout = getClass<VtxOut>(topNode,"VtxOut");
  PHGlobal    *global = getClass<PHGlobal>(topNode,"PHGlobal");
  //TrigLvl1    *trglvl1 = getClass<TrigLvl1>(topNode,"TrigLvl1");
  //BbcOut *bbc    = getClass<BbcOut>(topNode,"BbcOut");

  int RunNumber = (run!=NULL)    ? run->get_RunNumber() : -9999;
  int EventNum  = (evthdr!=NULL) ? evthdr->get_EvtSequence() : -9999;

  if(EventNumber==0) {
    cout << "run = "<<RunNumber<<" event = "<<EventNum<<endl;
  }


  PHCentralTrack *cntlist = getClass<PHCentralTrack>(topNode,"PHCentralTrack");
  int ntrk = (cntlist!=NULL) ? cntlist->get_npart() : 0;

/*
  SvxCentralTrackMap *svxcntmap  = findNode::getClass<SvxCentralTrackMap>(topNode, "SvxCentralTrack_comp");
  if(svxcntmap==NULL){
    cout<<"Error SvxanaPair::process_event no SvxCentralTrackMap"<<endl;
    return EVENT_OK;
  }
*/
  SvxCentralTrackRecalList *svxcntlist  = findNode::getClass<SvxCentralTrackRecalList>(topNode, "SvxCentralTrackRecalList");
  if(svxcntlist==NULL){
    cout<<"Error SvxanaPair::process_event no SvxCentralTrackRecalList"<<endl;
    return EVENT_OK;
  }

  McEvalSingleList *mctrklist  = findNode::getClass<McEvalSingleList>(topNode, "McSingle");
  if(mctrklist==NULL){
    cout<<"Error SvxanaPair::process_event no McEvalSingleList"<<endl;
    return EVENT_OK;
  }


  m_evt      = EventNumber;
  float cent = (global!=NULL) ? global->getCentrality() : -999.;
  m_bbcq     = (global!=NULL) ? (global->getBbcChargeN()+global->getBbcChargeS()) : -999.0;

  PHPoint vtxpos(-999., -999., -999.); 
  if(vtxout!=NULL)  {
    if(m_simmode==0){
      vtxpos = vtxout->get_Vertex();
    } 
    else {
      vtxpos = vtxout->get_Vertex("SIM");
    }
  }
  m_xvtx = vtxpos.getX();
  m_yvtx = vtxpos.getY();
  m_zvtx = vtxpos.getZ();

  PHPoint vsim = vtxout->get_Vertex("SIM");
  PHPoint vbbc = vtxout->get_Vertex("BBC");
  cout<<"Vertex : "<<vtxout->which_Vtx()<<endl;
  cout<<"  SIM  "<<vsim.getX()<<" "<<vsim.getY()<<" "<<vsim.getZ()<<endl;
  cout<<"  BBC  "<<vbbc.getX()<<" "<<vbbc.getY()<<" "<<vbbc.getZ()<<endl;
  cout<<"  Use  "<<m_xvtx<<" "<<m_yvtx<<" "<<m_zvtx<<endl;

  if(fabs(m_zvtx)>8) return EVENT_OK; // zvtx<8cm is OK, if not, skip this event


  // track ana
  analyze_track(cntlist, svxcntlist, mctrklist);
  
  // pair ana
  analyze_pair();

  if(init_ana==0) {
    cout<<"init ana"<<endl;
    topNode->print(); 
    init_ana=1;
 
    cout<<endl<<endl;
    cout<<"init ana ends"<<endl;
  }

  if(verbosity>1){
    if(EventNumber%1000==0) {
      cout << "------- Event # " << EventNumber 
           << " zvtx= "  << m_zvtx
           << " ntrack= " << ntrk 
           << " centrality= " << cent
           << " bbcq= " << m_bbcq
           << " vtx-x-y-z= "<<m_xvtx<<" "<<m_yvtx<<" "<<m_zvtx
           << endl;
    }
  }

  EventNumber++;
  return 0;
}

void SvxanaPair::analyze_track(PHCentralTrack* cntlist, SvxCentralTrackRecalList *svxcntlist, McEvalSingleList *mctrklist){
  if(cntlist==NULL||svxcntlist==NULL){
    if(cntlist==NULL)    cout<<"PHCentralTrack object is NULL"<<endl;
    if(svxcntlist==NULL) cout<<"SvxCentralTrackRecalList object is NULL"<<endl;
    return;
  }

  // clear electron vector
  //clear_VAnaTrk(m_vep);
  //clear_VAnaTrk(m_vem);
  m_vep.clear();
  m_vem.clear();

  m_vhp.clear();
  m_vhm.clear();


  int ntrk = cntlist->get_npart();

  // initialize Map 
  SvxCentralTrackRecal** svxCntAry = new SvxCentralTrackRecal*[ntrk];
  for(int itrk=0; itrk<ntrk; itrk++){ svxCntAry[itrk] = NULL; }

  int nsvxcnt = svxcntlist->GetNentries();
  for(int itrk=0; itrk<nsvxcnt; itrk++){ 
    SvxCentralTrackRecal *svxcnt = svxcntlist->GetSvxCNTRecal(itrk);
    int  cntidx = svxcnt->get_CNTID();
    if(cntidx<0||ntrk<=cntidx){ cerr<<" DchIndex is out of range"<<endl; }
    svxCntAry[cntidx] = svxcnt; 
  }

  // --- simmode ---
  int* mcTrkIdx = new int[ntrk];
  for(int itrk=0; itrk<ntrk; itrk++){ mcTrkIdx[itrk] = -999; }
  if(m_simmode>0){
    for(unsigned int itrk=0; itrk<mctrklist->get_McEvalSingleTrackN(); itrk++){
  
      int nreco = mctrklist->get_Nreco(itrk);
      for(int ireco=0; ireco<nreco; ireco++){
        int recoid  = mctrklist->get_recoid(itrk, ireco);
        if(0<=recoid&&recoid<ntrk){ mcTrkIdx[recoid] = itrk; }
        else { cout<<"Error :: out of range : "<<recoid<<" "<<itrk<<" "<<ireco<<endl; }
      }
    }
  }



  // CNT track loop
  for(int itrk=0; itrk<ntrk; itrk++){
    PHSnglCentralTrack *trk = cntlist->get_track(itrk);

    // overall momentum cut
    if(trk->get_mom()*sin(trk->get_the0())<=0.3) continue;

    ////////////////////////////////
    // goodTrack & electron ID
    if( isGoodTrack(trk) && 
        isElectron(trk))
      { 
        AnaTrk *etrk = new AnaTrk();
        fill_AnaTrk(etrk, trk, svxCntAry[itrk]);
        if(m_simmode>0) fill_AnaTrk_sim(etrk, mcTrkIdx[itrk], mctrklist);

        //etrk->print();
   
        if(etrk->c>0){ m_vep.push_back(etrk);}
        else         { m_vem.push_back(etrk); }
      }

    ////////////////////////////////
    // goodTrack for hadron
    if( isGoodTrack(trk) )
      {
        AnaTrk *etrk = new AnaTrk();
        fill_AnaTrk(etrk, trk, svxCntAry[itrk]);
        if(m_simmode>0) fill_AnaTrk_sim(etrk, mcTrkIdx[itrk], mctrklist);

        if(etrk->c>0){ m_vhp.push_back(etrk);}
        else         { m_vhm.push_back(etrk); }
      }

  }



  // delete
  delete [] svxCntAry;
  delete [] mcTrkIdx;
}

void SvxanaPair::analyze_pair(){

  ///////////////////////////////////
  // foreground
  fill_pair(m_vep, m_vem, M_e, M_e, m_eepair);

  fill_single(m_vep);
  fill_single(m_vem);

  fill_pair(m_vhp, m_vhm, M_pi, M_pi, m_pipipair);
  fill_pair(m_vhp, m_vhm, M_pi, M_K,  m_pikpair);
  fill_pair(m_vhp, m_vhm, M_K,  M_pi, m_pikpair);
  fill_pair(m_vhp, m_vhm, M_K,  M_K,  m_kkpair);

  ///////////////////////////////////
  // event mixing
  int iz=0, icent=0;
  calc_iz_icent(m_zvtx, m_bbcq, iz, icent);

  for(int imix=0; imix<NBUF; imix++){
    vector<AnaTrk*>& bg_vep = vepbuf[imix][iz][icent];
    vector<AnaTrk*>& bg_vem = vembuf[imix][iz][icent];

    fill_pair(m_vep,  bg_vem, M_e, M_e, m_eepair_bg);  // +- mixing
    fill_pair(bg_vep, m_vem,  M_e, M_e, m_eepair_bg);  // -+ mixing
  }

  ///////////////////////////////////
  //  update the buffer for the event mixing
 
  // copy the current list into the buf
  int ibuf = idxbuf[iz][icent];
  // clear the buffer 
  clear_VAnaTrk(vepbuf[ibuf][iz][icent]);
  clear_VAnaTrk(vembuf[ibuf][iz][icent]);

  // copy 
  vepbuf[ibuf][iz][icent] = m_vep;
  vembuf[ibuf][iz][icent] = m_vem;

  idxbuf[iz][icent]++;
  if(idxbuf[iz][icent]>=NBUF) {
    idxbuf[iz][icent] = 0;
    //cout<<"idxbuf cleared : "<<iz<<" "<<icent<<endl;
  }

}

bool SvxanaPair::isGoodTrack(PHSnglCentralTrack* trk){
  return( (trk->get_quality()==31||trk->get_quality()==63) &&
          fabs(trk->get_emcdphi())<0.03 &&
          fabs(trk->get_emcdz())<12) ;
}

bool SvxanaPair::isElectron(PHSnglCentralTrack *trk){
  return isElectron(trk->get_n0(), 
                    trk->get_chi2()/trk->get_npe0(), 
                    trk->get_disp(), 
                    trk->get_ecore()/trk->get_mom(), 
                    trk->get_prob());
}

bool SvxanaPair::isElectron(int n0, float ch2npe0, float disp, float Ep, float prob){
  return (n0>=3&&Ep>0.5);
}

void SvxanaPair::fill_AnaTrk(AnaTrk* etrk, PHSnglCentralTrack* trk, SvxCentralTrackRecal* svxtrk)
{
  if(etrk==NULL) return;

  etrk->fillTrackInfo(trk, svxtrk, m_fieldScale);
}

void SvxanaPair::fill_AnaTrk_sim(AnaTrk* etrk, int itrk, McEvalSingleList *mctrklist)
{
  if(etrk==NULL) return;

  if(itrk<0) return; //should be itrk>=0

  if(m_simmode>0){
    etrk->simpaid = mctrklist->get_parentid(itrk);
    //cout<<"simpaid : "<<simpidpa<<endl;
    etrk->simvx =  mctrklist->get_vertexx(itrk);
    etrk->simvy =  mctrklist->get_vertexy(itrk);
    etrk->simvz =  mctrklist->get_vertexz(itrk);

    etrk->simpx = mctrklist->get_momentumx(itrk);
    etrk->simpy = mctrklist->get_momentumy(itrk);
    etrk->simpz = mctrklist->get_momentumz(itrk);

    float pxpr = mctrklist->get_primarymomentumx(itrk);
    float pypr = mctrklist->get_primarymomentumy(itrk);

    float& vx = etrk->simvx;
    float& vy = etrk->simvy;
    etrk->simvr = sqrt(vx*vx+vy*vy);

    etrk->simptpr = sqrt(pxpr*pxpr + pypr*pypr);

    etrk->evalid  = itrk;
  }
}

void SvxanaPair::fill_pair(vector<AnaTrk*>& vep, vector<AnaTrk*>& vem, float Mep, float Mem, TTree *eepair){
  vector<AnaTrk*>::iterator itr_ep, itr_em;
  for(itr_ep=vep.begin(); itr_ep!=vep.end(); ++itr_ep){
    for(itr_em=vem.begin(); itr_em!=vem.end(); ++itr_em){
      AnaTrk *ep = *itr_ep;
      AnaTrk *em = *itr_em;

      fill_pair(ep, em, Mep, Mem, eepair);
    }
  }
}

void SvxanaPair::fill_pair(AnaTrk* ep, AnaTrk *em, float Mep, float Mem, TTree *eepair){
  calc_pair(ep, em, Mep, Mem,
    m_mee, m_px, m_py, m_pz, m_pt, m_thev, m_ptep, m_ptem, m_phiv);

  if(ep->convtag==0){
    if     (0   <=m_mee&&m_mee<0.06) ep->convtag = 1;
    else if(0.06<=m_mee&&m_mee<0.12) ep->convtag = 2;
  }
  if(em->convtag==0){
    if     (0   <=m_mee&&m_mee<0.06) em->convtag = 1;
    else if(0.06<=m_mee&&m_mee<0.12) em->convtag = 2;
  }

  { // using updated mom
    float a_px, a_py, a_pz, a_pt, a_thev, a_ptep, a_ptem;
    calc_pair(ep->svxpx(), ep->svxpy(), ep->svxpz(), Mep,
              em->svxpx(), em->svxpy(), em->svxpz(), Mem,
      m_svxmee, a_px, a_py, a_pz, a_pt, a_thev, a_ptep, a_ptem, m_svxphiv);

    // calc crossing position
    int ncross=0;
    vector<float> v_crossR, v_crossX, v_crossY;
    calc_crosspos(ep, em, ncross, v_crossR, v_crossX, v_crossY);

    m_ncross = ncross;
    m_crossr[0] = m_crossr[1] = -999;
    m_xmee = m_xpx = m_xpy = m_xpz = -999;
    if(ncross>0){
      for(int ic=0; ic<ncross; ic++){ m_crossr[ic] = v_crossR[ic]; }

      // calc mom_vector at crossing point (closest)
      float ep_cpx, ep_cpy, ep_cpz, ep_cdphi;
      calc_vector_at_cross(
            ep->svxpx(), ep->svxpy(), ep->svxpz(), 
            ep->dcapos[0], ep->dcapos[1], ep->dcapos[2],
            ep->c, ep->R,
            v_crossX[0], v_crossY[0],
            &ep_cpx, &ep_cpy, &ep_cpz, &ep_cdphi
          );

      float em_cpx, em_cpy, em_cpz, em_dphi;
      calc_vector_at_cross(
            em->svxpx(), em->svxpy(), em->svxpz(), 
            em->dcapos[0], em->dcapos[1], em->dcapos[2],
            em->c, em->R,
            v_crossX[0], v_crossY[0],
            &em_cpx, &em_cpy, &em_cpz, &em_dphi
          );

      // calc pair using vector at cross
      float x_phiv, x_thev, x_ptep, x_ptem;
      calc_pair(ep_cpx, ep_cpy, ep_cpz, Mep,
                em_cpx, em_cpy, em_cpz, Mem,
             m_xmee, m_xpx, m_xpy, m_xpz, m_xpt, x_thev, x_ptep, x_ptem, x_phiv);

      // calc dca using x_vector
      float xvtx, yvtx, zvtx;
      getPrimVertex(&xvtx, &yvtx, &zvtx);
      calc_dca_with_xvector(
        m_xpx, m_xpy, m_xpz, m_xpt,
        v_crossX[0], v_crossY[0],
        xvtx, yvtx, zvtx,
        &m_xdca2d, &m_xlength, &m_lxy
      );
    }
  }

  // ep
  fill_trkvalue(ep,
     m_momp, m_phi0p, m_the0p, m_n0p,  m_ch2npe0p, m_dispp,
     m_ecorep, m_epp,
     m_emcdpp, m_emcdzp,
     m_nhitp, m_chi2ndfp, m_d2dcap, m_zdcap,
     m_dphip, m_dzp,
     m_simpaidp, m_simvrp, m_simvzp, m_simptprp);

  // em
  fill_trkvalue(em,
     m_momm, m_phi0m, m_the0m, m_n0m,  m_ch2npe0m, m_dispm,
     m_ecorem, m_epm,
     m_emcdpm, m_emcdzm,
     m_nhitm, m_chi2ndfm, m_d2dcam, m_zdcam,
     m_dphim, m_dzm,
     m_simpaidm, m_simvrm, m_simvzm, m_simptprm);

  // sim
  if(m_simmode>0){
    if(fabs(ep->simvx - em->simvx)<0.0001 &&
       fabs(ep->simvy - em->simvy)<0.0001 &&
       fabs(ep->simvz - em->simvz)<0.0001 )
    {
      float t_px, t_py, t_thev, t_ptep, t_ptem, t_phiv;
      calc_pair(ep->simpx, ep->simpy, ep->simpz, Mep,
                em->simpx, em->simpy, em->simpz, Mem,
        m_simmee, t_px, t_py, m_simxpz, m_simxpt, t_thev, t_ptep, t_ptem, t_phiv);

        float xvtx, yvtx, zvtx;
        getPrimVertex(&xvtx, &yvtx, &zvtx);

      calc_dca_with_xvector(
        t_px, t_py, m_simxpz, m_simxpt,
        ep->simvx, ep->simvy,
        xvtx, yvtx, zvtx,
        &m_simxdca2d, &m_simxlength, &m_simlxy
      );
    } else {
      m_simmee     = -999;
      m_simxpz     = -999;
      m_simxpt     = -999;
      m_simxdca2d  = -999;
      m_simxlength = -999;
      m_simlxy     = -999;
    }
  }


  eepair->Fill();
}

void SvxanaPair::fill_single(std::vector<AnaTrk*>& ve)  // include loop
{
  vector<AnaTrk*>::iterator itr;
  for(itr=ve.begin(); itr!=ve.end(); ++itr){
    AnaTrk *etrk = *itr;

    fill_single(etrk);
  }
}

void SvxanaPair::fill_single(AnaTrk *etrk)
{
  m_mom     = etrk->mom;
  m_phi0    = etrk->phi0;
  m_the0    = etrk->the0;
  m_c       = etrk->c;
  m_dcqual  = etrk->dcqual;
  m_emcdphi = etrk->emcdphi;
  m_emcdz   = etrk->emcdz;
  m_n0      = etrk->n0;
  m_ch2npe0 = etrk->cch2npe0;
  m_disp    = etrk->disp;
  m_ecore   = etrk->ecore;
  m_ep      = etrk->ecore/etrk->mom;

  m_nhit    = etrk->nhit;
  m_chi2ndf = etrk->svxchi2ndf;
  m_d2dca   = etrk->d2dca;
  m_zdca    = etrk->zdca;
  m_convtag = etrk->convtag;
  for(int i=0; i<4; i++){
    m_dphi[i] = etrk->svxdphi[i];
    m_dz[i]   = etrk->svxdz[i];
  }
  m_svxmom    = etrk->svxmom;
  m_svxphi0   = etrk->svxphi0;
  m_svxthe0   = etrk->svxthe0;
  m_svxc      = etrk->svxc;

  if(m_simmode>0){
    m_simpaid = etrk->simpaid;
    m_simpx   = etrk->simpx;
    m_simpy   = etrk->simpy;
    m_simpz   = etrk->simpz;
    m_simvr   = etrk->simvr;
    m_simvz   = etrk->simvz;
    m_simptpr = etrk->simptpr;
  }

  m_etree->Fill();
}


void SvxanaPair::fill_trkvalue(AnaTrk *etrk,
     float& mom, float& phi0, float& the0,
     float& n0,  float& ch2npe0, float& disp,
     float& ecore, float& ep,
     float& emcdp, float& emcdz,
     int&   nhit,
     float& chi2ndf, float& d2dca, float& zdca,
     float* dphi_ptr, float* dz_ptr,
     int& simpaid, float& simvr, float& simvz, float& simptpr)
{
  mom     = etrk->mom;
  phi0    = etrk->phi0;
  the0    = etrk->the0;
  n0      = etrk->n0;
  ch2npe0 = etrk->cch2npe0;
  disp    = etrk->disp;
  ecore   = etrk->ecore;
  ep      = ecore/mom;
  emcdp   = etrk->emcdphi;
  emcdz   = etrk->emcdz;

  nhit    = etrk->nhit;
  chi2ndf = etrk->svxchi2ndf;
  d2dca   = etrk->d2dca;
  zdca    = etrk->zdca;
  for(int i=0; i<4; i++){
    dphi_ptr[i] = etrk->svxdphi[i];
    dz_ptr[i]   = etrk->svxdz[i];
  }

  if(m_simmode>0){
    simpaid = etrk->simpaid;
    simvr   = etrk->simvr;
    simvz   = etrk->simvz;
    simptpr = etrk->simptpr;
  }
}



void SvxanaPair::clear_VAnaTrk(vector<AnaTrk*>& vetrk){
  vector<AnaTrk*>::iterator itr;

  for(itr=vetrk.begin(); itr!=vetrk.end(); ++itr){
    AnaTrk* etrk = *itr;
    delete etrk;
  }
  vetrk.clear();
}



int SvxanaPair::End(PHCompositeNode *topNode) {
  cout << "SvxanaPair::End:  Writing out..." << endl;
  m_OutputNtupleFile->Write();
  cout << "SvxanaPair::End:  Closing output file..." << endl;
  m_OutputNtupleFile->Close();
  delete m_OutputNtupleFile;
  m_OutputNtupleFile=0;
  return 0;
}


void SvxanaPair::initPairTree(){
  m_eepair = new TTree("eepair", "E+E- pair tree");
  m_eepair_bg = new TTree("eepair_bg", "E+E- pair tree(mix)");

  m_pipipair = new TTree("pipipair", "Pi+Pi- pair tree");
  m_pikpair  = new TTree("pikpair",  "Pi+K- and Pi-K+ pair tree");
  m_kkpair   = new TTree("kkpair",   "K+K- pair tree");

  TTree *pair[5] = {m_eepair, m_eepair_bg, m_pipipair, m_pikpair, m_kkpair};

  // event
  for(int i=0; i<5; i++){
    // event
    pair[i]->Branch("evt",     &m_evt,     "evt/I");
    pair[i]->Branch("bbcq",    &m_bbcq,    "bbcq/F");
    pair[i]->Branch("xvtx",    &m_xvtx,    "xvtx/F");
    pair[i]->Branch("yvtx",    &m_yvtx,    "yvtx/F");
    pair[i]->Branch("zvtx",    &m_zvtx,    "zvtx/F");
    // pair
    pair[i]->Branch("mee",     &m_mee,   "mee/F");
    pair[i]->Branch("px",      &m_px,    "px/F");
    pair[i]->Branch("py",      &m_py,    "py/F");
    pair[i]->Branch("pz",      &m_pz,    "pz/F");
    pair[i]->Branch("pt",      &m_pt,    "pt/F");
    pair[i]->Branch("thev",    &m_thev,  "thev/F");
    pair[i]->Branch("phiv",    &m_phiv,  "phiv/F");
    pair[i]->Branch("ptep",    &m_ptep,  "ptep/F");
    pair[i]->Branch("ptem",    &m_ptem,  "ptem/F");
    pair[i]->Branch("svxmee",  &m_svxmee, "svxmee/F");
    pair[i]->Branch("svxphiv", &m_svxphiv,"svxphiv/F");
    pair[i]->Branch("ncross",  &m_ncross, "ncross/I");
    pair[i]->Branch("crossr",   m_crossr, "crossr[2]/F");
    pair[i]->Branch("xmee",    &m_xmee,  "xmee/F");
    pair[i]->Branch("xpx",     &m_xpx,   "xpx/F");
    pair[i]->Branch("xpy",     &m_xpy,   "xpy/F");
    pair[i]->Branch("xpz",     &m_xpz,   "xpz/F");
    pair[i]->Branch("xpt",     &m_xpt,   "xpt/F");
    pair[i]->Branch("xdca2d",  &m_xdca2d, "xdca2d/F");
    pair[i]->Branch("xlength", &m_xlength,"xlength/F");
    pair[i]->Branch("lxy",     &m_lxy,    "lxy/F");
    // ep
    pair[i]->Branch("momp",     &m_momp,     "momp/F");
    pair[i]->Branch("phi0p",    &m_phi0p,    "phi0p/F");
    pair[i]->Branch("the0p",    &m_the0p,    "the0p/F");
    pair[i]->Branch("n0p",      &m_n0p,      "n0p/F");
    pair[i]->Branch("ch2npe0p", &m_ch2npe0p, "ch2npe0p/F");
    pair[i]->Branch("dispp",    &m_dispp,    "dispp/F");
    pair[i]->Branch("ecorep",   &m_ecorep,   "ecorep/F");
    pair[i]->Branch("epp",      &m_epp,      "epp/F");
    pair[i]->Branch("emcdpp",   &m_emcdpp,   "emcdpp/F");
    pair[i]->Branch("emcdzp",   &m_emcdzp,   "emcdzp/F");
    pair[i]->Branch("dphip",    m_dphip,     "dphip[4]/F");
    pair[i]->Branch("dzp",      m_dzp,       "dzp[4]/F");
    pair[i]->Branch("chi2ndfp", &m_chi2ndfp, "chi2ndfp/F");
    pair[i]->Branch("nhitp",    &m_nhitp,    "nhitp/I");
    pair[i]->Branch("d2dcap",   &m_d2dcap,   "d2dcap/F");
    pair[i]->Branch("zdcap",    &m_zdcap,    "zdcap/F");
    // em
    pair[i]->Branch("momm",     &m_momm,     "momm/F");
    pair[i]->Branch("phi0m",    &m_phi0m,    "phi0m/F");
    pair[i]->Branch("the0m",    &m_the0m,    "the0m/F");
    pair[i]->Branch("n0m",      &m_n0m,      "n0m/F");
    pair[i]->Branch("ch2npe0m", &m_ch2npe0m, "ch2npe0m/F");
    pair[i]->Branch("dispm",    &m_dispm,    "dispm/F");
    pair[i]->Branch("ecorem",   &m_ecorem,   "ecorem/F");
    pair[i]->Branch("epm",      &m_epm,      "epm/F");
    pair[i]->Branch("emcdpm",   &m_emcdpm,   "emcdpm/F");
    pair[i]->Branch("emcdzm",   &m_emcdzm,   "emcdzm/F");
    pair[i]->Branch("dphim",    m_dphim,     "dphim[4]/F");
    pair[i]->Branch("dzm",      m_dzm,       "dzm[4]/F");
    pair[i]->Branch("chi2ndfm", &m_chi2ndfm, "chi2ndfm/F");
    pair[i]->Branch("nhitm",    &m_nhitm,    "nhitm/I");
    pair[i]->Branch("d2dcam",   &m_d2dcam,   "d2dcam/F");
    pair[i]->Branch("zdcam",    &m_zdcam,    "zdcam/F");

    if(m_simmode>0){
      pair[i]->Branch("simmee",     &m_simmee,     "simmee/F");
      pair[i]->Branch("simxpt",     &m_simxpt,     "simxpt/F");
      pair[i]->Branch("simxpz",     &m_simxpz,     "simxpz/F");
      pair[i]->Branch("simxdca2d",  &m_simxdca2d,  "simxdca2d/F");
      pair[i]->Branch("simxlength", &m_simxlength, "simxlength/F");
      pair[i]->Branch("simlxy",     &m_simlxy,     "simlxy/F");
      pair[i]->Branch("simpaidp", &m_simpaidp, "simpaidp/I");
      pair[i]->Branch("simpaidm", &m_simpaidm, "simpaidm/I");
      pair[i]->Branch("simvrp",   &m_simvrp,   "simvrp/F");
      pair[i]->Branch("simvrm",   &m_simvrm,   "simvrm/F");
      pair[i]->Branch("simvzp",   &m_simvzp,   "simvzp/F");
      pair[i]->Branch("simvzm",   &m_simvzm,   "simvzm/F");
      pair[i]->Branch("simptprp", &m_simptprp, "simptprp/F");
      pair[i]->Branch("simptprm", &m_simptprm, "simptprm/F");
    }
  }

  ///////////////
  // e trk
  m_etree = new TTree("etree",    "electron tree");

//  m_etree->Branch("bbcz",    &m_bbcz,    "bbcz/F");
  m_etree->Branch("bbcq",    &m_bbcq,    "bbcq/F");
  m_etree->Branch("xvtx",    &m_xvtx,    "xvtx/F");
  m_etree->Branch("yvtx",    &m_yvtx,    "yvtx/F");
  m_etree->Branch("zvtx",    &m_zvtx,    "zvtx/F");
  m_etree->Branch("mom",     &m_mom,     "mom/F");
  m_etree->Branch("phi0",    &m_phi0,    "phi0/F");
  m_etree->Branch("the0",    &m_the0,    "the0/F");
  m_etree->Branch("c",       &m_c,       "c/F");
  m_etree->Branch("dcqual",  &m_dcqual,  "dcqual/I");
  m_etree->Branch("emcdphi", &m_emcdphi, "emcdphi/F");
  m_etree->Branch("emcdz",   &m_emcdz,   "emcdz/F");
  m_etree->Branch("n0",      &m_n0,      "n0/F");
  m_etree->Branch("ch2npe0", &m_ch2npe0, "ch2npe0/F");
  m_etree->Branch("disp",    &m_disp,    "disp/F");
  m_etree->Branch("ecore",   &m_ecore,   "ecore/F");
  m_etree->Branch("ep",      &m_ep,      "ep/F");
  m_etree->Branch("dphi",    m_dphi,     "dphi[4]/F");
  m_etree->Branch("dz",      m_dz,       "dz[4]/F");
  m_etree->Branch("chi2ndf", &m_chi2ndf, "chi2ndf/F");
  m_etree->Branch("nhit",    &m_nhit,    "nhit/I");
  m_etree->Branch("d2dca",   &m_d2dca,   "d2dca/F");
  m_etree->Branch("zdca",    &m_zdca,    "zdca/F");
  m_etree->Branch("convtag", &m_convtag, "convtag/I");
  m_etree->Branch("svxmom",  &m_svxmom,  "svxmom/F");
  m_etree->Branch("svxphi0", &m_svxphi0, "svxphi0/F");
  m_etree->Branch("svxthe0", &m_svxthe0, "svxthe0/F");
  m_etree->Branch("svxc",    &m_svxc,    "svxc/F");

  if(m_simmode==1){
    m_etree->Branch("simpaid", &m_simpaid, "simpaid/I");
    m_etree->Branch("simpx",   &m_simpx,   "simpx/F");
    m_etree->Branch("simpy",   &m_simpy,   "simpy/F");
    m_etree->Branch("simpz",   &m_simpz,   "simpz/F");
    m_etree->Branch("simvr",   &m_simvr,   "simvr/F");
    m_etree->Branch("simvz",   &m_simvz,   "simvz/F");
    m_etree->Branch("simptpr", &m_simptpr, "simptpr/F");
  }


}
