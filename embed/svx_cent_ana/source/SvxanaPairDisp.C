#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "phool.h"
#include "PHTypedNodeIterator.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "Fun4AllReturnCodes.h"

#include "SvxanaPairDisp.h"

#include "VtxOut.h"
#include "PHCentralTrack.h"
#include "PHSnglCentralTrack.h"
#include "SvxCentralTrack.h"
#include "SvxCentralTrackList.h"
#include "McEvalSingleList.h"

#include "getClass.h"
#include "PHPoint.h"

#include "TCanvas.h"
#include "TH2.h"
#include "TPolyLine.h"
#include "TText.h"
#include "TMarker.h"
#include "TArc.h"
#include "TArrow.h"
#include "TLine.h"

using namespace std;
using namespace findNode;

static const float M_e  = 0.000511; // GeV/c^2
static const float M_pi = 0.1395;   // GeV/c^2
static const float M_K  = 0.4936;   // GeV/c^2



//==============================================================
class TrkCont {
  public:
    TrkCont(int cntid=-1, PHSnglCentralTrack* cnt=NULL, SvxCentralTrack* svxcnt=NULL, int mcid=-1) :
      m_cntid(cntid), m_cnt(cnt), m_svxcnt(svxcnt), m_mcid(mcid)
    { fillCntTrack(); }

    void fillCntTrack(){
      if(m_cnt==NULL) return;

      mom  = m_cnt->get_mom();
      phi0 = m_cnt->get_phi0();
      the0 = m_cnt->get_the0();
      charge = m_cnt->get_charge();

      pt = mom*sin(the0);
      px = pt*cos(phi0);
      py = pt*sin(phi0);
      pz = mom*cos(phi0);
    }

    void fillSvxCntTrack(SvxCentralTrack *svxcnt){
      if(svxcnt==NULL) return;

      m_svxcnt = svxcnt;
      svxpx = m_svxcnt->get3MomentumAtPrimaryVertex(0);
      svxpy = m_svxcnt->get3MomentumAtPrimaryVertex(1);
      svxpz = m_svxcnt->get3MomentumAtPrimaryVertex(2);
      dcapos[0] = m_svxcnt->getClosestApproach(0);
      dcapos[1] = m_svxcnt->getClosestApproach(1);
      dcapos[2] = m_svxcnt->getClosestApproach(2);

      svxpt =sqrt(svxpx*svxpx + svxpy*svxpy);

      calc_circleinfo(1.0);
    }

    void fillMcTrack(int mcid, float px, float py, float pz, 
                               float vx, float vy, float vz){
      m_mcid = mcid;
      simpx = px; simpy = py; simpz = pz;
      simvx = vx; simvy = vy; simvz = vz;
      simpt = sqrt(simpx*simpx + simpy*simpy);
    }

    void calc_circleinfo(float fieldScale){
      static const float B = 0.9;
      static const float b = 0.003*B;
    
      float b_sign = (fieldScale>0) ? 1.0 : -1.0;
    
      R   = svxpt/b;
      phi = atan2(svxpy, svxpx);
      cx  = dcapos[0] + (b_sign)*charge*R*(-sin(phi));
      cy  = dcapos[1] + (b_sign)*charge*R*( cos(phi));
      //cx  = dcapos[0] + (b_sign)*charge*R*( sin(phi));
      //cy  = dcapos[1] + (b_sign)*charge*R*(-cos(phi));
    
      //cout<<"AnaTrk : "<<R<<" "<<phi<<" "<<cx<<" "<<cy<<endl;
    };


    void fillCrossInfo(float px, float py, float pz){
      xpx = px; xpy = py; xpz = pz;
      xpt = sqrt(xpx*xpx + xpy*xpy);
    }

  public:
    int                 m_cntid;
    PHSnglCentralTrack *m_cnt;
    SvxCentralTrack    *m_svxcnt; 
    int                 m_mcid;

    // PHCentralTrack
    float mom, phi0, the0, charge;
    float px, py, pz, pt;

    // SvxCentralTrack
    float svxpx, svxpy, svxpz, svxpt; // mom vector at DCA pos
    float dcapos[3];
    float cx, cy, R, phi;

    float xpx, xpy, xpz, xpt; // mom at crossing position

    // McEvalSingle
    float simpx, simpy, simpz, simpt;
    float simvx, simvy, simvz;
};

class PairCont {
  public:
    PairCont()
    {
      trk[0]=trk[1]=NULL;
      ncross=0;
      cross[0] = cross[1] = 0;
      px=py=pz=0;
      kind=0;
    };

  public:
    TrkCont *trk[2];
    int   ncross;
    float cross[2];
    float px, py, pz, pt;
    int   kind;
};

//==============================================================

SvxanaPairDisp::SvxanaPairDisp() :
  SubsysReco("SvxanaPairDisp"),
  m_init_ana(0),
  m_EventNumber(0),
  m_isInitCanvas(false)
{
  m_c1 = NULL;
  m_frame = NULL;
  m_mvtx = NULL;

}

//==============================================================

SvxanaPairDisp::~SvxanaPairDisp() {
  if(m_c1!=NULL)    delete m_c1;
  if(m_frame!=NULL) delete m_frame;
  if(m_mvtx!=NULL)  delete m_mvtx;
}

//==============================================================

int SvxanaPairDisp::Init(PHCompositeNode *topNode) {

  cout << "SvxanaPairDisp::Init started..." << endl;
  cout << "SvxanaPairDisp::Init ended." << endl;
  return 0;
}

//==============================================================
  
int SvxanaPairDisp::InitRun(PHCompositeNode *topNode) {
  cout << "SvxanaPairDisp::InitRun started..." << endl;

  initCanvas();

  cout << "SvxanaPairDisp::InitRun ended." << endl;
  return 0;
}

//==============================================================

int SvxanaPairDisp::process_event(PHCompositeNode *topNode) {

  VtxOut *vtxout = getClass<VtxOut>(topNode,"VtxOut");

  PHCentralTrack      *trk           = getClass<PHCentralTrack>(topNode,"PHCentralTrack");
  SvxCentralTrackList *svxcnttrklist = getClass<SvxCentralTrackList>(topNode,"SvxCentralTrackList");
  McEvalSingleList    *mctrklist     = getClass<McEvalSingleList>(topNode, "McSingle");

  int ntrk = (trk!=NULL) ? trk->get_npart() : 0;

  if(m_init_ana==0) {
    cout<<"init ana"<<endl;
    topNode->print(); 
    m_init_ana=1;
 
    cout<<endl<<endl;
    cout<<"init ana ends"<<endl;
  }

  if(m_EventNumber%1==0) {
    cout << "------- Event # " << m_EventNumber << " ntrack= " << ntrk << endl;
  }


//------------------------ VTX information --------------------------
  //PHPoint vtxpos = vtxout->get_Vertex();
  //PHPoint vtxpos = vtxout->get_Vertex("BBC");
  PHPoint vtxpos = vtxout->get_Vertex("SIM");
  m_vtx[0] = vtxpos.getX();
  m_vtx[1] = vtxpos.getY();
  m_vtx[2] = vtxpos.getZ();
  cout<<"Vtx : "<<m_vtx[0]<<" "<<m_vtx[1]<<endl;


  //fillCntTrack(trk);

  fillCentralTrack(0, m_vtx[2], svxcnttrklist, trk, mctrklist);

  fillpair();


  m_EventNumber++;
  return 0;
}


void SvxanaPairDisp::fillCentralTrack(float bbcq, float zvtx, 
    SvxCentralTrackList* svxcnttrklist, PHCentralTrack *cntlist,
    McEvalSingleList *mctrklist
  )
{
  if(svxcnttrklist==NULL || cntlist==NULL){
    return;
  }

  // clear TrkCont vector
  vector<TrkCont*>::iterator itr;
  for(itr=m_vTrkCont.begin(); itr!=m_vTrkCont.end(); ++itr){
    TrkCont *tmp = *itr;
    delete tmp;
  }
  m_vTrkCont.clear();

 
  // fill PHCentralTrack
  int ncnt = (cntlist!=NULL) ? cntlist->get_npart() : 0;
  for(int itrk=0; itrk<ncnt; itrk++){
    PHSnglCentralTrack *cnt = cntlist->get_track(itrk);
    m_vTrkCont.push_back(new TrkCont(itrk, cnt));
  }

  // fill SvxCentralTrack
  int ntrk = svxcnttrklist->get_nCentralTracks();
  for(int itrk=0; itrk<ntrk; itrk++){
    SvxCentralTrack *svxtrk = svxcnttrklist->getCentralTrack(itrk);

    int cntidx=svxtrk->getDchIndex();
    if(cntidx>=(int)m_vTrkCont.size()){
      cerr<<"CntIdx is exceeded "<<cntidx<<" "<<m_vTrkCont.size()<<endl;
      continue;
    }

    TrkCont *trkcont = m_vTrkCont[cntidx];
    if(trkcont->m_cntid != cntidx){
      cout<<"Id is different : "<<cntidx<<" "<<trkcont->m_cntid<<endl;
      continue;
    }

    trkcont->fillSvxCntTrack(svxtrk);
  }

  // fill mc info
  if(mctrklist!=NULL){
    cout<<"N mctrack : "<<mctrklist->get_McEvalSingleTrackN()<<endl;
    for(unsigned int itrk=0; itrk<mctrklist->get_McEvalSingleTrackN(); itrk++){
  
      int nreco = mctrklist->get_Nreco(itrk);
      cout<<"   N mc reco : "<<nreco<<endl;
      for(int ireco=0; ireco<nreco; ireco++){
        int recoid  = mctrklist->get_recoid(itrk, ireco);
        if(0<=recoid&&recoid<ncnt){
          TrkCont *trkcont = m_vTrkCont[recoid];
          trkcont->fillMcTrack(itrk,
            mctrklist->get_momentumx(itrk),
            mctrklist->get_momentumy(itrk),
            mctrklist->get_momentumz(itrk),
            mctrklist->get_vertexx(itrk),
            mctrklist->get_vertexy(itrk),
            mctrklist->get_vertexz(itrk) 
          );
        }
        else { cout<<"Error :: out of range : "<<recoid<<" "<<itrk<<" "<<ireco<<endl; }
      }
    }
  }
  else {
    cout<<"No McEvalSingleList"<<endl;
  }

}

void SvxanaPairDisp::fillpair(){
  // clear PairCont vector
  vector<PairCont*>::iterator itr_pair;
  for(itr_pair=m_vPairCont.begin(); itr_pair!=m_vPairCont.end(); ++itr_pair){
    PairCont *tmp = *itr_pair;
    delete tmp;
  }
  m_vPairCont.clear();


  // charge by charge
  vector<TrkCont*> vTCp, vTCm;

  vector<TrkCont*>::iterator itr;
  for(itr=m_vTrkCont.begin(); itr!=m_vTrkCont.end(); ++itr){
    TrkCont *trkcont = *itr;
    if(trkcont->charge>0){ vTCp.push_back(trkcont); }
    else                 { vTCm.push_back(trkcont); }
  }

  // calc pair product
  vector<TrkCont*>::iterator itrp, itrm;
  for(itrp=vTCp.begin(); itrp!=vTCp.end(); ++itrp){
    for(itrm=vTCm.begin(); itrm!=vTCm.end(); ++itrm){
      TrkCont *trkp = *itrp;
      TrkCont *trkm = *itrm;

      int ncross;
      vector<float> v_crossR, v_crossX, v_crossY;
      int kind = calc_crosspos(trkp, trkm, ncross, v_crossR, v_crossX, v_crossY);

      PairCont *pair = new PairCont();
      m_vPairCont.push_back(pair);
 
      if(ncross>0){
        pair->ncross = ncross;
        pair->kind   = kind;
        pair->cross[0] = v_crossX[0];
        pair->cross[1] = v_crossY[0];
        pair->trk[0] = trkp;
        pair->trk[1] = trkm;

        // calc mom_vector at crossing point (closest)
        float ep_cpx, ep_cpy, ep_cpz, ep_cdphi;
        calc_vector_at_cross(
              trkp->svxpx, trkp->svxpy, trkp->svxpz, 
              trkp->dcapos[0], trkp->dcapos[1], trkp->dcapos[2],
              trkp->charge, trkp->R,
              v_crossX[0], v_crossY[0],
              &ep_cpx, &ep_cpy, &ep_cpz, &ep_cdphi
            );
        trkp->fillCrossInfo(ep_cpx, ep_cpy, ep_cpz);

        float em_cpx, em_cpy, em_cpz, em_dphi;
        calc_vector_at_cross(
              trkm->svxpx, trkm->svxpy, trkm->svxpz, 
              trkm->dcapos[0], trkm->dcapos[1], trkm->dcapos[2],
              trkm->charge, trkm->R,
              v_crossX[0], v_crossY[0],
              &em_cpx, &em_cpy, &em_cpz, &em_dphi
            );
        trkm->fillCrossInfo(em_cpx, em_cpy, em_cpz);

        // calc pair using vector at cross
        float x_mee, xpx, xpy, xpz, xpt, x_phiv, x_thev, x_ptep, x_ptem;
        calc_pair(ep_cpx, ep_cpy, ep_cpz, M_e,
                  em_cpx, em_cpy, em_cpz, M_e,
               x_mee, xpx, xpy, xpz, xpt, x_thev, x_ptep, x_ptem, x_phiv);
        pair->px = xpx;
        pair->py = xpy;
        pair->pz = xpz;
        pair->pt = sqrt(xpx*xpx+xpy*xpy);

/*

        // calc dca using x_vector
        float xvtx, yvtx, zvtx;
        getPrimVertex(&xvtx, &yvtx, &zvtx);
        calc_dca_with_xvector(
          m_xpx, m_xpy, m_xpz,
          v_crossX[0], v_crossY[0],
          xvtx, yvtx, zvtx,
          m_xdca2d 
        );
*/
      }


    }
  }
  
}


//==============================================================

int SvxanaPairDisp::End(PHCompositeNode *topNode) {
  cout << "SvxanaPairDisp::End" << endl;
  return 0;
}

void SvxanaPairDisp::initCanvas(){

  m_c1 = new TCanvas("c1","c1", 600, 600);

  //m_frame = new TH2F("frame", "", 100, -20, 20, 100, -20, 20);
  m_frame = new TH2F("frame", "", 100, -0.05, 0.05, 100, -0.05, 0.05);
  //m_frame = new TH2F("frame", "", 100, -0.05, 0.05, 100, -0.15, 0.15);

  m_isInitCanvas = true;
}

void SvxanaPairDisp::drawCanvas(){
  if(!m_isInitCanvas){
    cout<<"Canvas is Not initialized"<<endl;
    return;
  }

  m_c1->cd();
  m_frame->Draw();
}

void SvxanaPairDisp::drawEvent(bool drawcross, bool drawsim){
  if(!m_isInitCanvas){
    cout<<"Canvas is Not initialized"<<endl;
    return;
  }

  // reset TArc vector
  vector<TArc*>::iterator itrArc;
  for(itrArc=m_vArc.begin(); itrArc!=m_vArc.end(); ++itrArc){
    TArc *tmp = *itrArc;
    delete tmp;
  }
  m_vArc.clear();

  vector<TArrow*>::iterator itrArrow;
  for(itrArrow=m_vArrow.begin(); itrArrow!=m_vArrow.end(); ++itrArrow){
    TArrow *tmp = *itrArrow;
    delete tmp;
  }
  m_vArrow.clear();

  vector<TMarker*>::iterator itrSp;
  for(itrSp=m_vMarker.begin(); itrSp!=m_vMarker.end(); ++itrSp){
    TMarker *tmp = *itrSp;
    delete tmp;
  }
  m_vMarker.clear();

  vector<TArrow*>::iterator itrArrowCross;
  for(itrArrowCross=m_vArrowCross.begin(); 
      itrArrowCross!=m_vArrowCross.end(); ++itrArrowCross){
    TArrow *tmp = *itrArrowCross;
    delete tmp;
  }
  m_vArrowCross.clear();

  vector<TArrow*>::iterator itrArrowMc;
  for(itrArrowMc=m_vArrowMc.begin(); itrArrowMc!=m_vArrowMc.end(); ++itrArrowMc){
    TArrow *tmp = *itrArrowMc;
    delete tmp;
  }
  m_vArrowMc.clear();

  
  // draw
  m_c1->cd();

  
  if(m_mvtx!=NULL) { delete m_mvtx; m_mvtx = NULL; }
  m_mvtx = new TMarker(m_vtx[0], m_vtx[1], 20);
  m_mvtx->SetMarkerSize(0.5);
  m_mvtx->SetMarkerColor(4);
  m_mvtx->Draw("same");


  int idx=0;
  vector<TrkCont*>::iterator itr;
  for(itr=m_vTrkCont.begin(); itr!=m_vTrkCont.end(); ++itr){
    cout<<"draw "<<idx<<endl;
    TrkCont *trkcont = *itr;

    if(trkcont->m_svxcnt==NULL) {
      cout<<"no svx track "<<idx<<endl;
      continue;
    }

    //TArc *tmp = makeLine(trkcont->svxpx, trkcont->svxpy, trkcont->svxpz, trkcont->charge, 
    //                     trkcont->dcapos[0], trkcont->dcapos[1]);
    TArc *tmp = makeLine(trkcont->cx, trkcont->cy, trkcont->R);

    tmp->Draw("sameonly");
    m_vArc.push_back(tmp);

    cout<<trkcont->dcapos[0]<<" "<<trkcont->dcapos[1]<<" "<<endl;

    // arrow
    double ax[2] = {trkcont->dcapos[0], 0.02*trkcont->svxpx/trkcont->svxpt + trkcont->dcapos[0]};
    double ay[2] = {trkcont->dcapos[1], 0.02*trkcont->svxpy/trkcont->svxpt + trkcont->dcapos[1]};
    //TArrow *ar = new TArrow(0, 0, 0.02*trkcont->px/trkcont->pt, 0.02*trkcont->py/trkcont->pt);
    TArrow *ar = new TArrow(ax[0], ay[0], ax[1], ay[1], 0.03);
    ar->Draw("same>");
    m_vArrow.push_back(ar);

    // draw sim info
    if(drawsim){
      if(trkcont->m_mcid>=0){
        double simax[2] = {trkcont->simvx, 0.01*trkcont->simpx/trkcont->simpt + trkcont->simvx};
        double simay[2] = {trkcont->simvy, 0.01*trkcont->simpy/trkcont->simpt + trkcont->simvy};
        TArrow *simar = new TArrow(simax[0], simay[0], simax[1], simay[1], 0.03);
        simar->SetLineColor(6);
        simar->Draw("same>");
        m_vArrow.push_back(simar);
      }
      else {
        cout<<"No McTrackInfo "<<idx<<endl;
      }
    }

    idx++;
  }

  idx=0;
  // draw pair
  vector<PairCont*>::iterator itr_pair;
  for(itr_pair=m_vPairCont.begin(); itr_pair!=m_vPairCont.end(); ++itr_pair){
    PairCont *pair = *itr_pair;
    cout<<"Pair kind="<<pair->kind<<" : "<<pair->cross[0]<<" "<<pair->cross[1]<<endl;


    TMarker *m = new TMarker(pair->cross[0], pair->cross[1], 22);
    m->SetMarkerColor(7);
    m->Draw("same");

    m_vMarker.push_back(m);
   
    double ax[2] = {pair->cross[0], 0.01*pair->px/pair->pt + pair->cross[0]};
    double ay[2] = {pair->cross[1], 0.01*pair->py/pair->pt + pair->cross[1]};
    TArrow *ar = new TArrow(ax[0], ay[0], ax[1], ay[1], 0.03);
    ar->SetLineColor(3);
    ar->Draw("same>");
    m_vArrowCross.push_back(ar);

    // arrow at crossing point
    if(drawcross){
      for(int i=0; i<2; i++){
        if(pair->trk[i]!=NULL){
          double xax[2] = {pair->cross[0], 0.01*pair->trk[i]->xpx/pair->trk[i]->xpt + pair->cross[0]};
          double xay[2] = {pair->cross[1], 0.01*pair->trk[i]->xpy/pair->trk[i]->xpt + pair->cross[1]};
          TArrow *xar = new TArrow(xax[0], xay[0], xax[1], xay[1], 0.03);
          xar->SetLineColor(7);
          xar->Draw("same>");
          m_vArrowCross.push_back(xar);
        }
      }
    }

    idx++;
  }

}


TArc* SvxanaPairDisp::makeLine(
    float cx, float cy, float R
  )
{
  float draw_phi = atan2(-cy, -cx);

  //cout<<"  "<<cx<<" "<<cy<<" "<<R<<" "<<draw_phi*180./3.1415<<endl;

  TArc *a = new TArc(cx, cy, R, 0, draw_phi*180./3.1415);
  //TArc *a = new TArc(cx, cy, R);
  //TArc *a = new TArc(cx, cy, R, 0, 180);
  a->SetLineStyle(3);
  a->SetLineColor(2);
  a->SetFillStyle(0);
  a->SetFillColor(1);
  a->SetNoEdges();

  return a;
}

int SvxanaPairDisp::calc_crosspos(TrkCont *ep, TrkCont *em,
                               int& ncross, vector<float>& v_crossR, 
                   vector<float>& v_crossX, vector<float>& v_crossY)
{

  float dx = em->cx - ep->cx;
  float dy = em->cy - ep->cy;
  float L  = sqrt(dx*dx + dy*dy);

  float Rdiff = fabs(ep->R - em->R);
  float Radd  = ep->R+em->R;
  if(L<Rdiff) { // small circle are wrapped by big circle (small is in the large circle)
    //cout<<" No crossing point"<<endl;
    ncross = 0;
    return ncross;
  }

  // if 2 circle are crossing  (or touching)
  if(Rdiff<=L && L<=Radd) {
    float phi      = atan2(dy, dx);
    float cosAlpha = (L*L + ep->R*ep->R - em->R*em->R)/(2*L*ep->R);
    float alpha    = acos(cosAlpha);

    float cross_x[2], cross_y[2], cross_r[2];
    for(int i=0; i<2; i++){
      float sign = (i==0) ? 1:-1;
      float x_x = ep->cx + ep->R*cos(phi+sign*alpha);
      float x_y = ep->cy + ep->R*sin(phi+sign*alpha);

      cross_x[i] = x_x;
      cross_y[i] = x_y;
      cross_r[i] = sqrt(x_x*x_x + x_y*x_y);
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

void SvxanaPairDisp::calc_vector_at_cross(
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

  //float sign = (charge>0) ? -1 : 1;
  float sign = (charge>0) ? 1 : -1;
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

float SvxanaPairDisp::calc_pair(float pxp, float pyp, float pzp, float Mep,
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

