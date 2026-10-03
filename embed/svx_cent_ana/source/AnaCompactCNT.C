#include <iostream>
#include <string>

#include "AnaCompactCNT.h"

#include "phool.h"
#include "PHTypedNodeIterator.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "Fun4AllReturnCodes.h"


#include "RunHeader.h"
#include "EventHeader.h"
#include "PHGlobal.h"
#include "VtxOut.h"
#include "PHPoint.h"
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
static const float M_e = 0.000511; // GeV/c^2

static const int NBUF  = 20;
static const int NZ    = 16;
static const int NCENT = 18;

static int           idxbuf[NZ][NCENT];
static vector<eTrk*> vepbuf[NBUF][NZ][NCENT];
static vector<eTrk*> vembuf[NBUF][NZ][NCENT];


class eTrk { // only dcqual==31 or 63 + n0>2 + E/p>0.7
  public:
    eTrk() : mom(-9999.), phi0(-9999.), the0(-9999.), c(-9999.),
             dcqual(-999),
             emcdphi(-9999.), emcdz(-9999.), 
             n0(-9999), cch2npe0(-9999.), disp(-9999.),
             ecore(-9999.), ep(-9999.),
             svxchi2ndf(-9999.), 
             nhit(-9999), 
             d2dca(-9999.), zdca(-9999.)
      {
        for(int i=0; i<4; i++){
          svxdphi[i] = -9999.;
          svxdz[i]   = -9999.;
        }
        convtag=0;
      }
    virtual ~eTrk(){}

    float pt() { return mom*sin(the0); }
    float px() { return pt()*cos(phi0); }
    float py() { return pt()*sin(phi0); }
    float pz() { return mom*cos(the0); }

    void print(){
      cout<<"eTrk : "<<mom<<" "<<phi0<<" "<<the0<<" "<<c<<" : ";
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
 
    float svxdphi[4], svxdz[4];
    float svxchi2ndf;
    int   nhit;
    float d2dca, zdca;

    int convtag;

    int   simpaid;
    float simvr, simvz;

};

void AnaCompactCNT::init_idxbuf(){
 for(int iz=0; iz<NZ; iz++){
   for(int icent=0; icent<NCENT; icent++){
     idxbuf[iz][icent] = 0;
   }
 }
}

void AnaCompactCNT::calc_iz_icent(float zvtx, float bbcq, int& iz, int& icent){
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

//  cout<<"iz : icent = "<<zvtx<<" "<<iz<<" : "<<bbcq<<" : "<<icent<<endl;
}


float AnaCompactCNT::calc_pair(eTrk *ep, eTrk *em,
                float& Mee,  float& px,   float& py,   float& pz, float& pt,
                float& thev, float& ptep, float& ptem, float& phiv)
{
  float pxep = ep->px();
  float pyep = ep->py();
  float pzep = ep->pz();
  float Eep  = sqrt(ep->mom*ep->mom+M_e*M_e);

  float pxem = em->px();
  float pyem = em->py();
  float pzem = em->pz();
  float Eem  = sqrt(em->mom*em->mom+M_e*M_e);

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


//==============================================================

AnaCompactCNT::AnaCompactCNT(string filename) : m_OutputFileName(filename)
{
  ThisName = "AnaCompactCNT";
  init_ana=0;
  EventNumber=0;

  m_simmode = 0;

  init_idxbuf();
}

//==============================================================

AnaCompactCNT::~AnaCompactCNT() {
}

//==============================================================

int AnaCompactCNT::Init(PHCompositeNode *topNode) {

  cout << "AnaCompactCNT::Init started..." << endl;
  m_OutputNtupleFile = new TFile(m_OutputFileName.c_str(),"RECREATE");
  cout << "AnaCompactCNT::Init: output file " << m_OutputFileName << " opened." << endl;


  initPairTree();

  cout << "AnaCompactCNT::Init ended." << endl;
  return 0;
}

//==============================================================
  
int AnaCompactCNT::InitRun(PHCompositeNode *topNode) {
  cout << "AnaCompactCNT::InitRun started..." << endl;

  cout << "AnaCompactCNT::InitRun ended." << endl;
  return 0;
}

//==============================================================


int AnaCompactCNT::process_event(PHCompositeNode *topNode) {

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
    cout<<"Error AnaCompactCNT::process_event no SvxCentralTrackMap"<<endl;
    return EVENT_OK;
  }
*/
  SvxCentralTrackRecalList *svxcntlist  = findNode::getClass<SvxCentralTrackRecalList>(topNode, "SvxCentralTrackRecalList");
  if(svxcntlist==NULL){
    cout<<"Error AnaCompactCNT::process_event no SvxCentralTrackRecalList"<<endl;
    return EVENT_OK;
  }

  McEvalSingleList *mctrklist  = findNode::getClass<McEvalSingleList>(topNode, "McEvalSingle");
  if(m_simmode>0){
    if(mctrklist==NULL){
      cout<<"Error AnaCompactCNT::process_event no McEvalSingleList"<<endl;
      return EVENT_OK;
    }
  }


  m_evt      = EventNumber;
  float cent = (global!=NULL) ? global->getCentrality() : -999.;
  m_bbcq     = (global!=NULL) ? (global->getBbcChargeN()+global->getBbcChargeS()) : -999.0;

  PHPoint vtxpos(-999., -999., -999.); 
  if(vtxout!=NULL)  vtxpos = vtxout->get_Vertex();
  m_xvtx = vtxpos.getX();
  m_yvtx = vtxpos.getY();
  m_zvtx = vtxpos.getZ();

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
      cout << "------- Event # " << EventNumber << " zvtx="<<m_zvtx<< " ntrack= " << ntrk << " centrality="<<cent<<" bbcq="<<m_bbcq<<" vtx-x-y-z="<<m_xvtx<<" "<<m_yvtx<<" "<<m_zvtx<<endl;
    }
  }




  EventNumber++;
  return 0;
}

void AnaCompactCNT::analyze_track(PHCentralTrack* cntlist, SvxCentralTrackRecalList *svxcntlist, McEvalSingleList *mctrklist){
  if(cntlist==NULL||svxcntlist==NULL){
    if(cntlist==NULL)    cout<<"PHCentralTrack object is NULL"<<endl;
    if(svxcntlist==NULL) cout<<"SvxCentralTrackRecalList object is NULL"<<endl;
    return;
  }

  // clear electron vector
  //clear_VeTrk(m_vep);
  //clear_VeTrk(m_vem);
  m_vep.clear();
  m_vem.clear();


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


  // CNT track loop
  for(int itrk=0; itrk<ntrk; itrk++){
    PHSnglCentralTrack *trk = cntlist->get_track(itrk);

    // goodTrack & electron ID
    if( isGoodTrack(trk) &&
        isElectron(trk->get_n0(), 
                   trk->get_chi2()/trk->get_npe0(), 
                   trk->get_disp(), 
                   trk->get_ecore()/trk->get_mom(), 
                   trk->get_prob()) &&
        trk->get_mom()*sin(trk->get_the0())>0.3  // momentum cut
      )
      { 
        eTrk *etrk = new eTrk();
        fill_eTrk(etrk, trk, svxCntAry[itrk]);
        if(m_simmode>0) fill_eTrk_sim(etrk, itrk, mctrklist);

        //etrk->print();
   
        if(trk->get_charge()>0){
          m_vep.push_back(etrk);
        } else {
          m_vem.push_back(etrk);
        }
      }
  }


  // delete
  delete [] svxCntAry;
}

void AnaCompactCNT::analyze_pair(){

  ///////////////////////////////////
  // foreground
  fill_pair(m_vep, m_vem, m_eepair);

  fill_single(m_vep);
  fill_single(m_vem);


  ///////////////////////////////////
  // event mixing
  int iz=0, icent=0;
  calc_iz_icent(m_zvtx, m_bbcq, iz, icent);

  for(int imix=0; imix<NBUF; imix++){
    vector<eTrk*>& bg_vep = vepbuf[imix][iz][icent];
    vector<eTrk*>& bg_vem = vembuf[imix][iz][icent];

    fill_pair(m_vep,  bg_vem, m_eepair_bg);  // +- mixing
    fill_pair(bg_vep, m_vem,  m_eepair_bg);  // -+ mixing
  }

  ///////////////////////////////////
  //  update the buffer for the event mixing
 
  // copy the current list into the buf
  int ibuf = idxbuf[iz][icent];
  // clear the buffer 
  clear_VeTrk(vepbuf[ibuf][iz][icent]);
  clear_VeTrk(vembuf[ibuf][iz][icent]);

  // copy 
  vepbuf[ibuf][iz][icent] = m_vep;
  vembuf[ibuf][iz][icent] = m_vem;

  idxbuf[iz][icent]++;
  if(idxbuf[iz][icent]>=NBUF) {
    idxbuf[iz][icent] = 0;
    //cout<<"idxbuf cleared : "<<iz<<" "<<icent<<endl;
  }

}

bool AnaCompactCNT::isGoodTrack(PHSnglCentralTrack* trk){
  return( (trk->get_quality()==31||trk->get_quality()==63) &&
          fabs(trk->get_emcdphi())<0.03 &&
          fabs(trk->get_emcdz())<12) ;
}

bool AnaCompactCNT::isElectron(int n0, float ch2npe0, float disp, float Ep, float prob){
  return (n0>=3&&Ep>0.5);
}

void AnaCompactCNT::fill_eTrk(eTrk* etrk, PHSnglCentralTrack* trk, SvxCentralTrackRecal* svxtrk)
{
  if(etrk==NULL) return;

  if(trk!=NULL){
    etrk->mom      = trk->get_mom();
    etrk->phi0     = trk->get_phi0()     ;
    etrk->the0     = trk->get_the0()     ;
    etrk->c        = trk->get_charge()   ;
    etrk->dcqual   = trk->get_quality()   ;
    etrk->emcdphi  = trk->get_emcdphi()  ;
    etrk->emcdz    = trk->get_emcdz()    ;
    etrk->n0       = trk->get_n0()       ;
    etrk->cch2npe0 = trk->get_chi2()/trk->get_npe0();
    etrk->disp     = trk->get_disp()     ;
    etrk->ecore    = trk->get_ecore()    ;
    etrk->ep       = trk->get_ecore()/trk->get_mom();
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
      etrk->svxdphi[i] = dphi[i];
      etrk->svxdz[i]   = dz[i];
    }
    int nhitall = svxtrk->get_nhit();
    int ndf = (nhitall==0) ? 0 : (nhitall-1)*2;
    etrk->svxchi2ndf = svxtrk->get_Chisquare()/ndf;
    etrk->nhit       = nhitall;
    etrk->d2dca      = svxtrk->get_DCA2D();
    etrk->zdca       = svxtrk->get_DCAZ();
  }

}

void AnaCompactCNT::fill_eTrk_sim(eTrk* etrk, int itrk, McEvalSingleList *mctrklist)
{
  if(etrk==NULL) return;

  if(m_simmode>0){
    etrk->simpaid = mctrklist->get_parentid(itrk);
    //cout<<"simpaid : "<<simpidpa<<endl;
    float vx =  mctrklist->get_vertexx(itrk);
    float vy =  mctrklist->get_vertexy(itrk);
    float vz =  mctrklist->get_vertexz(itrk);

    etrk->simvr = sqrt(vx*vx+vy*vy);
    etrk->simvz = vz;
  }
}

void AnaCompactCNT::fill_pair(vector<eTrk*>& vep, vector<eTrk*>& vem, TTree *eepair){
  vector<eTrk*>::iterator itr_ep, itr_em;
  for(itr_ep=vep.begin(); itr_ep!=vep.end(); ++itr_ep){
    for(itr_em=vem.begin(); itr_em!=vem.end(); ++itr_em){
      eTrk *ep = *itr_ep;
      eTrk *em = *itr_em;

      fill_pair(ep, em, eepair);
    }
  }
}

void AnaCompactCNT::fill_pair(eTrk* ep, eTrk *em, TTree *eepair){
  calc_pair(ep, em,
    m_mee, m_px, m_py, m_pz, m_pt, m_thev, m_ptep, m_ptem, m_phiv);

  if(ep->convtag==0){
    if     (0   <=m_mee&&m_mee<0.06) ep->convtag = 1;
    else if(0.06<=m_mee&&m_mee<0.12) ep->convtag = 2;
  }
  if(em->convtag==0){
    if     (0   <=m_mee&&m_mee<0.06) em->convtag = 1;
    else if(0.06<=m_mee&&m_mee<0.12) em->convtag = 2;
  }


  // ep
  fill_trkvalue(ep,
     m_momp, m_phi0p, m_the0p, m_n0p,  m_ch2npe0p, m_dispp,
     m_ecorep, m_epp,
     m_nhitp, m_chi2ndfp, m_d2dcap, m_zdcap,
     m_dphip, m_dzp,
     m_simpaidp, m_simvrp, m_simvzp);

  // em
  fill_trkvalue(em,
     m_momm, m_phi0m, m_the0m, m_n0m,  m_ch2npe0m, m_dispm,
     m_ecorem, m_epm,
     m_nhitm, m_chi2ndfm, m_d2dcam, m_zdcam,
     m_dphim, m_dzm,
     m_simpaidm, m_simvrm, m_simvzm);


  eepair->Fill();
}

void AnaCompactCNT::fill_single(std::vector<eTrk*>& ve)  // include loop
{
  vector<eTrk*>::iterator itr;
  for(itr=ve.begin(); itr!=ve.end(); ++itr){
    eTrk *etrk = *itr;

    fill_single(etrk);
  }
}

void AnaCompactCNT::fill_single(eTrk *etrk)
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

  if(m_simmode==1){
    m_simpaid = etrk->simpaid;
    m_simvr   = etrk->simvr;
    m_simvz   = etrk->simvz;
  }

  m_etree->Fill();
}


void AnaCompactCNT::fill_trkvalue(eTrk *etrk,
     float& mom, float& phi0, float& the0,
     float& n0,  float& ch2npe0, float& disp,
     float& ecore, float& ep,
     int&   nhit,
     float& chi2ndf, float& d2dca, float& zdca,
     float* dphi_ptr, float* dz_ptr,
     int& simpaid, float& simvr, float& simvz)
{
  mom     = etrk->mom;
  phi0    = etrk->phi0;
  the0    = etrk->the0;
  n0      = etrk->n0;
  ch2npe0 = etrk->cch2npe0;
  disp    = etrk->disp;
  ecore   = etrk->ecore;
  ep      = etrk->ecore;

  nhit    = etrk->nhit;
  chi2ndf = etrk->svxchi2ndf;
  d2dca   = etrk->d2dca;
  zdca    = etrk->zdca;
  for(int i=0; i<4; i++){
    dphi_ptr[i] = etrk->svxdphi[i];
    dz_ptr[i]   = etrk->svxdz[i];
  }
}



void AnaCompactCNT::clear_VeTrk(vector<eTrk*>& vetrk){
  vector<eTrk*>::iterator itr;

  for(itr=vetrk.begin(); itr!=vetrk.end(); ++itr){
    eTrk* etrk = *itr;
    delete etrk;
  }
  vetrk.clear();
}



int AnaCompactCNT::End(PHCompositeNode *topNode) {
  cout << "AnaCompactCNT::End:  Writing out..." << endl;
  m_OutputNtupleFile->Write();
  cout << "AnaCompactCNT::End:  Closing output file..." << endl;
  m_OutputNtupleFile->Close();
  delete m_OutputNtupleFile;
  m_OutputNtupleFile=0;
  return 0;
}


void AnaCompactCNT::initPairTree(){
  m_eepair = new TTree("eepair", "E+E- pair tree");
  m_eepair_bg = new TTree("eepair_bg", "E+E- pair tree(mix)");

  TTree *pair[2] = {m_eepair, m_eepair_bg};

  // event
  for(int i=0; i<2; i++){
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
    // ep
    pair[i]->Branch("momp",     &m_momp,     "momp/F");
    pair[i]->Branch("phi0p",    &m_phi0p,    "phi0p/F");
    pair[i]->Branch("the0p",    &m_the0p,    "the0p/F");
    pair[i]->Branch("n0p",      &m_n0p,      "n0p/F");
    pair[i]->Branch("ch2npe0p", &m_ch2npe0p, "ch2npe0p/F");
    pair[i]->Branch("dispp",    &m_dispp,    "dispp/F");
    pair[i]->Branch("ecorep",   &m_ecorep,   "ecorep/F");
    pair[i]->Branch("epp",      &m_epp,      "epp/F");
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
    pair[i]->Branch("dphim",    m_dphim,     "dphim[4]/F");
    pair[i]->Branch("dzm",      m_dzm,       "dzm[4]/F");
    pair[i]->Branch("chi2ndfm", &m_chi2ndfm, "chi2ndfm/F");
    pair[i]->Branch("nhitm",    &m_nhitm,    "nhitm/I");
    pair[i]->Branch("d2dcam",   &m_d2dcam,   "d2dcam/F");
    pair[i]->Branch("zdcam",    &m_zdcam,    "zdcam/F");

    if(m_simmode>0){
      pair[i]->Branch("simpaidp", &m_simpaidp, "simpaidp/I");
      pair[i]->Branch("simpaidm", &m_simpaidm, "simpaidm/I");
      pair[i]->Branch("simvrp",   &m_simvrp,   "simvrp/F");
      pair[i]->Branch("simvrm",   &m_simvrm,   "simvrm/F");
      pair[i]->Branch("simvzp",   &m_simvzp,   "simvzp/F");
      pair[i]->Branch("simvzm",   &m_simvzm,   "simvzm/F");
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

  if(m_simmode==1){
    m_etree->Branch("simpaid", &m_simpaid, "simpaid/I");
    m_etree->Branch("simvr",   &m_simvr,   "simvr/F");
    m_etree->Branch("simvz",   &m_simvz,   "simvz/F");
  }


}
