#define AnaNtpTrkPair_cxx
#include "AnaNtpTrkPair.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

class eTrk;

static const float Me = 0.000511; // GeV/c^2

static const int NBUF  = 20;
static const int NZ    = 16;
static const int NCENT = 18;

int           idxbuf[NZ][NCENT];
vector<eTrk*> vepbuf[NBUF][NZ][NCENT];
vector<eTrk*> vembuf[NBUF][NZ][NCENT];

void init_idxbuf(){
 for(int iz=0; iz<NZ; iz++){
   for(int icent=0; icent<NCENT; icent++){
     idxbuf[iz][icent] = 0;
   }
 }
}

void calc_iz_icent(float zvtx, float bbcq, int& iz, int& icent){
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


class eTrk { // only dcqual==31 or 63 + n0>2 + E/p>0.7
  public:
    eTrk(){ convtag=0;}
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
    float n0, cch2npe0, disp;
    float sn0, scch2npe0, sdisp;
    float ecore, ep;
 
    float svxdphi[4], svxdz[4];
    float svxfitdphi[4], svxfitdz[4];
    float svxchi2ndf;
    int   nhit;
    float d2dca, zdca;

    int   svxexn[4];
    float svxexdp[4][3], svxexdz[4][3];

    int   convtag;
    
    int   simpaid;
    float simvr, simvz;

};

float calc_pair(eTrk *ep, eTrk *em,
                float& Mee,  float& px,   float& py,   float& pz, float& pt,
                float& thev, float& ptep, float& ptem, float& phiv) 
{
  float pxep = ep->px();
  float pyep = ep->py();
  float pzep = ep->pz();
  float pep  = sqrt(pxep*pxep + pyep*pyep + pzep*pzep);
  float Eep  = sqrt(pep*pep + Me*Me);

  float pxem = em->px();
  float pyem = em->py();
  float pzem = em->pz();
  float pem  = sqrt(pxem*pxem + pyem*pyem + pzem*pzem);
  float Eem  = sqrt(pem*pem + Me*Me);

  px = pxep + pxem;
  py = pyep + pyem;
  pz = pzep + pzem;
  float E  = Eep + Eem;

  pt = sqrt(px*px + py*py);
  Mee = sqrt(E*E - px*px - py*py - pz*pz);

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

void AnaNtpTrkPair::fill_eTrk(eTrk* etrk){
  etrk->mom      = mom;
  etrk->phi0     = phi0     ;
  etrk->the0     = the0     ;
  etrk->c        = c        ;
  etrk->dcqual   = dcqual   ;
  etrk->emcdphi  = emcdphi  ;
  etrk->emcdz    = emcdz    ;
  etrk->n0       = n0       ;
  etrk->cch2npe0 = cch2/npe0;
  etrk->disp     = disp     ;
  etrk->ecore    = ecore    ;
  etrk->ep       = ecore/mom;

  etrk->sn0       = sn0       ;
  etrk->scch2npe0 = scch2/snpe0;
  etrk->sdisp     = sdisp     ;

  etrk->svxdphi[0] = (dproj0>-100) ? dproj0+bend0 : -999.;
  etrk->svxdphi[1] = (dproj1>-100) ? dproj1+bend1 : -999.;
  etrk->svxdphi[2] = (dproj2>-100) ? dproj2+bend2 : -999.;
  etrk->svxdphi[3] = (dproj3>-100) ? dproj3+bend3 : -999.;
  etrk->svxdz[0]   = (dproj0>-100) ? zproj0-zv0   : -999.;
  etrk->svxdz[1]   = (dproj1>-100) ? zproj1-zv1   : -999.;
  etrk->svxdz[2]   = (dproj2>-100) ? zproj2-zv2   : -999.;
  etrk->svxdz[3]   = (dproj3>-100) ? zproj3-zv3   : -999.;
  etrk->svxfitdphi[0] = (dproj0>-100) ? fitdp0 : -999.;
  etrk->svxfitdphi[1] = (dproj1>-100) ? fitdp1 : -999.;
  etrk->svxfitdphi[2] = (dproj2>-100) ? fitdp2 : -999.;
  etrk->svxfitdphi[3] = (dproj3>-100) ? fitdp3 : -999.;
  etrk->svxfitdz[0]   = (dproj0>-100) ? fitdz0 : -999.;
  etrk->svxfitdz[1]   = (dproj1>-100) ? fitdz1 : -999.;
  etrk->svxfitdz[2]   = (dproj2>-100) ? fitdz2 : -999.;
  etrk->svxfitdz[3]   = (dproj3>-100) ? fitdz3 : -999.;

  etrk->svxchi2ndf = chi2/ndf;
  etrk->nhit       = nhit;
  etrk->d2dca      = d2dca;
  etrk->zdca       = zdca;

  int     vexn[4]  = {exn0, exn1, exn2, exn3};
  float* vexdp[4] = {exdp0, exdp1, exdp2, exdp3};
  float* vexdz[4] = {exdz0, exdz1, exdz2, exdz3};
  for(int ilay=0; ilay<4; ilay++){
    etrk->svxexn[ilay] = vexn[ilay];
    for(int i=0; i<3; i++){
      etrk->svxexdp[ilay][i] = vexdp[ilay][i];
      etrk->svxexdz[ilay][i] = vexdz[ilay][i];
    }
  }

  if(m_simmode==1){
    etrk->simpaid = simpidpa;
    //cout<<"simpaid : "<<simpidpa<<endl;
    etrk->simvr = sqrt(simvx*simvx+simvy*simvy);
    etrk->simvz = simvz;
  }
}

void AnaNtpTrkPair::clear_VeTrk(vector<eTrk*>& vetrk){
  vector<eTrk*>::iterator itr;

  for(itr=vetrk.begin(); itr!=vetrk.end(); ++itr){
    eTrk* etrk = *itr;
    delete etrk;
  }
  vetrk.clear();
}

void AnaNtpTrkPair::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L AnaNtpTrkPair.C
//      Root > AnaNtpTrkPair t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   float pre_zvtx = -2000.;

   float pre_bbcq = -999, pre_bbcz = -999;
   float pre_xvtx=-999, pre_yvtx=-999;

   int ievt=0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if(fabs(pre_zvtx-zvtx)>0.0001){ // event based analysis
        if(ievt%1000==0) cout<<ievt<<" "<<zvtx<<" "<<pre_zvtx<<endl;

        m_bbcz = pre_bbcz;
        m_bbcq = pre_bbcq;
        m_xvtx = pre_xvtx;
        m_yvtx = pre_yvtx;
        m_zvtx = pre_zvtx;

        analyze_event();

  
        // clear the vector
        //clear_VeTrk(m_vep);
        //clear_VeTrk(m_vem);
        m_vep.clear(); // just clear the vector, not delete the object in the vector
        m_vem.clear();

        pre_bbcz = bbcz;
        pre_bbcq = bbcq;
        pre_xvtx = xvtx;
        pre_yvtx = yvtx;
        pre_zvtx = zvtx;
        ievt++;
      }

      // fill track to the pool
      if(n0>2&&ecore/mom>0.7){
        eTrk *etrk = new eTrk();
        fill_eTrk(etrk);

        h_pt->Fill(etrk->pt());
        if(c>0){ m_vep.push_back(etrk); } 
        else   { m_vem.push_back(etrk); }
      }
   }
}

void AnaNtpTrkPair::analyze_event(){
  //cout<<"ntrk : "<<m_vep.size()<<" "<<m_vem.size()<<endl;

  //float Mee, px, py, pz, pt, thev, ptep, ptem, phiv; 

  //unsigned int n_ep = m_vep.size();
  //unsigned int n_em = m_vem.size();
  //cout<<"ntrk : "<<n_ep+n_em<<" "<<n_ep<<" "<<n_em<<endl;

  // pair analysis
  fill_pair(m_vep, m_vem, m_eepair);  // include loop

  // single
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

void AnaNtpTrkPair::fill_pair(std::vector<eTrk*>& vep, std::vector<eTrk*>& vem, TTree *eepair)  // include loop
{
  vector<eTrk*>::iterator itr_ep, itr_em;
  for(itr_ep=vep.begin(); itr_ep!=vep.end(); ++itr_ep){
    for(itr_em=vem.begin(); itr_em!=vem.end(); ++itr_em){
      eTrk *ep = *itr_ep;
      eTrk *em = *itr_em;

      //cout<<"ep "; ep->print();
      //cout<<"em "; em->print();
     
      fill_pair(ep, em, eepair);
    }
  }
}

void AnaNtpTrkPair::fill_pair(eTrk* ep, eTrk *em, TTree *eepair){
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
     m_simpaidm, m_simvrp, m_simvzp);
  

  eepair->Fill();
}

void AnaNtpTrkPair::fill_trkvalue(eTrk *etrk, 
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

  if(m_simmode==1){
    simpaid = etrk->simpaid;
    simvr   = etrk->simvr;
    simvz   = etrk->simvz;
  }
}

void AnaNtpTrkPair::fill_single(std::vector<eTrk*>& ve)  // include loop
{
  vector<eTrk*>::iterator itr;
  for(itr=ve.begin(); itr!=ve.end(); ++itr){
    eTrk *etrk = *itr;

    fill_single(etrk);
  }
}

void AnaNtpTrkPair::fill_single(eTrk *etrk)
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



////////////
void AnaNtpTrkPair::init_ana(const char *outname){
  m_outfile = TFile::Open(outname, "recreate");

  // eepair
  m_eepair    = new TTree("eepair",    "E+E- pair tree");
  m_eepair_bg = new TTree("eepair_bg", "E+E- pair tree(mix)");

  TTree *pair[2] = {m_eepair, m_eepair_bg};

  // event
  for(int i=0; i<2; i++){
    pair[i]->Branch("bbcz",    &m_bbcz,    "bbcz/F");
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

    if(m_simmode==1){
      pair[i]->Branch("simpaidp", &m_simpaidp, "simpaidp/I");
      pair[i]->Branch("simpaidm", &m_simpaidm, "simpaidm/I");
      pair[i]->Branch("simvrp",   &m_simvrp,   "simvrp/F");
      pair[i]->Branch("simvrm",   &m_simvrm,   "simvrm/F");
      pair[i]->Branch("simvzp",   &m_simvzp,   "simvzp/F");
      pair[i]->Branch("simvzm",   &m_simvzm,   "simvzm/F");
    }
  }

  h_pt = new TH1F("h_pt", "pT", 100, 0, 10);

  ///////////////
  // e trk
  m_etree = new TTree("etree",    "electron tree");

  m_etree->Branch("bbcz",    &m_bbcz,    "bbcz/F");
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



  init_idxbuf();
}

void AnaNtpTrkPair::loop_a_file(const char *listfile){
  ifstream fin(listfile);

  if(!fin||fin.bad()){
    cout<<"listfile is bad: "<<listfile<<endl;
    return;
  }

  string fname;
  while( getline(fin, fname) ){
    cout<<fname.c_str()<<" "<<fname.size()<<endl;
    if(fname.size()==0) {
      //cout<<"File size is zero : "<<fname.c_str()<<endl;
      continue;
    }

    TFile *fntp = TFile::Open(fname.c_str());
    if(fntp==NULL){
      cout<<"File does not exist : "<<fname.c_str()<<endl;
      continue;
    }


    // cnttrk
    TTree *ntp_cnttrk = (TTree*)fntp->Get("ntp_cnttrk");
    cout<<"start cnttrk"<<endl;
    Init(ntp_cnttrk);
    Loop();


    fntp->Close();
    delete fntp;
    fChain = NULL;
  }
  fin.close();

  cout<<"end of loop_a_file"<<endl;
}

