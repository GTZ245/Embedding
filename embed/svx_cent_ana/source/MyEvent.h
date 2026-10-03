#ifndef __MYEVENT_H__
#define __MYEVENT_H__

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MyEvent                                                              //
//                                                                      //
// Description of the event and track parameters                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TMath.h"
#include <climits>
#include <map>
#include <vector>
#include <deque>
#include <iostream>
#include "MyEventConstants.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TTree.h"
#include "stdio.h"
#include "TGraph.h"
#include "TF1.h"

class TDirectory;

namespace DileptonAnalysis
{

      class MyTrack : public TObject
      {

      private:
            short trkid;       // DC track index
            short trkinfo;     // arm: 0 east, 1 west
            short trk_quality; // side 0 south, 1 north
                             // sector 0-7 with 0 being lowest west sector

            // basic parameters from tracking
            short q;        // charge of the track
            float pt;     // transverse momentum
            float phi_DC; // phi of track at 220 cm
            float z_DC;   // z poition of track at 220 cm (measured by PC1)
            float phi0;   // phi at interaction point - assuming track is primary
            float the0;   // dito but theta
            float alpha;  // alpha inclination at DC (total field bent = alpha + (phi_DC-phi0))
            short arm;
            short dcside;
            short sect;
            //float qpt;
            short n_hits;

            // Updated track parameters

            short q_prime;
            float phi0_prime;
            float the0_prime;
            float pt_prime;
            float alpha_prime;

            // EMCal variables
            short emcid;     // EMC cluster
            int emctower;  // emctower = 1000*iy + 10*iz + sector
                           // EMC sector 0-4 from bottom of arm
                           // iy tower coordinates in rphi
                           // iz tower coordinates in z
            float ecore;   // energy of cluster - assuning its a photon cluster
            float prob;    // similar to chi2 (is possibly redundant)
float x;
float y;
float z;
            float emcdz;   // delta z of cluster and track projection
            float emcdphi; // delta phi of cluster and track projection
            float emctof;  // TOF of EMCal (PbSc only)
            float dep;     // sigmalized E/p

            // TOF variables (should include TOF east and TOF west)
            float tofe; // TOF west??
            float tofdphi;
            float tofdz;
            // pc3 info
            float pc3sdphi;
            float pc3sdz;

            float crkphi; // phi position of ring
            float crkz;   // z position of ring

            // MC variables
            short mcid; // associated MC track
            

      public:
            MyTrack()
            {
                  trkid = -8888; // max value of unsigned short
                  trkinfo = -8888;
                  trk_quality = -8888;
                  arm = -8888;
                  dcside = -8888;
                  sect = -8888;
                  emcid = -8888;
                  emctower = -8888;
                  mcid = -8888;
                  q = -8888;
                  q_prime = -8888;
                  n_hits = -8888;
                  pt = -8888;
                  pt_prime = -8888;
                  phi_DC = -8888;
                  z_DC = -8888;
                  phi0 = -8888;
                  phi0_prime = -8888;
                  the0 = -8888;
                  the0_prime = -8888;
                  alpha = -8888;
                  alpha_prime = -8888;
                  ecore = -8888;
                  dep = -8888;
                  prob = -8888;
x = -8888;
y = -8888;
z = -8888;
                  emcdz = -8888;
                  emcdphi = -8888;
                  emctof = -8888;
                  tofe = -8888;
                  tofdphi = -8888;
                  tofdz = -8888;
                  pc3sdphi = -8888;
                  pc3sdz = -8888;
                  crkphi = -8888;
                  crkz = -8888;
            };

            virtual ~MyTrack(){};

            int GetTrkId() const { return trkid; };
            int GetTrkInfo() const { return trkinfo; };
            int GetTrkQuality() const { return trk_quality; };
            int GetArm() const { return arm; };
            int GetDCside() const { return dcside; };
            int GetSect() const { return sect; };
            float GetPx() const
            {
                  if (pt < -100)
                  {
                        return pt;
                  }
                  else
                  {
                        return pt_prime * cos(phi0);
                  }
            };
            float GetPy() const
            {
                  if (pt < -100)
                  {
                        return pt;
                  }
                  else
                  {
                        return pt_prime * sin(phi0);
                  }
            };
            float GetPz() const
            {
                  if (pt < -100)
                  {
                        return pt;
                  }
                  else
                  {
                        return pt_prime / tan(the0_prime);
                  }
            };
            float GetPtot() const
            {
                  if (pt < -100)
                  {
                        return pt;
                  }
                  else
                  {
                        return pt_prime * sqrt( (float) ( 1.+1./tan(the0_prime)/tan(the0_prime) ) );
                  }
            };
            float GetPt() const { return pt; };
            float GetPtPrime() const { return pt_prime; };
            float GetPhiDC() const { return phi_DC; };
            float GetZDC() const { return z_DC; };
            float GetPhi0() const { return phi0; };
            int GetNHits() const { return n_hits; };
            float GetPhi0Prime() const { return phi0_prime; };
            float GetThe0() const { return the0; };
            float GetThe0Prime() const { return the0_prime; };
            float GetAlpha() const { return alpha; };
            float GetAlphaPrime() const { return alpha_prime; };
            int GetCharge() const { return q; };
            int GetChargePrime() const { return q_prime; };
            int GetEmcId() const { return emcid; };
            float GetEcore() const { return ecore; };
            float GetDep() const { return dep; };
            float GetProb() const { return prob; };
float GetX() const { return x; };
float GetY() const { return y; };
float GetZ() const { return z; };
            float GetEmcdz() const { return emcdz; };
            float GetEmcdphi() const { return emcdphi; };
            int GetEsect() const { return int(emctower - 10 * ((emctower - 1000 * (emctower / 1000)) / 10) - 1000 * (emctower / 1000)); };
            int GetYsect() const { return int(emctower / 1000); };
            int GetZsect() const { return int((emctower - 1000 * (emctower / 1000)) / 10); };
            float GetEmcTOF() const { return emctof; };
            float GetTOFE() const { return tofe; };
            float GetTOFDPHI() const { return tofdphi; };
            float GetTOFDZ() const { return tofdz; };
            float GetPC3SDPHI() const { return pc3sdphi; };
            float GetPC3SDZ() const { return pc3sdz; };
            int GetMcId() const { return mcid; };
            float GetCrkphi() const { return crkphi; };
            float GetCrkz() const { return crkz; };

            void SetTrkId(int strkid) { trkid = strkid; };
            void SetTrkQuality(int strk_quality) { trk_quality = strk_quality; };
            void SetArm(int sarm) { arm = sarm; };
            void SetTrkInfo(int sarm, int sdcside, int ssect) { trkinfo = short(ssect + 10 * sdcside + 100 * sarm); };
            void SetDCSide(int sdcside) { dcside = sdcside; };
            void SetSect(int ssect) { sect = ssect; };
            void SetPt(float spt) { pt = spt; pt_prime = spt; };
            void SetPtPrime(float spt_prime) { pt_prime = spt_prime; };
            void SetQ(int sq) { q = sq; };
            void SetQPrime(int sq_prime) { q_prime = sq_prime; };
            void SetPhiDC(float sphi_DC) { phi_DC = sphi_DC; };
            void SetZDC(float sz_DC) { z_DC = sz_DC; };
            void SetPhi0(float sphi0) { phi0 = sphi0; };
            void SetPhi0Prime(float sphi0_prime) { phi0_prime = sphi0_prime; };
            void SetThe0(float sthe0) { the0 = sthe0; };
            void SetThe0Prime(float sthe0_prime) { the0_prime = sthe0_prime; };
            void SetAlpha(float salpha) { alpha = salpha; };
            void SetAlphaPrime(float salpha_prime) { alpha_prime = salpha_prime; };
            void SetEmcId(int sid) { emcid = sid; };
            void SetEcore(float secore) { ecore = secore; };
            void SetDep(float sdep) { dep = sdep; };
            void SetProb(float sprob) { prob = sprob; };
void SetX(float sx) { x = sx; };
void SetY(float sy) { y = sy; };
void SetZ(float sz) { z = sz; };
            void SetEmcdz(float semcdz) { emcdz = semcdz; };
            void SetEmcdphi(float semcdphi) { emcdphi = semcdphi; }
            void SetEmcTower(float sesect, float sysect, float szsect) { emctower = 1000 * sysect + 10 * szsect + sesect; };
            void SetEmcTOF(float semctof) { emctof = semctof; };
            void SetTOFE(float stofe) { tofe = stofe; };
            void SetTOFDPHI(float stofdphi) { tofdphi = stofdphi; };
            void SetTOFDZ(float stofdz) { tofdz = stofdz; };
            void SetPC3SDPHI(float spc3sdphi) { pc3sdphi = spc3sdphi; };
            void SetPC3SDZ(float spc3sdz) { pc3sdz = spc3sdz; };
            void SetMcId(int smcid) { mcid = smcid; };
            void SetCrkphi(float scrkphi) { crkphi = scrkphi; };
            void SetCrkz(float scrkz) { crkz = scrkz; };

            void SetNHits(int sn_hits) { n_hits = sn_hits; };

            void SetPrimes(const float bbcz = 0, const float svxz = 0, const int rg_beamoffset = 0, const int rungroup = 0);

            float get_mean_phi_data(int rungroup, int centr_bin, int layer) { return mean_phi_pars[rungroup][centr_bin][layer]; };
            float get_mean_theta_data(int rungroup, int centr_bin, int layer) { return mean_theta_pars[rungroup][centr_bin][layer]; };
            float get_sigma_phi_data(int rungroup, int centr_bin, int layer)
            {
                  layer = arm + (1-q_prime) + 4 * layer;
                  return sigma_phi_pars[rungroup][centr_bin][layer][0] + sigma_phi_pars[rungroup][centr_bin][layer][1] * exp(sigma_phi_pars[rungroup][centr_bin][layer][2] * pt_prime)
                  +sigma_phi_pars[rungroup][centr_bin][layer][3] * pt_prime;
            };
            float get_sigma_theta_data(int rungroup, int centr_bin, int layer)
            {
                  layer = arm + (1-q_prime) + 4 * layer;
                  return sigma_theta_pars[rungroup][centr_bin][layer][0] + sigma_theta_pars[rungroup][centr_bin][layer][1] * exp(sigma_theta_pars[rungroup][centr_bin][layer][2] * pt_prime)
                  +sigma_theta_pars[rungroup][centr_bin][layer][3] * pt_prime;
            };

            float get_dynamic_mean_phi_data(int ilay, float phi_prev, int rungroup = 0 ) 
            { 
                  const int arg0 = 4*ilay + (1-q_prime) + arm;
                  return phi_mean_phi_params[arg0][0]+phi_mean_phi_params[arg0][1]*phi_prev; 
            };
            float get_dynamic_mean_theta_data(int ilay, float phi_prev, int rungroup = 0 )
            { 
                  const int arg0 = 4*ilay + (1-q_prime) + arm;
                  return (the_mean_the_params[arg0][0]+the_mean_the_params[arg0][1]*phi_prev)*(1+exp( (ilay<4?-0.4:-1) * pt_prime)); 
            };
            float get_dynamic_sigma_phi_data(int ilay, float phi_prev, int rungroup = 0 )
            {
                  const int arg0 = 4*ilay + (1-q_prime) + arm;
                  const float sigma_pt = phi_sigma_pt_params[arg0][0] + phi_sigma_pt_params[arg0][1] * exp(phi_sigma_pt_params[arg0][2] * pt_prime);
                  return sigma_pt;//*(phi_sigma_phi_params[rungroup][arg0][0]+phi_sigma_phi_params[rungroup][arg0][1]* phi_prev+phi_sigma_phi_params[rungroup][arg0][2]* phi_prev* phi_prev);
            };
            float get_dynamic_sigma_theta_data(int ilay, float phi_prev, int rungroup = 0 )
            {
                  const int arg0 = 4*ilay + (1-q_prime) + arm;
                  return the_sigma_pt_params[arg0][0] + the_sigma_pt_params[arg0][1] * exp(the_sigma_pt_params[arg0][2] * pt_prime);
            };
            float get_dynamic_smean_phi_data(int ilay, float phi_prev, int rungroup = 0 )
            { 
                  const int arg0 = 2*ilay + (1-q_prime)/2;
                  return phi_sMean_pt_params[rungroup][arg0][0] + phi_sMean_pt_params[rungroup][arg0][1] * exp(phi_sMean_pt_params[rungroup][arg0][2] * pt_prime); 
            };
            void ResetPrimes(const float bbcz = 0, const float svxz = 0, const int rg_beamoffset = 0);

            /// below defined virtual functions for ide to suggest based on inherentence in electron and Cluster classes

            virtual float GetChi2() const { return 0; };
            virtual int GetN0() const { return 0; };
            virtual int GetNpe0() const { return 0; };
            virtual float GetDisp() const { return 0; };
            virtual float GetEmcdz_e() const { return 0; };
            virtual float GetEmcdphi_e() const { return 0; };

            virtual void SetChi2(float schi2) { };
            virtual void SetN0(int sn0) { };
            virtual void SetNPE0(int snpe0) { };
            virtual void SetDISP(float sdisp) { };
            virtual void SetEmcdz_e(float semcdz) { };
            virtual void SetEmcdphi_e(float semcdphi) { }

            virtual void SetHitIndex(int iindex, int iilayer) { };
            virtual void AddHitCounter(int iilayer) {  };
            virtual void RemoveHitCounter(int iilayer) {  };
            virtual void SetMinDist(float imindist, int iilayer) { };
            virtual void SetGhost(int value) { };

            virtual int GetHitIndex(int iilayer) { return 0; };
            virtual int GetHitCounter(int iilayer) { return 0; };
            virtual float GetMinDist(int iilayer) { return 0; };
            virtual int GetGhost() { return 0; } 

            virtual void SetRConv(float val) {  };
            virtual void SetPhiConv(float val) {  };
            virtual void SetTheConv(float val) {  };
            virtual void SetdZedConv(float val) {  };
            virtual void SetPhiVConv(float val) {  };

            ClassDef(MyTrack, 1)
      };

      class MyElectron : public MyTrack
      {

      private:
            int n0;       // number of photo tubes in mask
            float disp;
            float npe0;
            float chi2;

            float emcdz_e;   // delta z of cluster and track projection
            float emcdphi_e; // delta phi of cluster and track projection

            int hit_index[4];
            int hit_counter[4];
            int isGhost;
            float min_dist[4];
            float min_dphi[4];
            float min_dthe[4];
            float min_sdphi[4];
            float min_sdthe[4];
            float DCA, DCA2;
            float sDCA;
            float DCA_X,DCA_X2;
            float DCA_Y,DCA_Y2;
            float true_pt;
            // Wenching information
            int isConv; // indicates is there Solution for Wenqing algoritm and how tigt is it
            float r_conv, phi_conv, the_conv, phiv_conv, dzed_conv;

            std::vector<float> dangle0,dangle1,dangle2,dangle3;

      public:
            MyElectron() : MyTrack()
            {
                  n0 = -8888; // max value for byte
                  disp = -8888;
                  chi2 = -8888;
                  npe0 = -8888;
                  emcdz_e = -8888;
                  emcdphi_e = -8888;
                  DCA = -8888;
                  sDCA = -8888;
                  DCA_X = -8888;
                  DCA_Y = -8888;
                  DCA2 = -8888;
                  DCA_X2 = -8888;
                  DCA_Y2 = -8888;
                  true_pt = -8888;
                  isConv=0;
                  r_conv = -8888;
                  phi_conv = -8888;
                  the_conv = -8888;
                  phiv_conv = -8888;
                  dzed_conv = -8888;

                  for (int iteri = 0; iteri < 4; iteri++)
                  {
                        hit_index[iteri] = -8888;
                        hit_counter[iteri] = 0;
                        min_dist[iteri] = 100.;
                        min_dphi[iteri] = 100;
                        min_dthe[iteri] = 100;
                        min_sdphi[iteri] = 100;
                        min_sdthe[iteri] = 100;
                  }
                  isGhost = 0;
            };

            virtual ~MyElectron(){};

            float GetChi2() const { return chi2; };
            int GetN0() const { return n0; };
            int GetNpe0() const { return npe0; };
            float GetDisp() const { return disp; };
            float GetEmcdz_e() const { return emcdz_e; };
            float GetEmcdphi_e() const { return emcdphi_e; };

            void SetChi2(float schi2) { chi2 = schi2; };
            void SetN0(int sn0) { n0 = sn0; };
            void SetNPE0(int snpe0) { npe0 = snpe0; };
            void SetDISP(float sdisp) { disp = sdisp; };
            void SetEmcdz_e(float semcdz) { emcdz_e = semcdz; };
            void SetEmcdphi_e(float semcdphi) { emcdphi_e = semcdphi; }

            void SetHitIndex(int iindex, int iilayer) { hit_index[iilayer] = iindex; };
            void AddHitCounter(int iilayer) { hit_counter[iilayer]++; };
            void SetHitCounter(int iilayer, int value) { hit_counter[iilayer] = value; };
            void RemoveHitCounter(int iilayer) { hit_counter[iilayer]--; };
            void ZeroHitCounters() { for (int i = 0; i < 4; i++) hit_counter[i]=0; };
            void SetMinDist(float imindist, int iilayer) { min_dist[iilayer] = imindist; };
            void SetMinDphi(float imindist, int iilayer) { min_dphi[iilayer] = imindist; };
            void SetMinDthe(float imindist, int iilayer) { min_dthe[iilayer] = imindist; };
            void SetGhost(int value) { isGhost=value; };
            void SetDCA(float value) {DCA = value;};
            void SetsDCA(float value) {sDCA = value;};
            void SetDCAX(float value) {DCA_X = value;};
            void SetDCAY(float value) {DCA_Y = value;};
            void SetDCA2(float value) {DCA2 = value;};
            void SetDCAX2(float value) {DCA_X2 = value;};
            void SetDCAY2(float value) {DCA_Y2 = value;};
            void SetReconPT(float value) {true_pt = value;};
            void SetMinsDphi(float imindist, int iilayer) { min_sdphi[iilayer] = imindist; };
            void SetMinsDthe(float imindist, int iilayer) { min_sdthe[iilayer] = imindist; };
            void SetIsConv(int value) {isConv = value;};
            void SetRConv(float val) { r_conv = val; };
            void SetPhiConv(float val) { phi_conv = val; };
            void SetTheConv(float val) { the_conv = val; };
            void SetdZedConv(float val) { dzed_conv = val; };
            void SetPhiVConv(float val) { phiv_conv = val; };

            int GetHitIndex(int iilayer) { return hit_index[iilayer]; };
            int GetHitCounter(int iilayer) { return hit_counter[iilayer]; };
            float GetMinDist(int iilayer) { return min_dist[iilayer]; };
            float GetMinDphi(int iilayer) { return min_dphi[iilayer]; };
            float GetMinDthe(int iilayer) { return min_dthe[iilayer]; };
            int GetGhost() { return isGhost; }
            float GetDCA() {return DCA;};
            float GetsDCA() {return sDCA;};
            float GetDCAX() {return DCA_X;};
            float GetDCAY() {return DCA_Y;};
            float GetDCA2() {return DCA2;};
            float GetDCAX2() {return DCA_X2;};
            float GetDCAY2() {return DCA_Y2;};
            float GetReconPT() {return true_pt;};
            float GetMinsDphi(int iilayer) { return min_sdphi[iilayer]; };
            float GetMinsDthe(int iilayer) { return min_sdthe[iilayer]; };
            int GetIsConv() { return isConv; };
            float GetRConv() { return r_conv; };
            float GetPhiConv() { return   phi_conv; };
            float GetTheConv() { return   the_conv; };
            float GetdZedConv() { return dzed_conv; };
            float GetPhiVConv() { return phiv_conv; };

            void SetdPhidThe(int layer, float val1, float val2, float val3, float val4, float val5, float val6) 
            {
                  if(layer==0) {dangle0.push_back(val1);dangle0.push_back(val2);dangle0.push_back(val3);dangle0.push_back(val4);dangle0.push_back(val5);dangle0.push_back(val6);}
                  if(layer==1) {dangle1.push_back(val1);dangle1.push_back(val2);dangle1.push_back(val3);dangle1.push_back(val4);dangle1.push_back(val5);dangle1.push_back(val6);}
                  if(layer==2) {dangle2.push_back(val1);dangle2.push_back(val2);dangle2.push_back(val3);dangle2.push_back(val4);dangle2.push_back(val5);dangle2.push_back(val6);}
                  if(layer==3) {dangle3.push_back(val1);dangle3.push_back(val2);dangle3.push_back(val3);dangle3.push_back(val4);dangle3.push_back(val5);dangle3.push_back(val6);}
            };

            float GetdPhi(int layer, unsigned int iter) 
            {
                  if     (layer==0&& 6*iter<dangle0.size()) return dangle0[iter*6];
                  else if(layer==1&& 6*iter<dangle1.size()) return dangle1[iter*6];
                  else if(layer==2&& 6*iter<dangle2.size()) return dangle2[iter*6];
                  else if(layer==3&& 6*iter<dangle3.size()) return dangle3[iter*6];
                  std::cout<<"GetdPhi "<<layer<<" "<<6*iter+0<<" "<<dangle0.size()<<std::endl;
                  return -999;
            };
            float GetdThe(int layer, unsigned int iter) 
            {
                  if     (layer==0 && 6*iter+1<dangle0.size()) return dangle0[iter*6+1];
                  else if(layer==1 && 6*iter+1<dangle1.size()) return dangle1[iter*6+1];
                  else if(layer==2 && 6*iter+1<dangle2.size()) return dangle2[iter*6+1];
                  else if(layer==3 && 6*iter+1<dangle3.size()) return dangle3[iter*6+1];
                  std::cout<<"GetdThe "<<layer<<" "<<6*iter+1<<" "<<dangle0.size()<<std::endl;
                  return -999;
            };
            float GetsdPhi(int layer, unsigned int iter) 
            {
                  if     (layer==0&& 6*iter+2<dangle0.size()) return dangle0[iter*6+2];
                  else if(layer==1&& 6*iter+2<dangle1.size()) return dangle1[iter*6+2];
                  else if(layer==2&& 6*iter+2<dangle2.size()) return dangle2[iter*6+2];
                  else if(layer==3&& 6*iter+2<dangle3.size()) return dangle3[iter*6+2];
                  std::cout<<"GetsdPhi "<<layer<<" "<<6*iter+2<<" "<<dangle0.size()<<std::endl;
                  return -999;
            };
            float GetsdThe(int layer, unsigned int iter) 
            {
                  if     (layer==0 && 6*iter+3<dangle0.size()) return dangle0[iter*6+3];
                  else if(layer==1 && 6*iter+3<dangle1.size()) return dangle1[iter*6+3];
                  else if(layer==2 && 6*iter+3<dangle2.size()) return dangle2[iter*6+3];
                  else if(layer==3 && 6*iter+3<dangle3.size()) return dangle3[iter*6+3];
                  std::cout<<"GetsdThe "<<layer<<" "<<6*iter+3<<" "<<dangle0.size()<<std::endl;
                  return -999;
            };
            float GetDist(int layer, unsigned int iter) 
            {
                  if     (layer==0 && 6*iter+4<dangle0.size()) return dangle0[iter*6+4];
                  else if(layer==1 && 6*iter+4<dangle1.size()) return dangle1[iter*6+4];
                  else if(layer==2 && 6*iter+4<dangle2.size()) return dangle2[iter*6+4];
                  else if(layer==3 && 6*iter+4<dangle3.size()) return dangle3[iter*6+4];
                  std::cout<<"GetDist "<<layer<<" "<<6*iter+4<<" "<<dangle0.size()<<std::endl;
                  return -999;
            };
            int GetHits(int layer, unsigned int iter) 
            {
                  if     (layer==0 && 6*iter+5<dangle0.size()) return (int) dangle0[iter*6+5];
                  else if(layer==1 && 6*iter+5<dangle1.size()) return (int) dangle1[iter*6+5];
                  else if(layer==2 && 6*iter+5<dangle2.size()) return (int) dangle2[iter*6+5];
                  else if(layer==3 && 6*iter+5<dangle3.size()) return (int) dangle3[iter*6+5];
                  std::cout<<"GetHits "<<layer<<" "<<6*iter+5<<" "<<dangle0.size()
                  <<" "<<dangle1.size()<<" "<<dangle2.size()<<" "<<dangle3.size()<<std::endl;
                  return -999;
            };
            void ClearNumberVectors(){dangle0.clear();dangle1.clear();dangle2.clear();dangle3.clear();};


            ClassDef(MyElectron, 1)
      };

      class MyCluster : public MyTrack
      {
      private:
            short hit_index[4];
            short hit_counter[4];
            short isGhost;
            float min_dist[4];
int tower;
float chi2;

      public:
            MyCluster() : MyTrack()
            {
                  for (int iteri = 0; iteri < 4; iteri++)
                  {
                        hit_index[iteri] = -8888;
                        hit_counter[iteri] = 0;
                        min_dist[iteri] = 100.;
                  }
                  isGhost = 0;
tower = -8888;
chi2 = -8888;
            };

            virtual ~MyCluster(){};

            void SetHitIndex(int iindex, int iilayer) { hit_index[iilayer] = iindex; };
            void AddHitCounter(int iilayer) { hit_counter[iilayer]++; };
            void RemoveHitCounter(int iilayer) { hit_counter[iilayer]--; };
            void SetMinDist(float imindist, int iilayer) { min_dist[iilayer] = imindist; };
            void SetGhost(int value) { isGhost=value; };
void SetChi2(float schi2) { chi2 = schi2; };
void SetTower(float sesect, float sysect, float szsect) { tower = (1000*sysect + 10*szsect + sesect); };
            int GetHitIndex(int iilayer) { return hit_index[iilayer]; };
            int GetHitCounter(int iilayer) { return hit_counter[iilayer]; };
            float GetMinDist(int iilayer) { return min_dist[iilayer]; };
            int GetGhost() { return isGhost; }
float GetChi2() const { return chi2; };
int GetIY() const { int iy; if(tower==-8888) {iy=tower;} else {iy = tower/1000;} return iy;};
int GetIZ() const { int iz; if(tower==-8888) {iz=tower;} else {iz = (tower-1000*(tower/1000))/10;} return iz; };
            ClassDef(MyCluster, 1)
      };

      class MyVTXHit : public TObject
      {

      private:
            unsigned int clustid; // = layer + ilayr*10 + sensor*100 + ladder * 1000 + clustid*100000 /// max  = 4,294,967,295 => clustid_max = 42949
            /*int layer;
            int ladder;
            int sensor;*/

            float xhit; // xyz position of cluster [cm]
            float yhit;
            float zhit;

            /*float rhit;
            float phihit;
            float thetahit;*/

            //short ilayer;

            std::map<int, float> associated_tracks;
            //std::vector<double> dangle1;

      public:
            MyVTXHit()
            {
                  clustid = 0;

                  xhit = -8888;
                  yhit = -8888;
                  zhit = -8888;

            };
            virtual ~MyVTXHit(){};

            int GetClustId() const { return clustid/100000; };
            int GetLayer() const { return clustid%10; };
            int GetiLayer() const { return clustid/10%10; };
            int GetSensor() const { return clustid/100%10; };
            int GetLadder() const { return clustid/1000%100; };

            float GetXHit() const { return xhit; };
            float GetYHit() const { return yhit; };
            float GetZHit() const { return zhit; };

            float GetPhiHit(const float xvtx=0, const float yvtx=0, const float zvtx=0);
            float GetTheHit(const float xvtx=0, const float yvtx=0, const float zvtx=0);


            float GetPhi() const
            {
                  float phi = atan2(yhit, xhit);
                  if (phi < -TMath::ACos(-1) / 2)
                        phi += 2 * TMath::ACos(-1);
                  return phi;
            };

            void SetClustId(int sclustid) { clustid = clustid - (clustid/100000)*100000 + sclustid*100000; };
            void SetLayer (int slayer) { clustid = clustid - clustid%10 + slayer; };
            void SetiLayer (int silayer) { clustid = clustid - (clustid/10%10)*10 + silayer*10;};
            void SetSensor(int ssensor) { clustid = clustid - (clustid/100%10)*100 + ssensor*100; };
            void SetLadder(int sladder) { clustid = clustid - (clustid/1000%100)*1000 + sladder*1000; };

            void SetXHit(float value) { xhit=value; };
            void SetYHit(float value) { yhit=value; };
            void SetZHit(float value) { zhit=value; };


            /*void SetPolars(const float xhit = 0, const float yhit = 0, const float zhit = 0, 
                           const float xvtx = 0, const float yvtx = 0, const float zvtx = 0);*/
            void SetiLayerFromR();

            void AddAssociatedTrack(int trackID, float dist) { associated_tracks[trackID] = dist; };

            int N_AssociatedTracks() { return associated_tracks.size(); };

            int GetAssociatedTrack(int itrack)
            {
                  std::map<int, float>::const_iterator iterator = associated_tracks.begin();
                  for (int iass = 0; iass < itrack; iass++)
                        ++iterator;
                  return iterator->first;
            };
            float GetAssociatedTrackDistance(int itrackid) { return associated_tracks[itrackid]; };
            void DeleteAssociatedTrack(int itrackid) { associated_tracks.erase(itrackid); };

            ClassDef(MyVTXHit, 1)
      };


      class MyGenTrack : public TObject
      {
      private:
            float px;
            float py;
            float pz;
            float vx;
            float vy;
            float vz;

            int id;

      public:
            MyGenTrack()
            {
                  px = -8888;
                  py = -8888;
                  pz = -8888;
                  vx = -8888;
                  vy = -8888;
                  vz = -8888;
                  id = -8888;
            };
            virtual ~MyGenTrack(){};

            float GetPx() const { return px; }
            float GetPy() const { return py; }
            float GetPz() const { return pz; }
            float GetPt() const { return sqrt(SQR(px)+SQR(py));}
            float GetVx() const { return vx; }
            float GetVy() const { return vy; }
            float GetVz() const { return vz; }
            float GetID() const { return id; }

            void SetPx(float val) { px = val; }
            void SetPy(float val) { py = val; }
            void SetPz(float val) { pz = val; }
            void SetVx(float val) { vx = val; }
            void SetVy(float val) { vy = val; }
            void SetVz(float val) { vz = val; }
            void SetID(int val) { id = val; }

            ClassDef(MyGenTrack, 1)
      };

      class MyBDTHit : public TObject
      {
      private:
            float sdphi[4];
            float sdthe[4];
            int isTrue[4];
            float secondhitphileft[4][2];
            float secondhitphiright[4][2];
            float secondhittheleft[4][2];
            float secondhittheright[4][2];
            int layer[2];
            int layersecondhitphileft[2][2];
            int layersecondhitphiright[2][2];
            int layersecondhittheleft[2][2];
            int layersecondhittheright[2][2];

      public:
            MyBDTHit()
            {
                  for (int iteri = 0; iteri < 4; iteri++)
                  {
                        sdphi[iteri] = -10;
                        sdthe[iteri] = -10;
                        isTrue[iteri] = 0;
                        for (int jteri = 0; jteri < 2; jteri++) 
                        {
                              secondhitphileft[iteri][jteri]  = -99;
                              secondhitphiright[iteri][jteri] = -99;
                              secondhittheleft[iteri][jteri]  = -99;
                              secondhittheright[iteri][jteri] = -99;
                        }
                  }
                  for (int iteri = 0; iteri < 2; iteri++)
                  {
                        layer[iteri] = -10;
                        layer[iteri] = -10;
                        for (int jteri = 0; jteri < 2; jteri++) 
                        {
                              layersecondhitphileft[iteri][jteri]  = -99;
                              layersecondhitphiright[iteri][jteri] = -99;
                              layersecondhittheleft[iteri][jteri]  = -99;
                              layersecondhittheright[iteri][jteri] = -99;
                        }
                  }
            };
            virtual ~MyBDTHit(){};

            float Getsdphi(int ilayer = 0) const { return sdphi[ilayer]; }
            float Getsdthe(int ilayer = 0) const { return sdthe[ilayer]; }
            float GetIsTrue(int ilayer = 0) const { return isTrue[ilayer]; }
            float GetSecondHitPhiR(int ilayer = 0, int ihit = 0) const { return secondhitphileft[ilayer][ihit]; }
            float GetSecondHitPhiL(int ilayer = 0, int ihit = 0) const { return secondhitphiright[ilayer][ihit]; }
            float GetSecondHitTheR(int ilayer = 0, int ihit = 0) const { return secondhittheleft[ilayer][ihit]; }
            float GetSecondHitTheL(int ilayer = 0, int ihit = 0) const { return secondhittheright[ilayer][ihit]; }
            int GetLayerSecondHitPhiR(int ilayer = 0, int ihit = 0) const { return layersecondhitphileft[ilayer][ihit]; }
            int GetLayerSecondHitPhiL(int ilayer = 0, int ihit = 0) const { return layersecondhitphiright[ilayer][ihit]; }
            int GetLayerSecondHitTheR(int ilayer = 0, int ihit = 0) const { return layersecondhittheleft[ilayer][ihit]; }
            int GetLayerSecondHitTheL(int ilayer = 0, int ihit = 0) const { return layersecondhittheright[ilayer][ihit]; }
            int GetOutiLayer(int ilayer = 2) const { return layer[ilayer-2];}
            
            void Setsdphi(int ilayer = 0, float val = 0) { sdphi[ilayer] = val; }
            void Setsdthe(int ilayer = 0, float val = 0) { sdthe[ilayer] = val; }
            void SetIsTrue(int ilayer = 0, float val = 0) { isTrue[ilayer] = val; }
            void SetOutiLayer(int ilayer = 2, int val = 0) { layer[ilayer-2] = val;}
            void SetSecondHitPhiR(int ilayer = 0, int ihit = 0, float val = 0) { secondhitphileft[ilayer][ihit] = val; }
            void SetSecondHitPhiL(int ilayer = 0, int ihit = 0, float val = 0) { secondhitphiright[ilayer][ihit] = val; }
            void SetSecondHitTheR(int ilayer = 0, int ihit = 0, float val = 0) { secondhittheleft[ilayer][ihit] = val; }
            void SetSecondHitTheL(int ilayer = 0, int ihit = 0, float val = 0) { secondhittheright[ilayer][ihit] = val; }
            void SetLayerSecondHitPhiR(int ilayer = 2, int ihit = 0, int val = 0) { layersecondhitphileft[ilayer-2][ihit] = val; }
            void SetLayerSecondHitPhiL(int ilayer = 2, int ihit = 0, int val = 0) { layersecondhitphiright[ilayer-2][ihit] = val; }
            void SetLayerSecondHitTheR(int ilayer = 2, int ihit = 0, int val = 0) { layersecondhittheleft[ilayer-2][ihit] = val; }
            void SetLayerSecondHitTheL(int ilayer = 2, int ihit = 0, int val = 0) { layersecondhittheright[ilayer-2][ihit] = val; }

            ClassDef(MyBDTHit, 1)
      };

      class MyBDTrack : public TObject
      {
      private:

            std::vector<DileptonAnalysis::MyBDTHit> BDTHITlist;
            float pt;
            float phi0;
            float the0;
            float phiDC;
            float zDC;
            float alpha;
            float ecore;
            int centrality;
            int charge;
            int arm;
float x;
float y;
float z;

      public:
            MyBDTrack()
            {
                  BDTHITlist.clear();
                  pt = -999;
                  phi0 = -999;
                  the0 = -999;
                  phiDC = -999;
                  zDC = -999;
                  alpha = -999;
                  ecore = -999;
                  centrality = -999;
                  charge = -999;
                  arm = -999;
x = -999;
y = -999;
z = -999;
            };
            virtual ~MyBDTrack(){};

            void AddBDTHit(const MyBDTHit *newBDTHit) { BDTHITlist.push_back(*newBDTHit); };
            Long64_t GetNBDThit() { return BDTHITlist.size(); };
            DileptonAnalysis::MyBDTHit *GetBDTHitEntry(unsigned int i) const { return const_cast<DileptonAnalysis::MyBDTHit *>(&(BDTHITlist[i])); };
            void RemoveBDTHitEntry(const unsigned int i){ BDTHITlist.erase(BDTHITlist.begin() + i); };
            
            void SetPt(float val) { pt = val; };
            void SetPhi0(float val) { phi0 = val; };
            void SetThe0(float val) { the0 = val; };
            void SetPhiDC(float val) { phiDC = val; };
            void SetZDC(float val) { zDC = val; };
            void SetAlpha(float val) { alpha = val; };
            void SetEcore(float val) { ecore = val; };
void SetX(float val) { x = val; };
void SetY(float val) { y = val; };
void SetZ(float val) { z = val; };
            void SetCentrality(int val) { centrality = val; };
            void SetCharge(int val) { charge = val; };
            void SetArm(int val) { arm = val; };

            float GetPt() const { return pt; };
            float GetPhi0() const { return phi0; };
            float GetThe0() const { return the0; };
            float GetPhiDC() const { return phiDC; };
            float GetZDC() const { return zDC; };
            float GetAlpha() const { return alpha; };
            float GetEcore() const { return ecore; };
            int GetCentrality() const { return centrality; };
            int GetCharge() const { return charge; };
            int GetArm() const { return arm; };

            ClassDef(MyBDTrack, 1)
      };

      class MyPair : public TObject
      {
      private:
            float phi_e;
            float phi_p;
            float theta_e;
            float theta_p;
            float r_pair;

            UShort_t id_e;
            UShort_t id_p;

      public:
            MyPair()
            {
                  phi_e = -8888;
                  phi_p = -8888;
                  theta_e = -8888;
                  theta_p = -8888;
                  r_pair = -8888;
                  id_e = -8888;
                  id_p = -8888;
            };
            virtual ~MyPair(){};

            float GetPhiElectron() const { return phi_e; }
            float GetPhiPositron() const { return phi_p; }
            float GetThetaElectron() const { return theta_e; }
            float GetThetaPositron() const { return theta_p; }
            float GetRPair() const { return r_pair; }
            float GetIDElectron() const { return id_e; }
            float GetIDPositron() const { return id_p; }

            void SetPhiElectron(float sphi_e) { phi_e = sphi_e; }
            void SetPhiPositron(float sphi_p) { phi_p = sphi_p; }
            void SetThetaElectron(float stheta_e) { theta_e = stheta_e; }
            void SetThetaPositron(float stheta_p) { theta_p = stheta_p; }
            void SetRPair(float sr_pair) { r_pair = sr_pair; }
            void SetIDElectron(int sid_e) { id_e = sid_e; }
            void SetIDPositron(int sid_p) { id_p = sid_p; }

            ClassDef(MyPair, 1)
      };

      class MyEvent : public TObject
      {

      private:
            int evtno;
            float BBCcharge;
            float BBCchargeN;
            float BBCchargeS;
            float BBCtimeN;
            float BBCtimeS;
            float zvertex;
            float centrality;
            float preciseX;
            float preciseY;
            float preciseZ;
            float run_number;
            float psi2_BBC; // BBC (2nd and 3rd)
            float psi3_BBC;
            float psi2_FVTXA0; // FVTXA0
            float psi3_FVTXA0;

            std::vector<DileptonAnalysis::MyElectron> TrackList;
            std::vector<DileptonAnalysis::MyCluster> ClusterList;
            std::vector<DileptonAnalysis::MyVTXHit> VTXHitList;
            std::vector<DileptonAnalysis::MyElectron> ElecCandList;
            std::vector<DileptonAnalysis::MyGenTrack> GenTrackList;

      public:
            MyEvent()
            {
                  evtno = -999;
                  BBCcharge = -999;
                  BBCchargeN = -999;
                  BBCchargeS = -999;
                  BBCtimeN = -999;
                  BBCtimeS = -999;
                  zvertex = -999;
                  centrality = -999;

                  preciseX = -999;
                  preciseY = -999;
                  preciseZ = -999;
                  run_number = -999;

                  psi2_BBC    = -999;
                  psi3_BBC    = -999;
                  psi2_FVTXA0    = -999;
                  psi3_FVTXA0    = -999;
            };

            virtual ~MyEvent(){};

            void ClearEvent();

            void SetEvtNo(int sevtno) { evtno = sevtno; };
            int GetEvtNo() { return evtno; };

            void SetRunNumber(int srun_number) { run_number = srun_number; };
            int GetRunNumber() { return run_number; };

            void SetBBCcharge(float sBBCcharge) { BBCcharge = sBBCcharge; };
            float GetBBCcharge() { return BBCcharge; };

            void SetBBCchargeN(float sBBCchargeN) { BBCchargeN = sBBCchargeN; };
            float GetBBCchargeN() { return BBCchargeN; };

            void SetBBCchargeS(float sBBCchargeS) { BBCchargeS = sBBCchargeS; };
            float GetBBCchargeS() { return BBCchargeS; };

            void SetBBCtimeN(float sBBCtimeN) { BBCtimeN = sBBCtimeN; };
            float GetBBCtimeN() { return BBCtimeN; };

            void SetBBCtimeS(float sBBCtimeS) { BBCtimeS = sBBCtimeS; };
            float GetBBCtimeS() { return BBCtimeS; };

            void SetVtxZ(float szvertex) { zvertex = szvertex; };
            float GetVtxZ() { return zvertex; };

            void SetCentrality(float scentrality) { centrality = scentrality; };
            float GetCentrality() { return centrality; };

            void SetPreciseX(float sPreciseX) { preciseX = sPreciseX; };
            float GetPreciseX() { return preciseX; };

            void SetPreciseY(float sPreciseY) { preciseY = sPreciseY; };
            float GetPreciseY() { return preciseY; };

            void SetPreciseZ(float sPreciseZ) { preciseZ = sPreciseZ; };
            float GetPreciseZ() { return preciseZ; };

            void SetPsi2BBC(float spsi2_BBC) { psi2_BBC = spsi2_BBC; };
            float GetPsi2BBC() { return psi2_BBC; };

            void SetPsi3BBC(float spsi3_BBC) { psi3_BBC = spsi3_BBC; };
            float GetPsi3BBC() { return psi3_BBC; };

            void SetPsi2FVTXA0(float spsi2_FVTXA0) { psi2_FVTXA0 = spsi2_FVTXA0; };
            float GetPsi2FVTXA0() { return psi2_FVTXA0; };

            void SetPsi3FVTXA0(float spsi3_FVTXA0) { psi3_FVTXA0 = spsi3_FVTXA0; };
            float GetPsi3FVTXA0() { return psi3_FVTXA0; };

            int GetRunGroup(int in_run_number = 0) const;

            void AddTrack(const MyElectron *newTrack) { TrackList.push_back(*newTrack); };
            void AddCluster(const MyCluster *newTrack) { ClusterList.push_back(*newTrack); };
            void AddVTXHit(const MyVTXHit *newVTXHit) { VTXHitList.push_back(*newVTXHit); };
            void AddElecCand(const MyElectron *newTrack) { ElecCandList.push_back(*newTrack); };
            void AddGenTrack(const MyGenTrack *newTrack) { GenTrackList.push_back(*newTrack); }

            Long64_t GetNtrack() { return TrackList.size(); };
            Long64_t GetNcluster() { return ClusterList.size(); };
            Long64_t GetNVTXhit() { return VTXHitList.size(); };
            Long64_t GetNeleccand() { return ElecCandList.size(); };
            Long64_t GetNgentrack() { return GenTrackList.size(); };

            MyElectron *GetEntry(unsigned int i) const { return const_cast<MyElectron *>(&(TrackList[i])); };
            MyCluster *GetClusterEntry(unsigned int i) const { return const_cast<MyCluster *>(&(ClusterList[i])); };
            MyVTXHit *GetVTXHitEntry(unsigned int i) const { return const_cast<MyVTXHit *>(&(VTXHitList[i])); };
            MyElectron *GetElecCand(unsigned int i) const { return const_cast<MyElectron *>(&(ElecCandList[i])); };
            MyGenTrack *GetGenTrack(unsigned int i) const { return const_cast<MyGenTrack *>(&(GenTrackList[i])); };

            void RemoveTrackEntry(const unsigned int i){ TrackList.erase(TrackList.begin() + i); };
            void RemoveClusterEntry(const unsigned int i){ ClusterList.erase(ClusterList.begin() + i); };
            void RemoveVTXHitEntry(const unsigned int i){ VTXHitList.erase(VTXHitList.begin() + i); };

            std::vector<DileptonAnalysis::MyElectron> GetTracks() { return TrackList; };
            std::vector<DileptonAnalysis::MyCluster> GetClusters() { return ClusterList; };
            std::vector<DileptonAnalysis::MyVTXHit> GetVTXHits() { return VTXHitList; };

            void SetDCA(const unsigned int itrk = 0, const int layer2 = 1);
            void SetDCA2(const unsigned int itrk = 0, const int layer3 = 2);
            void ReshuffleElectrons();

            ClassDef(MyEvent, 1) // MyEvent structure
      };

      class MyEventContainer : public TObject
      {

      private:
            MyEvent *event;
            std::vector<DileptonAnalysis::MyEvent> EventList;
            std::deque<DileptonAnalysis::MyEvent> evtbuff_list[MIX_CENTBIN][MIX_ZVTXBIN][MIX_RP2BIN];
            std::vector<DileptonAnalysis::MyBDTrack> BDTracklist;
            TFile *infile, *outfile;
            TH2D *hist_br, *hist_bz;
            TTree *tree;
            TH1D *event_hist, *centr_hist;
            TH3D *el_pt_hist;
            TH3D *dphi_hist[N_centr], *dthe_hist[N_centr], *sdphi_hist[N_centr], *sdthe_hist[N_centr];
            TH3D *chi2_ndf[N_centr], *truehithist, *truehitsigmahist, *charge_recover_hist, *ilayerhitshist[N_centr];
            TH3D *dphi_hist_el[N_centr], *dthe_hist_el[N_centr], *sdphi_hist_el[N_centr], *sdthe_hist_el[N_centr];
            TH3D *dphi_hist_el_dynamic[N_dynamic], *dthe_hist_el_dynamic[N_dynamic], *sdphi_hist_el_dynamic[N_dynamic], *sdthe_hist_el_dynamic[N_dynamic];
            TH3D *dphi_phi0_init_hist[nvtx_layers], *dthe_the0_init_hist[nvtx_layers];
            TH3D *dphi_phi0_corr_hist[nvtx_layers], *dthe_the0_corr_hist[nvtx_layers]; 
            TH3D *dphi_the0_init_hist[nvtx_layers], *dthe_phi0_init_hist[nvtx_layers];
            TH3D *dphi_the0_corr_hist[nvtx_layers], *dthe_phi0_corr_hist[nvtx_layers]; 
            TH3D *d_dphi_hist[N_centr], *d_dthe_hist[N_centr], *DCA_hist[N_centr];
            TH3D *sd_dphi_hist[N_centr], *sd_dthe_hist[N_centr], *sDCA_hist[N_centr];
            TH3D *temc, *ttof, *n0_hist, *ep_hist, *prob_hist, *disp_hist, *chi2npe0_hist;
            TH3D *DCA_2D_hist[N_centr], *sDCA_2D_hist[N_centr];
            TH3D *emc_dphi_el, *emc_dz_el, *n0_hist_el, *ep_hist_el, *prob_hist_el, *disp_hist_el, *chi2npe0_hist_el, *rich_prob1, *rich_prob2, *rich_prob3;
            TH3D *el_had_dphi, *el_had_dz, *el_had_dr;
            TH3D *DCA2_hist[N_centr], *sDCA2_hist[N_centr], *DCA2_2D_hist[N_centr], *sDCA2_2D_hist[N_centr], *DCA12_hist[N_centr];
            TH3D *DCPT_ReconPT, *sDCPT_ReconPT,  *charge_hist, *phi_hist, *the_hist, *pt_corr;
            TH3D *veto_phi_hist[N_centr], *veto_the_hist[N_centr], *veto_phi_phi_hist[N_centr], *veto_the_the_hist[N_centr];
            TH3D *veto_sphi_phi_hist[N_centr*4], *veto_sthe_the_hist[N_centr*4], *veto_sphi_sphi_hist[N_centr*5], *veto_sthe_sthe_hist[N_centr*5];
            TH3D *couter_veto_hist, *counter_assoc_eff_hist, *counter_assoc_ghost_hist, *veto_type_hist;
            TH3D *inv_mass_dca_fg0[N_centr*3],*inv_mass_dca_fg1[N_centr*3],*inv_mass_dca_fg2[N_centr*3],*inv_mass_dca_fg3[N_centr*3],*inv_mass_dca_fg4[N_centr*3];
            TH3D *inv_mass_dca_bg0[N_centr*3],*inv_mass_dca_bg1[N_centr*3],*inv_mass_dca_bg2[N_centr*3],*inv_mass_dca_bg3[N_centr*3],*inv_mass_dca_bg4[N_centr*3];
            TH3D *delt_phi_dca_fg0[N_centr*3],*delt_phi_dca_fg1[N_centr*3],*delt_phi_dca_fg2[N_centr*3],*delt_phi_dca_fg3[N_centr*3],*delt_phi_dca_fg4[N_centr*3];
            TH3D *delt_phi_dca_bg0[N_centr*3],*delt_phi_dca_bg1[N_centr*3],*delt_phi_dca_bg2[N_centr*3],*delt_phi_dca_bg3[N_centr*3],*delt_phi_dca_bg4[N_centr*3];
            TH3D *inv_mass_dca_gen[N_centr*3];
            TH3D* myvtx_hist;
            TH3D *BBC_psi_hist, *FVTX_psi_hist, *cos_BBC_hist, *cos_FVTX_hist, *v2_BBC_hist, *v2_FVTX_hist; 
            int is_fill_hsits, is_fill_cluster_hsits, is_fill_tree, is_fill_dphi_hist, is_fill_DCA_hist, is_fill_track_QA, 
            is_fill_flow, is_fill_DCA2_hist, is_check_veto, is_fill_inv_mass;
           
      public:
            MyEventContainer()
            {
                  event = nullptr;
                  infile = nullptr;
                  outfile = nullptr;
                  hist_br = nullptr;
                  hist_bz = nullptr;
                  tree = nullptr;
                  event_hist = nullptr; centr_hist = nullptr; el_pt_hist = nullptr;
                  temc = nullptr; ttof = nullptr; n0_hist = nullptr; ep_hist=nullptr; prob_hist=nullptr; disp_hist=nullptr; chi2npe0_hist=nullptr;
                  emc_dphi_el = nullptr; emc_dz_el = nullptr; n0_hist_el = nullptr; ep_hist_el=nullptr; prob_hist_el=nullptr; disp_hist_el=nullptr; chi2npe0_hist_el=nullptr;
                  rich_prob1 = nullptr; rich_prob2 = nullptr; rich_prob3 = nullptr;
                  el_had_dphi = nullptr, el_had_dz = nullptr, el_had_dr = nullptr;
                  DCPT_ReconPT = nullptr, sDCPT_ReconPT = nullptr, charge_hist = nullptr, phi_hist = nullptr, the_hist = nullptr, pt_corr = nullptr;
                  couter_veto_hist = nullptr; counter_assoc_eff_hist = nullptr;counter_assoc_ghost_hist=nullptr, veto_type_hist = nullptr;
                  myvtx_hist = nullptr;
                  truehithist = nullptr; truehitsigmahist = nullptr;charge_recover_hist=nullptr;
                  v2_BBC_hist = nullptr; v2_FVTX_hist = nullptr; 
                  BBC_psi_hist = nullptr; FVTX_psi_hist = nullptr; cos_BBC_hist = nullptr; cos_FVTX_hist = nullptr;
                  for (int i = 0; i < N_dynamic; i++)
                  {
                        dphi_hist_el_dynamic[i] = nullptr;
                        dthe_hist_el_dynamic[i] = nullptr;
                        sdphi_hist_el_dynamic[i] = nullptr;
                        sdthe_hist_el_dynamic[i] = nullptr;
                  }
                  for (int i = 0; i < nvtx_layers; i++)
                  {
                        dphi_phi0_init_hist[i] = nullptr;
                        dthe_the0_init_hist[i] = nullptr;
                        dphi_phi0_corr_hist[i] = nullptr;
                        dthe_the0_corr_hist[i] = nullptr;
                        dphi_the0_init_hist[i] = nullptr;
                        dthe_phi0_init_hist[i] = nullptr;
                        dphi_the0_corr_hist[i] = nullptr;
                        dthe_phi0_corr_hist[i] = nullptr;
                  }
                  for (int i = 0; i < N_centr; i++)
                  {
                        dphi_hist[i] = nullptr;
                        dthe_hist[i] = nullptr;
                        sdphi_hist[i] = nullptr;
                        sdthe_hist[i] = nullptr;
                        dphi_hist_el[i] = nullptr;
                        dthe_hist_el[i] = nullptr;
                        sdphi_hist_el[i] = nullptr;
                        sdthe_hist_el[i] = nullptr;
                        chi2_ndf[i] = nullptr;
                        ilayerhitshist[i] = nullptr;
                        d_dphi_hist[i] = nullptr;
                        d_dthe_hist[i] = nullptr;
                        DCA_hist[i] = nullptr;
                        sd_dphi_hist[i] = nullptr;
                        sd_dthe_hist[i] = nullptr;
                        sDCA_hist[i] = nullptr;
                        DCA_2D_hist[i] = nullptr;
                        sDCA_2D_hist[i] = nullptr;
                        DCA2_hist[i] = nullptr;
                        sDCA2_hist[i] = nullptr;
                        DCA2_2D_hist[i] = nullptr;
                        sDCA2_2D_hist[i] = nullptr;
                        DCA12_hist[i]  = nullptr;
                        veto_phi_hist[i] = nullptr, veto_the_hist[i] = nullptr, veto_phi_phi_hist[i] = nullptr, veto_the_the_hist[i] = nullptr;
                  }
                  for (int i = 0; i < N_centr*3; i++)
                  { 
                        inv_mass_dca_fg0[i]=nullptr;inv_mass_dca_fg1[i]=nullptr;inv_mass_dca_fg2[i]=nullptr;inv_mass_dca_fg3[i]=nullptr;inv_mass_dca_fg4[i]=nullptr;
                        inv_mass_dca_bg0[i]=nullptr;inv_mass_dca_bg1[i]=nullptr;inv_mass_dca_bg2[i]=nullptr;inv_mass_dca_bg3[i]=nullptr;inv_mass_dca_bg4[i]=nullptr;
                        delt_phi_dca_fg0[i]=nullptr;delt_phi_dca_fg1[i]=nullptr;delt_phi_dca_fg2[i]=nullptr;delt_phi_dca_fg3[i]=nullptr;delt_phi_dca_fg4[i]=nullptr;
                        delt_phi_dca_bg0[i]=nullptr;delt_phi_dca_bg1[i]=nullptr;delt_phi_dca_bg2[i]=nullptr;delt_phi_dca_bg3[i]=nullptr;delt_phi_dca_bg4[i]=nullptr;
                        inv_mass_dca_gen[i]=nullptr;
                  }
                  for (int i = 0; i < N_centr*4; i++)
                  { 
                        veto_sphi_phi_hist[i] = nullptr, veto_sthe_the_hist[i] = nullptr;
                  }
                  for (int i = 0; i < N_centr*5; i++)
                  { 
                        veto_sphi_sphi_hist[i] = nullptr, veto_sthe_sthe_hist[i] = nullptr;
                  }
                  is_fill_hsits = 0;
                  is_fill_cluster_hsits = 0;
                  is_fill_tree = 0;
                  is_fill_dphi_hist = 0;
                  is_fill_DCA_hist = 0;
                  is_fill_track_QA = 0;
                  is_fill_flow = 0;
                  is_fill_DCA2_hist = 0;
                  is_check_veto = 0;
                  is_fill_inv_mass = 0;
            };

            virtual ~MyEventContainer()
            {
                  std::cout << "Starting to Delete Event in my Container" << std::endl;
                  if(is_fill_tree||is_fill_cluster_hsits||is_fill_hsits||is_fill_dphi_hist||is_fill_DCA_hist||is_fill_track_QA
                  ||is_fill_flow||is_fill_DCA2_hist||is_check_veto||is_fill_inv_mass) {delete outfile;}
                  delete event;
                  std::cout << "Event in the container was deleted" << std::endl;
            };

            void InitEvent() { event = new MyEvent; };
            void ClearEvent() { event->ClearEvent(); };
            void GetHistsFromFile(const std::string &loc);
            void CreateOutFileAndInitHists(std::string outfilename, const int fill_el = 0, const int fill_had = 0, const int fill_tree = 0, const int fill_dphi = 0, 
            const int fill_DCA = 0, const int fill_track_QA = 0, const int fill_flow = 0, const int fill_true_DCA= 0, const int check_veto= 0, const int fill_inv_mas = 0);
            void ResetTree() {tree->Reset();};
            void FillTree() {tree->Fill();};
            void WriteOutFile();
            void Associate_Hits_to_Leptons(float sigma = 2, float sigma_veto = 2, float sigma_inner = 2, int not_fill = 0);
            void Associate_Hits_to_Leptons_OLD(float sigma = 2, float sigma_veto = 2, float sigma_inner = 2, int not_fill = 0);
            void Associate_Hits_to_Clusters(float sigma = 2);

            void Reveal_Cluster();

            void FillDphiHists();
            void FillTrueDCA();

            MyEvent *GetEvent() { return event; };
            void SetEvent(DileptonAnalysis::MyEvent *ev) { event = ev; };

            void CheckVeto();
            void FillEventHist(const float varX){ event_hist->Fill(varX); };
            void FillCentrHist(const float centr){ centr_hist->Fill(centr); };
            int GetNGoodElectrons();

            void fill_evtbuff_list(const unsigned int pool_depth = 10);
            void fill_inv_mass(const unsigned int pool_depth = 10);
            void fill_inv_mass_sim();
            void correct_beam_offset();
            void CleanUpHitList();
            void FillQAHist(const int mc_id = -999);
            void FillFlow(const float psi_BBCS=-999, const float psi_BBCN=-999, const float psi_FVTXS=-999, const float psi_FVTXN=-999);

            void AddBDTHit(const MyBDTrack *newBDTrack) { BDTracklist.push_back(*newBDTrack); };
            Long64_t GetNBDThit() { return BDTracklist.size(); };
            DileptonAnalysis::MyBDTrack *GetBDTHitEntry(unsigned int i) const { return const_cast<DileptonAnalysis::MyBDTrack *>(&(BDTracklist[i])); };
            void RemoveBDTHitEntry(const unsigned int i){ BDTracklist.erase(BDTracklist.begin() + i); };    

            ClassDef(MyEventContainer, 1) // MyEvent structure
      };
}

#endif /*__MYEVENT_H__*/
