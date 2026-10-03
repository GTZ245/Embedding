#include <iostream>
#include <math.h>

#include "SvxCheckCompactCNT.h"

#include "phool.h"
#include "PHTypedNodeIterator.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "Fun4AllReturnCodes.h"

#include "SvxClusterList.h"
#include "SvxCluster.h"
#include "SvxClusterv4.h"
#include "SvxSegmentList.h"
#include "SvxSegmentv6.h"
#include "SvxCentralTrackList.h"
#include "SvxCentralTrack.h"
#include "SvxCentralTrackv5.h"

#include "compactCNT/SvxHitMap.h"
#include "compactCNT/SvxHitMapEntry.h"
#include "compactCNT/SvxTrackMap.h"
#include "compactCNT/SvxTrackMapEntry.h"
#include "compactCNT/SvxCentralTrackMap.h"
#include "compactCNT/SvxCentralTrackMapEntry.h"

#include "getClass.h"

#include "TString.h"

using namespace std;
using namespace findNode;

//==============================================================

SvxCheckCompactCNT::SvxCheckCompactCNT(const string& name) : SubsysReco(name) 
{
  init_ana=0;
  EventNumber=0;
  m_print=false; // 
}

//==============================================================

SvxCheckCompactCNT::~SvxCheckCompactCNT() {
}

//==============================================================

int SvxCheckCompactCNT::Init(PHCompositeNode *topNode) {

  cout << "SvxCheckCompactCNT::Init started..." << endl;



  cout << "SvxCheckCompactCNT::Init ended." << endl;
  return 0;
}

//==============================================================
  
int SvxCheckCompactCNT::InitRun(PHCompositeNode *topNode) 
{
  cout << "SvxCheckCompactCNT::InitRun started..." << endl;


 // int i = CreateNodeTree(topNode);
 // if(verbosity>0) cout << "SvxCheckCompactCNT::InitRun-I: CreateNodeTree returned " << i << endl;
 // if(!(i==EVENT_OK)) {return EVENT_OK;}


  topNode->print();

  cout << "SvxCheckCompactCNT::InitRun ended." << endl;
  return 0;
}

//==============================================================


int SvxCheckCompactCNT::process_event(PHCompositeNode *topNode) 
{
  cout<<"EventNumber : "<<EventNumber<<endl;

  cout<<"compareClusterList"<<endl;
  compareClusterList(topNode);

  cout<<"compareSegmentList"<<endl;
  compareSegmentList(topNode);

  cout<<"compareCentralTrackList"<<endl;
  compareCentralTrackList(topNode);
 
  EventNumber++;
  return 0;
}

bool SvxCheckCompactCNT::compareClusterList(PHCompositeNode *topNode){


  SvxClusterList *clslist = getClass<SvxClusterList>(topNode, "SvxSelectedClusterList");
  if(clslist==NULL){
    cout<<"Error SvxCheckCompactCNT::compareClusterList no SvxSelectedClusterList"<<endl;
    return false;
  }

  SvxHitMap *svxmap  = findNode::getClass<SvxHitMap>(topNode, "SvxHit_comp");
  if(svxmap==NULL){
    cout<<"Error SvxCheckCompactCNT::compareClusterList no SvxHitMap"<<endl;
    return false;
  }

  int n_org = clslist->get_nClusters();
  int n_new = svxmap->GetNentries();

  bool result = (n_org==n_new);
  if(!result){
    cout<<"Inconsistent SvxCheckCompactCNT::compareClusterList Ncluster : "<<n_org<<"!="<<n_new<<endl;
  }
  else {
    
    for(int iclus=0; iclus<n_org; iclus++){
      SvxCluster      *cls_org   = clslist->get_Cluster(iclus);

      float r = 0.0;
      SvxClusterv4 cls_new;
      {
        const SvxHitMapEntry *svxmapent = svxmap->GetHit(iclus);
        cls_new.set_hitID(svxmapent->get_id());
        float x_new = svxmapent->get_x() * 1e-3; // 10mkm
        float y_new = svxmapent->get_y() * 1e-3; // 10mkm
        float z_new = svxmapent->get_z() * 1e-3; // 10mkm
        cls_new.set_xyz_global(0, x_new);
        cls_new.set_xyz_global(1, y_new);
        cls_new.set_xyz_global(2, z_new);

        //float r = sqrt(x_new*x_new+y_new*y_new);
        r = sqrt(x_new*x_new+y_new*y_new);
        int layer = 0;

        if     (r<4.0)  layer=0;
        else if(r<7.0)  layer=1;
        else if(r<14.0) layer=2;
        else            layer=3;

        cls_new.set_layer(layer);
        if(layer<2){
          cls_new.set_size(svxmapent->get_adcandsize());
        } else {
          int x_adc = svxmapent->get_adcandsize()&0xFF;
          int z_adc = (svxmapent->get_adcandsize()>>8)&0xFF;
          cls_new.set_adc(0, x_adc);
          cls_new.set_adc(1, z_adc);
        }
      }

      result &= compareCluster(cls_org, &cls_new);
      if(!result){
        cout<<"Inconsistent SvxCheckCompactCNT::compareClusterList ClusterID : "<<iclus<<endl;
        cout<<"  r : "<<r<<endl;
        break; 
      }
    }
  }

  cout<<"Checked SexCheckCompactCNT::compareClusterList Ncluster : "<<n_org<<" "<<flush;
  cout<<(result ? "OK" : "FAILED")<<endl;
  cout<<endl;
  

  return result;
}

bool SvxCheckCompactCNT::compareSegmentList(PHCompositeNode *topNode){
  SvxSegmentList *seglist = getClass<SvxSegmentList>(topNode, "SvxSegmentList");
  if(seglist==NULL){
    cout<<"Error SvxCheckCompactCNT::process_event no SvxSegmentList"<<endl;
    return false;
  }

  SvxTrackMap *svxmap  = findNode::getClass<SvxTrackMap>(topNode, "SvxTrack_comp");
  if(svxmap==NULL){
    cout<<"Error SvxCheckCompactCNT::compareSegmentList no SvxTrackMap"<<endl;
    return false;
  }

  int n_org = seglist->get_nSegments();
  int n_new = svxmap->GetNentries();
  bool result = (n_org==n_new);
  if(!result){
    cout<<"Inconsistent SvxCheckCompactCNT::compareSegmentList Nsegment : "<<n_org<<"!="<<n_new<<endl;
  }
  else {
    
    for(int iseg=0; iseg<n_org; iseg++){
      SvxSegment      *seg_org   = seglist->get_segment(iseg);

      SvxSegmentv6 seg_new;
      {
        const SvxTrackMapEntry *svxmapent = svxmap->GetHit(iseg);

        seg_new.setIsPositive(svxmapent->get_charge()==1);
//        seg_new.setNhits(     svxmapent->get_nhits()); // nhit is not recoverable
        seg_new.setQuality(   svxmapent->get_quality() * 1e-2);
        seg_new.setDCA2D(     svxmapent->get_dca2d() *1e-4);
        seg_new.setDCA(       svxmapent->get_dca3d() *1e-4);
        seg_new.setClosestApproach(svxmapent->get_x()*1e-4, 
                                   svxmapent->get_y()*1e-4,
                                   svxmapent->get_z()*1e-3);
        seg_new.set3MomentumAtPrimaryVertex(svxmapent->get_px()*1e-3, 
                                            svxmapent->get_py()*1e-3,
                                            svxmapent->get_pz()*1e-3);
        seg_new.setLivePercentage(0, svxmapent->get_livePercentage(0)*1e-2);
        seg_new.setLivePercentage(1, svxmapent->get_livePercentage(1)*1e-2);
        seg_new.setLivePercentage(2, svxmapent->get_livePercentage(2)*1e-2);
        seg_new.setLivePercentage(3, svxmapent->get_livePercentage(3)*1e-2);
        seg_new.setSegmentQuality(svxmapent->get_segmentQuality()*1e-4);
        seg_new.setSegmentScore(  svxmapent->get_segmentScore()*1e-2);
        seg_new.set_dEdX1(  svxmapent->get_dEdX1()*1e-1 );
        seg_new.set_dEdX2(  svxmapent->get_dEdX2()*1e-1 );
        for(int ilay=0; ilay<4; ilay++) { 
          for(int ihit=0; ihit<2; ihit++) { 
            seg_new.setClusterGoodFraction(ilay, ihit, svxmapent->get_clusterGoodFraction(ilay, ihit)*1e-4 );
          }
        }
      }

      result &= compareSegment(seg_org, &seg_new);
      if(!result){
        cout<<"Inconsistent SvxCheckCompactCNT::compareSegmentList SegmentID : "<<iseg<<endl;
        break; 
      }
    }
  }

  cout<<"Checked SexCheckCompactCNT::compareSegmentList Nsegment : "<<n_org<<" "<<flush;
  cout<<(result ? "OK" : "FAILED")<<endl;
  cout<<endl;

  return result;
}


bool SvxCheckCompactCNT::compareCentralTrackList(PHCompositeNode *topNode){
  SvxCentralTrackList *cntlist = getClass<SvxCentralTrackList>(topNode, "SvxCentralTrackList");
  if(cntlist==NULL){
    cout<<"Error SvxCheckCompactCNT::process_event no SvxCentralTrackList"<<endl;
    return false;
  }

  SvxCentralTrackMap *svxmap  = findNode::getClass<SvxCentralTrackMap>(topNode, "SvxCentralTrack_comp");
  if(svxmap==NULL){
    cout<<"Error SvxCheckCompactCNT::compareCentralTrackList no SvxCentralTrackMap"<<endl;
    return false;
  }

  SvxHitMap *svxhitmap  = findNode::getClass<SvxHitMap>(topNode, "SvxHit_comp");
  if(svxhitmap==NULL){
    cout<<"Error SvxCheckCompactCNT::compareCentralTrackList no SvxHitMap"<<endl;
    return false;
  }

  int n_org = cntlist->get_nCentralTracks();
  int n_new = svxmap->GetNentries();
  bool result = (n_org==n_new);
  if(!result){
    cout<<"Inconsistent SvxCheckCompactCNT::compareCentralTrackList Ntrack : "<<n_org<<"!="<<n_new<<endl;
  }
  else {
    
    for(int itrk=0; itrk<n_org; itrk++){
      SvxCentralTrack      *scnt_org   = cntlist->getCentralTrack(itrk);

      SvxCentralTrackv5 scnt_new;
      {
        const SvxCentralTrackMapEntry *svxmapent = svxmap->GetHit(itrk);

        scnt_new.setDchIndex(svxmapent->get_DchIndex());
        scnt_new.setUnique(  svxmapent->get_Unique());
        scnt_new.setDCA2D(   svxmapent->get_DCA2D()*1e-4);
        scnt_new.setDCAZ(    svxmapent->get_DCAZ()*1e-4);
        scnt_new.setClosestApproach(svxmapent->get_ClosestApproach(0)*1e-4,
                                    svxmapent->get_ClosestApproach(1)*1e-4,
                                    svxmapent->get_ClosestApproach(2)*1e-3);
        scnt_new.setChiSquareDPHI(  svxmapent->get_ChisquarePhi()*1e-2);
        scnt_new.setChiSquareDZ(    svxmapent->get_ChisquareZ()*1e-2);
        scnt_new.setChiSquare(      svxmapent->get_Chisquare()*1e-2);
        scnt_new.setChiSquare2(     svxmapent->get_Chisquare2()*1e-2);

        short pattern = svxmapent->get_HitPattern();
        

        int sublay[8] = {-1,-1,-1,-1,-1,-1,-1,-1};
        int nhitpt = 0;
        for(int ibit=0; ibit<8; ibit++){
          if( ((pattern>>(7-ibit))&0x1)==0x1 ){
            sublay[nhitpt] = 7-ibit;
            //sublay[ibit] = 7-ibit;
            nhitpt++;
          }
        }
/*
        cout<<"hitpattern : 0x"<<hex<<pattern<<dec<<" ";
        for(int ibit=0; ibit<8; ibit++){
            cout<<sublay[ibit]<<" ";
        }
        cout<<endl;
*/
       
        int nhit = svxmapent->get_NClusters();
 //       if(nhit!=nhitpt){
 //         cout<<" hit is inconsistent : "<<nhit<<"!="<<nhitpt<<endl;
 //       }

        scnt_new.setNDF(2*nhit);

        for(int icls=0; icls<nhit; icls++){
          SvxClusterInfov3 info;
          info.setClusterId(svxmapent->get_ClusterID(icls));
          info.setdproj(    svxmapent->get_ClusterDPhi(icls)*1e-4);
          info.setbend(     0.0);
          info.setzproj(    svxmapent->get_ClusterDZ(icls) *1e-4);
          info.setPosition( 0,0,0); // need to be z=0 for dz calculation 

          // dummy layer and ladder
          int layer = 0, ladder=0;
          if(sublay[icls]<2)      { layer = sublay[icls];}
          else if(sublay[icls]<5) { layer = 2; ladder = 7-(sublay[icls]-2); }
          else                    { layer = 3; ladder = sublay[icls]-5; }
          info.setLayer(layer);
          info.setLadder(ladder);

          scnt_new.addClusterInfo(&info);
        }
       
      }

      result &= compareCentralTrack(scnt_org, &scnt_new);
      if(!result){
        cout<<"Inconsistent SvxCheckCompactCNT::compareCentralTrackList CentralTrackID : "<<itrk<<endl;
        break; 
      }
    }
  }

  cout<<"Checked SexCheckCompactCNT::compareCentralTrackList Ntrack : "<<n_org<<" "<<flush;
  cout<<(result ? "OK" : "FAILED")<<endl;
  cout<<endl;

  return result;
}

//==============================================================

int SvxCheckCompactCNT::End(PHCompositeNode *topNode) {
  cout << "SvxCheckCompactCNT::End:  Writing out..." << endl;
  return 0;
}

//==============================================================
int SvxCheckCompactCNT::CreateNodeTree(PHCompositeNode* topNode)
{

/*
  PHNodeIterator iter(topNode);

  // Find SVX node.
  PHCompositeNode* svxNode = dynamic_cast<PHCompositeNode*> (iter.findFirst("PHCompositeNode", "SVX"));
  if (! svxNode)
  {
    cerr << PHWHERE << "SVX node missing, doing nothing." << endl;
    return ABORTEVENT;
  }


  PHIODataNode<PHObject> *SvxClusterListOrgNode = 
     (PHIODataNode<PHObject>*)iter.findFirst("PHIODataNode", "SvxSelectedClusterListOrg");
  if(!SvxClusterListOrgNode)
  {
    SvxClusterList* svxcluster = new SvxClusterListv4();
    SvxClusterListOrgNode = new PHIODataNode<PHObject>(svxcluster, "SvxSelectedClusterListOrg", "PHObject");
    svxNode->addNode(SvxClusterListOrgNode);
  }
*/

  return EVENT_OK;
}

bool SvxCheckCompactCNT::compareCluster(SvxCluster* corg, SvxCluster* cnew){
  if(corg==NULL||cnew==NULL){
    cerr<<"SvxCheckCompactCNT::compareCluster"<<endl;
    return false;
  }

  bool result = true;

  result &= compareInt(corg->get_hitID(),cnew->get_hitID(), "hitID");
//  result &= compareInt(corg->get_svxSection(), cnew->get_svxSection(), "svxSection");
  result &= compareInt(corg->get_layer()     , cnew->get_layer()     , "layer");
//  result &= compareInt(corg->get_ladder()    , cnew->get_ladder()    , "ladder");
//  result &= compareInt(corg->get_sensor()    , cnew->get_sensor()    , "sensor");

//  result &= compareInt(corg->get_sensorType(), cnew->get_sensorType(), "sensorType");
//  result &= compareInt(corg->get_edgeflag  (), cnew->get_edgeflag  (), "edgeFlag");

//  result &= compareInt(corg->get_ambiguous ()          , cnew->get_ambiguous(), "ambiguous");
//  result &= compareInt(corg->get_AssociatedCGL()       , cnew->get_AssociatedCGL(), "AssociatedCLG");
//  result &= compareInt(corg->get_AssociatedStandalone(), cnew->get_AssociatedStandalone(), "AssociatedCLG");
//  result &= compareInt(corg->get_circumference()       , cnew->get_circumference(), "circumference");


  if(corg->get_layer()<2){ // size check for layer 0,1
    result &= compareInt(corg->get_size()                , cnew->get_size(),      "size");
  }
  else { // adc check for layer 2, 3
    for(int i=0; i<2; i++){
      float adc_org = corg->get_adc(i)<256 ? corg->get_adc(i) : 255;
      result &= compareInt(adc_org, cnew->get_adc     (i), "adc");
    //  result &= compareInt(corg->get_xz_size (i), cnew->get_xz_size (i), "xz_size");
    }
  }
  for(int i=0; i<3; i++){ 
//    result &= compareFloat(corg->get_xyz_local (i), cnew->get_xyz_local (i), "XYZ_Local");
    result &= compareFloat(corg->get_xyz_global(i), cnew->get_xyz_global(i), 0.001, "XYZ_Global");
  }

  return result;
}

bool SvxCheckCompactCNT::compareSegment(SvxSegment* corg, SvxSegment* cnew){
  if(corg==NULL||cnew==NULL){
    cerr<<"SvxCheckCompactCNT::compareSegment"<<endl;
    return false;
  }

  bool result = true;

  result &= compareInt(corg->IsPositive(),   cnew->IsPositive(),        "IsPositive");
//  result &= compareInt(corg->getNHits()  ,   cnew->getNHits(),          "NHits");

  result &= compareFloat(corg->getQuality(), cnew->getQuality(),0.01,   "Quality");
  result &= compareFloat(corg->getDCA2D(),   cnew->getDCA2D(),  0.0001, "DCA2D");
  result &= compareFloat(corg->getDCA(),     cnew->getDCA(),    0.0001, "DCA");

  for(int i=0; i<3; i++){ 
    float limit = (i<2) ? 0.0001 : 0.001;
    result &= compareFloat(corg->getClosestApproach(i), cnew->getClosestApproach(i), 
                           limit, TString::Format("ClosestApproach %d", i));
  }
  for(int i=0; i<3; i++){ 
    result &= compareFloat(corg->get3MomentumAtPrimaryVertex(i), cnew->get3MomentumAtPrimaryVertex(i),
                           0.001, TString::Format("MomentumAtPrimaryVertex %d", i));
  }

  for(int i=0; i<3; i++){ 
    if(corg->getLivePercentage(i)<-9900){ // uninit value
      result &= compareFloat(cnew->getLivePercentage(i), -320.,
                             0.01, TString::Format("LivePercentage %d", i));
    } 
    else {
      result &= compareFloat(corg->getLivePercentage(i), cnew->getLivePercentage(i),
                             0.01, TString::Format("LivePercentage %d", i));
    }
  }
  result &= compareFloat(corg->getSegmentQuality(), cnew->getSegmentQuality(), 0.01, "SegmentQuality");
  result &= compareFloat(corg->getSegmentScore(),   cnew->getSegmentScore(),   0.01, "SegmentScore");

  result &= compareFloat(corg->get_dEdX1(),   cnew->get_dEdX1(),   0.1, "dEdX1");
  result &= compareFloat(corg->get_dEdX2(),   cnew->get_dEdX2(),   0.1, "dEdX2");
  for(int ilay=0; ilay<4; ilay++){ 
    for(int ihit=0; ihit<2; ihit++){ 
      if(corg->getClusterGoodFraction(ilay, ihit)<-9900){ // uninit value
        result &= compareFloat(cnew->getClusterGoodFraction(ilay, ihit), -3.2,
                               0.0001, TString::Format("ClusterGoodFrac %d %d", ilay, ihit));
      } 
      else {
        result &= compareFloat(corg->getClusterGoodFraction(ilay, ihit), cnew->getClusterGoodFraction(ilay, ihit),
                               0.0001, TString::Format("ClusterGoodFrac %d %d", ilay, ihit));
      }
    }
  }

  return result;
}

bool SvxCheckCompactCNT::compareCentralTrack(SvxCentralTrack* strk_org, SvxCentralTrack* strk_new){
  if(strk_org==NULL||strk_new==NULL){
    cerr<<"SvxCheckCompactCNT::compareSegment"<<endl;
    return false;
  }

  bool result = true;

  result &= compareInt(  strk_org->getDchIndex(),      strk_new->getDchIndex(),            "DchIndex");
  result &= compareInt(  strk_org->getUnique(),        strk_new->getUnique(),              "Unique");
  result &= compareFloat(strk_org->getDCA2D(),         strk_new->getDCA2D(),         1e-4, "DCA2D");
  result &= compareFloat(strk_org->getDCAZ(),          strk_new->getDCAZ(),          1e-4, "DCAZ");
  result &= compareFloat(strk_org->getChiSquareDPHI(), strk_new->getChiSquareDPHI(), 1e-2, "ChiSquareDPHI");
  result &= compareFloat(strk_org->getChiSquareDZ(),   strk_new->getChiSquareDZ(),   1e-2, "ChiSquareDZ");
  result &= compareFloat(strk_org->getChiSquare(),     strk_new->getChiSquare(),     1e-2, "ChiSquare");
  result &= compareFloat(strk_org->getChiSquare2(),    strk_new->getChiSquare2(),    1e-2, "ChiSquare2");

  for(int i=0; i<3; i++){ 
    float limit = (i<2) ? 0.0001 : 0.001;
    result &= compareFloat(strk_org->getClosestApproach(i), strk_new->getClosestApproach(i), 
                           limit, TString::Format("ClosestApproach %d", i));
  }
 
  int nhit_org = strk_org->getNhits();
  int nhit_new = strk_new->getNhits();
  result &= compareInt(nhit_org, nhit_new, "Nhits");

  for(int icls=0; icls<nhit_org; icls++){ 
    SvxClusterInfo *cls_org = strk_org->getClusterInfo(icls);
    SvxClusterInfo *cls_new = strk_new->getClusterInfo(icls);

    result &= compareInt(cls_org->getClusterId(), cls_new->getClusterId(), 
                         TString::Format("ClusterId %d", icls));

    result &= compareInt(cls_org->get_sublayer(), cls_new->get_sublayer(), 
                         TString::Format("sublayer %d", icls));

    result &= compareFloat(cls_org->getdphi(), cls_new->getdphi(), 
                           1e-4, TString::Format("dphi %d", icls));
    result &= compareFloat(cls_org->getdz(), cls_new->getdz(), 
                           1e-4, TString::Format("dz %d", icls));

  }

  return result;
}

bool SvxCheckCompactCNT::compareInt(int orgval, int newval, const char* err){
  bool result = (orgval==newval);
  if(!result) {cout<<"Faild : "<<err<<" "<<orgval<<"!="<<newval<<endl;}
  if(m_print) {cout<<err<<"  "<<orgval<<"=="<<newval<<"  "<<(result ? "OK" : "FAIL")<<endl;}
  return result;
}

bool SvxCheckCompactCNT::compareFloat(float orgval, float newval, float judge, const char* err){
  bool result = (fabs(orgval-newval)<judge);
  if(!result) {cout<<"Faild : "<<err<<" "<<orgval<<"!="<<newval<<endl;}
  if(m_print) {cout<<err<<"  "<<orgval<<"=="<<newval<<" "<<judge<<"  "<<(result ? "OK" : "FAIL")<<endl;}
  return result;
}

