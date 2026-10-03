#include <iostream>
#include <string>

#include "phool.h"
#include "PHTypedNodeIterator.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "Fun4AllReturnCodes.h"

#include "getClass.h"

#include "DebugTest.h"

#include "SvxClusterList.h"
#include "SvxCluster.h"
#include "SvxCentralTrackList.h"
#include "SvxCentralTrack.h"

using namespace std;
using namespace findNode;
//==============================================================

DebugTest::DebugTest() {
  ThisName = "DebugTest";
}

//==============================================================

DebugTest::~DebugTest() {
}

int DebugTest::Init(PHCompositeNode *topNode) {

  cout << "DebugTest::Init started..." << endl;

  cout << "DebugTest::Init ended." << endl;
  return 0;
}

//==============================================================
  
int DebugTest::InitRun(PHCompositeNode *topNode) {
  cout << "DebugTest::InitRun started..." << endl;

  cout << "DebugTest::InitRun ended." << endl;
  return 0;
}

//==============================================================


int DebugTest::process_event(PHCompositeNode *topNode) {

  SvxClusterList *svxsellist = getClass<SvxClusterList>(topNode, "SvxSelectedClusterList");
  if(svxsellist==NULL){
    cout<<"SvxSelectedClusterList is NULL"<<endl;
    return 0;
  }
  int nclus = svxsellist->get_nClusters();

  SvxCentralTrackList *svxcntlist = getClass<SvxCentralTrackList>(topNode, "SvxCentralTrackList");
  if(svxcntlist==NULL){
    cout<<"SvxCentralTrackList is NULL"<<endl;
    return 0;
  }
  SvxCentralTrackList *svxcntbglist = getClass<SvxCentralTrackList>(topNode, "SvxCentralTrackBackList");
  if(svxcntbglist==NULL){
    cout<<"SvxCentralTrackBackList is NULL"<<endl;
    return 0;
  }
  
  cout<<"Nclusters : "<<nclus<<endl;
  if(nclus>0){
    SvxCluster *cls = svxsellist->get_Cluster(0);
    cout<<"   x:y:z = ";
    cout<<cls->get_xyz_global(0)<<" ";
    cout<<cls->get_xyz_global(1)<<" ";
    cout<<cls->get_xyz_global(2)<<endl;
  }

  /////////
  testSelectedCluster(svxsellist, svxcntlist, svxcntbglist);


/*
  m_svxrawhitclus = getClass<SvxRawhitClusterList>(topNode,"SvxRawhitClusterList");

  int nsvx       = m_svx->get_nClusters();
  int nsvxrawhit = m_svxrawhit->get_nRawhits();

  if(init_ana==0) {
    cout<<"init ana"<<endl;
    topNode->print(); 
    init_ana=1;
 
    if(m_svx!=NULL)          m_svx->identify();          else cout<<"no SvxCluster object"<<endl; 
    if(m_svxrawhitclus!=NULL)m_svxrawhitclus->identify();else cout<<"no SvxRawhitClus object"<<endl;
    if(m_svxrawhit!=NULL)    m_svxrawhit->identify();    else cout<<"no SvxRawhit object"<<endl;

    cout<<endl<<endl;
    cout<<"init ana ends"<<endl;
  }

  if(EventNumber%1==0) {
    cout << "------- Event # " << EventNumber << " nsvxrawhit= " << nsvxrawhit << " nsvx= " << nsvx << endl;
  }
*/

  return 0;
}

bool DebugTest::testSelectedCluster(SvxClusterList      *svxsellist, 
                                    SvxCentralTrackList *svxcntlist,
                                    SvxCentralTrackList *svxcntbglist)
{
  if(svxsellist==NULL ||svxcntlist==NULL ){
    cout<<"No object : SvxSelectedClusterList, SvxCentralTrackList"<<endl;
    return false;
  }

  bool result=true;

  // SvxCentralTrack
  int nsvxcnt=svxcntlist->get_nCentralTracks();
  for(int isvxcnt=0; isvxcnt<nsvxcnt; isvxcnt++){
    SvxCentralTrack *svxcnt = svxcntlist->getCentralTrack(isvxcnt);

    int nhits = svxcnt->getNhits();
    int hitptn = svxcnt->getHitPattern();
    int nhits2=0;
    for(int i=0; i<8; i++){
      if(((hitptn>>i)&0x1)==0x1) nhits2++;
    }
    //cout<<"itrk : "<<isvxcnt<<", "<<nhits<<" "<<nhits2<<endl;
    int nfound=0;
    for(int ihit=0; ihit<nhits; ihit++){
      SvxClusterInfo *cls = svxcnt->getClusterInfo(ihit);
      if(cls!=NULL){
        int id = cls->getClusterId();

        // check if the id is included
        bool found = false;
        for(int i=0; i<svxsellist->get_nClusters();  i++){
          SvxCluster *cluster = svxsellist->get_Cluster(i);
          int id_cls = (cluster!=NULL) ? cluster->get_hitID() : -1;
          if(id==id_cls) { found=true; nfound++; break; }
        }
        if(!found){
          cout<<"Not included : id="<<id<<endl;
          result=false;
          break;
        }
        else {
        //  cout<<"Cluster found : id= "<<id<<endl;
        }
      }
    }
    if(!result){
      cout<<" result is false : "<<isvxcnt<<endl;
      break;
    }
    else {
      cout<<"ncluster : nfound = "<<nhits<<" : "<<nfound<<" ";
      cout<< ((nhits==nfound) ? "Pass" : "Fail") <<endl;
    }
  }

  // SvxCentralTrackBack
  cout<<"SvxCentralTrackBack"<<endl;
  bool result2=true;
  int nsvxcntbg=svxcntbglist->get_nCentralTracks();
  for(int isvxcnt=0; isvxcnt<nsvxcntbg; isvxcnt++){
    SvxCentralTrack *svxcnt = svxcntbglist->getCentralTrack(isvxcnt);

    int nhits = svxcnt->getNhits();
    int hitptn = svxcnt->getHitPattern();
    int nhits2=0;
    for(int i=0; i<8; i++){
      if(((hitptn>>i)&0x1)==0x1) nhits2++;
    }
    //cout<<"itrk : "<<isvxcnt<<", "<<nhits<<" "<<nhits2<<endl;
    int nfound=0;
    for(int ihit=0; ihit<nhits; ihit++){
      SvxClusterInfo *cls = svxcnt->getClusterInfo(ihit);
      if(cls!=NULL){
        int id = cls->getClusterId();

        // check if the id is included
        bool found = false;
        for(int i=0; i<svxsellist->get_nClusters();  i++){
          SvxCluster *cluster = svxsellist->get_Cluster(i);
          int id_cls = (cluster!=NULL) ? cluster->get_hitID() : -1;
          if(id==id_cls) { found=true; nfound++; break; }
        }
        if(!found){
          cout<<"Not included : id="<<id<<endl;
          result2=false;
          break;
        }
        else {
        //  cout<<"Cluster found : id= "<<id<<endl;
        }
      }
    }
    if(!result2){
      cout<<" result is false : "<<isvxcnt<<endl;
      break;
    }
    else {
      cout<<"ncluster : nfound = "<<nhits<<" : "<<nfound<<" ";
      cout<< ((nhits==nfound) ? "Pass" : "Fail") <<endl;
    }
  }

  return result&result2;
}

//==============================================================

int DebugTest::End(PHCompositeNode *topNode) {
  cout << "DebugTest::End:  Writing out..." << endl;
  cout << "DebugTest::End:  Closing output file..." << endl;
  return 0;
}

