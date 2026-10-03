// $Id: embed_IOManager.C,v 1.2 2012/06/01 15:36:25 lebedev Exp $
/*!
   \file pisaToDST_IOManager.C
   \brief simulation output managers
   \author <a href="mailto:pereira@hep.saclay.cea.fr">Hugo Pereira</a>
   \version $Revision: 1.2 $
   \date $Date: 2012/06/01 15:36:25 $
*/

#include <stdio.h> 
#include <time.h> 

//_______________________________________________________________
void DST_IOManager(const char* dstout, Fun4AllServer* se)
{

  // Control for simulated DST
  Fun4AllDstOutputManager *dst_manager  = new Fun4AllDstOutputManager("EMBEDDST", dstout);
 
  recoConsts* rc = recoConsts::instance();

/*
  // run and event header
  dst_manager->AddNode("RunHeader");	
  dst_manager->AddNode("EventHeader");
  
  // pisa nodes
  dst_manager->AddNode("fkin");
  dst_manager->AddNode("pythia");
  dst_manager->AddNode("primary");
  dst_manager->AddNode("header");

  // global vertex and t0
  dst_manager->AddNode("T0Out");
  dst_manager->AddNode("VtxOut");
  
  // BBC
  dst_manager->AddNode("bbcghit");
  dst_manager->AddNode("BbcOut");
  dst_manager->AddNode("BbcRaw");

  // ZDC
  dst_manager->AddNode("ZdcOut");
  dst_manager->AddNode("ZdcRaw");
  
  
  // DCH
  dst_manager->AddNode("dcghit");
  dst_manager->AddNode("DchHitLineTablev1");
  dst_manager->AddNode("DchHitLineTable");
  dst_manager->AddNode("DchTrack");
  dst_manager->AddNode("dDchTracksPerf");
  dst_manager->AddNode("dDchTracksExtPerf");
  dst_manager->AddNode("dDchGhitHits");

  // PC
  dst_manager->AddNode("pc1ghit");
  dst_manager->AddNode("pc2ghit");
  dst_manager->AddNode("pc3ghit");
  dst_manager->AddNode("dPc1GhitClus");
  dst_manager->AddNode("dPc2GhitClus");
  dst_manager->AddNode("dPc3GhitClus");
  dst_manager->AddNode("Pc1Cluster");
  dst_manager->AddNode("Pc2Cluster");
  dst_manager->AddNode("Pc3Cluster");
  dst_manager->AddNode("Pc1Raw");
  dst_manager->AddNode("Pc2Raw");
  dst_manager->AddNode("Pc3Raw");

  // RICH
  dst_manager->AddNode("crkghit");
  dst_manager->AddNode("CrkHit");
  
  // these nodes, when added to the output
  // makes the read-back of the DSTs crash
  // it is removed for the moment
  // dst_manager->AddNode("CrkRing");
  // dst_manager->AddNode("CrkRingBack");

  // TEC
  //dst_manager->AddNode("tecghit");
  //dst_manager->AddNode("dTecGhitRaw");
  //dst_manager->AddNode("TecOut");
  //dst_manager->AddNode("TecHitOut");

  // TOF
  //dst_manager->AddNode("tofghit");
  //dst_manager->AddNode("TofOut");
  //dst_manager->AddNode("dTofGdigi");
  //dst_manager->AddNode("dTofGhitGdigi");
  //dst_manager->AddNode("dTofGdigiRec");
  
  // TOF West
  //dst_manager->AddNode("TofwRaw");
  //dst_manager->AddNode("TofwHit");

  // Aerogel
  //dst_manager->AddNode("AerGeaHits");
  //dst_manager->AddNode("AccCluster");
  //dst_manager->AddNode("AccRaw");
  //dst_manager->AddNode("AccHit");

  // EMCal
  dst_manager->AddNode("emcghit");
  dst_manager->AddNode("emcClusterContainer");
  dst_manager->AddNode("emcTowerContainer");
  
  // additional EMCal evaluation nodes
  if( rc->FlagExist("EVALUATIONFLAG") && rc->get_IntFlag("EVALUATIONFLAG")==1 )
  { 

    // Evaluation output from EMCal
    dst_manager->AddNode("dEmcGeaClusterTrack");
    dst_manager->AddNode("dEmcGeaTrack");
    dst_manager->AddNode("dEmcGeaTrackCluster");
  
  } 

  // CGL
  dst_manager->AddNode("CglTrack");
  dst_manager->AddNode("CglTrackBack");
  dst_manager->AddNode("PHDchTrackOut");
  dst_manager->AddNode("PHTrackOut");
  dst_manager->AddNode("PHTrackOutBack");
  
  // copied from CNT node
  dst_manager->AddNode("PHCentralTrack");
  dst_manager->AddNode("PHGlobal");
  dst_manager->AddNode("PHGlobal_CENTRAL");
  dst_manager->AddNode("PHGlobal_MUON");

  dst_manager->AddNode("McSingle");
  
  // muon nodes
  //dst_manager->AddNode("TMCPrimary");
  //dst_manager->AddNode("TMuiMCHitO");
  //dst_manager->AddNode("TMutMCHit");
  //dst_manager->AddNode("TMutMCTrk");
*/
 
  dst_manager->AddNode("fkin");
/*
  dst_manager->AddNode("eepair");
  dst_manager->AddNode("eepair_bg");
  dst_manager->AddNode("GammaObj");
  dst_manager->AddNode("Pi0Obj");
*/
  dst_manager->AddNode("Sync");
  dst_manager->AddNode("emcClusterContainer");
  dst_manager->AddNode("PHCentralTrack");
  dst_manager->AddNode("SvxRawhitList");
  dst_manager->AddNode("SvxClusterList");
  dst_manager->AddNode("SvxGhitClusterList");

  dst_manager->AddNode("PHGlobal");
  dst_manager->AddNode("McSingle");
  dst_manager->AddNode("PHEmbedMcRecoTrackList");
  dst_manager->AddNode("PHEmbedMcRecoTrack");
  dst_manager->AddNode("PHEmbededEvent");
  dst_manager->AddNode("McEvalSingleList");
  dst_manager->AddNode("dEmcGeaClusterTrack");
  dst_manager->AddNode("dEmcGeaTrack");
  dst_manager->AddNode("dEmcGeaTrackCluster");
  
  dst_manager->AddNode("CrkHit");
  dst_manager->AddNode("CrkRing");
  dst_manager->AddNode("CrkRingBack");
  dst_manager->AddNode("Pc1Cluster");
  dst_manager->AddNode("Pc2Cluster");
  dst_manager->AddNode("Pc3Cluster");

//  dst_manager->AddEventSelector("PairReco_ee");

  se->registerOutputManager(dst_manager);

  dst_manager->Print();

  return;

}

//____________________________________________________
void UDST_IOManager(const char* udstout, Fun4AllServer* se)
{

  // Control for simulated microDST (this is not actually used)
  Fun4AllDstOutputManager *usimDST  = new Fun4AllDstOutputManager("USIMDST", udstout);
 
  usimDST->AddNode("TecHitOut");
  usimDST->AddNode("emcTowerContainer");
  usimDST->AddNode("EventHeader");

  se->registerOutputManager(usimDST);

  usimDST->Print();

  return;

}

//____________________________________________________
void CNT_IOManager(const char* ndstout, Fun4AllServer* se)
{

  // Control for simulated CNT
  Fun4AllDstOutputManager *manager  = new Fun4AllDstOutputManager("SIMCNT", ndstout);
 
  manager->AddNode("PHCentralTrack");
  manager->AddNode("McSingle");
  manager->AddNode("PHGlobal");
  manager->AddNode("PHGlobal_CENTRAL");
  manager->AddNode("EventHeader");
  manager->AddNode("AccCluster"); 
  se->registerOutputManager(manager);

  manager->Print();

  return;

}

//____________________________________________________
void CWG_IOManager(const char* ndstout, Fun4AllServer* se)
{
  Fun4AllDstOutputManager *manager  = new Fun4AllDstOutputManager("SIMCWG", ndstout);
  manager->AddNode("PhCglList");
  manager->AddNode("PhPhotonList");
  manager->AddNode("HadronPhCglList");
  manager->AddNode("PHGlobal");
  manager->AddNode("PHGlobal_CENTRAL");
  manager->AddNode("EventHeader");
  manager->Print();

  se->registerOutputManager(manager);

  return;
}

//____________________________________________________
void HWG_IOManager(const char* ndstout, Fun4AllServer* se)
{

  // Control for simulated HWG
  Fun4AllDstOutputManager *manager  = new Fun4AllDstOutputManager("SIMHWG", ndstout);
 
  manager->AddNode("HWGCentralTrack");
  manager->AddNode("McSingle");
  manager->AddNode("PHGlobal");
  manager->AddNode("PHGlobal_CENTRAL");
  manager->AddNode("EventHeader");
  manager->Print();

  se->registerOutputManager(manager);
  return;

}

//____________________________________________________
void EWG_IOManager(const char* ndstout, Fun4AllServer* se)
{

  // Control for simulated HWG
  Fun4AllDstOutputManager *manager  = new Fun4AllDstOutputManager("SIMEWG", ndstout);
 
  manager->AddNode("EWGCentralTrack");
  manager->AddNode("McSingle");
  manager->AddNode("PHGlobal");
  manager->AddNode("PHGlobal_CENTRAL");
  manager->AddNode("RpSumXYObject");
  se->registerOutputManager(manager);

  manager->Print();

  return;

}

//____________________________________________________
void MWG_IOManager(const char* ndstout, Fun4AllServer* se)
{

  // Control for simulated HWG
  Fun4AllDstOutputManager *manager  = new Fun4AllDstOutputManager("SIMMWG", ndstout);
 
  manager->AddNode("RunHeader");
  manager->AddNode("PHGlobal");
  manager->AddNode("PHGlobal_MUON");
  manager->AddNode("PHMuoTracksOO");
  manager->AddNode("EventHeader");
  manager->AddNode("TMuiPseudoBLTO");
  manager->AddNode("header");
  manager->AddNode("TrigLvl1");
  
  se->registerOutputManager(manager);

  manager->Print();

  return;

}
