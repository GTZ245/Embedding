#include <stdio.h>
#include <time.h>

//_______________________________________________________________
void DST_IOManager(const char* dstout, Fun4AllServer* se)
{

  // Control for simulated DST
  Fun4AllDstOutputManager *manager  = new Fun4AllDstOutputManager("SIMDST", dstout);

  recoConsts* rc = recoConsts::instance();

  // run and event header
  manager->AddNode("RunHeader");
  manager->AddNode("EventHeader");
  manager->AddNode("TrigLvl1");

  // pisa nodes

  manager->AddNode("fkin");
  manager->AddNode("pythia");
  manager->AddNode("primary");
  manager->AddNode("header");

  // global vertex and t0
  manager->AddNode("T0Out");
  manager->AddNode("VtxOut");

  // BBC
  manager->AddNode("bbcghit");
  manager->AddNode("BbcOut");
  manager->AddNode("BbcRaw");

  // ZDC
  manager->AddNode("ZdcOut");
  manager->AddNode("ZdcRaw");
  
  // DCH
  manager->AddNode("dcghit");
  manager->AddNode("DchHitLineTablev1");
  manager->AddNode("DchHitLineTable");
  manager->AddNode("DchTrack");
  manager->AddNode("dDchTracksPerf");
  manager->AddNode("dDchTracksExtPerf");
  manager->AddNode("dDchGhitHits");

  // PC
  manager->AddNode("pc1ghit");
  manager->AddNode("pc2ghit");
  manager->AddNode("pc3ghit");
  manager->AddNode("dPc1GhitClus");
  manager->AddNode("dPc2GhitClus");
  manager->AddNode("dPc3GhitClus");
  manager->AddNode("Pc1Cluster");
  manager->AddNode("Pc2Cluster");
  manager->AddNode("Pc3Cluster");
  manager->AddNode("Pc1Raw");
  manager->AddNode("Pc2Raw");
  manager->AddNode("Pc3Raw");

  // RICH
  manager->AddNode("crkghit");
  manager->AddNode("CrkHit");

  // these nodes, when added to the output
  // makes the read-back of the DSTs crash
  // it is removed for the moment
  // they can be added back starting from pro.80
  // manager->AddNode("CrkRing");
  // manager->AddNode("CrkRingBack");
  
  // HBD
  // for now these are commented out
  // on decision from the HBD experts
  // as they are still working on the reconstruction algorithm
  // manager->AddNode("HbdGhitList");
  // manager->AddNode("HbdBlobList");
  // manager->AddNode("HbdCellList");
  // manager->AddNode("HbdHitList");

//--  // TEC
//--  manager->AddNode("tecghit");
//--  manager->AddNode("dTecGhitRaw");
//--  manager->AddNode("TecOut");
//--  manager->AddNode("TecHitOut");

  // TOF
  manager->AddNode("tofghit");
  manager->AddNode("TofOut");
  manager->AddNode("dTofGdigi");
  manager->AddNode("dTofGhitGdigi");
  manager->AddNode("dTofGdigiRec");

  // TOF West
  manager->AddNode("TofwRaw");
  manager->AddNode("TofwHit");

//--  // Aerogel
//--  manager->AddNode("AerGeaHits");
//--  manager->AddNode("AccCluster");
//--  manager->AddNode("AccRaw");
//--  manager->AddNode("AccHit");
  




  // EMCal
  manager->AddNode("emcghit");
  manager->AddNode("emcClusterContainer");
  manager->AddNode("emcTowerContainer");

  /*
  // additional EMCal evaluation nodes
  if( rc->FlagExist("EVALUATIONFLAG") && rc->get_IntFlag("EVALUATIONFLAG")==1 )
  {

    // Evaluation output from EMCal
    manager->AddNode("dEmcGeaClusterTrack");
    manager->AddNode("dEmcGeaTrack");
    manager->AddNode("dEmcGeaTrackCluster");

  }
  */

  // CGL
  manager->AddNode("CglTrack");
  manager->AddNode("CglTrackBack");
  manager->AddNode("PHDchTrackOut");
  manager->AddNode("PHTrackOut");
  manager->AddNode("PHTrackOutBack");

  // SVX
  manager->AddNode("SvxPisaHit");
  manager->AddNode("SvxGhitList");
  manager->AddNode("SvxRawhitList");
  manager->AddNode("SvxGhitRawhitList");
  manager->AddNode("SvxRawhitClusterList");
  manager->AddNode("SvxGhitClusterList");
  manager->AddNode("SvxClusterList");
  manager->AddNode("SvxSegmentList");
  manager->AddNode("SvxCentralTrackList");
  manager->AddNode("SvxCentralTrackBackList");


  // copied from CNT node
  manager->AddNode("PHGlobal");
  manager->AddNode("PHGlobal_CENTRAL");
  manager->AddNode("PHCentralTrack");
//  manager->AddNode("PHGlobal_MUON");


  manager->AddNode("McSingle");


/*
  // compact nodes
  manager->AddNode("DchHit_VarArray");
  manager->AddNode("EmcHit_VarArray");
  manager->AddNode("Pc1Hit_VarArray");
  manager->AddNode("Pc2Hit_VarArray");
  manager->AddNode("Pc3Hit_VarArray");
  manager->AddNode("TofeHit_VarArray");
  manager->AddNode("TofwHit_VarArray");
  manager->AddNode("CrkHit_VarArray");
  manager->AddNode("AccHit_VarArray");
  manager->AddNode("CglTrackHits_VarArray");
  manager->AddNode("CglTrackBackHits_VarArray");
  manager->AddNode("TrackProjection_VarArray");
  manager->AddNode("TrackLineProjection_VarArray");
  manager->AddNode("TrackPathLength_VarArray");
  manager->AddNode("emcHitContainer");

  manager->AddNode("SvxHit_VarArray");
  manager->AddNode("SvxTrack_VarArray");
  manager->AddNode("SvxCentralTrack_VarArray");
  manager->AddNode("SvxCentralTrackBG_VarArray");
*/


  se->registerOutputManager(manager);

  manager->Print();

  return;

}

//____________________________________________________
void UDST_IOManager(const char* udstout, Fun4AllServer* se)
{

  // Control for simulated microDST (this is not actually used)
  Fun4AllDstOutputManager *manager  = new Fun4AllDstOutputManager("USIMDST", udstout);

  manager->AddNode("TrigLvl1");

  // global vertex and t0
  manager->AddNode("T0Out");
  manager->AddNode("VtxOut");

  manager->AddNode("TecHitOut");
  manager->AddNode("emcTowerContainer");
  manager->AddNode("EventHeader");

  se->registerOutputManager(manager);

  manager->Print();

  return;

}

//____________________________________________________
void CNT_IOManager(const char* ndstout, Fun4AllServer* se)
{

  // Control for simulated CNT
  Fun4AllDstOutputManager *manager  = new Fun4AllDstOutputManager("SIMCNT", ndstout);

  // global vertex and t0
  manager->AddNode("T0Out");
  manager->AddNode("VtxOut");
  manager->AddNode("TrigLvl1");

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

  // global vertex and t0
  manager->AddNode("T0Out");
  manager->AddNode("VtxOut");
  manager->AddNode("TrigLvl1");

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

  // global vertex and t0
  manager->AddNode("T0Out");
  manager->AddNode("VtxOut");
  manager->AddNode("TrigLvl1");

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

  // global vertex and t0
  manager->AddNode("T0Out");
  manager->AddNode("VtxOut");
  manager->AddNode("TrigLvl1");

  manager->AddNode("EWGCentralTrack");
  manager->AddNode("McSingle");
  manager->AddNode("PHGlobal");
  manager->AddNode("PHGlobal_CENTRAL");
  manager->AddNode("RpSumXYObject");
  manager->AddNode("EventHeader");
  se->registerOutputManager(manager);

  manager->Print();

  return;

}

//____________________________________________________
void MWG_IOManager(const char* ndstout, Fun4AllServer* se)
{

  // Control for simulated HWG
  Fun4AllDstOutputManager *manager  = new Fun4AllDstOutputManager("SIMMWG", ndstout);

  // global vertex and t0
  manager->AddNode("T0Out");
  manager->AddNode("VtxOut");

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
