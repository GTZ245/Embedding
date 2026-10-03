#include <stdio.h> 
#include <time.h> 
#include <string>
#include <iostream>

char output[200];
vector<string> outfiles;

string MakeOutput(int runnumber, int segment, const char *file ="")
{
  sprintf(output, "%s_%s-%010d-%04d.root", file, gSystem->Getenv("PRODTAG"), runnumber, segment);
}

string MakePRDFOutput(int runnumber, int segment, const char *file ="")
{
  sprintf(output, "%s_%s-%010d-%04d.PRDFF", file, gSystem->Getenv("PRODTAG"), runnumber, segment);
}

void addCommon(Fun4AllDstOutputManager *manager)
{
  manager->AddNode("EventHeader");
  manager->AddNode("Sync");
  manager->AddNode("TrigLvl1");
  manager->AddNode("PreviousEvent");
  manager->AddNode("ErtOut");
  manager->AddNode("VtxOut");
}

void DST_EVE(const int runnumber,const int segment,const char *file,const char *trgsel = 0)
{

  MakeOutput(runnumber,segment,file);
  Fun4AllServer *se = Fun4AllServer::instance();
  Fun4AllDstOutputManager *DSTEVE_Manager  = new Fun4AllDstOutputManager(file,output);

  if (trgsel)
    {
      DSTEVE_Manager->AddEventSelector(trgsel);
    }
  addCommon(DSTEVE_Manager);
  DSTEVE_Manager->AddNode("BbcOut");
  DSTEVE_Manager->AddNode("BbcRaw");
  DSTEVE_Manager->AddNode("ZdcOut");
  DSTEVE_Manager->AddNode("ZdcRaw");
  DSTEVE_Manager->AddNode("SmdOut");
  DSTEVE_Manager->AddNode("PHGlobal");
  DSTEVE_Manager->AddNode("PHGlobal_CENTRAL");
  DSTEVE_Manager->AddNode("SpinDataEventOut");
  //  DSTEVE_Manager->AddNode("RpSumXYObject");

  se->registerOutputManager(DSTEVE_Manager);

}

void DST_SVX(const int runnumber, const int segment,const char *file,const char *trgsel = 0)
{

  MakeOutput(runnumber,segment,file);
  Fun4AllServer *se = Fun4AllServer::instance();
  Fun4AllDstOutputManager *OutManager  = new Fun4AllDstOutputManager(file,output);

  if (trgsel)
    {
      OutManager->AddEventSelector(trgsel);
    }
  addCommon(OutManager);
  OutManager->AddNode("SvxEventInfo");
  OutManager->AddNode("SvxHit_VarArray");
  OutManager->AddNode("SvxTrack_VarArray");
  OutManager->AddNode("SvxCentralTrack_VarArray");
  OutManager->AddNode("SvxCentralTrackBG_VarArray");

//  OutManager->AddNode("SvxClusterList");
//  OutManager->AddNode("SvxSegmentList");
//  OutManager->AddNode("SvxCentralTrackList");
//  OutManager->AddNode("SvxCentralTrackBackList");

  se->registerOutputManager(OutManager);
}


void CNT_Compact(const int runnumber,const int segment,const char *file,const char *trgsel = NULL)
{

  MakeOutput(runnumber,segment,file);
  Fun4AllServer *se = Fun4AllServer::instance();
  Fun4AllDstOutputManager *CNTCOMPACT_Manager  = new Fun4AllDstOutputManager(file,output);

  if (trgsel)
    {
      CNTCOMPACT_Manager->AddEventSelector(trgsel);
    }

  addCommon(CNTCOMPACT_Manager);
  CNTCOMPACT_Manager->AddNode("PHGlobal");
  CNTCOMPACT_Manager->AddNode("PHGlobal_CENTRAL");
  //CNTCOMPACT_Manager->AddNode("RpSumXYObject");
  CNTCOMPACT_Manager->AddNode("DchHit_VarArray");
  CNTCOMPACT_Manager->AddNode("EmcHit_VarArray");
  CNTCOMPACT_Manager->AddNode("Pc1Hit_VarArray");
  CNTCOMPACT_Manager->AddNode("Pc2Hit_VarArray");
  CNTCOMPACT_Manager->AddNode("Pc3Hit_VarArray");
  CNTCOMPACT_Manager->AddNode("TofeHit_VarArray");
  CNTCOMPACT_Manager->AddNode("TofwHit_VarArray");
  CNTCOMPACT_Manager->AddNode("CrkHit_VarArray");
  CNTCOMPACT_Manager->AddNode("AccHit_VarArray");
  CNTCOMPACT_Manager->AddNode("CglTrackHits_VarArray");
  CNTCOMPACT_Manager->AddNode("CglTrackBackHits_VarArray");
  CNTCOMPACT_Manager->AddNode("TrackProjection_VarArray");
  CNTCOMPACT_Manager->AddNode("TrackLineProjection_VarArray");
  CNTCOMPACT_Manager->AddNode("TrackPathLength_VarArray");
  CNTCOMPACT_Manager->AddNode("emcHitContainer");
  CNTCOMPACT_Manager->AddNode("TFvtxCompactTrk");

  se->registerOutputManager(CNTCOMPACT_Manager);

}

void CNT_merge(const int runnumber,const int segment,const char *file,const char *trgsel = NULL)
{

  MakeOutput(runnumber,segment,file);
  Fun4AllServer *se = Fun4AllServer::instance();
  Fun4AllDstOutputManager *CNTmerge_Manager  = new Fun4AllDstOutputManager(file,output);

  if (trgsel)
    {
      CNTmerge_Manager->AddEventSelector(trgsel);
    }

  addCommon(CNTmerge_Manager);
  CNTmerge_Manager->AddNode("PHGlobal");
//CNTmerge_Manager->AddNode("PHGlobal_CENTRAL");
//CNTmerge_Manager->AddNode("RpSumXYObject");
  CNTmerge_Manager->AddNode("BbcOut");
  CNTmerge_Manager->AddNode("ZdcOut");

  CNTmerge_Manager->AddNode("Pc1Cluster");
  CNTmerge_Manager->AddNode("Pc2Cluster");
  CNTmerge_Manager->AddNode("Pc3Cluster");
  CNTmerge_Manager->AddNode("Pc1Raw");
  CNTmerge_Manager->AddNode("Pc2Raw");
  CNTmerge_Manager->AddNode("Pc3Raw");

  CNTmerge_Manager->AddNode("DchHitLineTable");
  CNTmerge_Manager->AddNode("DchTrack");

  CNTmerge_Manager->AddNode("CrkHit");
  CNTmerge_Manager->AddNode("CrkRing");
  CNTmerge_Manager->AddNode("CrkRingBack");

  CNTmerge_Manager->AddNode("emcClusterContainer");
  CNTmerge_Manager->AddNode("emcTowerContainer");

  CNTmerge_Manager->AddNode("CglTrack");
  CNTmerge_Manager->AddNode("CglTrackBack");

  CNTmerge_Manager->AddNode("PHDchTrackOut");
  CNTmerge_Manager->AddNode("PHTrackOut");
  CNTmerge_Manager->AddNode("PHTrackOutBack");

  CNTmerge_Manager->AddNode("TofOut");
  CNTmerge_Manager->AddNode("TofwRaw");

  CNTmerge_Manager->AddNode("AccRaw");
  CNTmerge_Manager->AddNode("AccCluster");


  //CNTmerge_Manager->AddNode("AccHit_VarArray");
  //CNTmerge_Manager->AddNode("emcHitContainer");
  //CNTmerge_Manager->AddNode("TFvtxCompactTrk");

  CNTmerge_Manager->AddNode("SvxRawhitList");
  CNTmerge_Manager->AddNode("SvxClusterList");
  CNTmerge_Manager->AddNode("SvxGhitClusterList");
  CNTmerge_Manager->AddNode("McSingle");
  CNTmerge_Manager->AddNode("McEvalSingleList");
  CNTmerge_Manager->AddNode("PHEmbedMcRecoTrackList");
  CNTmerge_Manager->AddNode("PHEmbedMcRecoTrack");
  CNTmerge_Manager->AddNode("PHEmbededEvent");
  //CNTmerge_Manager->AddNode("PHCentralTrack");

  se->registerOutputManager(CNTmerge_Manager);
}

void JPSI_IOManager(const int runnumber,const int segment,const char *file,const char *trgsel = 0)
{
  MakeOutput(runnumber,segment,file);

  Fun4AllServer *se = Fun4AllServer::instance();
  Fun4AllDstOutputManager *MWG_Manager  = new Fun4AllDstOutputManager(file,output);

  if (trgsel)
    {
      MWG_Manager->AddEventSelector(trgsel);
    }
  MWG_Manager->AddNode("PHMuoTracksOO");
  MWG_Manager->AddNode("PHGlobal");
  MWG_Manager->AddNode("PHGlobal_MUON");
  addCommon(MWG_Manager);
  se->registerOutputManager(MWG_Manager);
}

void MWG_IOManager(const int runnumber,const int segment,const char *file,const char *trgsel = 0)
{
  MakeOutput(runnumber,segment,file);

  Fun4AllServer *se = Fun4AllServer::instance();
  Fun4AllDstOutputManager *MWG_Manager  = new Fun4AllDstOutputManager(file,output);

  if (trgsel)
    {
      MWG_Manager->AddEventSelector(trgsel);
    }
  addCommon(MWG_Manager);
	
	// Important group
  MWG_Manager->AddNode("PHMuoTracksOO");
  MWG_Manager->AddNode("PHGlobal");
  MWG_Manager->AddNode("PHGlobal_MUON");

  se->registerOutputManager(MWG_Manager);
}


void MuonDST_IOManager(const int runnumber, const int segment, const char *file, const char *trgsel=0)
{
  MakeOutput(runnumber,segment,file);

  Fun4AllServer *se = Fun4AllServer::instance();
  Fun4AllDstOutputManager *MuonDST_Manager  = new Fun4AllDstOutputManager(file,output);

  if (trgsel)
    {
      MuonDST_Manager->AddEventSelector(trgsel);
    }
  addCommon(MuonDST_Manager);

  //  MuonDST_Manager->AddNode("PHMuoTracksAdc");
  MuonDST_Manager->AddNode("TMuiPseudoLL1");
  MuonDST_Manager->AddNode("TMuiRoadO");
  MuonDST_Manager->AddNode("TMutTrk");
  MuonDST_Manager->AddNode("TMutTrkMap");
  MuonDST_Manager->AddNode("TFvtxCompactTrk");

  se->registerOutputManager(MuonDST_Manager);
}


void MuTRnDST_IOManager(const int runnumber,const int segment,const char *file,const char *trgsel = 0)
{
  MakeOutput(runnumber,segment,file);

  Fun4AllServer *se = Fun4AllServer::instance();
  Fun4AllDstOutputManager *MuTRnDST_Manager  = new Fun4AllDstOutputManager(file,output);

  if (trgsel)
  {
    MuTRnDST_Manager->AddEventSelector(trgsel);
  }
  addCommon(MuTRnDST_Manager);
  //  MuTRnDST_Manager->AddNode("SpinDataEventOut");
  MuTRnDST_Manager->AddNode("BbcOut");
  MuTRnDST_Manager->AddNode("PHMuoTracksOO");
  MuTRnDST_Manager->AddNode("PHGlobal");
  MuTRnDST_Manager->AddNode("PHGlobal_MUON");
  MuTRnDST_Manager->AddNode("TMuiPseudoLL1");
  //  MuTRnDST_Manager->AddNode("RpSumXYObject");
  se->registerOutputManager(MuTRnDST_Manager);
}

char* muid_eff_IOManager(const int runnumber,const int segment)
{
  MakeOutput(runnumber,segment,"MUID_EFF");
  return output;
}

char* pDST_IOManager(const int runnumber,const int segment)
{
  MakeOutput(runnumber,segment,"pDST");
  return output; // output char filled within MakeOutput function.
}

std::string LPanaOut(const int runnumber,const int segment, const char *file)
{
  std::string out = MakeOutput(runnumber,segment,file);
  cout << "full fname: " << out << endl;
  return out;
}


void QA_IOManager(const int runnumber,const int segment)
{
  MakeOutput(runnumber,segment,"qaRoot");
  Fun4AllServer *se = Fun4AllServer::instance();
  Fun4AllHistoManager *hm = se->getHistoManager("QA");
  if(hm) hm->dumpHistos(output);
}

void SVXQA_IOManager(const int runnumber,const int segment)
{
  MakeOutput(runnumber,segment,"SVXQA");
  Fun4AllServer *se = Fun4AllServer::instance();
  Fun4AllHistoManager *hm = se->getHistoManager("SVXQA");
  if(hm) hm->dumpHistos(output);
}

void PRDF_IOManager(const int runnumber,const int segment, const char *file, char *trgsel=0, char *trgsel2=0)
{
  MakePRDFOutput(runnumber,segment,file);
  Fun4AllServer *se = Fun4AllServer::instance();
  Fun4AllEventOutputManager *prdf_manager       = new Fun4AllEventOutputManager(file, output);

  if (trgsel)
    {
      prdf_manager->AddEventSelector(trgsel);
    }
  if (trgsel2)
    {
      prdf_manager->AddEventSelector(trgsel2);
    }
  se->registerOutputManager(prdf_manager);
}

void DST_MPC(const int runnumber,const int segment,const char *file, char *trgsel=0)
{

  MakeOutput(runnumber,segment,file);
  Fun4AllServer *se = Fun4AllServer::instance();
  Fun4AllDstOutputManager *DST_MPC_Manager  = new Fun4AllDstOutputManager(file,output);

  if (trgsel)
    {
      DST_MPC_Manager->AddEventSelector(trgsel);
    }

  addCommon(DST_MPC_Manager);
  DST_MPC_Manager->AddNode("PHGlobal");
  DST_MPC_Manager->AddNode("MpcRaw2");

  se->registerOutputManager(DST_MPC_Manager);
}

void FileSummary(const char *file = "FileSummary.txt")
{
  ofstream fp_out;
  fp_out.open(file, ios::app);
  for (unsigned int i=0; i<outfiles.size(); i++)
    {
      fp_out << gSystem->Getenv("PRODTAG") << " " << outfiles[i] << endl;
    }
  fp_out.close();
}
