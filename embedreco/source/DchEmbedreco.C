#include "TObject.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "PHTypedNodeIterator.h"
#include "PHNodeIterator.h"
#include "PHNodeReset.h"
#include "recoConsts.h"
#include "phool.h"

#include "RunHeader.h"

// Modules
#include "mNewDchCandidatory.hh"
#include "mNewDchInitializer.hh"
#include "mNewDchNoFieldCandidatory.hh"
#include "PHDchGeometryObject.h"
#include "DchEmbedreco.h"

// i/o tables (staff and new classes)
#include "DchHitLineTablev2.hh"
#include "DchHitLineTablev1.hh"
#include "DchTrackv1.h"
#include "McEvalSingleList_v1.h"
#include "dDchReconstructionParWrapper.h"
#include "dDchTracksExtWrapper.h"
#include "mNewDchNoiseAnalyzer.hh"
#include "dDchTracksWrapper.h"
#include "dDchHitWrapper.h"

#include "DchPISAHit.h"
#include "DchMixer.hh"
#include "Fun4AllServer.h"

#include <cmath>
#include <cstdlib>

using namespace std;

typedef PHIODataNode<PHObject> PHObjectNode_t;
typedef PHIODataNode<PHTable> PHTableNode_t;
typedef PHIODataNode<TObject> TObjectNode_t;
typedef PHIODataNode<RunHeader> RunHeaderNode_t;
typedef PHIODataNode <DchHitLineTable> DchHitLineTableNode_t;


DchEmbedreco::DchEmbedreco(const std::string &name)
{
  ThisName = name;
}

int DchEmbedreco::Init(PHCompositeNode *topNode)
{
  if(verbosity>0) { cout << PHWHERE << "DchEmbedreco::Init" << endl; }

  dDchTrack                = NULL;
  dDchTrackPerf            = NULL;
  dDchTrackExt             = NULL;
  dDchTrackExtPerf         = NULL;
  mDchInitializer          = NULL;
  mDchCandidatory          = NULL;
  dchmixer                 = NULL;
  return 0;
}

DchEmbedreco::~DchEmbedreco()
{
  
  if (dchmixer )
    {
      delete dchmixer;
    }
  if (mDchInitializer)
    {
      delete mDchInitializer;
    }

  if (mDchCandidatory)
    {
      delete mDchCandidatory;
    }

}

int DchEmbedreco::InitRun(PHCompositeNode *topNode)
{
  if(verbosity>0) { cout << PHWHERE << "DchEmbedreco::InitRun started..." << endl; }

  recoConsts* rc = recoConsts::instance();
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode;
  dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","DST"));
  if (!dstNode)
    {
      cout << PHWHERE << "DST Node missing doing nothing" << endl;
      return -1;
    }
  PHCompositeNode *parNode;
  parNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","PAR"));

  PHCompositeNode* dchNode = new PHCompositeNode("DCH");
  topNode->addNode(dchNode);

  //set up common tables
  PHIODataNode<PHObject>* DchHitLineTableNode =
    (PHIODataNode<PHObject>*)iter.findFirst("PHIODataNode", "DchHitLineTable");

  if (!DchHitLineTableNode)
    {
      DchHitLineTable* dchhitlinetablev2 = new DchHitLineTablev2();
      DchHitLineTableNode =
        new PHIODataNode<PHObject>(dchhitlinetablev2, "DchHitLineTable", "PHObject"); // contain PHObject
      dstNode->addNode(DchHitLineTableNode);
    }

  PHIODataNode<PHObject>* DchHitLineTablev1Node =
    (PHIODataNode<PHObject>*)iter.findFirst("PHIODataNode", "DchHitLineTablev1");

  if (!DchHitLineTablev1Node)
    {
      DchHitLineTable* dchhitlinetablev1 = new DchHitLineTablev1();
      DchHitLineTablev1Node =
        new PHIODataNode<PHObject>(dchhitlinetablev1, "DchHitLineTablev1", "PHObject"); // contain PHObject
      dstNode->addNode(DchHitLineTablev1Node);
    }

  PHIODataNode<PHObject>* DchTrackNode =
    (PHIODataNode<PHObject>*)iter.findFirst("PHIODataNode", "DchTrack");
  if (!DchTrackNode)
    {
      DchTrack* dchtrackv1 = new DchTrackv1();
      DchTrackNode =
        new PHIODataNode<PHObject>(dchtrackv1, "DchTrack", "PHObject"); // contain PHObject
      dstNode->addNode(DchTrackNode);
    }

  PHIODataNode<PHTable>* dDchHitNode =
    (PHIODataNode<PHTable>*)iter.findFirst("PHIODataNode", "dDchHit");
  if (!dDchHitNode)
    {
      dDchHitWrapper* dDchHit = new dDchHitWrapper("dDchHit", 60000);
      dDchHitNode = new PHIODataNode<PHTable>(dDchHit, "dDchHit");
      dchNode->addNode(dDchHitNode);
    }
  PHIODataNode<PHTable>* dDchTrackNode =
    (PHIODataNode<PHTable>*)iter.findFirst("PHIODataNode", "dDchTracks");
  if (!dDchTrackNode)
    {
      dDchTrack = new dDchTracksWrapper("dDchTracks", 4000);
      dDchTrackNode = new PHIODataNode<PHTable>(dDchTrack, "dDchTracks");
      dchNode->addNode(dDchTrackNode);
    }
  PHIODataNode<PHTable>* dDchTrackExtNode =
    (PHIODataNode<PHTable>*)iter.findFirst("PHIODataNode",
					   "dDchTracksExt");
  if (!dDchTrackExtNode)
    {
      dDchTrackExt = new dDchTracksExtWrapper("dDchTracksExt", 4000);
      dDchTrackExtNode =
	new PHIODataNode<PHTable>(dDchTrackExt, "dDchTracksExt");
      dchNode->addNode(dDchTrackExtNode);
    }
  PHIODataNode<PHTable>*dDchTrackPerfNode =
    (PHIODataNode<PHTable>*)iter.findFirst("PHIODataNode",
					   "dDchTracksPerf");
  if (!dDchTrackPerfNode)
    {
      dDchTrackPerf = new dDchTracksWrapper("dDchTracksPerf", 4000);
      dDchTrackPerfNode =
	new PHIODataNode<PHTable>(dDchTrackPerf, "dDchTracksPerf");
      dchNode->addNode(dDchTrackPerfNode);
    }
  PHIODataNode<PHTable>*dDchTrackExtPerfNode =
    (PHIODataNode<PHTable>*)iter.findFirst("PHIODataNode",
					   "dDchTracksExtPerf");
  if (!dDchTrackExtPerfNode)
    {
      dDchTrackExtPerf = new dDchTracksExtWrapper("dDchTracksExtPerf",
						  4000);
      dDchTrackExtPerfNode = 
	new PHIODataNode<PHTable>(dDchTrackExtPerf,
				  "dDchTracksExtPerf");
      dchNode->addNode(dDchTrackExtPerfNode);
    }
    
  dDchReconstructionParWrapper* dDchRecoPar = new dDchReconstructionParWrapper("dDchRecoPar", 1);
  PHIODataNode<PHTable>* dDchRecoParNode = new PHIODataNode<PHTable>(dDchRecoPar, "dDchRecoPar");
  parNode->addNode(dDchRecoParNode);

  dDchRecoPar->SetRowCount(1);
  dDchRecoPar->set_houghThresholdOnXCell(0, 10);  //
  dDchRecoPar->set_houghThresholdOnXMask(0, 15);  //
  dDchRecoPar->set_houghThresholdOnUVCell(0, 3);  //
  dDchRecoPar->set_houghThresholdOnUVMask(0, 6); //
  dDchRecoPar->set_purgeCandidateThreshold(0, 15); //
  dDchRecoPar->set_firstXHoughThreshold(0, 15);  // mc is 15
  dDchRecoPar->set_secondXHoughThreshold(0, 15); // mc is 15
  dDchRecoPar->set_minimumNumberOfXHits(0, 8);
  dDchRecoPar->set_minimumNumberOfUVHits(0, 0);
  dDchRecoPar->set_XHitsThreshold(0, 10);
  dDchRecoPar->set_cellDifferenceCut(0, 8);
  dDchRecoPar->set_delBetaCut(0, 0.2);
  dDchRecoPar->set_deltaBetaCut(0, 0.2);
  dDchRecoPar->set_wireResolution(0, 0.015); //
  dDchRecoPar->set_initUVChi2(0, 10); //
  dDchRecoPar->set_initXChi2(0, 5);  //
  dDchRecoPar->set_deltaBetaVertexCut(0, 0.5); //
  dDchRecoPar->set_numberOfAlphaBins(0, 300);
  dDchRecoPar->set_numberOfPhiBins(0, 6000);
  dDchRecoPar->set_maxAlpha(0, 0.8);
  dDchRecoPar->set_minAlpha(0, -0.8);
  dDchRecoPar->set_maxPhi(0, 1.);
  dDchRecoPar->set_minPhi(0, -0.65);
  dDchRecoPar->set_numberOfBetaBins(0, 60);
  dDchRecoPar->set_numberOfZedBins(0, 200);
  dDchRecoPar->set_maxBeta(0, 2.5);
  dDchRecoPar->set_minBeta(0, 0.6);
  dDchRecoPar->set_maxZed(0, 100.);
  dDchRecoPar->set_minZed(0, -100.);
  dDchRecoPar->set_mirrorHitAnalysis(0, 1); // 0 for no stereo. 1 for stereo
  
  int tmprunnumber = -1;
  if (rc->FlagExist("RUNNUMBER")) {
    tmprunnumber = rc->get_IntFlag("RUNNUMBER"); 
  }
  cout << PHWHERE << "DchEmbedreco::InitRun : run number for initialization = " << tmprunnumber << endl;

  // Call the initializer which does the most internal initialization work
  mDchInitializer = new mNewDchInitializer(1, 0, 0, tmprunnumber);
  mDchInitializer->setGeometryFileNames
    ("/afs/rhic.bnl.gov/phenix/software/calibration/sim00/DchGeometry.info",
     "/afs/rhic.bnl.gov/phenix/software/calibration/sim00/DchGeometry.wireMc",
     "/afs/rhic.bnl.gov/phenix/software/calibration/sim00/DchGeometry.frame00NoRetracted");

  int DchSet = 0;
  if (rc->FlagExist("PPFLAG") && rc->get_IntFlag("PPFLAG"))
    {
      DchSet = 1;
      cout << "\n  DchEmbedreco <I>: Using p+p dead channel and efficiency maps" << endl;
      mDchInitializer->setNoiseFileName
	("/afs/rhic.bnl.gov/phenix/PHENIX_LIB/calibration/new/AlwaysDeadCh_pp.dat",
	 "/afs/rhic.bnl.gov/phenix/PHENIX_LIB/calibration/new/DchEfficiency0_run2_pp-average.Real");        // Average efficiency map for pp run
    }
  else if (rc->FlagExist("RUN2AUAU") && rc->get_IntFlag("RUN2AUAU"))
    {
      DchSet = 1;
      cout << "\n  DchEmbedreco <I>: Using Run2 Au+Au dead channel and efficiency maps" << endl;
      mDchInitializer->setNoiseFileName
	("/afs/rhic.bnl.gov/phenix/PHENIX_LIB/calibration/new/AlwaysDeadCh.dat",
	 "/afs/rhic.bnl.gov/phenix/PHENIX_LIB/calibration/new/DchEfficiency0_run2_v03-average.Real");  // Average efficiency map used.
    }
  else if (rc->FlagExist("RUN3DAU") && rc->get_IntFlag("RUN3DAU"))
    {
      DchSet = 1;
      cout << "\n  DchEmbedreco <I>: Using Run3 dead channel and efficiency maps" << endl;
      mDchInitializer->setNoiseFileName
	("/afs/rhic.bnl.gov/phenix/PHENIX_LIB/calibration/run2003/AlwaysDeadCh_dAu_72096.dat",
	 "/afs/rhic.bnl.gov/phenix/PHENIX_LIB/calibration/run2003/DchEfficiency_dAu.Real");   // Average efficiency map used.
    }
  else if (rc->FlagExist("RUN4AUAU63GEV") && rc->get_IntFlag("RUN4AUAU63GEV")) 
    {
      DchSet = 1;
      cout << "\n  DchEmbedreco <I>: Using Run4 Au+Au 62.4 GeV dead channel and efficiency maps" << endl;
      mDchInitializer->setNoiseFileName
	("/afs/rhic.bnl.gov/phenix/PHENIX_LIB/calibration/run2004/AlwaysDeadCh_AuAu63_122929.dat",
	 "/afs/rhic.bnl.gov/phenix/PHENIX_LIB/calibration/run2004/DchEfficiency_AuAu63.Real");   // Average efficiency map used.
    } 
  else if (rc->FlagExist("RUN4AUAU200GEV") && rc->get_IntFlag("RUN4AUAU200GEV")) 
    {
      DchSet = 1;
      cout << "\n  DchEmbedreco <I>: Using Run4 Au+Au 200 GeV dead channel and efficiency maps" << endl;
      mDchInitializer->setNoiseFileName
	("/afs/rhic.bnl.gov/phenix/PHENIX_LIB/calibration/run2004/AlwaysDeadCh_AuAu_Run04.dat",
	 "/afs/rhic.bnl.gov/phenix/PHENIX_LIB/calibration/run2004/DchEfficiency_Constant.Real");   // Average efficiency map used.
    }
  else if (rc->FlagExist("RUN7AUAU200GEV") && rc->get_IntFlag("RUN7AUAU200GEV"))
    {
      DchSet = 1;
      cout << "\n  DchEmbedreco <I>: Using Run7 AuAu 200 GeV dead channel and efficiency maps (Run 231428" << endl;
      mDchInitializer->setNoiseFileName
	("/afs/rhic.bnl.gov/phenix/PHENIX_LIB/calibration/run2007/AlwaysDeadCh_231428.dat",
	 "/afs/rhic.bnl.gov/phenix/PHENIX_LIB/calibration/run2007/DchEfficiency_AuAu_231428.Real");  
    }
  
  if(DchSet==0) {
    cerr << endl << PHWHERE << endl;
    cerr << "   DchEmbedreco initialization fails; RHIC Run not specified" << endl;
    exit(1);
  } // Safety check
  
  mDchInitializer->setCalibrationFileName("/afs/rhic.bnl.gov/phenix/software/calibration/sim00/DchCalibration.Mc");
  
  mDchCandidatory = new mNewDchCandidatory;
  mDchCandidatory->setMaxEventForUnsuppressedHits(100000);
  mDchInitializer->event(topNode);
  mDchInitializer->getDGO()->rotateAndTranslate();
  
  dchmixer = new DchMixer;
  if(rc->FlagExist("VERBOSITY"))
    dchmixer ->setVerbose( rc->get_IntFlag("VERBOSITY"));
  if(rc->FlagExist("EMBED_KickOutHitsToSpeedupReco"))
    dchmixer ->setkickOutHitsToSpeedupReco(rc->get_IntFlag("EMBED_KickOutHitsToSpeedupReco"));
  if(rc->FlagExist("EMBED_DCEASTT0"))
    dchmixer ->setdcT0East(rc->get_FloatFlag("EMBED_DCEASTT0"));
  if(rc->FlagExist("EMBED_DCWESTT0"))
    dchmixer ->setdcT0West(rc->get_FloatFlag("EMBED_DCWESTT0"));
  
  
  Fun4AllServer* se = Fun4AllServer::instance();

  PHCompositeNode* mcnode;
  if ( rc->FlagExist("EMBED_MC_TOPNODE") )
    mcnode = se->topNode(rc->get_CharFlag("EMBED_MC_TOPNODE"));
  else
    {
      cout << PHWHERE << "EMBED_MC_TOPNODE char flag should be set" << endl;
      exit(0);
    }
  
  PHCompositeNode* realnode;
  if ( rc->FlagExist("EMBED_REAL_TOPNODE") )
    realnode = se->topNode(rc->get_CharFlag("EMBED_REAL_TOPNODE"));
  else
    {
      cout << PHWHERE << "EMBED_REAL_TOPNODE char flag should be set" << endl;
      exit(0);
    }
  
  PHCompositeNode* mergednode = se->topNode("TOP");
  
  dchmixer->InitRun(mcnode,realnode,mergednode);
  if(verbosity>0) { cout << PHWHERE << "DchEmbedreco::InitRun finished." << endl; }
  return 0;

}

int DchEmbedreco::process_event(PHCompositeNode *topNode)
{
  if(verbosity>0) { cout << PHWHERE << "DchEmbedreco::process_event started..." << endl; }
  
  dchmixer->merge();
  mDchCandidatory->event(topNode);
  return 0;

}

int DchEmbedreco::ResetEvent(PHCompositeNode *topNode)
{  
  mDchCandidatory->ResetEvent(topNode);
  PHNodeIterator mainIter(topNode);
  PHNodeReset reset;
  if (mainIter.cd("DCH"))
    {
      mainIter.forEach(reset);
      mainIter.cd();
    }
  
  return 0;
}
