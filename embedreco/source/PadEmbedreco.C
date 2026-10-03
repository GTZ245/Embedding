#include <sstream>
#include <cmath>

#include "PHObject.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "PHTypedNodeIterator.h"
#include "PHNodeIterator.h"
#include "PHNodeReset.h"
#include "recoConsts.h"
#include "phool.h"
#include "PadEmbedreco.h"

#include "padEvtToRaw.hh"
#include "padInclBad.hh"
#include "PadRecModule.hh"
#include "dPad23ParWrapper.h"  /* is this still necessary in simulation? */
#include "dPadClusterWrapper.h"
#include "dPadGeomWrapper.h"
#include "dPadRawClusWrapper.h"
#include "dPadRawWrapper.h"
#include "dPadRecParWrapper.h"  /* is this still necessary in simulation? */
#include "PadRawv1.h"
#include "PadClusterv2.h"
#include "PadMixer.hh"
#include "Fun4AllServer.h"
#include "RunToTime.hh"

using namespace std;

typedef PHIODataNode<PHObject>  PHObjectNode_t;
typedef PHIODataNode<PHTable>   PHTableNode_t;
typedef PHIODataNode<dPadClusterWrapper> dPadClusterNode_t;
typedef PHIODataNode<dPadRawWrapper> dPadRawNode_t;
typedef PHIODataNode<PadCluster> PadClusterNode_t;
typedef PHIODataNode<PadRaw> PadRawNode_t;


PadEmbedreco::PadEmbedreco(const char *name)
{
  ThisName = name;
  
  PadInclBad        = NULL;
  mPadDetGeo        = NULL;
  PadEvtToRaw       = NULL;
  Pc1Rec            = NULL;
  Pc2Rec            = NULL;
  Pc3Rec            = NULL;
  padmixer          = NULL;
  padSplitMode = 0;  // default 0 means no splitting of clusters (Run2 v03 DSTs choice)
}

PadEmbedreco::~PadEmbedreco()
{
  delete PadInclBad;
  delete PadEvtToRaw;
  delete Pc1Rec;
  delete Pc2Rec;
  delete Pc3Rec;
  delete padmixer;
  return ;
}

int PadEmbedreco::Init(PHCompositeNode *topNode)
{
  return 0;
}

int PadEmbedreco::InitRun(PHCompositeNode *topNode)
{

  recoConsts *rc = recoConsts::instance();
  size_t mr, nrc;

  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode;
  dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      cout << PHWHERE << "DST Node missing doing nothing" << endl;
      return -1;
    }
  PHCompositeNode *parNode;
  parNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "PAR"));

 if (!parNode)
    {
      cout << PHWHERE << "PAR Node missing; PadEmbedreco exiting" << endl;
      exit(1);
    }

  PHCompositeNode* padNode = new PHCompositeNode("PAD");
  topNode->addNode(padNode);

  mr = 1;
  dPadGeomWrapper* dPadGeom = new dPadGeomWrapper("dPadGeom", mr);
  PHIODataNode<PHTable>* dPadGeomNode = new PHIODataNode<PHTable>(dPadGeom, "dPadGeom");
  parNode->addNode(dPadGeomNode);

  //!!! Should these pars be hard-coded here?  
  // However, I verified that they are identical to those in v1.35 in
  // offline/framework/simreco/PadReco.C
  // Andrew 11/12/2008

  // Setting dPadGeom Parameters
  nrc = 1;
  dPadGeom->SetRowCount(nrc);
  dPadGeom->set_pdxoff(0, 0, -24.31);
  dPadGeom->set_pdxoff(1, 0, -81.2);
  dPadGeom->set_pdxoff(2, 0, -95.7);
  dPadGeom->set_pdzoff(0, 0, -89.5575);
  dPadGeom->set_pdzoff(1, 0, -152.475);
  dPadGeom->set_pdzoff(2, 0, -178.69);
  dPadGeom->set_pdgas(0, 0, 0.60);
  dPadGeom->set_pdgas(1, 0, 1.00);
  dPadGeom->set_pdgas(2, 0, 1.20);
  dPadGeom->set_aasep(0, 0, 0.84);
  dPadGeom->set_aasep(1, 0, 1.40);
  dPadGeom->set_aasep(2, 0, 1.65);
  dPadGeom->set_pxlen(0, 0, 0.82);
  dPadGeom->set_pxlen(1, 0, 1.375);
  dPadGeom->set_pxlen(2, 0, 1.622);
  dPadGeom->set_wside(0, 0, 0.27);
  dPadGeom->set_wside(1, 0, 0.47);
  dPadGeom->set_wside(2, 0, 0.55);
  dPadGeom->set_wcent(0, 0, 0.15);
  dPadGeom->set_wcent(1, 0, 0.26);
  dPadGeom->set_wcent(2, 0, 0.31);
  dPadGeom->set_pxsep(0, 0, 0.025);
  dPadGeom->set_pxsep(1, 0, 0.025);
  dPadGeom->set_pxsep(2, 0, 0.025);
  dPadGeom->set_clsep(0, 0, 0.1);
  dPadGeom->set_clsep(1, 0, 0.15);
  dPadGeom->set_clsep(2, 0, 0.2);
  dPadGeom->set_npdsec(0, 0, 16);
  dPadGeom->set_npdsec(1, 0, 8);
  dPadGeom->set_npdsec(2, 0, 8);
  dPadGeom->set_npdwr(0, 0, 58);
  dPadGeom->set_npdwr(1, 0, 116);
  dPadGeom->set_npdwr(2, 0, 116);
  dPadGeom->set_npdx(0, 0, 20);
  dPadGeom->set_npdx(1, 0, 40);
  dPadGeom->set_npdx(2, 0, 40);
  dPadGeom->set_npdz(0, 0, 216);
  dPadGeom->set_npdz(1, 0, 216);
  dPadGeom->set_npdz(2, 0, 216);
  dPadGeom->set_sectperarm(0, 0, 8);
  dPadGeom->set_sectperarm(1, 0, 4);
  dPadGeom->set_sectperarm(2, 0, 4);
  dPadGeom->set_inradius(0, 0, 248.891);
  dPadGeom->set_inradius(1, 0, 419.173);
  dPadGeom->set_inradius(2, 0, 492.012);
  dPadGeom->set_zgap(0, 0, 0.0);
  dPadGeom->set_zgap(1, 0, 8.106);
  dPadGeom->set_zgap(2, 0, 8.106);
  dPadGeom->set_phibote(0, 213.75);
  dPadGeom->set_phitope(0, 123.75);
  dPadGeom->set_phibotw(0, -33.75);
  dPadGeom->set_phitopw(0, 56.25);

  mr = 1;
  dPad23ParWrapper* dPad23Par = new dPad23ParWrapper("dPad23Par", mr);
  PHIODataNode<PHTable>* dPad23ParNode = new PHIODataNode<PHTable>(dPad23Par, "dPad23Par");
  parNode->addNode(dPad23ParNode);
  
  mr = 1;
  dPadRecParWrapper* dPadRecPar = new dPadRecParWrapper("dPadRecPar", mr);
  PHIODataNode<PHTable>* dPadRecParNode = new PHIODataNode<PHTable>(dPadRecPar, "dPadRecPar");
  parNode->addNode(dPadRecParNode);
  
  // Setting dPad23Par Parameters
  nrc = 1;
  dPad23Par->SetRowCount(nrc);
  dPad23Par->set_idatePC23(0, 19981225);
  
  // Setting dPadRecPar Parameters
  nrc = 1;
  dPadRecPar->SetRowCount(nrc);
  dPadRecPar->set_verbose(0, 0);
  dPadRecPar->set_method(0, 0, 0);
  dPadRecPar->set_method(1, 0, 0);
  dPadRecPar->set_method(2, 0, 0);
  
  PadEvtToRaw = new padEvtToRaw();
  PadInclBad = new padInclBad();
  mPadDetGeo = new padDetectorGeo();
  
  Pc1Rec = new PadRecModule(0);
  Pc2Rec = new PadRecModule(-1);
  Pc3Rec = new PadRecModule(-2);

  mr = 60000;
  PHTableNode_t  *PHTableIONode;

  dPadRawWrapper *dpadraw;
  const char *NodeNamedRaw[3] = {"dPc1Raw","dPc2Raw","dPc3Raw"};
  PHTypedNodeIterator<dPadRawWrapper> dPadRawIter(topNode);
  for (short int j = 1; j <= 3 ; j++)
    {
      if (!dPadRawIter.find(NodeNamedRaw[j-1]))
        {
          dpadraw = new dPadRawWrapper(NodeNamedRaw[j-1],mr);
	  PHTableIONode = new PHTableNode_t(dpadraw, NodeNamedRaw[j-1]);
          //parNode->addNode(PHTableIONode);
          dstNode->addNode(PHTableIONode);
        }
    }

  dPadClusterWrapper *dpadcluster;
  const char *NodeNamedCluster[3] = {"dPc1Cluster","dPc2Cluster","dPc3Cluster"};
  PHTypedNodeIterator<dPadClusterWrapper> dPadClusterIter(topNode);
  for (short int j = 1; j <= 3 ; j++)
    {
      if (!dPadClusterIter.find(NodeNamedCluster[j-1]))
        {
          dpadcluster = new dPadClusterWrapper(NodeNamedCluster[j-1],mr);
	  PHTableIONode = new PHTableNode_t(dpadcluster, NodeNamedCluster[j-1]);
          dstNode->addNode(PHTableIONode);
        }
    }

  dPadRawClusWrapper *dpadrawclus;
  const char *NodeNamedRawClus[3] = {"dPc1RawClus","dPc2RawClus","dPc3RawClus"};
  PHTypedNodeIterator<dPadRawClusWrapper> dPadRawClusIter(topNode);
  for (short int j = 1; j <= 3 ; j++)
    {
      if (!dPadRawClusIter.find(NodeNamedRawClus[j-1]))
        {
          dpadrawclus = new dPadRawClusWrapper(NodeNamedRawClus[j-1],mr);
	  PHTableIONode = new PHTableNode_t(dpadrawclus, NodeNamedRawClus[j-1]);
          padNode->addNode(PHTableIONode);
        }
    }
  
  PHObjectNode_t *PHObjectIONode;
  PadRaw *padraw = 0;
  const char *NodeNameRaw[3] = {"Pc1Raw", "Pc2Raw", "Pc3Raw"};
  PHTypedNodeIterator<PadRaw> PadRawIter(topNode);
  for (short int j = 1; j <= 3 ; j++)
    {
      if (!PadRawIter.find(NodeNameRaw[j-1]))
        {
          padraw = new PadRawv1();
          PHObjectIONode = new PHObjectNode_t(padraw, NodeNameRaw[j-1], "PHObject");
          dstNode->addNode(PHObjectIONode);
        }
    }

  PadCluster *padcluster = 0;
  const char *NodeNameCluster[3] = {"Pc1Cluster", "Pc2Cluster", "Pc3Cluster"};
  PHTypedNodeIterator<PadCluster> PadClusterIter(topNode);
  for (short int j = 1; j <= 3 ; j++)
    {
      if (!PadClusterIter.find(NodeNameCluster[j-1]))
        {
          padcluster = new PadClusterv2();
          PHObjectIONode = new PHObjectNode_t(padcluster, NodeNameCluster[j-1], "PHObject");
          dstNode->addNode(PHObjectIONode);
        }
    }

  PHDataNode<padDetectorGeo> *mPadDetGeoNode = new PHDataNode<padDetectorGeo>(mPadDetGeo, "mPadDetGeo");
  parNode->addNode(mPadDetGeoNode);
  
  // Setting padInclBad parameters
  PadInclBad->setDebugLevel(0);
  PadInclBad->doNotInclBadChs();
  PadInclBad->doInclBadROCs();
  PadInclBad->RemoveHotROCs();
  PadInclBad->RemoveUnSynchROCs();
  PadInclBad->doNotAddHotROCs();
  PadInclBad->doNotAddInactiveROCs();
  PadInclBad->doNotAddUnSynchROCs();

  // Setting mPadDetGeo Parameters
  double pc1Radius = 248.891;
  double pc2Radius = 419.173;
  double pc3Radius = 492.012;
  double ThetaArm[2] = { -0.589049, 2.15984};
  mPadDetGeo->set_pc1Radius(pc1Radius);
  mPadDetGeo->set_pc2Radius(pc2Radius);
  mPadDetGeo->set_pc3Radius(pc3Radius);
  mPadDetGeo->set_Theta0(ThetaArm);

  // Setting PadRecModule parameters
  short level = 0;              // Debug level 0: no debug information
  Pc1Rec->setDebugLevel(level); // Debug level 1: print method calls
  Pc2Rec->setDebugLevel(level); // Debug level 2: show reconstruction
  //                statistics
  Pc3Rec->setDebugLevel(level); // Debug level 3: show dPadCluster fillings
  // Debug level 4: Save hit information in a
  //                text file -> ./padrec.dat
  // Debug level 5: Save hit information for
  //                Birdseye event display
  //                -> ./birdseye.dat
  Pc1Rec->setSplitMode(padSplitMode); // Split mode 0: Do not split large clusters
  Pc2Rec->setSplitMode(padSplitMode); // Split mode 1: Split large clusters
  Pc3Rec->setSplitMode(padSplitMode);

  // The maximum allowed cluster size (number of cells)
  short splitMax = 100;
  Pc1Rec->setSplitMax(splitMax);
  Pc2Rec->setSplitMax(splitMax);
  Pc3Rec->setSplitMax(splitMax);

  short limit = 2; // Demand two (default) out of three good pads per cell
  Pc1Rec->setPadTypeLimit(limit);
  Pc2Rec->setPadTypeLimit(limit);
  Pc3Rec->setPadTypeLimit(limit);

  // The idea behind the reconstruction code is to place an imaginative box
  // around each cluster and to measure its size. These values set the maximum
  // size of the box that covers the cluster created by a single particle.
  short oneW = 4; // Box height
  short oneZ = 4; // Box width
  short oneL = 11; // Maximum number of allowed cells within a one particle box
  short twoW = 8;
  short twoZ = 8;
  short twoL = 16;
  short threeW = 100;
  short threeZ = 100;
  short threeL = 100;
  Pc1Rec->setClusterSizes(oneW, oneZ, oneL, twoW, twoZ, twoL, threeW, threeZ, threeL);
  Pc2Rec->setClusterSizes(oneW, oneZ, oneL, twoW, twoZ, twoL, threeW, threeZ, threeL);
  Pc3Rec->setClusterSizes(oneW, oneZ, oneL, twoW, twoZ, twoL, threeW, threeZ, threeL);


  RunToTime *rt = RunToTime::instance();
  int runnumber = rc->get_IntFlag("RUNNUMBER");
  PHTimeStamp *ts = rt->getBeginTime(runnumber);

  PHTimeStamp TS = *ts;
  delete ts;
  
  // if one cares about the real data tracks, then you cannot fetch from the sim database
  //  run with a different flag if you care about the tracks in the real event

  PHBoolean sim;

  if( rc->get_IntFlag("LOOKATSIMFROMEMBED",1) == 1)    // 1 means sim data, 0 means real data
    {
      cout << "reading from the sim database for pc..." << endl;
      sim = mPadDetGeo->FetchFromSimDatabase(TS); //reading from sim DB
    }
  else if( rc->get_IntFlag("LOOKATSIMFROMEMBED",1) == 0 )
    {
      cout << "reading from the real database for pc..." << endl;
      sim = mPadDetGeo->FetchFromDatabase(TS);
    }
  else
    {
      cout << PHWHERE << "you picked an undefined number...you should stop and rerun, but I will treat it like looking at sim..." << endl;
      sim = mPadDetGeo->FetchFromSimDatabase(TS); //reading from sim DB
    }
  
  if (!sim) {
      cout << "\n PadEmbedreco <I>: FetchFromSimDatabase failed, using ASCII file instead\n";
      mPadDetGeo->FetchFromFile();  //reading from ASCII file
  }else {
      cout << "\n PadEmbedreco <I>: FetchFromSimDatabase succeeded, used TimeStamp December 31, 2001\n";
  }
  mPadDetGeo->Fetch_dPadGeom(topNode);
  
  //TimeStp =  PHTimeStamp(2004,3,1,23,3,0); // value suggested by Tim Miller for Run4 
  //Now obtain timestamp through get_TimeStamp() method
  PHTimeStamp TimeStp = rc->get_TimeStamp();
  
  cout << "\n  For bad channels PadEmbedreco is using TimeStamp = " << TS;  
  PHBoolean padInclStat = PadInclBad->FetchCalDataFromObjy(TimeStp);
  cout << ",  PadInclBad status return = " << padInclStat << endl << endl;


  padmixer = new PadMixer;
  padmixer ->setVerbose(rc->get_IntFlag("VERBOSITY"));
  Fun4AllServer* se = Fun4AllServer::instance();
  
  PHCompositeNode* mcnode     = se->topNode(rc->get_CharFlag("EMBED_MC_TOPNODE"));
  PHCompositeNode* realnode   = se->topNode(rc->get_CharFlag("EMBED_REAL_TOPNODE"));
  PHCompositeNode* mergednode = se->topNode("TOP");
  
  padmixer->InitRun(mcnode,realnode,mergednode);
  return 0;
}

int PadEmbedreco::CreateNodeTree(PHCompositeNode *topNode)
{
  return 0;
}

int PadEmbedreco::process_event(PHCompositeNode *topNode)
{

  padmixer->merge();
  PadInclBad->event(topNode);
  Pc1Rec->event(0,mPadDetGeo,topNode);
  Pc2Rec->event(1,mPadDetGeo,topNode);
  Pc3Rec->event(2,mPadDetGeo,topNode);
  copyWrapper(topNode);
  //  cout<<"pcdone"<<endl;
  return 0;
  //copyWrapper(topNode);
}

int PadEmbedreco::ResetEvent(PHCompositeNode *topNode)
{
  PHNodeIterator mainIter(topNode);
  PHNodeReset reset;
  if (mainIter.cd("PAD"))
    {
      mainIter.forEach(reset);
      mainIter.cd();
    }

  return 0;
}

void PadEmbedreco::setPadSplitMode(short padMode)
{
  switch (padMode)
    {
    case 0:
      cout << "\n Pad Chamber cluster splitting is disabled\n" << endl;
      padSplitMode = padMode;
      break;
    case 1:
      cout << "\n Pad Chamber cluster splitting is enabled\n" << endl;
      padSplitMode = padMode;
      break;
    default:
      cout << "\n Unrecognized Pad Chamber cluster splitting choice " << padMode << endl;
    }
  return ;  // setPadSplitMode()
}

short PadEmbedreco::getPadSplitMode()
{
  switch (padSplitMode)
    {
    case 0:
      cout << "\n Pad Chamber cluster splitting is disabled\n" << endl;
      break;
    case 1:
      cout << "\n Pad Chamber cluster splitting is enabled\n" << endl;
      break;
    default:
      cout << "\n Unrecognized Pad Chamber cluster splitting choice" << padSplitMode << endl;
    }
  return padSplitMode;
}


int PadEmbedreco::copyWrapper(PHCompositeNode *topNode)
{
  int iret = 0;
  ostringstream WrapperNodeName, PhobjectNodeName;

  dPadClusterWrapper *dpadcluster = NULL;
  for (short int j = 1; j <= 3; j++)
    {
      WrapperNodeName << "dPc" << j << "Cluster";
      PHTypedNodeIterator<dPadClusterWrapper> padcluster_iter(topNode);
      dPadClusterNode_t *dPadClusterNode = padcluster_iter.find(WrapperNodeName.str().c_str());
      if (dPadClusterNode)
        {
          dpadcluster = dPadClusterNode->getData();
        }
      else // no wrapper means we have a microDST
        {
	  cout<<"no Wrapper"<<endl;
          dpadcluster = NULL;
        }
      PhobjectNodeName << "Pc" << j << "Cluster";
      PadCluster *padcluster = NULL;
      PHTypedNodeIterator<PadCluster> padclusteriter(topNode);
      PadClusterNode_t *PadClusterNode = padclusteriter.find(PhobjectNodeName.str().c_str());
      if (PadClusterNode)
        {
          padcluster = PadClusterNode->getData();
        }
      else
        {
          padcluster = NULL;
        }
      if (dpadcluster && padcluster)
        {
          if (!padcluster->isValid() && dpadcluster->RowCount())
            {
              padcluster->FillFromWrapper(dpadcluster);
            }
        }
      WrapperNodeName.str("");
      PhobjectNodeName.str("");
    }

  dPadRawWrapper *dpadraw = NULL;
  for (short int j = 1; j <= 3; j++)
    {
      WrapperNodeName << "dPc" << j << "Raw";
      PHTypedNodeIterator<dPadRawWrapper> padraw_iter(topNode);
      dPadRawNode_t *dPadRawNode = padraw_iter.find(WrapperNodeName.str().c_str());
      if (dPadRawNode)
        {
          dpadraw = dPadRawNode->getData();
        }
      else // no wrapper means we have a microDST
        {
          dpadraw = NULL;
        }
      PhobjectNodeName << "Pc" << j << "Raw";
      PadRaw *padraw = NULL;
      PHTypedNodeIterator<PadRaw> padrawiter(topNode);
      PadRawNode_t *PadRawNode = padrawiter.find(PhobjectNodeName.str().c_str());
      if (PadRawNode)
        {
          padraw = PadRawNode->getData();
        }
      else
        {
          padraw = NULL;
        }
      if (dpadraw && padraw)
        {
          if (!padraw->isValid() && dpadraw->RowCount())
            {
              padraw->FillFromWrapper(dpadraw);
            }
        }
      WrapperNodeName.str("");
      PhobjectNodeName.str("");
    }
  return iret;
}
