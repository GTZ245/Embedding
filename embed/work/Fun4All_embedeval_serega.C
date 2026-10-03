#include "Riostream.h"
//This macro merges simulated DST defined by "mcdst" with real dst. 
////The simulated DST and real DST have been already sorted to have matched BBC 
////vertexes(vertex bin specified by input variable vtxbin)
////The output evaluation ntuple can be specified by variable "ntname"

void Fun4All_embedeval(
const int nevent = 100,
const char *rdinputname = "/phenix/hhj/lebedev/chi_c/embed/skip251500_0/DST_ERT_run8dAu_Central_200GeV_pro82-0000251500-0001.root", 
const char *mcinputname = "/phenix/hhj/lebedev/chi_c/simulation/simdst/simDST_chisig_100_shvtx251500.root", 
const char *dstout      = "kaon.root", 
const char *ntname      = "embed.root",
const int   runnum      = 231429)
{

	cout << "Run number : " << runnum << endl;
	cout << "Number of events : " << nevent << endl;
	cout << "Input RD dst : " << rdinputname << endl;
	cout << "Input MC dst : " << mcinputname << endl;
	cout << "Output dst name          : " << dstout << endl;
	cout << "Output evaluation ntuple : " << ntname << endl;

	float vtxmatch = 1.0;
	float vtxmax = 30.0;

	gSystem->Load("libfun4all.so");
	gSystem->Load("libfun4allfuncs.so");

	gSystem->Load("libcgl.so");
	gSystem->Load("libCrkPID.so");
	
	gSystem->Load("libembed.so");
	gSystem->Load("libembedreco.so");

	gSystem->Load("libemcEmbed4all.so");
	
	//gROOT->ProcessLine(".L embed_IOManager.C");

	gSystem->ListLibraries();

	///////////////////////////////////////////
	// recoConsts setup
	//////////////////////////////////////////
	recoConsts *rc = recoConsts::instance();
	rc->set_IntFlag("RUN7AUAU200GEV",1);    // flag for Run7 200 GeV Au+Au
	rc->set_IntFlag("RUNNUMBER", runnum);  // for 200GEV run8 dAu
	rc->set_IntFlag("EVALUATIONFLAG", 1);  // Requested by EMCal
	rc->set_IntFlag("DCHREQFLAG", 0);
	rc->set_IntFlag("VERBOSITY", 0);
	
	//rc->set_IntFlag("AFSABSENT",1);  // use local files

	Fun4AllServer *se = Fun4AllServer::instance(); 
	se->Verbosity(0);

	/*
	SubsysRecoStack *chroot = new SubsysRecoStack("emcTowerContainerDSTImp", 
		Fun4AllServer::instance()->topNode(rdinputname));
	chroot->x_push_back(new EmcTowerContainerResurrector());
	se->RegisterSubsystem(chroot);
	*/
	
	SubsysReco *mixrec       = new MixEmbedreco("MIX");
	SubsysReco *bbcrec       = new BbcEmbedreco("BBC");
	SubsysReco *vtxrec       = new VtxReco("VTX");
	SubsysReco *padrec       = new PadEmbedreco("PAD");
	SubsysReco *dchrec       = new DchEmbedreco("DCH");
	SubsysReco *tecrec       = new TecEmbedreco("TEC");
	SubsysReco *tofrec       = new TofEmbedreco("TOF");
	SubsysReco *tofwrec      = new TofwEmbedreco("TFW");
	SubsysReco *crkrec       = new CrkEmbedreco("CRK");
	SubsysReco *emcrec       = new EmcEmbedreco("EMC");
	SubsysReco *accrec       = new AccEmbedreco("ACC");
	
	SubsysReco *cglrec       = new CglEmbedreco("CGL");
	SubsysReco *ringrec      = new RingReco("RING");
	SubsysReco *evarec       = new EvaEmbedreco("ChargedEVA");
	
	mixrec->Verbosity(1);
	evarec->Verbosity(1);
	dchrec->Verbosity(0);


	//EmbedVertexSelect enforce the matching of the bbc vertex between real DST and single DST. 
	//The range of matching can be specified by the SetVertexRange function.
	EmbedVertexSelect *vtxsel = new EmbedVertexSelect("VTXSEL","REAL");
	vtxsel->SetVertexRange(vtxmatch); //match vertex with in vtxmatch cm
	vtxsel->Verbosity(1);


	//MC TopNode name
	rc->set_CharFlag("EMBED_MC_TOPNODE","SINGLE");
	// real event TopNode name
	rc->set_CharFlag("EMBED_REAL_TOPNODE","REAL");

	//if one arm has no MC hits, then kick out the hits from real DST, the reconstruction will be much faster.
	rc->set_IntFlag ("EMBED_KickOutHitsToSpeedupReco",0);
	//T0 information for DC East and DC West
	rc->set_FloatFlag("EMBED_DCEASTT0",40);
	rc->set_FloatFlag("EMBED_DCWESTT0",39);

	//The output evaluation ntuples for charged tracks.
	//you can add your own evaluation modules to EvaEmbedreco class
	rc->set_CharFlag("EMBED_CHARGED_EVAOUT",ntname);
	std::cout << "Ntuple name: " << ntname << std::endl;
	// or you can set the output file name directly: 
	//PHEmbedHistogrammer::instance()->setFileName(ntname);

	//sevarl subsystem create node tree in Init method, 
	//this prevent the InputManager from reading the table from DST files
	//these tables includes:
	//DetectorGeometry, VtxOut,CglTrack,CglTrackBack,PHTrackOut,PHTrackOutBack,PHDchTrackOut,AccRaw;
	//this was already fixed by modify the subsystem reco modules.

	//Mix module is suppose to do some initializations, right now is a dummy module
	se->registerSubsystem(mixrec);  
	se->registerSubsystem(bbcrec);
	se->registerSubsystem(vtxrec);
	se->registerSubsystem(padrec);
	se->registerSubsystem(dchrec);
	se->registerSubsystem(tecrec);
	se->registerSubsystem(tofrec);
	se->registerSubsystem(tofwrec);
	se->registerSubsystem(crkrec);
	se->registerSubsystem(accrec);
	se->registerSubsystem(emcrec);
	se->registerSubsystem(cglrec);
	se->registerSubsystem(ringrec);
	se->registerSubsystem(evarec);
	
	Fun4AllInputManager *in1 = new Fun4AllNoSyncDstInputManager("DSTin1","DST","SINGLE"); 
	Fun4AllInputManager *in2 = new Fun4AllDstInputManager("DSTin2","DST","REAL");//real data tree

	in2->registerSubsystem(vtxsel);

	in1->AddFile(mcinputname);   //read into "SINGLE" Node  
	in2->AddFile(rdinputname);   //read into "REAL" Node

	se->registerInputManager(in1);
	se->registerInputManager(in2);

	cout << "running ..." << endl;
	gBenchmark->Start("eventLoop");
	se->run(nevent);
	se->End();
	gBenchmark->Show("eventLoop");
	cout<<"deleting se"<<endl;
	delete se;

	cout<<"finished"<<endl;
	gSystem -> Exit(0);
}
