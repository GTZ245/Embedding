#include "Riostream.h"
//This macro merges simulated DST defined by "mcdst" with real dst. 
//The simulated DST and real DST have been already sorted to have matched BBC 
//vertexes(vertex bin specified by input variable vtxbin)
//The output evaluation ntuple can be specified by variable "ntname"
void Fun4All_embedeval(
   const int nevent = 100,
   const char *mcinputname = "/phenix/hhj/lebedev/chi_c/simulation/simdst/simDST_chisig_100_shvtx251500.root", 
   const char *rdinputname = "/phenix/hhj/lebedev/chi_c/embed/skip251500_0/DST_ERT_run8dAu_Central_200GeV_pro82-0000251500-0001.root", 
   const char *dstout      = "/phenix/hhj/lebedev/chi_c/simulation/pairobj/pairobject_chisigembed_251500_100.root", 
   const char *ntname      = "embed.root",
   const char *ntananame   = "embedana.root",
   const int   runnum      = 409471
                  )
{
  
  cout << "Input RD dst : " << rdinputname << endl;
  cout << "Input MC dst : " << mcinputname << endl;
  cout << "Output dst name          : " << dstout << endl;
  cout << "Output evaluation ntuple : " << ntname << endl;
  cout << "Output evaluation anantuple : " << ntananame << endl;
   
  float vtxmatch = 1000.0;
  float vtxmax = 30.0;
  
/*
  Int_t magField=3;
  switch(magField)
    {
    case 1:  // Run1/Run2 (3D01)
      gSystem->Exec("ln -fs /afs/rhic/phenix/software/calibration/run2001/fieldIntegral.dat fieldIntegral.dat");
      cout << "\n Magnetic field map fieldIntegral.dat file set for Run1/Run2 (3D01)" << endl;
      break;
    case 2:  // Run3 (3D03)
      gSystem->Exec("ln -fs /afs/rhic/phenix/software/calibration/run2003/fieldIntegral.dat fieldIntegral.dat");
      cout << "\n Magnetic field map fieldIntegral.dat file set for Run3 (3D03), same as for Run4 (3D+0)" << endl;
      break;
    case 3:  // Run4 (3D++)
      gSystem->Exec("ln -fs  /afs/rhic/phenix/software/calibration/run2004/fieldIntegral++.dat.run04 fieldIntegral.dat");
      cout << "\n Magnetic field map fieldIntegral.dat file set for Run4 (3D++), same polarity in both coils" << endl;
      break;
    case 4: // Run4 (3D+-) not yet implemented
      gSystem->Exec("ln -fs /afs/rhic/phenix/software/calibration/run2004/fieldIntegral+-.dat.run04 fieldIntegral.dat");
      cout << "\n Magnetic field map fieldIntegral.dat file set for 3D+-, fully reversed polarity in two coils" << endl;
      break;
    case 5: // Run7 (3D+-)
      gSystem->Exec("ln -fs /afs/rhic/phenix/software/calibration/run2007/fieldIntegral+-.dat.run07 fieldIntegral.dat");
      cout << "\n Magnetic field map fieldIntegral.dat file set for 3D+-, fully reversed polarity in two coils as in run7" << endl;
      break;
    default:
      cout << "\n magField value " << magField << " is not recognized; job is aborting" << endl;
      return;
    }
*/
  
  bool local = true;
  bool checkin = true;
  
  gSystem->Load("libfun4all.so");
  gSystem->Load("libfun4allfuncs.so");
  
  gSystem->Load("libembed.so");
  gSystem->Load("libembedreco.so");

  gSystem->Load("libcompactCNT.so");

  gSystem->Load("libsvxcentana.so");

	gSystem->Load("libcgl.so");
	gSystem->Load("libCrkPID.so");
	gSystem->Load("libemcEmbed4all.so");

  gROOT->ProcessLine(".L embed_IOManager.C");

  gSystem->ListLibraries();

  ///////////////////////////////////////////
  // recoConsts setup
  //////////////////////////////////////////
  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUN7AUAU200GEV",1);    // flag for Run4 200 GeV Au+Au
  rc->set_IntFlag("RUNNUMBER", runnum);  // for 200GEV run8 dAu
  rc->set_IntFlag("EVALUATIONFLAG", 1);  // Requested by EMCal
//  rc->set_IntFlag("DCHREQFLAG", 0);
  rc->set_IntFlag("DCHREQFLAG", 3);
  rc->set_IntFlag("PC1REQFLAG", 0);
  rc->set_IntFlag("PC2REQFLAG", 0);
  rc->set_IntFlag("PC3REQFLAG", 0);
  rc->set_IntFlag("TOFREQFLAG", 0);
  // not yet operational ???
  rc->set_IntFlag("EMCREQFLAG", 0);

  rc->set_IntFlag("VERBOSITY", 0); 

  rc->set_IntFlag("AFSABSENT",1);  // use local files
  
  // 2 means PISA-To-DST
  rc->set_IntFlag("SIMULATIONFLAG",2);

  // this should be moved to the Init method of TofSimreco
  rc->set_FloatFlag("TOFTIMINGRESOLUTION", 0.100);



  Fun4AllServer *se = Fun4AllServer::instance(); 
  se->Verbosity(40);
  
  SubsysReco *mixrec       = new MixEmbedreco("MIX");
  mixrec->Verbosity(40);
  SubsysReco *bbcrec       = new BbcEmbedreco("BBC");
  SubsysReco *vtxrec       = new VtxReco("VTX");
  SubsysReco *padrec       = new PadEmbedreco("PAD");
  SubsysReco *dchrec       = new DchEmbedreco("DCH");
  //dchrec->Verbosity(0);
  
  //SubsysReco *tecrec       = new TecEmbedreco("TEC");
  SubsysReco *tofrec       = new TofEmbedreco("TOF");
  SubsysReco *crkrec       = new CrkEmbedreco("CRK");
  SubsysReco *emcrec       = new EmcEmbedreco("EMC");
  //SubsysReco *accrec       = new AccEmbedreco("ACC");
    SubsysReco *cglrec       = new CglEmbedreco("CGL");
  //SubsysReco *cglrec       = new CglReco("CGL");
  
  //SubsysReco *ringrec      = new RingEmbedreco("RING");
  SubsysReco *ringrec      = new RingReco("RING");


  SubsysReco *evarec       = new EvaEmbedreco("ChargedEVA");
  evarec->Verbosity(40);

  
  //EmbedVertexSelect enforce the matching of the bbc vertex between real DST and single DST. 
  //The range of matching can be specified by the SetVertexRange function.
  EmbedVertexSelect *vtxmatch1 = new EmbedVertexSelect("VTX1","REAL");
  vtxmatch1->SetVertexRange(vtxmatch); //match vertex with in vtxmatch cm
  vtxmatch1->Verbosity(40);
 

  //MC TopNode name
  rc->set_CharFlag("EMBED_MC_TOPNODE","SINGLE");
  // real event TopNode name
  rc->set_CharFlag("EMBED_REAL_TOPNODE","REAL");
  
  //if one arm has no MC hits, then kick out the hits from real DST, the reconstruction will be much faster.
  rc->set_IntFlag ("EMBED_KickOutHitsToSpeedupReco",0);
  //T0 information for DC East and DC West
  rc->set_FloatFlag("EMBED_DCEASTT0",40);
  rc->set_FloatFlag("EMBED_DCWESTT0",39);
  //  rc->set_IntFlag("VERBOSITY", 0);
  
  //The output evaluation ntuples for charged tracks.
  //you can add your own evaluation modules to EvaEmbedreco class
  rc->set_CharFlag("EMBED_CHARGED_EVAOUT",ntname);
  // or you can set the output file name directly: 
  //PHEmbedHistogrammer::instance()->setFileName(ntname);
  
  //sevarl subsystem create node tree in Init method, this prevent the InputManager from reading the table from DST files
  //these tables includes:
  //DetectorGeometry, VtxOut,CglTrack,CglTrackBack,PHTrackOut,PHTrackOutBack,PHDchTrackOut,AccRaw;
  //this was already fixed by modify the subsystem reco modules.
  
  //Mix module is suppose to do some initializations, right now is a dummy module
  se->registerSubsystem(mixrec);  
  se->registerSubsystem(bbcrec);
  se->registerSubsystem(vtxrec);
  se->registerSubsystem(padrec);
  se->registerSubsystem(dchrec);
  //se->registerSubsystem(tecrec);
  se->registerSubsystem(tofrec);
  //se->registerSubsystem(crkrec);
  //se->registerSubsystem(emcrec);
  // se->registerSubsystem(accrec);
  se->registerSubsystem(cglrec);
  //se->registerSubsystem(ringrec);
  
  //se->registerSubsystem(new CentraltrackEmbedReco( 22 ));

  se->registerSubsystem(new GlobalReco());

  se->registerSubsystem(new GlobalReco_central());
  //se->registerSubsystem(new CentraltrackReco( 24 ));
  //

  //=========================================
  // These fill the compactCNT storage nodes
  //=========================================
  SubsysReco *fillprojections = new FillTrackProjections();
  SubsysReco *filllineprojections = new FillTrackLineProjections();
  SubsysReco *fillpl = new FillTrackPathLengths();
  SubsysReco *filltrkhits = new FillTrackHits();
  SubsysReco *fillpadhits = new FillPadHits();
  SubsysReco *filldchits = new FillDchHits();
  SubsysReco *filltofehits = new FillTofeHits();
  SubsysReco *filltofwhits = new FillTofwHits();
  SubsysReco *fillcrkhits = new FillCrkHits();
  //SubsysReco *fillacchits = new FillAccHits();
  SubsysReco *fillemchits = new FillEmcHits();

  se->registerSubsystem(fillprojections);
  se->registerSubsystem(filllineprojections);
  se->registerSubsystem(fillpl);
  se->registerSubsystem(filltrkhits);
  se->registerSubsystem(filldchits);
  se->registerSubsystem(fillpadhits);
  se->registerSubsystem(filltofehits);

  se->registerSubsystem(filltofwhits);
  //se->registerSubsystem(fillcrkhits);
  //se->registerSubsystem(fillacchits);

  // This one requires that EmcClusterContainer is already on the node tree
  //se->registerSubsystem(fillemchits);

  //==============================================
  // These modules read the compactCNT nodes
  // and create hits objects for each subsystem
  //================================================


  se->registerSubsystem(new RecoverTrackProjections());
  se->registerSubsystem(new RecoverTrackLineProjections());
  se->registerSubsystem(new RecoverTrackPathLengths());
  se->registerSubsystem(new RecoverTrackHits());
  se->registerSubsystem(new RecoverDchHits());
  se->registerSubsystem(new RecoverPadHits());

  se->registerSubsystem(new RecoverTofeHits());
  se->registerSubsystem(new RecoverTofwHits());

  //se->registerSubsystem(new RecoverCrkHits());
  //se->registerSubsystem(new RecoverAccHits());
  //se->registerSubsystem(new RecoverEmcHits());

  //========================
  // Creates PHCentralTrack
  //========================

  se->registerSubsystem(new CreateCNT());

  //=================================================
  // These modules re-associate hits with tracks and
  // fill the PHCentralTrack fields
  //==================================================

  se->registerSubsystem(new FillCNT_TrackProjections());
  se->registerSubsystem(new FillCNT_TrackPathLengths());
  se->registerSubsystem(new FillCNT_TrackHits());
  se->registerSubsystem(new FillCNT_DchHits());
  se->registerSubsystem(new FillCNT_TofeHits());
  se->registerSubsystem(new FillCNT_TofwHits());

  se->registerSubsystem(new FillCNT_PadHits());
  //se->registerSubsystem(new FillCNT_CrkHits());
  //se->registerSubsystem(new FillCNT_AccHits());
  // This one needs EmcClusterContainer also
  //se->registerSubsystem(new FillCNT_EmcHits());
  
 

 
  //evaluation module should always be the last one, there could be multiple evaluation modules
  se->registerSubsystem(evarec);
  
  //  se->registerSubsystem(central);
  
  // my analysis module
  embedana *ana = new embedana(ntananame);
  se->registerSubsystem(ana);
  
  Fun4AllInputManager *in1 = new Fun4AllNoSyncDstInputManager("DSTin1","DST","SINGLE"); 
  Fun4AllInputManager *in2 = new Fun4AllDstInputManager("DSTin2","DST","REAL");//real data tree
  
  in2->registerSubsystem(vtxmatch1);
  
  in1->AddFile(mcinputname);   //read into "SINGLE" Node  
  in2->AddFile(rdinputname);   //read into "REAL" Node
  
  se->registerInputManager(in1);
  se->registerInputManager(in2);
  
  ///////////////////////////////////////////
  // OutputManagers Set up functions  
  ///////////////////////////////////////////
  if( dstout ) DST_IOManager(dstout, se);
  
  cout << "running ..." << endl;
  se->run(nevent);
  cout<<"finished"<<endl;
  se->End();
  cout<<"deleting se"<<endl;
  delete se;
  
  cout<<"done"<<endl;
  gSystem -> Exit(0);
  
}

