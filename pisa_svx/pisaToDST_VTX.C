/*!
   \file pisaToDST.C
   \brief pisa to DST reconstruction chain
   \author <a href="mailto:pereira@hep.saclay.cea.fr">Hugo Pereira</a>
   \version $Revision: 1.1 $
   \date $Date: 2015/11/16 21:12:08 $
*/
// Modified 26 Jan 2016 D. McGlinchey
//    Included VTX modules
//
// filein :
//    Input PISA file
// dstout :
//    Name of output DST
// evalout :
//    Name of SvxEvaluator output (if registered)
// run_number :
//    Run number. Used for hot dead maps in VTX, if SvxApplyHotDead is registered
// dataFlag : 
//    0=Single : SimVtx in VtxOut is automatically used (since only SimVtx is kept in the VtxOut)
//    1=Hijing : Precise Vtx ordered by VtxOut is used as same as real data.
//               You need to set the beam offset (svxpar->set_BeamCenter(x, y)) by hand. default is (0,0)
//    2=Pythia : ppflag=1 in SvxStandaloneReco
// B-field
//  Direction of magnetic field in simulation reconstruction is taken 
//  from the PISA input file, without using run number.

void pisaToDST_VTX(
  Int_t nEvents = 1000,
  char *filein="PISAEvent.root",
  char *dstout = "dst_out.root",
  char *evalout = "svxeval.root",
  int run_number = 430238,
  int dataFlag = 0
)
{

  // print output
  cout << "pisaToDST - nEvents: " << nEvents << endl;
  cout << "pisaToDST - filein: " << filein << endl;
  if( dstout ) cout << "pisaToDST - dstout: " << dstout << endl;
  cout << "pisaToDST - run_number: " << run_number << endl;
  cout << endl;

  // allow to disable muon simulations
  // they are enabled by default
  bool do_muon_arms = false;

  // load libraries
  gSystem->Load("libsvx");
  gSystem->Load("libfun4all");
  gSystem->Load("libmutoo_subsysreco");
  gSystem->Load("libfun4allfuncs");
  gSystem->Load("libsimreco");
  gSystem->Load("libcompactCNT.so");
  gSystem->Load("libSvxDstQA.so");
  gSystem->Load("libsvxeval");

  gROOT->ProcessLine(".L pisaToDST_IOManager.C");

  gSystem->ListLibraries();

  ///////////////////////////////////////////
  // recoConsts setup
  //////////////////////////////////////////
  recoConsts *rc = recoConsts::instance();

  // 2 means PISA-To-DST
  rc->set_IntFlag("SIMULATIONFLAG",2);

  // disable embedding
  rc->set_IntFlag("EMBEDFLAG",0);

  // Reference run number used in 2007 Au+Au 200 GeV
  rc->set_IntFlag("RUNNUMBER",run_number);

  // Requested by EMCal
  rc->set_IntFlag("EVALUATIONFLAG", 0);
  //rc->set_IntFlag("EMCSIMULATIONV2", 1); // enable if you want the new code
  //rc->set_IntFlag("EMCSIMULATIONV2NOQA", 1); // do not apply QA, disabled

  // this should be moved to the Init method of TofSimreco
  rc->set_FloatFlag("TOFTIMINGRESOLUTION", 0.100);

  /*
  Flags to abort event if required number of GEANT hits is not present in the subsystem
  Defaults are all 0 except for the Drift Chamber
  default setting is 3 Drift Chamber wire plane hits
  */
  rc->set_IntFlag("DCHREQFLAG", 3);
  rc->set_IntFlag("PC1REQFLAG", 0);
  rc->set_IntFlag("PC2REQFLAG", 0);
  rc->set_IntFlag("PC3REQFLAG", 0);
  rc->set_IntFlag("TOFREQFLAG", 0);

  // not yet operational
  rc->set_IntFlag("EMCREQFLAG", 0);

  // assume AFS is present as at RCF
  rc->set_IntFlag("AFSABSENT", 0);

  //--------------- added
  // Kalman Flags
  /*
  rc->set_FloatFlag("KALPMIN",0.400);
  rc->set_IntFlag("KALFILTERDCUV",1);
  rc->set_IntFlag("KALFIT",1);
  rc->set_IntFlag("KALTESTNTUPLE",0);
  rc->set_IntFlag("KALREGENDERIV",1);
  rc->set_IntFlag("KALUSEDCHX1X2",1);
  rc->set_IntFlag("KALUSEDCHUV",1);
  rc->set_IntFlag("KALUSEPC1",1);
  rc->set_IntFlag("KALUSEPC2",0);
  rc->set_IntFlag("KALUSEPC3",0);
  rc->set_IntFlag("KALUSETEC",0);
  rc->set_IntFlag("KALUSETOF",0);
  rc->set_IntFlag("KALUSEEMC",0);
  rc->set_IntFlag("KALUSESVX",1);
  rc->set_IntFlag("KALSVXASSOC",0); // set to 0 if you want to use cgl results for svx
  */
 
  // simVertexFlag = 0 (default) means that the BBC Z0 value will be used
  // simVertexFlag = 1 means that the same simZ0Vertex centroid value is used for all events
  // simVertexFlag = 2 means that the Z0 centroid is taken from the PISA event header for each event
  // The centroid values are modified by the Width values which are Gaussian sigma values
  Int_t simVertexFlag=2;
  Float_t simZ0Vertex=0.0, simT0Vertex=0.0;
  Float_t simZ0VertexWidth=2.0, simT0VertexWidth=0.05;
  
  rc->set_IntFlag("SIMVERTEXFLAG",simVertexFlag);
  rc->set_FloatFlag("SIMZ0VERTEX",simZ0Vertex);             // checked in BbcSimreco only when simVertexFlag = 1
  rc->set_FloatFlag("SIMZ0VERTEXWIDTH",simZ0VertexWidth);   // checked in BbcSimreco only when simVertexFlag = 1 or 2
  rc->set_FloatFlag("SIMT0VERTEX",simT0Vertex);             // checked in BbcSimreco only when simVertexFlag = 1
  rc->set_FloatFlag("SIMT0VERTEXWIDTH",simT0VertexWidth);   // checked in BbcSimreco only when simVertexFlag = 1 or 2

  //rc->set_IntFlag("SVXACTIVE",0);


  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  ///////////////////////////////////////////
  // Activate the subsystems
  //////////////////////////////////////////

  // run header and trigger setting
  SubsysReco *head = new HeadSimreco();
  se->registerSubsystem(head );
  se->registerSubsystem( new TrigSimreco() );

  // event counter, set to true if you really want it
  if( false ) se->registerSubsystem( new MuonCounter() );

  // BBC simReco
  se->registerSubsystem(new BbcSimreco("BBC"));

  // pisa is used as an input vertex.
  // it overwrites the contents of the BBC out node.
  // NOTE: This smears according to the BBC z vertex resolution
  //VtxSimreco* vtx_sim = new VtxSimreco();
  //vtx_sim->UseVtx( VtxSimreco::PISA );
  //vtx_sim->UseXY( true );   // default is false
  //vtx_sim->SmearZ( false ); // default is true
  //vtx_sim->OverwriteBBC( true );  // this is the default
  //vtx_sim->ZVertexSigma( 0.5 );   // default error on the simulated vertex
  //se->registerSubsystem( vtx_sim );
  
  //NOTE: This is how to set the values for the SVXPRECISE resolution
  //  current values are for Run 11 AuAu - Need to be updated for the 
  //  system being run


  // We have a BBCout node and a VtxOut node. The BBCOut(VtxOut) node stores the vertex reconstructed by the BBC(VTX).
  // OverWriteBBC() by default is true. It overwrites the BBC reconstructed vertex
  // by the vertex stored in the event header (oscar.particles.dat).
  
   VtxSimreco *vtx_sim = new VtxSimreco();
   vtx_sim->UseVtx( VtxSimreco::PISA ); // Store the event vertex as given in the oscar.paerticles.dat 
   vtx_sim->SmearX(false); 
   vtx_sim->SmearY(false);
   vtx_sim->SmearZ(false);
   vtx_sim->XVertexSigma(0.0); //0.00962
   vtx_sim->YVertexSigma(0.0); //0.00426
   vtx_sim->ZVertexSigma(0.0); //0.00751
   se->registerSubsystem( vtx_sim );
  

  // t0
  T0Simreco* t0_sim = new T0Simreco();
  t0_sim->T0Sigma(0.04);
  se->registerSubsystem( t0_sim );

  // pad chambers
  se->registerSubsystem(new PadSimreco("PAD"));

  // The VtxReco works unchanged for both real and simulation events
  se->registerSubsystem(new VtxReco("VTX"));

  // The T0Reco works unchanged for both real and simulation events
  se->registerSubsystem(new T0Reco());

  DchSimreco *dch = new DchSimreco("DCH");
  //  dch->setSeed(1234);  // is you want to use a fixed seed to make running reproducible
  //  dch->perfecttracker(); // turn on perfect tracker (fitting geant hits directly - CAVEAT: dch hits do not match to tracks)
  //  dch->setDeadMapFile("none"); // turn off dead map (or if file exists, use that file)
  //  dch->setEfficiencyFile("none"); // turn off efficiency (use default 0.95)  (or if file exists, use that file)
  se->registerSubsystem(dch );

  // Time of flight detector
  se->registerSubsystem(new TofSimreco("TOF"));

  // Tof west
  se->registerSubsystem(new TfwSimreco("TFW"));

  // RICH
  se->registerSubsystem(new CrkSimreco("CRK"));

  // Aerogel subsystem as per e-mail from Narumi Kurihara on May 13, 2005
  se->registerSubsystem(new AccSimreco("ACC"));
  se->registerSubsystem(new AccReco());

  // EMCal uses the real data class
  // rc->set_FloatFlag("EMCTOWERLOWGAIN", 0.0015625);
  // rc->set_FloatFlag("EMCTOWERHIGHGAIN", 0.0125);
  se->registerSubsystem( new EmcReco3() );

  //------------------------------
  // VTX simulation
  // register first 
  SvxParManager *svxpar = new SvxParManager();
  svxpar->set_ReadGeoParFromFile(1);  // read parameters from ascii file
  svxpar->set_OffsetVtxToCnt  (0.0, 0.0, 0.0);
  svxpar->set_OffsetEastToWest(0.0, 0.0, 0.0);
  svxpar->set_BeamCenter(0.0, 0.0);
  svxpar->Load_ThresholdFile("svx_threshold.dat");
  svxpar->set_UseStripThresholdDatbase(false);
  se->registerSubsystem(svxpar);

  // Register to include VTX dead map                                                                                                          
  SvxSimulator *svxsim     = new SvxSimulator();
  se->registerSubsystem(svxsim);

  SubsysReco *svxapplyhotdead  = new SvxApplyHotDead();
  svxapplyhotdead->Verbosity(0);
  se->registerSubsystem(svxapplyhotdead);

  SvxReco *svxreco     = new SvxReco();
  svxreco->Verbosity(0);
  svxreco->set_ThisIsSimulation();
  svxreco->set_StripixelAdcSumThreshold(0);
  se->registerSubsystem(svxreco);

  SubsysReco *svxvtxseedfinder = new SvxPriVertexSeedFinder();
  se->registerSubsystem(svxvtxseedfinder);
  
  SvxStandAloneReco *svxstandalone = new SvxStandAloneReco();
  svxstandalone->Verbosity(0);
  if(dataFlag==2) svxstandalone->setPPFlag(true); //Makes windows wider in p+p
  svxstandalone->setVertexRecoFlag(2); //2 gets best vertex
  se->registerSubsystem( svxstandalone );

  SubsysReco *svxprimvtxfinder = new SvxPrimVertexFinder();
  svxprimvtxfinder->Verbosity(0);
  se->registerSubsystem(svxprimvtxfinder);


  //---------------

  // The CglReco works unchanged for both real and simulation events
  CglReco *cgl = new CglReco("CGL");
  cgl->set_SvxUseAsciiFile(true);
  se->registerSubsystem(cgl);

  //Aerogel cluster  (Needs to be after cglRec)
  se->registerSubsystem(new AccclusterReco());

  //  This is the class which makes the RICH Ring data structure
  se->registerSubsystem( new RingReco() );

  // This is the class which makes the Central Tracks nanoDST output
  // 22 corresponds to the version used in pro.78 for Run7 Au+Au
  se->registerSubsystem(new CentraltrackReco( 22 ));

  //  This is the class which makes the GlobalEvent data on the nanoDST output
  se->registerSubsystem(new GlobalReco());

  // This is the class which checks for charged particles going into EMCal
  se->registerSubsystem(new ChargedvetoReco());

  //added the DC based global evaluation module
  se->registerSubsystem( new McEvalSimreco() );


  // muon arm reconstruction
  if( do_muon_arms )
  {

    // unfortunately the muon arm need the sign of the magnetic field
    // and does not have yet the logic to retrieve it from the pisa event header
    // it is hard-coded here
    mMfmMT::setMapFileScale( 1.0 );
    MuonUtil::set_check_mapfile_scale( false );

    // pisa unpacking
    // The muon reconstruction is not performed
    // and only the simulated objects are created.
    // the reconstruction itself is performed in the Fun4All_RecoDST_sim afterburner macro
    se->registerSubsystem( new MuonUnpackPisa() );

  }

/*

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
  SubsysReco *fillacchits = new FillAccHits();
  SubsysReco *fillemchits = new FillEmcHits();

  se->registerSubsystem(fillprojections);
  se->registerSubsystem(filllineprojections);
  se->registerSubsystem(fillpl);
  se->registerSubsystem(filltrkhits);
  se->registerSubsystem(filldchits);
  se->registerSubsystem(fillpadhits);
  se->registerSubsystem(filltofehits);
  se->registerSubsystem(filltofwhits);
  se->registerSubsystem(fillcrkhits);
  se->registerSubsystem(fillacchits);

  // This one requires that EmcClusterContainer is already on the node tree
  se->registerSubsystem(fillemchits);

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
  se->registerSubsystem(new RecoverCrkHits());
  se->registerSubsystem(new RecoverAccHits());
  se->registerSubsystem(new RecoverEmcHits());

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
  se->registerSubsystem(new FillCNT_CrkHits());
  se->registerSubsystem(new FillCNT_AccHits());
  // This one needs EmcClusterContainer also
  se->registerSubsystem(new FillCNT_EmcHits());
*/
  
  // SvxCentralTrack
  SvxCentralTrackReco *svxcnttrackreco = new SvxCentralTrackReco();
  svxcnttrackreco->setSearchWindowFlag(1); // default(0), for single particle sim, use(1) = tigher cut
  svxcnttrackreco->setVertexFlag(1); //1 is simulated vertex, 0 is precise if it exists or else bbcz, 2 is just svxprecise
  svxcnttrackreco->setPrintLinkInfo(true);
  se->registerSubsystem( svxcnttrackreco );

  /*

  // Select only clusters associated with SvxCentralTracks
  SvxSelectClusters* svxselect = new SvxSelectClusters();
  svxselect->Verbosity(0);
  se->registerSubsystem(svxselect);


  // svx compactCNTs
  FillSvxHits *fillsvxhits = new FillSvxHits();
  fillsvxhits->Verbosity(0);
  se->registerSubsystem(fillsvxhits);

*/
  // dumper if you want to look at ascii files (good for checking for 
  // changes by diffing before and after)
  if( false )
  {
    gSystem->Load( "libphnodedump" );
    se->registerSubsystem( new Dumper() );
  }

  ///////////////////////////////////////////
  // SVX Evaluator
  ///////////////////////////////////////////
  //SvxEvaluator *svxeval = new SvxEvaluator("SVXEVALUATOR", evalout);
  //svxeval->Verbosity(0);
  //se->registerSubsystem(svxeval);

  ///////////////////////////////////////////
  // InputManager
  ///////////////////////////////////////////
  Fun4AllInputManager *input_manager = new Fun4AllPisaInputManager("PisaIn","TOP");
  se->registerInputManager(input_manager);

  ///////////////////////////////////////////
  // OutputManagers Set up functions
  ///////////////////////////////////////////
  if( dstout ) DST_IOManager(dstout, se);

  ///////////////////////////////////////////
  // open input file
  se->fileopen(input_manager->Name(),filein);

  // process input events
  gBenchmark->Start("eventLoop");
  se->run(nEvents);
  se->End();
  gBenchmark->Show("eventLoop");

  delete se;
  // If you do not see this message, the job failed
  cout << "Completed reconstruction." << endl;
  gSystem->Exit(0);
}
