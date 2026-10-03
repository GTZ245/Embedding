#include <string>
#include <sstream>
#include <iostream>
#include <vector>

using namespace std;

void Fun4All_CA_merge(int nEvents=100000,
//const char *inputfile="/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/yuri_embed/real/RealDST/CNT_MB_run12CuAu_200GeV_CA_pro99-0000372403-9000.root") 
	     const char *inputfile="/phenix/crs/PRDFs/run12/200CuAu/EVENTDATA_P00-0000372586-0011.PRDFF"
	     //const char *inputfile="/phenix/lustre01/phenix/fromhpss/run12/eventdata/run_0000373000_0000374000/EVENTDATA_P00-0000373140-0000.PRDFF"
	     )
{
  //Tell root to really crash when something goes wrong not start it's
  //signal handling.
  for (int i = 0; i < kMAXSIGNALS; i++)
    {
      gSystem->IgnoreSignal((ESignals)i);
    }

  gSystem->Exec("/bin/env");

  char ifile[5000];
  strcpy(ifile, inputfile);
  strtok(ifile, "-");
  int runnumber = atoi(strtok(0, "-"));
  int segnumber = atoi(strtok(strtok(0, "-"), "."));

  ///////////////////////////////////////////
  // Load Libraries
  //////////////////////////////////////////
  gSystem->Load("libfun4all.so");
  gSystem->Load("libfun4allfuncs.so");
  gSystem->Load("libcompactCNT.so");
  gSystem->Load("libspin.so");
  gSystem->Load("libSvxDstQA.so");
  gSystem->Load("libDstMerger.so");
  gSystem->Load("librecal.so");

  gROOT->ProcessLine(".L OutputManager.C");
  gROOT->ProcessLine(".L rawdatacheck.C");
  gROOT->ProcessLine(".L TrigSelect.C");

  gSystem->ListLibraries();

  SetCvsTag();

  ///////////////////////////////////////////
  // recoConsts setup
  //////////////////////////////////////////
  recoConsts *rc = recoConsts::instance();
  //rc->set_IntFlag("PADVTX_PC3_MULT_CUT",10);  // padvtx cut
  rc->set_IntFlag("MPC_RECO_MODE",0x3);
  rc->set_FloatFlag("EASTMAXSAG", -0.017);
  rc->set_FloatFlag("WESTMAXSAG", -0.017);

  // set the cvstag from the build to store in the dst
  //SetCvsTag();

  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////
  Fun4AllServer *se = Fun4AllServer::instance(); 
  se->Verbosity(0);

  ///////////////////////////////////////////
  // Make and register the Raw Data Checker
  //////////////////////////////////////////
  RawDataCheck *raw = rawdatacheck();

  ///////////////////////////////////////////
  // Make the Synchronization Object
  ///////////////////////////////////////////
  SubsysReco *sync = new SyncReco();

  //////////////////////////////////////////
  // Central arms
  //////////////////////////////////////////
  HeadReco *head    = new HeadReco();
  //head->SetRawDataCheck(raw); // add the rawdatacheck pointer so a list of 
                              // bad packets is added to the EventHeader
  SubsysReco *trig    = new TrigReco();
  SubsysReco *spin    = new SpinReco();
  SubsysReco *peve    = new PreviousEventReco();

  //////////////////////////////////
  // Define the triggers
  //////////////////////////////////

  ////////////////////////////////// 
  // Accounting
  ////////////////////////////////// 
  SubsysReco *trigacc = new TriggerAccounting();
  SubsysReco *outacc = new OutputAccounting();

  //////////////////////////////////////////
  // Central arms
  //////////////////////////////////////////
  BbcReco *bbc     = new BbcReco();
  bbc->setBbcVtxError( 1.0 );

  SubsysReco *ert     = new ErtReco();
  SubsysReco *zdc     = new ZdcReco();
  SubsysReco *t0      = new T0Reco();
  SubsysReco *pad     = new PadReco();
  SubsysReco *vtx     = new VtxReco();
  SubsysReco *emc     = new EmcReco3();
  SubsysReco *emcres  = new EmcClusterContainerResurrector();
  //SubsysReco *padvtx  = new PadVtxReco();
  SubsysReco *tof     = new TofReco();
  SubsysReco *aero    = new AccReco();
  SubsysReco *dch     = new DchReco();
  SubsysReco *crk     = new CrkReco();
  CglReco *cgl     = new CglReco();
  SubsysReco *aerocl  = new AccclusterReco();
  SubsysReco *ring    = new RingReco();
  SubsysReco *tofw    = new TofwReco();
  SubsysReco *global  = new GlobalReco(); 
  SubsysReco *global_central  = new GlobalReco_central(); 
  SubsysReco *mpc     = new MpcReco();

  SvxParManager *svxpar = new SvxParManager();
  //svxpar->set_BeamCenter(-0.2102,0.1887);//

  SvxDecode *svxdecode = new SvxDecode();
  svxdecode->includePixel(true);
  svxdecode->includeStripixel(true);
  svxdecode->setAdcOffset(24);
  svxdecode->setAdcCutoff(-24);

  SubsysReco *svxhotdead                   = new SvxApplyHotDead();
  SubsysReco *svxrec                       = new SvxReco();
  SvxPriVertexSeedFinder *svxvtxseedfinder = new SvxPriVertexSeedFinder();
  SvxStandAloneReco *svxstandalone         = new SvxStandAloneReco();
  svxstandalone->setVertexRecoFlag(2);
  svxstandalone->setPPFlag(false); /// takashi set is a true, but it takes much longer....
  SvxPrimVertexFinder *svxprimvtxfinder    = new SvxPrimVertexFinder();
  

  ////////////////////////////////// 
  // Register SubSystems 
  ////////////////////////////////// 
  se->registerSubsystem(head);
  se->registerSubsystem(sync);
  se->registerSubsystem(trig);
  se->registerSubsystem(trigacc);

  ///////////////////////////////////////////
  /// Trigger dicing now in loaded TrigSelect.C macro
  ///////////////////////////////////////////
  TrigSelect();

  se->registerSubsystem(outacc);

  se->registerSubsystem(peve);
  se->registerSubsystem(bbc);
  se->registerSubsystem(zdc);
  se->registerSubsystem(t0);
  se->registerSubsystem(pad);
  se->registerSubsystem(vtx);
  se->registerSubsystem(emc);
  se->registerSubsystem(emcres);
  //se->registerSubsystem(padvtx);
  se->registerSubsystem(ert);
  se->registerSubsystem(tof);
  se->registerSubsystem(aero);
  se->registerSubsystem(dch);
  se->registerSubsystem(crk);
  se->registerSubsystem(tofw);

  se->registerSubsystem(cgl);

  se->registerSubsystem(aerocl);
  se->registerSubsystem(ring);
  se->registerSubsystem(mpc);

  se->registerSubsystem(spin);
 
  se->registerSubsystem(global);
  se->registerSubsystem(global_central);

//MasterRecalibratorManager *mrm = new MasterRecalibratorManager("MASTERRECALIBRATORMANAGER");
  //mrm->FillHistos(0);
  //MasterRecalibrator *mr = new MasterRecalibrator();
  //mr->Unlock(0);
  //Recalibrator *escale = mr->getRecalibrator("EmcGenericEScaleRecalReco");

  Recalibrator *ReadbackCompactPWG = new ReadbackCompactPWG();
  ReadbackCompactPWG->SetInputNode("emcClusterContainer");
  ReadbackCompactPWG->SetBaseClass("emcClusterContainer");
  se->registerSubsystem(ReadbackCompactPWG);
  Recalibrator *ReadbackCompactCNT = new ReadbackCompactCNT();
  ReadbackCompactCNT->SetInputNode("PHCentralTrack");
  ReadbackCompactCNT->SetBaseClass("PHCentralTrack");
  //se->registerSubsystem(ReadbackCompactCNT);
  Recalibrator *escale = new EmcGenericEScaleRecalReco();
  escale->SetInputNode("emcClusterContainer");
  escale->SetBaseClass("emcClusterContainer");
  se->registerSubsystem(escale);
  Recalibrator *centr = new Run12CuAu200GeVCentralityReco();
  centr->SetInputNode("PHGlobal");
  centr->SetBaseClass("PHGlobal");
  se->registerSubsystem(centr);


  /*se->registerSubsystem(svxpar);
  se->registerSubsystem(svxdecode);
  se->registerSubsystem(svxhotdead);
  se->registerSubsystem(svxrec);
  se->registerSubsystem(svxvtxseedfinder);
  se->registerSubsystem(svxstandalone); /// vtx determination is a little better with it on,  but works a decade longer 
  se->registerSubsystem(svxprimvtxfinder);
*/
  SkipEvents* skipevents = new SkipEvents(0); // select events in a vertex range
  //skipevents->Verbosity(1);
  skipevents->setZVertexRange(-10.0,10.0);
  //skipevents->setZVertexRange(4.0,6.0);
  se->registerSubsystem(skipevents);

//  DstMerger *sel = new DstMerger();
//  sel -> setVtxAsciiFileName(vtxout);
//  se->registerSubsystem(sel);

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

  /// SvxCentralTrackReco should be called after PHCentralTrack is reconstructed.
  SubsysReco* svxcentraltrack = new SvxCentralTrackReco();
  //se->registerSubsystem(svxcentraltrack);

  SubsysReco* svxcentraltrackbg = new SvxCentralTrackReco("SVXCENTRALTRACKRECOBACK");
  (dynamic_cast<SvxCentralTrackReco*>svxcentraltrackbg)->RndmAssocDchFlag(1);
  svxcentraltrackbg->Verbosity(0);
  //se->registerSubsystem(svxcentraltrackbg);

  // SvxSelectClusters* svxselect = new SvxSelectClusters();
  //  se->registerSubsystem(svxselect);

  // svx compactCNTs
  FillSvxHits *fillsvxhits = new FillSvxHits();
  //  fillsvxhits->Verbosity(0);
  fillsvxhits->setSaveOnlySelectedHits(false);
  //se->registerSubsystem(fillsvxhits);

  //////////////////////////////////
  // Central arm output
  //////////////////////////////////
  // trgsel is a global vector from TrigSelect.C which
  // contains the names of the trigger selectors
  // which are used to determine the filenames
  char dstname[100];
  for (int i = 0; i < trgsel.size(); i++)
    {
/*
      sprintf(dstname, "CNT_%s", trgsel[i].c_str());
      CNT_Compact(runnumber, segnumber, dstname, trgsel[i].c_str());
      sprintf(dstname, "DST_EVE_%s", trgsel[i].c_str());
      DST_EVE(runnumber, segnumber, dstname, trgsel[i].c_str());
      sprintf(dstname, "DST_SVX_%s", trgsel[i].c_str());
      DST_SVX(runnumber, segnumber, dstname, trgsel[i].c_str());
      //      sprintf(dstname, "DST_MPC_%s", trgsel[i].c_str());
      //      DST_MPC(runnumber, segnumber, dstname, trgsel[i].c_str());
*/
      sprintf(dstname, "CNTmerge_%s", trgsel[i].c_str());
      CNT_merge(runnumber, segnumber, dstname, trgsel[i].c_str());
    }

  //MakeOutput(runnumber,segnumber,"Histos");
  
/*
  SubsysReco *svxproductionqa  = new SvxQAForProduction();
  se->registerSubsystem(svxproductionqa);
  
  SubsysReco *svxalignmentqa  = new SvxAlignment_QA();
  se->registerSubsystem(svxalignmentqa);
*/

  ///////////////////////////////////////////
  // Analyze the Data.
  //////////////////////////////////////////
  gSystem->Exec("ps -o sid,ppid,pid,user,comm,vsize,rssize,time");
  /*
  DeathToMemoryHogs * dt = new DeathToMemoryHogs();
  dt->event_frequency(100); // check every hundred events
  dt->SaveHisto();  // save data in histogram
  se->registerSubsystem(dt);
  */
  Fun4AllInputManager *in = new Fun4AllPrdfInputManager("PRDFin");
//Fun4AllDstInputManager *in = new Fun4AllDstInputManager("DST");//real data tree  
in->fileopen(inputfile);
  se->registerInputManager(in);
  se->run(nEvents);

  se->End();
/*
  SVXQA_IOManager(runnumber,segnumber);
//  se->dumpHistos(histofile);
*/
  PrintTrigSelect();

  int evts = se->PrdfEvents();

  std::cout << "Total Events:  " << evts << std::endl;


  cout << "Successfully Completed Analysis." << endl;
  PHTimeServer::get()->print_stat();

  delete se;

/*
  char sql[1000];
  sprintf(sql,"echo \"update %s set nevents=%d,time='now' where runnumber=%d and sequence=%d and prodtag = \'%s\';\" | isql phnxreco_odbc",gSystem->Getenv("PRODDBTABLE"),evts,runnumber,segnumber,gSystem->Getenv("PRODTAG"));
  std::cout << "Update DB statement: " << sql << std::endl;
  gSystem->Exec(sql);
*/

  gSystem->Exec("ps -o sid,ppid,pid,user,comm,vsize,rssize,time");

  cout<<"Fun4All successfully completed "<<endl;
  gSystem->Exit(0);

}
