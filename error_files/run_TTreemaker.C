void run_TTreemaker(const char* inFile = "dst_00000.root", const char *outFile = "dummy.root", const int nevent = 0, const int mode = 1)
{
  // gSystem->Setenv("ODBCINI","/opt/phenix/core/etc/odbc.ini.test");
  
  gSystem->Load("liblvl2.so");
  gSystem->Load("libcompactCNT.so");
  gSystem->Load("libsimreco");
  gSystem->Load("libfun4all.so");
  gSystem->Load("libfun4allfuncs.so");

  gSystem->Load("libcteval.so");
  gSystem->Load("libemc.so");
  gSystem->Load("libemcEmbed4all.so");
  gSystem->Load("libreactionplane.so");

  gSystem->Load("librecal.so"); 
  gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libDileptonAnalysisEvent.so");
  gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libDileptonAnalysisTTreeMaker_noVTX.so");
  gSystem->Load("libTOAD");
    
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  //RC flags
  recoConsts *reco_consts =  recoConsts::instance();
  reco_consts->set_IntFlag("RD_TO_MC", 0);
  reco_consts->set_IntFlag("TRIGGER_BIT", 2);
  reco_consts->set_IntFlag("RD_RUN_SELECTION", 12);
  reco_consts->set_IntFlag("RD_SYSTEM_SELECTION", 0);
  reco_consts->set_DoubleFlag("RD_BBCZCUT_VAL", 10);
  reco_consts->set_IntFlag("QA_DEBUG", 0);
  reco_consts->set_IntFlag("RD_CENTCUT", 0);

  reco_consts->set_IntFlag("EMCGEOFLAG", 1);

  reco_consts->set_IntFlag("FILL_TRKS", 1);
  reco_consts->set_IntFlag("FILL_CLRS", 1);
  reco_consts->set_IntFlag("FILL_MCTRK", 0);
  reco_consts->set_IntFlag("FILL_MCCLR", 0);

  if (mode==1){
    // import simulated data (no embedding) and re-clusterize
    SubsysRecoStack *simimp = new EmcGeaContainerImporter();
    simimp->x_push_back( new EmcUnclusterizer() );
    EmcTowerScalerSmearer *emcsm = new EmcTowerScalerSmearer(1., 0.);

    //d+Au Run16

    //const double myPbScScale= 0.985;
    //const double myPbGlScale= 1.028;
    //const double myPbScSmear= 0.05;
    //const double myPbGlSmear= 0.12;
    
    //emcsm->SetScale(0.985, 0.985, 0.985, 0.985, 0.985, 0.985, 1.028, 1.028);
    //emcsm->SetSmear(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.12, 0.12);
    
    //p+p 200 Run 15                                                                                                                                                                             
    //emcsm->SetScale(0.998, 0.999, 0.997, 0.998, 1.007, 1.009, 1.099, 1.095); //RE                                                                                                              
    //emcsm->SetSmear(0.055, 0.054, 0.055, 0.055, 0.056, 0.055, 0.040, 0.037); //RE                                                                                                               
    //emcsm->SetSmear2(0.025, 0.025, 0.026, 0.024, 0.024, 0.026, 0.028, 0.028); //RE    
        
    //p+Au 200 Run 15
    //emcsm->SetScale(0.983, 0.978, 0.971, 0.9715, 0.968, 1.002, 1.037, 1.027); //RE                                        
    //emcsm->SetSmear(0.055, 0.066, 0.055, 0.059, 0.057, 0.056, 0.072, 0.072); //RE                                                 
    //emcsm->SetSmear2(0.011, 0.017, 0.011, 0.013, 0.012, 0.012, 0.063, 0.062); //RE
        
    //p+p 200 Run 15
    //emcsm->SetScale(0.968, 0.962, 0.969, 0.964, 0.963, 0.957, 1.001, 0.995);
    //emcsm->SetSmear(0.048, 0.066, 0.055, 0.050, 0.053, 0.056, 0.122, 0.122);
    //emcsm->SetScale(0.969, 0.966, 0.970, 0.966, 0.966, 0.964, 1.006, 1.006); 
    //emcsm->SetSmear(0.055, 0.066, 0.055, 0.059, 0.057, 0.056, 0.072, 0.072);
    //emcsm->SetSmear2(0.011, 0.017, 0.011, 0.013, 0.012, 0.012, 0.063, 0.062);

    //Au+Au 200 Run 14
    //emcsm->SetScale(0.988, 0.990, 0.985, 0.981, 0.980, 0.985, 1.025, 1.023);
    //emcsm->SetSmear(0.055, 0.066, 0.055, 0.059, 0.057, 0.056, 0.072, 0.072);
    //emcsm->SetSmear2(0.011, 0.017, 0.011, 0.013, 0.012, 0.012, 0.063, 0.062);

    simimp->x_push_back(emcsm);
    se->registerSubsystem(simimp);
    se->registerSubsystem( new EmcEmbedReclusterizer("TOP", "TOP", "TOP", "", 1) ); // 0 --kReal, 1 -- kPISA for mEmcGeometryModule
  }


   
  //  Reconstruction Modules...    
  //se->registerSubsystem( new RecalEMCalTOF() );
  se->registerSubsystem( new DileptonAnalysis::TTreeMaker_noVTX(outFile));

  Fun4AllDstInputManager *in1 = new Fun4AllDstInputManager("DSTin1","DST");
  in1->Verbosity(0);

  se->registerInputManager(in1);
  se->fileopen("DSTin1",inFile);

  //  OK...Do the analysis!!!!
  se->run(nevent);
  se->dumpHistos(outFile); 
  se->End();
  delete se;
}
