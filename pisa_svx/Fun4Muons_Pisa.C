void Fun4Muons_Pisa(
int nEvents=0,
const char *filein = "PISAEvent.root", // a regular file with J/psi from PYTHIA with FVTX hits
const char *evtgen = "pythia.root",
const char *dstout = "dst.root",
const char *ndstfile = "mwg.root",
const char *singlepdstout = "singlemuon_pdst.root",
const char *dimuonpdstout = "dimuon_pdst.root",
string lutfile = "FVTX_LUT.dat",
const char *qafile = "qaout.root",
int run_number = 365313 // run number should be larger than 365300 for FVTX
)
{
  char *l2_inputfile = 0;
  const char *prdfout = "data_out.prdf";
  const char *qa_file = "qa_out.root";

    if( filein ) cout << "Fun4Muons_Pisa - filein: " << filein << endl;
    if( l2_inputfile ) cout << "Fun4Muons_Pisa - l2_inputfile  : " << l2_inputfile << endl;
    if( dstout ) cout << "Fun4Muons_Pisa - dstout: " << dstout << endl;
    if( ndstfile ) cout << "Fun4Muons_Pisa - ndstout: " << ndstfile << endl;
    if( singlepdstout ) cout << "Fun4Muons_Pisa - pdstout: " << singlepdstout << endl;
    if( dimuonpdstout ) cout << "Fun4Muons_Pisa - pdstout: " << dimuonpdstout << endl;
    if( qa_file ) cout << "Fun4Muons_Pisa - qa_file: " << qa_file << endl;
    cout << "Fun4Muons_Pisa - run_number: " << run_number << endl;

    // if true, dst output gets written
    bool write_dst = false;
    bool write_pdst = true;
    bool write_ndst = true;
    bool write_qa = false;
    bool do_evaluation = false;

    // load libraries
    gSystem->Load("libPHGeant" );
    gSystem->Load("libfun4all");
    gSystem->Load("libsimreco_base");
    gSystem->Load("libmutoo_modules.so");
    gSystem->Load("libmutoo_classes.so");
    gSystem->Load("libmutoo_subsysreco.so");
    gSystem->Load("liblvl2");
    gSystem->Load("libPythia6.so");
    gSystem->Load("libPHPythia.so");
    gSystem->Load("libPHPythiaEventGen.so");

    TMutDatabaseCntrl::set_database_access("disable_HV",true);
    TMutDatabaseCntrl::set_database_access( "use_local_landau_parameters_file", true );
    TMutDatabaseCntrl::set_filename("use_local_landau_parameters_file",
				    "landau_parameters.txt" );

    bool use_local_hv_files = false;
    if( use_local_hv_files )
      {
	// generate hv file from run number
	char* hv_file = "mut.disabledAnodes.dat_run307084";
	TMutDatabaseCntrl::set_database_access("use_local_dead_HV_file",true);
	TMutDatabaseCntrl::set_filename("use_local_dead_HV_file", hv_file );
	
      } 
    else TMutDatabaseCntrl::set_database_access("use_local_dead_HV_file",false);
    
    TMutDatabaseCntrl::set_database_access("attenuated_channels",true);
    TMutDatabaseCntrl::set_database_access("dead_channels",true);
    //TMutDatabaseCntrl::set_database_access("dead_FEMs",true);
    
    ///////////////////////////////////////////
    // recoConsts setup
    //////////////////////////////////////////
    recoConsts *rc = recoConsts::instance();
    rc->set_IntFlag("SIMVERTEXFLAG",2);
    rc->set_FloatFlag("SIMZ0VERTEX",0);
    rc->set_FloatFlag("SIMT0VERTEX",0);
    rc->set_FloatFlag("SIMZ0VERTEXWIDTH",0.02);
    rc->set_FloatFlag("SIMT0VERTEXWIDTH",0);
    rc->set_IntFlag("SIMULATIONFLAG",2);
    rc->set_IntFlag("EMBEDFLAG",0);
    rc->set_IntFlag("RUNNUMBER", run_number );

    /*
    disable check of mapFileScale since its decided by pisa,
    not by the run number.
    */
    mMfmMT::setMapFileFlag( mMfmMT::MAP_3D_PLUS_PLUS );
    mMfmMT::setMapFileScale( 1.0 );
    MuonUtil::set_check_mapfile_scale( false );

    ///////////////////////////////////////////
    // Make the Server
    ///////////////////////////////////////////
    Fun4AllServer *se = Fun4AllServer::instance();

    ///////////////////////////////////////////
    // Subsystems
    //////////////////////////////////////////

    // run header and trigger setting
    se->registerSubsystem( new HeadSimreco() );
    //    se->registerSubsystem( new TrigSimreco() );

    // vertex simulation
    // puts the vertex from the pisa header node into vtxOut object
    gSystem->Load( "libsimreco_base.so" );
    VtxSimreco* vtx_reco = new VtxSimreco();
    vtx_reco->SmearZ( true );
    vtx_reco->UseXY( true );
    vtx_reco->OverwriteBBC( true );
    vtx_reco->XVertexSigma( 0.002 ); // VTX resolution of 20 microns
    vtx_reco->YVertexSigma( 0.002 );
    vtx_reco->ZVertexSigma( 0.002 );
    se->registerSubsystem( vtx_reco );

    // external vertex source
    //TMutExtVtx::get().set_verbosity( MUTOO::SOME );

    // counter
    MuonCounter* counter = new MuonCounter();
    counter->set_event_dump(1000);
    se->registerSubsystem( counter );

    // muon prdf unpacker
    MuonUnpackPisa* muon_unpack_pisa( new MuonUnpackPisa() );
    muon_unpack_pisa->Verbosity( MUTOO::SOME );
    muon_unpack_pisa->set_flag(MuonUnpackPRDF::SKIP_ZERO_SUPPRESSION,1);
    se->registerSubsystem( muon_unpack_pisa );

    // mutoo reconstruction
    se->registerSubsystem( new MuiooReco() );
    se->registerSubsystem( new MuonDev() );

    gSystem->Load("libfvtx_subsysreco");
    FvtxUnpackPisa *Fvtx_pisa = new FvtxUnpackPisa();
    se->registerSubsystem(Fvtx_pisa);

    FvtxReco* fvtxreco = new FvtxReco();
    se->registerSubsystem(fvtxreco);

    // create LUT for FVTXTrig
    //    gSystem->Load("libfvtx_rp.so");
    //    createFvtxLUT* fvtxlut = new createFvtxLUT(lutfile);
    //    se->registerSubsystem(fvtxlut);

    if (do_evaluation) se->registerSubsystem(new FvtxEval());
    
    // picoDST
    if( singlepdstout && write_pdst )
    {
	//	se->registerSubsystem( new MpcReco() );
	//	gSystem->Load("librxnp_subsysreco.so");
	//	se->registerSubsystem( new RxnpReco() );
	//	se->registerSubsystem( new RpSumXYReco() ); // recalibrator and rp doesn't work together! 

        // global Reco
        se->registerSubsystem( new GlobalReco() );

        // MWG
        gSystem->Load("libMWGOO.so");
        PHInclusiveNanoCuts *MWGcuts = new MWGInclusiveNanoCutsv2();
        se->registerSubsystem(new MWGFvtxReco(MWGcuts));

	gSystem->Load("libpicodst_object.so");
	mFillSingleMuonContainer* msngl = new mFillSingleMuonContainer();
	msngl->set_bbcz_cut(15.0);
	msngl->set_is_sim(true);
	se->registerSubsystem(msngl);

	mFillDiMuonContainer* mdi = new mFillDiMuonContainer(false); // do not make mixed events
	se->registerSubsystem(mdi);
	mdi->set_mass_cut(0.5);
	mdi->set_is_sim(true);
	
	if ( write_ndst )
	  {
	    Fun4AllOutputManager *outndst = new Fun4AllDstOutputManager("Outndst", ndstfile);
	    outndst->AddNode("Sync");
	    outndst->AddNode("TrigLvl1");
	    outndst->AddNode("VtxOut");
	    outndst->AddNode("PHGlobal");
	    outndst->AddNode("PHPythiaHeader");
	    outndst->AddNode("PHPythia");
	    outndst->AddNode("PHMuoTracksOO");
	    se->registerOutputManager(outndst);
	  }

	Fun4AllOutputManager *outsmu = new Fun4AllDstOutputManager("Outsmu",singlepdstout);
	outsmu->AddNode("Sync");
	outsmu->AddNode("SingleMuonContainer");
	outsmu->AddNode("VtxOut");
	outsmu->AddNode("PHGlobal");
	outsmu->AddNode("PHPythiaHeader");
	outsmu->AddNode("PHPythia");
	outsmu->AddEventSelector("mFillSingleMuonContainer");
	se->registerOutputManager(outsmu);

	Fun4AllOutputManager *outdimu = new Fun4AllDstOutputManager("Outdimu",dimuonpdstout);
	outdimu->AddNode("Sync");
	outdimu->AddNode("DiMuonContainer");
	outdimu->AddNode("PHPythiaHeader");
	outdimu->AddNode("PHPythia");
	outdimu->AddEventSelector("mFillDiMuonContainer");
	se->registerOutputManager(outdimu);
   }

    // qa
    if( write_qa )
    {
        gSystem->Load("libdstqa_muons.so");
	se->registerSubsystem( new QAMui() );
        se->registerSubsystem( new QAMut() );
    }

    if( false )
    {
        MuonDumpMap* dumper = new MuonDumpMap();
        dumper->add_map( "TMuiRoadMapO" );
        se->registerSubsystem( dumper );
    }

    ///////////////////////////////////////////
    // Input manager
    ///////////////////////////////////////////
    Fun4AllPisaInputManager *input_manager = new Fun4AllPisaInputManager("PisaIn");
    se->registerInputManager(input_manager);

    //    Fun4AllDstInputManager *in2 = new Fun4AllNoSyncDstInputManager("DSTin2","DST");
    //    se->registerInputManager(in2);
    
    // level2 input file
    if( l2_inputfile )
    {
        Fun4AllInputManager *l2_signal = new Fun4AllNoSyncDstInputManager("IML2DST","DST");
        se->registerInputManager(l2_signal);
        se->fileopen(l2_signal->Name(),l2_inputfile);
    }

    ///////////////////////////////////////////
    // Output manager
    ///////////////////////////////////////////
    if( dstout && write_dst )
    {

        Fun4AllDstOutputManager *dst_manager	= new Fun4AllDstOutputManager("DSTOUT",dstout);
        se->registerOutputManager(dst_manager);

        // Header and vertex nodes
        dst_manager->AddNode("RunHeader");
        dst_manager->AddNode("EventHeader");
        dst_manager->AddNode("VtxOut");
        dst_manager->AddNode("BbcOut");
	dst_manager->AddNode("BbcRaw");
	dst_manager->AddNode("ZdcOut");
	dst_manager->AddNode("ZdcRaw");
	dst_manager->AddNode("TrigLvl1");

        // mutoo/muioo nodes
        dst_manager->AddNode("TMCPrimary");
	dst_manager->AddNode("PHPythiaHeader");
	dst_manager->AddNode("PHPythia");
	dst_manager->AddNode("TMuiMCHitO");
	dst_manager->AddNode("TMutMCHit");
	dst_manager->AddNode("TMutMCTrk");
	dst_manager->AddNode("TMuiPseudoLL1");
	dst_manager->AddNode("PHMuoTracksOO");

        // From EVA node
        dst_manager->AddNode("header");
        dst_manager->AddNode("fkin");
        dst_manager->AddNode("primary");
        dst_manager->AddNode("pythia");

	// PHGlobal
	dst_manager->AddNode("PHGlobal");
	dst_manager->AddNode("PHGlobal_MUON");

	//	dst_manager->AddNode("RpSumXYObject");
    }

    ///////////////////////////////////////////
    // Analyze the Data.
    //////////////////////////////////////////

    // open first file
    se->fileopen(input_manager->Name(), filein );

    // if you want to save PHPYTHIA info point to the pythia file created and comment out the lines below
    //    if ( evtgen ) 
      //      se->fileopen(in2->Name(), evtgen );

    // run all events
    se->run(nEvents);
    se->End();

    if( write_qa )
    {
        Fun4AllHistoManager *hm = se->getHistoManager("QA");
        if( hm ) hm->setOutfileName( qa_file );
        else cout << "Fun4Muons_Pisa - unable to find QA histograms" << endl;
    }

    se->dumpHistos();

    delete se;

    cout << "Completed reconstruction." << endl;

}
