#include "embedana.h"
//==============================================================


embedana::embedana(string filename, string filepath, string oscarpath, const int InID) : m_outFileName(filename)
{
  ThisName = "embedana";
  EventNumber = 0;
  InPartNumber = 0;
  InIDread = InID;
  event_container = nullptr;
  fill_TTree = 0;
  remove_cluster_hits = 0;
  local_filepath = filepath;
  local_oscarpath = oscarpath;
  vertexes.clear();
  
  InitParams();

  memset(InData_read, 0, 22000 * sizeof(InData));
}

//==============================================================

embedana::~embedana()
{
  delete event_container;
  std::cout << "Event Container was deleted" << std::endl;
}

//==============================================================

int embedana::Init(PHCompositeNode *topNode)
{

  std::cout << "embedana::Init started..." << std::endl;
  // OutputNtupleFile = new TFile(OutputFileName.c_str(),"RECREATE");
  // std::cout << "embedana::Init: output file " << OutputFileName << " opened." << std::endl;
  recoConsts *rc = recoConsts::instance();
  remove_cluster_hits = rc->get_IntFlag("Remove_cluster_hits", 0);
  const int fill_QA_cluster_hists = rc->get_IntFlag("Fill_QA_cluster_hists", 0);
  const int fill_QA_lepton_hists = rc->get_IntFlag("Fill_QA_lepton_hists", 0);
  fill_TTree = rc->get_IntFlag("Fill_TTree", 0);
  const int fill_d_dphi_hists = rc->get_IntFlag("Fill_d_dphi_hists", 0);
  const int fill_DCA_hists = rc->get_IntFlag("Fill_DCA_hists", 0);
  const int use_iden = rc->get_IntFlag("Use_ident", 0);
  const int do_track_QA = rc->get_IntFlag("Do_track_QA", 0);
  const int do_reveal_cluster = rc->get_IntFlag("Do_reveal_cluster", 0);
  const int fill_true_DCA = rc->get_IntFlag("Fill_true_DCA", 0);
  const int check_veto = rc->get_IntFlag("Check_Veto", 0);

  std::cout << "Remove_cluster_hits: " << remove_cluster_hits << std::endl;
  std::cout << "fill_QA_cluster_hists: " << fill_QA_cluster_hists << std::endl;
  std::cout << "fill_QA_lepton_hists: " << fill_QA_lepton_hists << std::endl;
  std::cout << "fill_TTree: " << fill_TTree << std::endl;
  std::cout << "fill_d_dphi_hists: " << fill_d_dphi_hists << std::endl;
  std::cout << "fill_DCA_hists: " << fill_DCA_hists << std::endl;
  std::cout << "use_iden: " << use_iden << std::endl;
  std::cout << "Do_track_QA: " << do_track_QA << std::endl;
  std::cout << "do_reveal_cluster: " << do_reveal_cluster << std::endl;
  std::cout << "fill_true_DCA: " << fill_true_DCA << std::endl;
  std::cout << "check_veto: " << check_veto << std::endl;

  TOAD toad("Run14AuAuLeptonComby");

  const std::string loc = toad.location("field_map.root");

  event_container = new DileptonAnalysis::MyEventContainer();
  event_container->InitEvent();
  //event_container->GetHistsFromFile(loc);
  event_container->CreateOutFileAndInitHists(m_outFileName.c_str(), fill_QA_lepton_hists, fill_QA_cluster_hists, fill_TTree, fill_d_dphi_hists,
                                             fill_DCA_hists, do_track_QA, do_reveal_cluster, fill_true_DCA, check_veto);

  if (fill_TTree)
    event_container->ResetTree();

  std::cout << "embedana::Init ended." << std::endl;
  return 0;
}

//==============================================================

int embedana::InitRun(PHCompositeNode *topNode)
{
  std::cout << "embedana::InitRun started..." << std::endl;

  
  const int error = ReadOrigPartMoms();
  if (error)
  {
    printf("Something wrong!!!!\n");
    return 1;
  }

  std::ifstream myfile(local_filepath.c_str());
  string line;
  if (myfile.is_open())
  {
    while (getline(myfile, line))
    {
      string s;
      std::stringstream ss(line);
      while (getline(ss, s, ' '))
      {
        vertexes.push_back(atof(s.c_str()));
      }
    }
    myfile.close();
  }

  std::cout << "embedana::InitRun ended." << std::endl;
  return 0;
}

//==============================================================
int embedana::ResetEvent(PHCompositeNode *topNode)
{
  event_container->ClearEvent();
  return 0;
}
//==============================================================

int embedana::process_event(PHCompositeNode *topNode)
{

  event_container->ClearEvent();
  static int first_event = 0;
  //if (first_event == 0)
  if (true)
  {
    std::cout << "embedana::process_event firstevent" << std::endl;
    topNode->print();
    first_event++;
  }
  
  if (EventNumber > (int)vertexes.size() / 4)
    return 0;
  
  // event info
  //--EventHeader*  evthdr  = getClass<EventHeader>(topNode,"EventHeader");

  PHCentralTrack *trk = getClass<PHCentralTrack>(topNode, "PHCentralTrack");
  PHDchTrackOut *phdchtrk = getClass<PHDchTrackOut>(topNode, "PHDchTrackOut");

  recoConsts *rc = recoConsts::instance();
  Fun4AllServer *se = Fun4AllServer::instance();
  PHCompositeNode *mcnode = se->topNode(rc->get_CharFlag("EMBED_MC_TOPNODE"));
  PHCentralTrack *trk_mc = getClass<PHCentralTrack>(mcnode, "PHCentralTrack");

  VtxOut *vtxout   = getClass<VtxOut>(topNode, "VtxOut");
  VtxOut *vtxoutmc = getClass<VtxOut>( mcnode, "VtxOut");
  
  int write_ttree = 1;
  if( true || (TMath::Abs( (vtxout->get_Vertex()).getZ() - (vtxoutmc->get_Vertex()).getZ() ) > 1 ))
  {
    std::cout<<"\n\nVTX missmatch: data : sim "<<(vtxout->get_Vertex()).getZ()<< " : " <<(vtxoutmc->get_Vertex()).getZ()<<"\n"<<std::endl;
    write_ttree = 0;
    //return 0;
  }

  SvxClusterList *svx = getClass<SvxClusterList>(topNode, "SvxClusterList");
  SvxClusterList *svxsim = getClass<SvxClusterList>(mcnode, "SvxClusterList");
  SvxGhitClusterList *svxembed = getClass<SvxGhitClusterList>(topNode, "SvxGhitClusterList");
  const RunHeader *runHDR = findNode::getClass<RunHeader>(topNode, "RunHeader");
  
  //std::cout << "RN, real embed and sim Nhits and Ntraks: " << runHDR->get_RunNumber() << " "<< svx->get_nClusters() << " " << svxsim->get_nClusters() << " " << svxembed->get_nGhitClusters() << " " << trk->get_npart() << " " << trk_mc->get_npart() << std::endl;
  // std::cout<<"event : "<<EventNumber<<"  "<<((evthdr!=NULL) ? evthdr->get_EvtSequence() : -1 ) <<endl;
  //std::cout<<"after print"<<endl;
  if (EventNumber >= 10)
  {
    CglTrack *cgl = getClass<CglTrack>(topNode, "CglTrack");
    cgl->ShutUp(1);
    PHTrackOut *phtrk = getClass<PHTrackOut>(topNode, "PHTrackOut");
    phtrk->ShutUp(1);
  }

  //////////////////////////////////////////////////////////////////////
PHCompositeNode *realnode = se->topNode(rc->get_CharFlag("EMBED_REAL_TOPNODE"));
  PHGlobal *globalCNT =
      findNode::getClass<PHGlobal>(realnode, "PHGlobal");
  
  if (!globalCNT)
      std::cout << "NO GLOBAL!!!!!!!!!!!!!!!\n";    
  
  const float bbc_vertex = globalCNT->getBbcZVertex();
  const float centrality = globalCNT->getCentrality();


  PHGlobal *mcglobalCNT =
      findNode::getClass<PHGlobal>(mcnode, "PHGlobal"); //yuri
  
  if (!mcglobalCNT)
      std::cout << "NO GLOBAL!!!!!!!!!!!!!!!\n";   //yuri
      
  const float mc_bbc_vertex = mcglobalCNT->getBbcZVertex(); //yuri 

  PHPointerList<PHEmbedMcRecoTrack> *embedtrk = getClass<PHPointerList<PHEmbedMcRecoTrack> >(topNode, "PHEmbedMcRecoTrack");
  DileptonAnalysis::MyEvent *event = event_container->GetEvent();
  event->SetPreciseX(vertexes[4 * EventNumber + 0]);
  event->SetPreciseY(vertexes[4 * EventNumber + 1]);
  event->SetPreciseZ(bbc_vertex);
  event->SetVtxZ(mc_bbc_vertex);///yuri
  event->SetCentrality(centrality);
  event->SetRunNumber(runHDR->get_RunNumber());

  const int nGtracks = InData_read[InPartNumber].nn;
  for (int igtrack = 0; igtrack < nGtracks; igtrack++)
  {  
    DileptonAnalysis::MyGenTrack gTrack;

    gTrack.SetID(InData_read[InPartNumber].id);
    gTrack.SetPx(InData_read[InPartNumber].px);
    gTrack.SetPy(InData_read[InPartNumber].py);
    gTrack.SetPz(InData_read[InPartNumber].pz);
    gTrack.SetVx(InData_read[InPartNumber].vx);
    gTrack.SetVy(InData_read[InPartNumber].vy);
    gTrack.SetVz(InData_read[InPartNumber].vz);
    event->AddGenTrack(&gTrack);
    InPartNumber++;
    if(true)
    std::cout<<"read track id, px, py, pz, vx, vy, vz:  "<<gTrack.GetID()<<" "<<gTrack.GetPx()<<" "<<gTrack.GetPz()<<" "<<gTrack.GetPy()<<" "
    <<gTrack.GetVx()<<" "<<gTrack.GetVy()<<" "<<gTrack.GetVz()<<std::endl;
  }
  

  //std::cout<<event->GetCentrality()<<"\t"<<emccont_merged->size()<<"\t"<<trk->get_npart()<<endl;
  if (false)
    std::cout <<"\033[32m"<< "x,y,z,cent: " << vertexes[4 * EventNumber] << " " << vertexes[4 * EventNumber + 1] << " " << vertexes[4 * EventNumber + 2] << " " << vertexes[4 * EventNumber + 3] << " "  <<"\033[0m"<< std::endl;
 if (true)
    std::cout << "n gen traks txt : read : " << nGtracks << " : " << event->GetNgentrack()<< std::endl;
 if(nGtracks!= (int)  event->GetNgentrack() ) std::cout<<"\nBad read "<<nGtracks<<" "<<event->GetNgentrack()<<"\n"<<std::endl;

    const emcClusterContainer* emccont_merged = findNode::getClass<emcClusterContainer>(topNode, "emcClusterContainer");
    const emcClusterContainer* emccont_sim = findNode::getClass<emcClusterContainer>(mcnode, "emcClusterContainer");
    const emcTowerContainer* twrcont_merged = findNode::getClass<emcTowerContainer>(topNode, "emcTowerContainer");
    const emcTowerContainer* twrcont_sim = findNode::getClass<emcTowerContainer>(mcnode, "emcTowerContainer");
std::cout<<event->GetCentrality()<<"\t"<<emccont_merged->size()<<"\t"<<trk->get_npart()<<endl;
    for (int imergedcluster = 0; imergedcluster < (int) emccont_merged->size(); imergedcluster++){
      emcClusterContent* emc_merged = emccont_merged->getCluster(imergedcluster);
      if (!emc_merged) continue;
      float ecore_merged = emc_merged->e();
      float e_sim = 0;
      if (ecore_merged < 0.3) continue; 
      float ecore_sim = 0;
      int sector_merged = (emc_merged->arm() == 0) ? emc_merged->sector() : 7 - emc_merged->sector();
      for (int icluster = 0; icluster < (int) emccont_sim->size(); icluster++){
        emcClusterContent* emc = emccont_sim->getCluster(icluster);
        if (!emc) continue;
        if (emc->ecore() < 0.1) continue;
        int sector = (emc->arm() == 0) ? emc->sector() : 7 - emc->sector(); 
        e_sim = emc->e();
        if (sector != sector_merged) continue;
        float dphi = emc_merged->phi() - emc->phi();
        if (dphi > M_PI)  dphi -= 2*M_PI;
        if (dphi < -M_PI) dphi += 2*M_PI;
        float dz = emc_merged->z() - emc->z();
        if (TMath::Abs(dphi) >= 0.1) continue;
        if (TMath::Abs(dz)   >= 8.0)  continue;

        for (int itower_merged = 0; itower_merged < emc_merged->multiplicity(); itower_merged++){
          emcTowerContent* twr_merged = twrcont_merged->findTower(emc_merged->towerid(itower_merged));
          if(!twr_merged) continue;
          for (int itower_sim = 0; itower_sim < emc->multiplicity(); itower_sim++){
            emcTowerContent* twr_sim = twrcont_sim->findTower(emc->towerid(itower_sim));
            if(!twr_sim) continue;
            if(emc_merged->towerid(itower_merged)==emc->towerid(itower_sim)){
              //"merged is sim"//?
              ecore_sim += twr_sim->Energy();
            }
          }
        }
        if(ecore_sim>0) break;
      }
      if (ecore_sim/ecore_merged>0) std::cout<<"check=================="<<"ratio = "<<ecore_sim/ecore_merged<<"\t"<<"e_sim = "<<e_sim<<"\t"<<"ecore_sim = "<<ecore_sim<<"\t"<<"ecore_merged = "<<ecore_merged<<endl;
   if(ecore_sim<0.01) continue; 
      DileptonAnalysis::MyCluster newCluster;
      newCluster.SetEcore(emc_merged->ecore());
      newCluster.SetSect(sector_merged);
      newCluster.SetEmcId(imergedcluster);
      newCluster.SetTower(sector_merged, emc_merged->iypos(), emc_merged->izpos() );
      newCluster.SetChi2(emc_merged->chi2());
      //Without electron: define Dep relative to the sim cluster we matched to
      //newCluster.SetDep(emc_merged->ecore()/emc->ecore());

      newCluster.SetZ(emc_merged->z());
      newCluster.SetX(emc_merged->x());
      newCluster.SetY(emc_merged->y());
      //newCluster.SetProb(emc_merged->prob_photon());

      //Without track direction: store merged-vs-sim residuals
      //newCluster.SetEmcdphi(dphi);
      //newCluster.SetEmcdz(dz);

      newCluster.SetArm(emc_merged->arm());  // more sensible than centrality
      newCluster.SetMcId(1);
      
      event->AddCluster(&newCluster);
    }
    
for (int icluster = 0; icluster < (int) emccont_sim->size(); icluster++){
        emcClusterContent* emc = emccont_sim->getCluster(icluster);
        if (!emc) continue;
        if (emc->ecore() < 0.1) continue;
        int sector = (emc->arm() == 0) ? emc->sector() : 7 - emc->sector();
        DileptonAnalysis::MyCluster newCluster;
        newCluster.SetEcore(emc->ecore());
        newCluster.SetSect(sector);
        newCluster.SetEmcId(icluster);
        newCluster.SetTower(sector, emc->iypos(), emc->izpos() );
        newCluster.SetChi2(emc->chi2());
        newCluster.SetZ(emc->z());
        newCluster.SetX(emc->x());
        newCluster.SetY(emc->y());
        newCluster.SetArm(emc->arm());  // more sensible than centrality
        newCluster.SetMcId(0);
        event->AddCluster(&newCluster);
}


  if(false){
    const emcClusterContainer* emccont = findNode::getClass<emcClusterContainer>(mcnode, "emcClusterContainer");
    const emcTowerContainer* twrcont = findNode::getClass<emcTowerContainer>(mcnode, "emcTowerContainer");
    if(!emccont||!twrcont) return 0;
    for (unsigned int iemc = 0; iemc < emccont->size(); iemc++){
      emcClusterContent* emc = emccont->getCluster(iemc);
      if(!emc||emc->ecore()<0.4) return 0;
      std::vector<double> towerid, twrenergy;
      for (int itower=0; itower < emc->multiplicity(); itower++){
        towerid.push_back(emc->towerid(itower));
	      emcTowerContent* twr = twrcont->findTower(emc->towerid(itower));
        if(!twr) continue;
        twrenergy.push_back(twr->Energy());
      }
      double twre = 0, twrecore = 0;
      //std::cout<<towerid.size()<<" "<<twrenergy.size()<<" "<<twrcont->size()<<std::endl;
      for (unsigned int itower=0; itower < twrenergy.size(); itower++){
          twre += twrenergy[itower];
      }
      for (unsigned int itower=0; itower < twrenergy.size(); itower++){
          if(twrenergy[itower]>0.02*twre) twrecore+=twrenergy[itower];
      }
      std::cout<<"simul ecore, e, twr ecore, e: "<<emc->ecore()<<" "<<emc->e()<<" "<<twrecore/0.932<<" "<<twre<<" "<<emc->multiplicity()<<std::endl;
    }
  }
 
  

  if (embedtrk == nullptr)
  {
    std::cout << " no PHEmbedMcRecoTrack" << std::endl;
  }
  else
  {
    //--cout<<"PHEmbedMcRecoTrack: length = "<<embedtrk->length()<<endl;
    for (unsigned int i = 0; i < embedtrk->length(); i++)
    {
      PHEmbedMcRecoTrack *embed = (*embedtrk)[i];
      //--cout<<"   id_g, _s, _r = "<<embed->get_dctrkidG()<<" "<<embed->get_dctrkidS()<<" "<<embed->get_dctrkidR()<<" ";
      //--cout<<" mom="<<embed->get_momG()<<" "<<embed->get_momR()<<endl;

      float momr = -9999, momr1 = -9999;
      int dcidr = embed->get_dctrkidR();
      if (0 <= dcidr && dcidr < (int)trk->get_npart())
      {
        PHSnglCentralTrack *sngl = trk->get_track(dcidr);
        momr = sngl->get_mom();
        momr1 = phdchtrk->get_momentum(dcidr);
        //--cout<<"  reco mom= "<<momr<<" "<<momr1<<endl;
      }
      else
      {
        std::cout << "out of range : " << dcidr << momr << momr1 << std::endl;
      }
    }

    int Ncnt = trk->get_npart();

    // ntp fill
    unsigned int Nembed = (unsigned int)(0.5 * embedtrk->length());

    for (unsigned int itrk = 0; itrk < Nembed; itrk++)
    {
      PHEmbedMcRecoTrack *embed = (*embedtrk)[itrk];
      // int idG = embed->get_dctrkidG();
      int idR = embed->get_dctrkidR();
      // int idS = embed->get_dctrkidS();
      if (false)
        std::cout << embed->get_bbccent() << " " << embed->get_ptG() * TMath::Cos(embed->get_phi0G()) << " " << embed->get_ptG() * TMath::Sin(embed->get_phi0G()) << " " << idR << " " << Nembed << std::endl;

      if (idR < 0 || Ncnt <= idR)
      {
        std::cout << " RealtrackID is out of range = " << idR << " : " << Ncnt << std::endl;
        continue;
      }

      PHSnglCentralTrack *sngl = trk->get_track(idR);
      if (!sngl)
        continue;
      /// SvxCentralTrack*    svxcntsngl = m_vsvxcnt[idR];
      for (unsigned int imctrk = 0; imctrk < trk_mc->get_npart(); imctrk++)
      {
        if( TMath::Abs( embed->get_momS() - trk_mc->get_mom(imctrk) ) < 0.002 ) trk_mc->set_mcid(imctrk,embed->get_partidG());
      }

      if(false) 
        std::cout<<"\n\n    TRACK DEADMAP CHECK: "<<applySingleTrackCut(trk, idR, event->GetPreciseZ(), event->GetCentrality(), event->GetRunNumber())<<"\n\n"<<std::endl;
      
      if(applySingleTrackCut(trk, idR, event->GetPreciseZ(), event->GetCentrality(), event->GetRunNumber())<0) continue;

      /*float ntp[100];
      for (int i = 0; i < 100; i++)
        ntp[i] = -9999.;

      float charge = sngl->get_charge();
    
      ntp[0] = vtxout->get_ZVertex(); // zvtx
      ntp[1] = sngl->get_mom() * sin(sngl->get_the0());
      ntp[2] = charge;
      ntp[3] = sngl->get_mom();
      ntp[4] = sngl->get_phi0();
      ntp[5] = sngl->get_the0();
      ntp[6] = sngl->get_quality();
      ntp[7] = sngl->get_emcdz();
      ntp[8] = sngl->get_emcdphi();
      ntp[9] = sngl->get_pc3dz();
      ntp[10] = sngl->get_pc3dphi();
      ntp[11] = sngl->get_ecore();
      ntp[12] = sngl->get_prob();
      ntp[13] = sngl->get_n0();
      ntp[14] = sngl->get_npe0();
      ntp[15] = sngl->get_chi2() / sngl->get_npe0();
      ntp[16] = sngl->get_disp();

      // if(svxcntsngl!=NULL){
      //   float ndf        = svxcntsngl->getNDF();
      //   float chisqr     = svxcntsngl->getChiSquare();
      //   ntp[17] = svxcntsngl->getDCA2D(); // dca2d
      //   ntp[18] = svxcntsngl->getDCAZ();  // dcaz
      //   ntp[19] = (ndf>0) ? chisqr/ndf : -9999.; // chi2ndf
      //   ntp[20] = svxcntsngl->getHitPattern(); // hitptn
      // }

      ntp[21] = embed->get_ptS();
      ntp[22] = embed->get_momS();
      ntp[23] = embed->get_phi0S();
      ntp[24] = embed->get_the0S();
      ntp[25] = embed->get_emcecoreS();
      ntp[26] = embed->get_emcprobphotS();
      ntp[27] = embed->get_crknpmt0S();
      ntp[28] = embed->get_crknpe0S();
      ntp[29] = embed->get_crkchi2S() / embed->get_crknpe0S();
      ntp[30] = embed->get_crkdispS();

      // use gen info
      float genpt = embed->get_ptG();
      float genphi = embed->get_phi0G();
      float genthe = embed->get_the0G();
      float genvx = embed->get_xvtxG();
      float genvy = embed->get_xvtxG() * TMath::Tan(embed->get_phi0G());
      float genvz = embed->get_zvtxG();

      float beamxy[2] = {0., 0.};
      float fieldPol = 1.0;
      float gendca[3], gendca2d;

      calcDCA_BCbyCircleProjection(
          genpt, genphi, genthe,              // pt in xy-plane, phi, theta at inner most layer
          charge,                             // charge of the track
          genvx, genvy, genvz,                // hit position at inner most layer
          beamxy[0], beamxy[1],               // beam center
          fieldPol,                           //
          &gendca[0], &gendca[1], &gendca[2], // dca position
          &gendca2d                           // return
      );

      ntp[31] = embed->get_ptG();
      ntp[32] = embed->get_momG();
      ntp[33] = embed->get_phi0G();
      ntp[34] = embed->get_the0G();
      ntp[35] = embed->get_xvtxG();
      ntp[36] = embed->get_xvtxG() * TMath::Tan(embed->get_phi0G());
      ntp[37] = embed->get_zvtxG();
      ntp[38] = embed->get_partidG();
      ntp[39] = gendca2d;

      m_ntp_embed->Fill(ntp);*/

      DileptonAnalysis::MyElectron newElectron;
      newElectron.SetPt(embed->get_momS() * sin(embed->get_the0S()));
      newElectron.SetPtPrime(sngl->get_mom() * sin(sngl->get_the0()));
      newElectron.SetReconPT(embed->get_momG() * sin(embed->get_the0G()));
      newElectron.SetTrkId(itrk);
      newElectron.SetTrkQuality(sngl->get_quality());
      newElectron.SetArm(sngl->get_dcarm());
      newElectron.SetDCSide(sngl->get_dcside());
      newElectron.SetSect(sngl->get_sect());
      newElectron.SetQ(sngl->get_charge());
      newElectron.SetQPrime(sngl->get_charge());
      newElectron.SetPhiDC(sngl->get_phi());
      newElectron.SetPhi0(embed->get_phi0S());
      newElectron.SetThe0(embed->get_the0S());
      newElectron.SetPhi0Prime(sngl->get_phi0());
      newElectron.SetThe0Prime(sngl->get_the0());
      newElectron.SetZDC(sngl->get_zed());
      newElectron.SetAlpha(sngl->get_alpha());
      newElectron.SetAlphaPrime(sngl->get_alpha());
      newElectron.SetEmcId(sngl->get_emcid());
      newElectron.SetEcore(sngl->get_ecore());
      newElectron.SetDep(sngl->get_emce());
      newElectron.SetProb(sngl->get_prob());
      newElectron.SetEmcdz(sngl->get_emcdz());
      newElectron.SetEmcdphi(sngl->get_emcdphi());
      newElectron.SetEmcTower(sngl->get_sect(), sngl->get_ysect(), sngl->get_zsect());
      newElectron.SetTOFDPHI(sngl->get_pemcy());
      newElectron.SetTOFDZ(sngl->get_pemcz());
      newElectron.SetPC3SDPHI(sngl->get_pc3sdphi());
      newElectron.SetPC3SDZ(sngl->get_pc3sdz());
      newElectron.SetCrkphi(sngl->get_center_phi());
      newElectron.SetCrkz(sngl->get_center_z());
      //newElectron.SetTOFE(sngl->get_m2tof());
      newElectron.SetEmcTOF(sngl->get_temc());
      newElectron.SetChi2(sngl->get_chi2());
      newElectron.SetN0(sngl->get_n0());
      newElectron.SetNPE0(sngl->get_npe0());
      newElectron.SetDISP(sngl->get_disp());
      newElectron.SetMcId(embed->get_partidG());
      if(true) {
        //std::cout<< " tof m2 using two approaches"<<sngl->get_m2tof()<<" "<<SQR(sngl->get_mom())*( SQR(sngl->get_ttof())/SQR(sngl->get_pltof()) - 1)<<std::endl;
        if(embed->get_toftS()>-99 && embed->get_pltofS()>-99) 
        newElectron.SetTOFE(SQR(sngl->get_mom())*( SQR(embed->get_toftS()/embed->get_pltofS()*30.) - 1));
        else newElectron.SetTOFE(-9999);
      }
      if(false) std::cout<< "\033[32m"<<"embed trk pt n0 e/p cent: "<<newElectron.GetPt()<<" "<<newElectron.GetN0()<<" "<<newElectron.GetEcore()/newElectron.GetPtot()<<" "<<event->GetCentrality()<< "\033[0m"<<std::endl;
      if(false)
      {
        float min_dist = 99999;
        int n_had = 0;
        if (sngl->get_n0()>0)  n_had = (int) trk->get_npart();
        for (int icluster = 0; icluster < n_had; icluster++)
        {
          if(trk->get_center_phi(icluster)<-99) continue;
          if(applySingleTrackCut(trk, icluster, event->GetPreciseZ(), event->GetCentrality(), event->GetRunNumber())<-1) continue;
          if((int)icluster == idR) continue;
          const float dcenter_z = (sngl->get_center_z() - trk->get_center_z(icluster)) / 5.0;
          const float dcenter_phi = (sngl->get_center_phi() - trk->get_center_phi(icluster)) / 0.013;
          const float dist = sqrt(dcenter_z*dcenter_z + dcenter_phi*dcenter_phi);
          if(dist<min_dist)
          {
            min_dist = dist;
            newElectron.SetEmcdz_e(dcenter_z);
            newElectron.SetEmcdphi_e(dcenter_phi);
            newElectron.SetTOFDPHI(trk->get_n0(icluster)+10*((int)(10*trk->get_disp(icluster))));
            if(trk->get_npe0(icluster)>0)newElectron.SetTOFE(trk->get_chi2(icluster)/trk->get_npe0(icluster));
            newElectron.SetTOFDZ(trk->get_mom(icluster));
            if(false) std::cout<<newElectron.GetPtPrime()<<" "<<newElectron.GetN0()<<" "<<newElectron.GetDisp()<<" "<<newElectron.GetChi2()/newElectron.GetNpe0()<<" "
                              <<newElectron.GetEcore()/trk->get_mom(icluster)<<" "<<newElectron.GetTOFDPHI()<<" "<<newElectron.GetTOFE()<<std::endl;
          }
        }
      }
            
      event->AddTrack(&newElectron);
    }
  }
//////////////////////////single sim trk////////////////////////////
  for (unsigned int itrk = 0; itrk < trk_mc->get_npart(); itrk++)
  {
    
      DileptonAnalysis::MyElectron newElectron;
      if(applySingleTrackCut(trk_mc, itrk, event->GetPreciseZ(), event->GetCentrality(), event->GetRunNumber())<0) continue;

      newElectron.SetPt(trk_mc->get_pt(itrk));
      newElectron.SetPtPrime(trk_mc->get_pt(itrk));
      newElectron.SetReconPT(trk_mc->get_pt(itrk));
      newElectron.SetTrkId(itrk);
      newElectron.SetTrkQuality(trk_mc->get_quality(itrk));
      newElectron.SetArm(trk_mc->get_dcarm(itrk));
      newElectron.SetDCSide(trk_mc->get_dcside(itrk));
      newElectron.SetSect(trk_mc->get_sect(itrk));
      newElectron.SetQ(trk_mc->get_charge(itrk));
      newElectron.SetQPrime(trk_mc->get_charge(itrk));
      newElectron.SetPhiDC(trk_mc->get_phi(itrk));
      newElectron.SetPhi0(trk_mc->get_phi0(itrk));
      newElectron.SetThe0(trk_mc->get_the0(itrk));
      newElectron.SetPhi0Prime(trk_mc->get_phi0(itrk));
      newElectron.SetThe0Prime(trk_mc->get_the0(itrk));
      newElectron.SetZDC(trk_mc->get_zed(itrk));
      newElectron.SetAlpha(trk_mc->get_alpha(itrk));
      newElectron.SetAlphaPrime(trk_mc->get_alpha(itrk));
      newElectron.SetEmcId(trk_mc->get_emcid(itrk));
      newElectron.SetEcore(trk_mc->get_ecore(itrk));
      newElectron.SetDep(trk->get_emce(itrk));
      newElectron.SetProb(trk_mc->get_prob(itrk));
      newElectron.SetEmcdz(trk_mc->get_emcdz(itrk));
      newElectron.SetEmcdphi(trk_mc->get_emcdphi(itrk));
      newElectron.SetEmcTower(trk_mc->get_sect(itrk), trk_mc->get_ysect(itrk), trk_mc->get_zsect(itrk));
      newElectron.SetTOFDPHI(trk_mc->get_pemcy(itrk));
      newElectron.SetTOFDZ(trk_mc->get_pemcz(itrk));
      newElectron.SetPC3SDPHI(trk_mc->get_pc3sdphi(itrk));
      newElectron.SetPC3SDZ(trk_mc->get_pc3sdz(itrk));
      newElectron.SetCrkphi(trk_mc->get_center_phi(itrk));
      newElectron.SetCrkz(trk_mc->get_center_z(itrk));
      if (trk_mc->get_ttof(itrk)>-999 && trk_mc->get_pltof(itrk)>-999 )
        newElectron.SetTOFE(SQR(trk_mc->get_mom(itrk))*( SQR(trk_mc->get_ttof(itrk)*30./trk_mc->get_pltof(itrk)) - 1));
      else newElectron.SetTOFE(-9999);
      newElectron.SetEmcTOF(trk_mc->get_temc(itrk));
      newElectron.SetChi2(trk_mc->get_chi2(itrk));
      newElectron.SetN0(trk_mc->get_n0(itrk));
      newElectron.SetNPE0(trk_mc->get_npe0(itrk));
      newElectron.SetDISP(trk_mc->get_disp(itrk));
      newElectron.SetEmcdz_e(trk_mc->get_emcsdz_e(itrk));
      newElectron.SetEmcdphi_e(trk_mc->get_emcsdphi_e(itrk));
      newElectron.SetMcId(trk_mc->get_mcid(itrk));
      if(false) std::cout<< "\033[31m"<<"sim trk pt n0 e/p id cent: "<<newElectron.GetPt()<<" "<<newElectron.GetN0()<<" "<<newElectron.GetEcore()/newElectron.GetPtot()<<" "<<newElectron.GetMcId()<<" "<<event->GetCentrality()<< "\033[0m"<<std::endl;
      event->AddElecCand(&newElectron);
  }
  //////////////////////////end of tracks////////////////////////
  if (true) 
  //if ( event->GetNtrack() == 0 && event->GetNeleccand() == 0 )
    {////for conversion sim////
      event_container->FillTree();
      //if (fill_TTree && InIDread != 0) event_container->FillTree();
      EventNumber++;
      return 0;
    }
  /////////////////////////single sim hits///////////////////////
  if(false)std::cout<<(vtxout->get_Vertex()).getX()<<" "<<(vtxout->get_Vertex()).getY()<<" "<<(vtxout->get_Vertex()).getZ()<<std::endl;
  if(false)std::cout<<event->GetPreciseX()<<" "<<event->GetPreciseY()<<" "<<event->GetPreciseZ()<<" "<<event->GetCentrality()<<std::endl;
  for (int ihit = 0; ihit < svxsim->get_nClusters(); ihit++)
  {

    SvxCluster *svxhit = svxsim->get_Cluster(ihit);

    if (svxhit == nullptr)
    {
      std::cout << "cluster NULL : " << ihit << std::endl;
      continue;
    }

    DileptonAnalysis::MyVTXHit newHit;

    newHit.SetClustId(ihit);
    newHit.SetLayer(svxhit->get_layer());
    newHit.SetLadder(svxhit->get_ladder());
    newHit.SetSensor(1);
    newHit.SetXHit(svxhit->get_xyz_global(0));
    newHit.SetYHit(svxhit->get_xyz_global(1));
    newHit.SetZHit(svxhit->get_xyz_global(2));
    if (false)
      std::cout<<" 1 " << newHit.GetXHit() << " " << newHit.GetYHit() << " " << newHit.GetZHit()<<" "<<svxhit->get_adc(0)<<" "<<svxhit->get_adc(1) << std::endl;
    newHit.SetiLayerFromR();
    if (svxhit->get_layer() != newHit.GetLayer() || svxhit->get_ladder() != newHit.GetLadder() ||
        newHit.GetSensor() != 1)
    {
      std::cout << " smth is wrong " << std::endl;
      continue;
    }
    event->AddVTXHit(&newHit);
  }
  ///////////////embeded+sim hits////////////////
  std::vector<int> embed_ids;
  for (int ihit = 0; ihit < svxembed->get_nGhitClusters(); ihit++)
  {

    SvxGhitCluster *svxembeshit = svxembed->get_GhitCluster(ihit);
    SvxCluster *svxhit = svx->get_Cluster(svxembeshit->get_clusterID());
    if (svxhit == nullptr)
    {
      std::cout << "cluster NULL : " << ihit << std::endl;
      continue;
    }

    embed_ids.push_back(svxembeshit->get_clusterID());

    DileptonAnalysis::MyVTXHit newHit;

    newHit.SetClustId(ihit);
    newHit.SetLayer(svxhit->get_layer());
    newHit.SetLadder(svxhit->get_ladder());
    newHit.SetSensor(0);
    newHit.SetXHit(svxhit->get_xyz_global(0));
    newHit.SetYHit(svxhit->get_xyz_global(1));
    newHit.SetZHit(svxhit->get_xyz_global(2) - event->GetPreciseZ()*0);

    if(event->GetNVTXhit()>=1)
    {
      DileptonAnalysis::MyVTXHit *prevHit = event->GetVTXHitEntry(event->GetNVTXhit()-1);
      if (prevHit->GetXHit()==newHit.GetXHit() && prevHit->GetYHit()==newHit.GetYHit() && prevHit->GetZHit()==newHit.GetZHit())
      {
        if (false)std::cout<<"repeat"<<std::endl;
        continue;
      }
    }else std::cout<<"found bug"<<std::endl;
    if (false)
      std::cout << " 0 " << newHit.GetXHit() << " " << newHit.GetYHit() << " " << newHit.GetZHit() <<" "<<svxhit->get_adc(0)<<" "<<svxhit->get_adc(1) << std::endl;
    newHit.SetiLayerFromR();
    if (svxhit->get_layer() != newHit.GetLayer() || svxhit->get_ladder() != newHit.GetLadder() ||
        newHit.GetSensor() != 0)
    {
      std::cout << " smth is wrong " << std::endl;
      continue;
    }
    event->AddVTXHit(&newHit);
  }
  ////////////////all hits////////////////////
  int n_layers_hits[4]= {0, 0, 0, 0};
  int n_good_hits[4]= {0, 0, 0, 0};
  float average_z[4] = {0, 0, 0, 0};
  if (remove_cluster_hits)
  {
    for (int ihit = 0; ihit < svx->get_nClusters(); ihit++)
    {
      
      SvxCluster *svxhit = svx->get_Cluster(ihit);

      if (svxhit == nullptr)
      {
        std::cout << "cluster NULL : " << ihit << std::endl;
        continue;
      }

      bool already_used = false;
      for (unsigned int iid = 0; iid < embed_ids.size(); iid++)
      {
        if(ihit==embed_ids[iid]) already_used = true;
      }
      if(already_used) continue;
      n_layers_hits[svxhit->get_layer()]++;
     if (false) average_z[svxhit->get_layer()] += svxhit->get_xyz_global(2);

      DileptonAnalysis::MyVTXHit newHit;

      newHit.SetClustId(ihit);
      newHit.SetLayer(svxhit->get_layer());
      newHit.SetLadder(svxhit->get_ladder() + 50);
      newHit.SetSensor(2);
      newHit.SetXHit(svxhit->get_xyz_global(0));
      newHit.SetYHit(svxhit->get_xyz_global(1));
      newHit.SetZHit(svxhit->get_xyz_global(2));
      if (false)
        std::cout << newHit.GetXHit() << " " << newHit.GetYHit() << " " << newHit.GetZHit() << std::endl;
      newHit.SetiLayerFromR();
      if (svxhit->get_layer() != newHit.GetLayer() || svxhit->get_ladder() != newHit.GetLadder() - 50 ||
          newHit.GetSensor() != 2)
      {
        std::cout << " smth is wrong " << std::endl;
        continue;
      }
      event->AddVTXHit(&newHit);
    }
    if (false) std::cout << "\033[32m"
          << "average_z: " << average_z[0] / n_layers_hits[0] << " " << average_z[1] / n_layers_hits[1] << " "
          << average_z[2] / n_layers_hits[2] << " " << average_z[3] / n_layers_hits[3]
          << " at centrality and zvertex " << event->GetCentrality() << " " << event->GetPreciseZ() << "\033[0m" << std::endl;
    if(n_layers_hits[0]+n_layers_hits[1]<(n_layers_hits[2]+n_layers_hits[3])*0.1)
    {
      std::cout << "\033[31m" << "Too few hits in layers 0 and 1: " << n_layers_hits[0] + n_layers_hits[1] << " < "
                << (n_layers_hits[2] + n_layers_hits[3]) << " at centrality " << event->GetCentrality() << "\033[0m" << std::endl;
      EventNumber++;
      return 0;
    }
    event_container->Associate_Hits_to_Leptons();
    if(true)event_container->CleanUpHitList();
    for (int ihit = 0; ihit < event->GetNVTXhit(); ihit++)
    {
      DileptonAnalysis::MyVTXHit *svxhit = event->GetVTXHitEntry(ihit);

      if(svxhit->GetSensor() == 2) n_good_hits[svxhit->GetLayer()]++;
    }
    if(false)
    {
      std::cout << "\033[31m"
            << "n layers hits: " << n_layers_hits[0] << " " << n_layers_hits[1] << " "
            << n_layers_hits[2] << " " << n_layers_hits[3]
            << "\033[0m" << std::endl;
      std::cout << "\033[34m"
            << "n good hits: " << n_good_hits[0] << " " << n_good_hits[1] << " "
            << n_good_hits[2] << " "<< n_good_hits[3]
            << "\033[0m" << std::endl;
    }
  }

    



    //std::cout<< "\033[32m centrality "<< event->GetCentrality() <<" cluster size "<< emccont_merged->size()<<" \033[0m"<<std::endl;
    // for (int itrk =0; itrk < event->GetNtrack(); itrk++)
    // {
    //   DileptonAnalysis::MyElectron *elec = event->GetEntry(itrk);
    //   if(elec->GetEcore()<0.1)  continue;
    //   std::vector<float> usedclustersenergy;
    //   usedclustersenergy.push_back(elec->GetEcore());
    //   const float phi0_layer2 =     elec->GetPhi0() + (elec->GetPhiDC()-elec->GetPhi0())*0.1*2;///aprxomiation of phi of track in layers of vtx where most likely bremsstrahlung happens
    //   const float phi0_hit_layer2 = elec->GetPhi0() + (elec->GetPhiDC()-elec->GetPhi0())*0.1;///aprxomiation of phi of layers of vtx where most likely bremsstrahlung happens
    //   for ( int imergedcluster = 0; imergedcluster < (int) emccont_merged->size(); imergedcluster++)
    //   {
    //     emcClusterContent* emc_merged = emccont_merged->getCluster(imergedcluster);
    //     if(!emc_merged) continue;
    //     if(elec->GetEmcId() != imergedcluster && emc_merged->ecore()>0.1 )// && elec->GetEcore()*0.9 > emc_merged->ecore() ) 
    //     {
          
    //       float phi_clus = TMath::ATan2(emc_merged->y()-15*TMath::Sin(phi0_hit_layer2), emc_merged->x()-15*TMath::Cos(phi0_hit_layer2));
    //       phi_clus = phi_clus<-TMath::Pi()/2 ? phi_clus+2*TMath::Pi() : phi_clus;
    //       if(TMath::Abs(phi_clus-phi0_layer2)>0.03 || TMath::Abs(emc_merged->theta()-elec->GetThe0())>0.02)
    //       {
    //         continue;
    //       }
    //       bool already_used = false;
    //       for (unsigned int iid = 0; iid < usedclustersenergy.size(); iid++)
    //       {
    //         if(TMath::Abs(emc_merged->ecore()/usedclustersenergy[iid]-1)<0.01) already_used = true;
    //       }
    //       if(already_used) continue;
    //       usedclustersenergy.push_back(emc_merged->ecore());
    //       int sector_merged = emc_merged->arm() == 0 ? emc_merged->sector() : 7 - emc_merged->sector();
    //       DileptonAnalysis::MyCluster newCluster;
    //       newCluster.SetEcore(emc_merged->ecore());
    //       newCluster.SetSect(sector_merged);
    //       newCluster.SetEmcId(imergedcluster);
    //       newCluster.SetDep(emc_merged->ecore()/elec->GetEcore());
    //       newCluster.SetZDC(emc_merged->z());
    //       newCluster.SetPhi0(emc_merged->x());
    //       newCluster.SetThe0(emc_merged->y()); 
    //       newCluster.SetProb(emc_merged->prob_photon());
    //       newCluster.SetEmcdphi(phi_clus-phi0_layer2);
    //       newCluster.SetEmcdz(emc_merged->theta()-elec->GetThe0());
    //       newCluster.SetArm(event->GetCentrality()+100);
    //       newCluster.SetMcId(0);
    //       event->AddCluster(&newCluster);
    //     }
    //   }
    // } 

  if (fill_TTree && write_ttree )
    event_container->FillTree();

  EventNumber++;
  return 0;
}

int embedana::End(PHCompositeNode *topNode)
{
  std::cout << "Writing out..." << std::endl;
  event_container->WriteOutFile();
  std::cout << "Closing.." << std::endl;

  return 0;
}

int embedana::ReadOrigPartMoms()
{

  int i, k, nn, j;

  printf("Small oscar file is %s\n", local_oscarpath.c_str());

  fstream fin(local_oscarpath.c_str());
  string line;

  if (!fin)
  {
    printf("No such small file!!\n");
    return -1;
  }

  i = 0; k = 0; j = 0; nn = 0;
  while (getline(fin, line))
  {
    if (line.find("#") == 0)
      continue;
    int i1, i2, i3;
    double px, py, pz, vx, vy, vz, dummy;
    stringstream s(line);
    s >> i1 >> i2;
    if ( i2 > 0 && i2 < 6)
    {
       nn = i2;
       j = 0;
       k++;
    }
    if (TMath::Abs(i2) > 5 && TMath::Abs(i2)<9999)
    {
      s >> i3 >> px >> py >> pz >> dummy >> dummy >> vx >> vy >> vz;

      if(j>100) continue;
      //std::cout<<i2<<" "<<px<<" "<<py<<" "<<pz<<" "<<nn<<endl;
      if (i >= 100000)
      {
        printf("Too much particles in small file!\n");
        return -1;
      }

      InData_read[i].id = i2;
      InData_read[i].nn = nn;
      InData_read[i].px = px;
      InData_read[i].py = py;
      InData_read[i].pz = pz;
      InData_read[i].vx = vx;
      InData_read[i].vy = vy;
      InData_read[i].vz = vz;

      ++i;
      ++j;
    }
  }
  std::cout << "Nlines: " << i << ";  Nev: " << k << std::endl;
  if (k != 10000 )
  {
    printf("Wrong number of events");
    return -1;
  }

  return 0;
}

void embedana::calcDCA_BCbyCircleProjection(
    float pt, float phi, float the,  // pt in xy-plane, phi, theta at inner most layer
    int charge,                      // charge of the track
    float hx, float hy, float hz,    // hit position at inner most layer
    float vx, float vy,              // beam center
    float fieldPolarity,             //
    float *dx, float *dy, float *dz, // dca position
    float *d2dca_bc)                 // return
{
  static const float B = 0.90;
  static const float b = 0.003 * B;

  // vector to rotation center at hit1
  float R = pt / b;
  float pz = pt / tan(the);

  float b_sign = (fieldPolarity > 0) ? -1.0 : 1.0;
  float dir = ((b_sign)*charge > 0.) ? -1.0 : 1.0;
  float cx = hx - dir * R * sin(phi);
  float cy = hy + dir * R * cos(phi);

  // L is a distance btw the rotation center and primary vtx
  // psi is a angle of the vector starting from the rotation center to primary vertex
  float L = sqrt((vx - cx) * (vx - cx) + (vy - cy) * (vy - cy));
  float psi = atan2((vy - cy), (vx - cx));

  // DCA point
  *dx = cx + R * cos(psi);
  *dy = cy + R * sin(psi);

  float dphi = phi - dir * 0.5 * M_PI - psi;
  if (dphi > M_PI * 1.5)
    dphi -= M_PI * 2.;
  else if (dphi < -M_PI * 1.5)
    dphi += M_PI * 2.;
  float dzdphi = dir * pz * R / pt;
  *dz = hz + dphi * dzdphi;

  // DCA value
  *d2dca_bc = b_sign * (charge * (R - L));
}


int embedana::applySingleTrackCut(const PHCentralTrack *d_trk, const int itrk, const float vertex, const float centrality, const int run_number)
{
    const float px = d_trk->get_px(itrk);
    const float py = d_trk->get_py(itrk);
    const float pz = d_trk->get_pz(itrk);
    const float ecore = d_trk->get_ecore(itrk);
    const float p = sqrt(px * px + py * py + pz * pz);
    const float pT = sqrt(px * px + py * py);
    const float rich_phi = d_trk->get_center_phi(itrk);
    
    const float dep = d_trk->get_dep(itrk);
    //if (Z_GLOBAL > -99 && fabs(d_trk->get_zed(itrk)) >= Z_GLOBAL) std::cout<<"WTF:::: "<<d_trk->get_zed(itrk)<< " "<< Z_GLOBAL<<std::endl;
    if (Z_GLOBAL > -99 && TMath::Abs(d_trk->get_zed(itrk)) >= Z_GLOBAL)
        return -1;

    if( ( pT<0.09 || d_trk->get_quality(itrk)<0 || pT>25 )) return -2;
    //if(IsCentralSupportCut(d_trk->get_the0(itrk), vertex)) 
    //    return -1;

    if (DC_DEADMAP)
    {
        const int run_grp = get_rungroup(run_number);
        const int side = d_trk->get_dcside(itrk);
        const int arm = d_trk->get_dcarm(itrk);
        const float board = get_board(d_trk->get_phi(itrk));
        const float alpha = d_trk->get_alpha(itrk);
        for (int iarea = 0; iarea < MAX_DEAD_AREA; ++iarea)
        {
            if (dead_region(board, alpha, dcmap_xx1[run_grp][side][arm][iarea], dcmap_yy1[run_grp][side][arm][iarea], dcmap_xx2[run_grp][side][arm][iarea], dcmap_yy2[run_grp][side][arm][iarea], dcmap_xx3[run_grp][side][arm][iarea], dcmap_yy3[run_grp][side][arm][iarea], dcmap_xx4[run_grp][side][arm][iarea], dcmap_yy4[run_grp][side][arm][iarea]))
            {
                return -3;
            }
        }
    }
    return 0;
    bool pass_quality0 = true, pass_quality1 = true, pass_quality2 = true;
    if (QUALITY[0] > -99 && d_trk->get_quality(itrk) != QUALITY[0])
        pass_quality0 = false;
    if (QUALITY[1] > -99 && d_trk->get_quality(itrk) != QUALITY[1])
        pass_quality1 = false;
    if (QUALITY[2] > -99 && d_trk->get_quality(itrk) != QUALITY[2])
        pass_quality2 = false;
    if (!pass_quality0 && !pass_quality1 && !pass_quality2 && rich_phi<-99)
        return 0;
    
    if (E_PT > -99 && pT <= E_PT*0.9 &&  rich_phi<-99)
        return 0;

    if (MAX_PT > -99 && pT >= MAX_PT*1.05 &&  rich_phi<-99)
        return 0;
        
    if (dep < -2 && ecore < 0.3 && d_trk->get_n0(itrk)<0 /*&& rich_phi<-99*/)
        return 1;

    if(rich_phi>-99 && ((E_PT > -99 && pT <= E_PT*0.8) || (!pass_quality0 && !pass_quality1 && !pass_quality2) ||
        (MAX_PT > -99 && pT >= MAX_PT*1.2) )) return 3;
    //if (TEMC>-99 && d_trk->get_temc(itrk) > TEMC)
    //    return 1;
    if (N0 > -99 && d_trk->get_n0(itrk) <= N0)
        return 3;
    if (DISP > -99 && d_trk->get_disp(itrk) >= DISP)
        return 3;
    if (CHI2_NPE0 > -99 && (d_trk->get_chi2(itrk) / d_trk->get_npe0(itrk)) >= CHI2_NPE0)
        return 3;
    if (EMCDPHI > -99 && TMath::Abs(d_trk->get_emcdphi(itrk)) >= EMCDPHI)
        return 3;
    if (EMCDZ > -99 && TMath::Abs(d_trk->get_emcdz(itrk)) >= EMCDZ)
        return 3;
    if (EOVERP > -99 && ecore / p <= EOVERP)
        return 3;
    if (DEP[0] > -99 && dep <= DEP[0])
        return 3;
    if (DEP[1] > -99 && dep >= DEP[1])
        return 3;
    if(centrality<60)
    {
        if(pT > 0.40)
        {
            if (EMCSDPHI > -99 && TMath::Abs(d_trk->get_emcsdphi_e(itrk)) >= EMCSDPHI)
                return 3;
            if (EMCSDZ > -99 && TMath::Abs(d_trk->get_emcsdz_e(itrk)) >= EMCSDZ)
                return 3;
        }else{
            if (TMath::Abs(d_trk->get_emcdphi(itrk)) >= 0.02)
                return 3;
            if (TMath::Abs(d_trk->get_emcdz(itrk)-0.8) >= 8)
                return 3;
        }
    }
    if (PROB > -99 && d_trk->get_prob(itrk) <= PROB)
        return 3;
    return 2;
}
