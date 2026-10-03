void QAInit()
{
  gSystem->Load("libdstqa_muons");
  gSystem->Load("libTHmul.so");

  SubsysReco* qamut = new QAMut();
  // SubsysReco* qamui = new QAMui();

  Fun4AllServer *se = Fun4AllServer::instance();

  se->registerSubsystem(qamut);
  // se->registerSubsystem(qamui);  when this runs it generates errors every event, commented out by cmv 03/04/11, following CesarLdS's suggestion 
}
