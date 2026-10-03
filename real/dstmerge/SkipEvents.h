#ifndef __SKIPEVENTS_H__
#define __SKIPEVENTS_H__

#include "SubsysReco.h"
#include <iostream>
#include <fstream>

class PHCompositeNode;

class SkipEvents: public SubsysReco {
  
 public:
  
  SkipEvents(int nevents);
  virtual ~SkipEvents() {}
  
  int Init(PHCompositeNode *topNode) {return 0;}
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  
  int Reset(PHCompositeNode *topNode) 		{return 0;}
  int ResetEvent(PHCompositeNode *topNode) 	{return 0;}
  void Print(const char *what) const 		{return;}
  void setZVertexRange(float a, float b) {minZVertex=a; maxZVertex=b;}

 protected:
  int Events2Skip;
  int EventNumber;
  float minZVertex;
  float maxZVertex;

  std::ofstream myfile;

};

#endif

