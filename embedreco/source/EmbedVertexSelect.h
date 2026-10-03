#ifndef __EMBEDVERTEXSELECT_H__
#define __EMBEDVERTEXSELECT_H__

#include <string>
#include <vector>

#include "SubsysReco.h"

class EmbedVertexSelect: public SubsysReco
{
 public:
  EmbedVertexSelect(const std::string &name = "EMBEDVERTEXSELECT", 
		    const std::string &topnode="TOP");
  virtual ~EmbedVertexSelect() {}

  int process_event(PHCompositeNode *topNode);

  int SetReturnCode(const char *action = "DISCARD");
  void SetVertexRange(const float Range) {MatchRange = Range;}

 protected:
  float MatchRange;
  int RetCode;
  std::string topNodeForThisSelector;
};

#endif
