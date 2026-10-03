#include "TObject.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "PHTypedNodeIterator.h"
#include "PHNodeIterator.h"
#include "PHNodeReset.h"
#include "recoConsts.h"

#include "RunHeader.h"
#include <cmath>

#include "GenericSimreco.h"

using namespace std;

typedef PHIODataNode<PHObject> PHObjectNode_t;
typedef PHIODataNode<PHTable> PHTableNode_t;
typedef PHIODataNode<TObject> TObjectNode_t;
typedef PHIODataNode<RunHeader> RunHeaderNode_t;

GenericSimreco::GenericSimreco(const char *name)
{
  ThisName = name;
}

int GenericSimreco::Init(PHCompositeNode *topNode)
{
  return 0;
}

int GenericSimreco::InitRun(PHCompositeNode *topNode)
{
  return 0;
}

int GenericSimreco::CreateNodeTree(PHCompositeNode *topNode)
{
  return 0;
}

int GenericSimreco::process_event(PHCompositeNode *topNode)
{
  return 0;
}

int GenericSimreco::ResetEvent(PHCompositeNode *topNode)
{
  return 0;
}
