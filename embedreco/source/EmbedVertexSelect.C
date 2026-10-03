#include "EmbedVertexSelect.h"

#include <PHCompositeNode.h>
#include <PHIODataNode.h>
#include <PHNodeIterator.h>
#include <TriggerHelper.h>

#include <Fun4AllReturnCodes.h>
#include <Fun4AllServer.h>
#include <recoConsts.h>
#include <PHGlobal.h>
#include <getClass.h>
#include "VtxOut.h"
#include "PHPoint.h"

#include <cmath>
#include <string>
#include <vector>

#include "McEvalSingleList.h"

using namespace std;
using namespace findNode;

// Here we are comparing the zvertex of two events.  In a mismatch,
// DISCARDEVENT gets rid of just the real data event.

EmbedVertexSelect::EmbedVertexSelect(const string &name,const string &topnode): SubsysReco(name) {
  topNodeForThisSelector = topnode;
  RetCode = DISCARDEVENT;
  MatchRange = 3;//+-3cm
}

int EmbedVertexSelect::process_event(PHCompositeNode *topNode) {

  RetCode = EVENT_OK;
  recoConsts* rc = recoConsts::instance();
  Fun4AllServer* se = Fun4AllServer::instance();
  PHCompositeNode* mcnode = se->topNode(rc->get_CharFlag("EMBED_MC_TOPNODE"));
  PHCompositeNode* realnode = se->topNode(topNodeForThisSelector.c_str());

  /*  
  // Get McSingle
  McEvalSingleList* mcsingle = 0;
  PHTypedNodeIterator<McEvalSingleList> singleiter(mcnode);
  PHIODataNode <McEvalSingleList> *singlenode = singleiter.find("McSingle");
  if (singlenode) mcsingle = singlenode->getData();
  else cout << PHWHERE << "Not finding node McSingle" << endl;
  */

  if((!mcnode)||(!realnode)) {
    RetCode = ABORTEVENT;
  }
  else {
    PHGlobal* evt1 = getClass<PHGlobal>(mcnode,"PHGlobal");
    PHGlobal* evt2 = getClass<PHGlobal>(realnode,"PHGlobal");
    VtxOut *d_vtx = findNode::getClass<VtxOut>(realnode,"VtxOut");
    if ( !(evt1 ) || !(evt2 ) || !(d_vtx) ) {
      topNode->print();
      cout<<"missing global nodes "<<evt1<<" "<< evt2<<endl;
      RetCode = ABORTEVENT;
      //RetCode = DISCARDEVENT;
    }
    else {
      //      float ztrue = mcsingle->get_vertexz(0);
      float z1 = evt1 -> getBbcZVertex(); // MC from PHGlobal
      float z2 = (d_vtx->get_Vertex()).getZ(); // RD from VTXOUT
      //      cout << "VERTEXZ = " << ztrue
      if (verbosity > 0) cout << "VtxMC= " << z1<<", VtxRD= "<<z2<< ", bbcz_real= "<<evt2 -> getBbcZVertex()<<endl;

      // MatchRange is currently set at 1cm in the macro. Agreement
      // should be well within this because pisaToDST.C was set to
      // use the MC z-vertex and z-vertex smearing was set to zero.
      if(fabs(z1-z2) < MatchRange) {
	return EVENT_OK;
      } 
      else {
	RetCode = ABORTEVENT;
	//RetCode = DISCARDEVENT;
      }
    }
  }
  if (verbosity > 0) {
    if (RetCode == DISCARDEVENT) {
      cout << ThisName << ": Discarding Event" << endl;
    }
    else if (RetCode == ABORTEVENT) {
      cout << ThisName << ": Aborting Event" << endl;
    }
    else {
      cout << ThisName << ": Unknown bad event return code: " << RetCode << endl;
    }
  }
  return RetCode;
}

int EmbedVertexSelect::SetReturnCode(const char *action)
{
  if (!strcmp(action,"DISCARD"))
    {
      cout << "EmbedVertexSelect action: Discard Event for specific IOManager" << endl;
      RetCode = DISCARDEVENT;
    }
  else if (!strcmp(action,"ABORT"))
    {
      cout << "EmbedVertexSelect action: Abort Event in reconstruction" << endl;
      RetCode = ABORTEVENT;
    }
  else if (!strcmp(action,"OKAY"))
    {
      cout << "EmbedVertexSelect action: Allways write Event (kind of pointless)" << endl;
      RetCode = EVENT_OK;
    }
  else
    {
      cout << PHWHERE << " Unknown option: " << action << endl;
      return -1;
    }
  return 0;
}
