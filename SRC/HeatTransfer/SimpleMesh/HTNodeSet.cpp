#include <HTNodeSet.h>

HTNodeSet::HTNodeSet(int tag):TaggedObject(tag), NodeID(0)
{
	
}

HTNodeSet::~HTNodeSet()
{
	//
}

int
HTNodeSet::addNodeID(ID& newNodeID)
{
	if (NodeID==0) {
    NodeID=newNodeID;
  } else {
    opserr<< "HTNodeSet::AddNodeID is not able to merge node id"<<newNodeID<<endln;
  }
	return 0;
}

const ID&
HTNodeSet::getNodeID()
{
	return NodeID;
}


