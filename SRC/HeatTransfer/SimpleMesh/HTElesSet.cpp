#include <HTEleSet.h>

HTEleSet::HTEleSet(int tag):TaggedObject(tag), EleID(0)
{
	
}

HTEleSet::~HTEleSet()
{
	//
}

int
HTEleSet::addEleID(ID& newEleID)
{
	if (EleID==0) {
    EleID=newEleID;
  } else {
    opserr<< "HTEleSet::AddEleID is not able to merge ele id"<<newEleID<<endln;
  }
	return 0;
}

const ID&
HTEleSet::getEleID()
{
	return EleID;
}


