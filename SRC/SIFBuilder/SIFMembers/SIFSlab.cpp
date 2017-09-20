#include <SIFSlab.h>

SIFSlab::SIFSlab(int tag, int typeTag, int jt1, int jt2, int jt3, int jt4,const ID& theMemInfo):
SIFMember(tag,10,typeTag, theMemInfo),theLoad(0)
{  
  ConnectedJoints = new ID(4);
  //ConnectedJoints->resize(4);
  (*ConnectedJoints)(0) = jt1;
  (*ConnectedJoints)(1) = jt2;
  (*ConnectedJoints)(2) = jt3;
  (*ConnectedJoints)(3) = jt4;
}


SIFSlab::SIFSlab(int tag, int typeTag, int jt1, int jt2, const ID& theMemInfo):
SIFMember(tag,11,typeTag, theMemInfo),theLoad(0)
{  
  ConnectedJoints = new ID(2);
  //ConnectedJoints->resize(2);
  (*ConnectedJoints)(0) = jt1;
  (*ConnectedJoints)(1) = jt2;
}


SIFSlab::~SIFSlab()
{
  //
}

void 
SIFSlab::addLoad(const Vector& theLoadVec)
{
if(theLoad==0){
		theLoad = theLoadVec;
	}

}
	
const Vector&
SIFSlab::getLoad()
{

	return theLoad;

}

/*
SIFMember*
SIFSlab::getConnectedBeam(void)
{
  if (TypeTag==2) {
    //return connectedSlab;
  }
  else
    return 0;
}
*/

