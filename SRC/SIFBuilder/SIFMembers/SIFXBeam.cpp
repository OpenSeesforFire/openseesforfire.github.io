#include <SIFXBeam.h>

SIFXBeam::SIFXBeam(int tag, int typeTag, int jt1, int jt2,const ID& theMemInfo ,double gamma):
SIFMember(tag,1,typeTag,theMemInfo),Gamma(gamma),theLoad(0),ConnectedSlab(0)
{
  //typeTag indicates composite beam or single beam
  //typeTag==1 : single beam;
  //typeTag==2 : composite beam
	if(SIFSectionPtr!=0){
		SIFSectionPtr = 0;
	}
  ConnectedJoints = new ID(2);
  ConnectedJoints->resize(2);
  (*ConnectedJoints)(0) = jt1;
  (*ConnectedJoints)(1) = jt2;
}




SIFXBeam::~SIFXBeam()
{
	//
}

void 
SIFXBeam::addLoad(const Vector& theLoadVec)
{
if(theLoad==0){
		theLoad = theLoadVec;
	}

}
	
const Vector&
SIFXBeam::getLoad()
{

	return theLoad;

}


int
SIFXBeam::setConnectedSlab(int slabID)
{
	ConnectedSlab = slabID;
	return 0;
}
int
SIFXBeam::getConnectedSlab(void)
{
    return ConnectedSlab;
}



double
SIFXBeam::getGamma(void)
{
  
  return Gamma;
  
}
