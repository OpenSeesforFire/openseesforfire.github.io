#include <SIFYBeam.h>

SIFYBeam::SIFYBeam(int tag, int typeTag, int jt1, int jt2,const ID& theMemInfo, double gamma):
SIFMember(tag,2,typeTag, theMemInfo),Gamma(gamma),theLoad(0),ConnectedSlab(0)
{
  //typeTag indicates composite beam or single beam
  //typeTag==1 : single beam;
  //typeTag==2 : composite beam
	ConnectedJoints = new ID(2);
  ConnectedJoints->resize(2);
  (*ConnectedJoints)(0) = jt1;
  (*ConnectedJoints)(1) = jt2;
}




SIFYBeam::~SIFYBeam()
{
  //
}

void 
SIFYBeam::addLoad(const Vector& theLoadVec)
{
if(theLoad==0){
		theLoad = theLoadVec;
	}

}
	
const Vector&
SIFYBeam::getLoad()
{

	return theLoad;

}



int
SIFYBeam::setConnectedSlab(int slabID)
{
	ConnectedSlab = slabID;
	return 0;
}

int 
SIFYBeam::getConnectedSlab(void)
{
  
   return ConnectedSlab;

}



double
SIFYBeam::getGamma(void)
{
  
  return Gamma;
  
}
