#include <HTConstants.h>

HTConstants::HTConstants(int tag, Vector& constants):TaggedObject(tag){

	if(constants.Size()!=6){
		opserr<<"WARNING: There should have 6 constants stored for HTConstants " <<this->getTag()<<endln;
		}
	else 
		Constants=constants;
}

HTConstants::~HTConstants()
{
	//
}

const Vector&
HTConstants::getConstants()
{
	if(Constants.Size()!=6)
		Constants=0;

	return Constants;
}


