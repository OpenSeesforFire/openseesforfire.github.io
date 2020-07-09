#include <HTConstants.h>

HTConstants::HTConstants(int tag, Vector& constants):TaggedObject(tag){

	if(constants.Size()!=5){
		opserr<<"WARNING: There should have 4 constants input from users " <<this->getTag()<<endln;
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
	
	return Constants;
}


