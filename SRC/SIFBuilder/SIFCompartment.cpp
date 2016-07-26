/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                   
/**********************************************************************
** This project is aiming to provide a Tcl interface to define models  **
** for simulating structural behaviours under fire action.           **
** Developed by:                  `                                   **
**   Liming Jiang (liming.jiang@ed.ac.uk)                            **
**   Praven Kamath(Praveen.Kamath@ed.ac.uk)                          **
**   Xu Dai(X.Dai@ed.ac.uk)                                          **
**   Asif Usmani(asif.usmani@ed.ac.uk)                               **
**********************************************************************/
// $Revision: 2.4.0.1 $
// This file constructs the SIFCompartment to store the corresponding information.
// Created by Liming Jiang @UoE

#include <SIFCompartment.h>

SIFCompartment::SIFCompartment(int tag, const ID& CompInfo, const Vector& CompOrigin):TaggedObject(tag),theCompInfo(0),
               theXBeamTags(0), theYBeamTags(0), theSecXBeamTags(0), theSecYBeamTags(0),theColumnTags(0),theSlabTags(0),theCompOrigin(0)
{
		theCompInfo = ID(3);
		theCompInfo= CompInfo;

	    theXBeamTags = new ID();
		theYBeamTags = new ID();
		theSecXBeamTags = new ID();
		theSecYBeamTags = new ID();

		theColumnTags = new ID();
		theSlabTags = new ID();

		theCompOrigin = Vector(3);
		theCompOrigin = CompOrigin;
	
}
		

SIFCompartment::~SIFCompartment()
{
	if(theXBeamTags!=0){
		delete theXBeamTags;
		theXBeamTags=0;
	}
	if(theYBeamTags!=0){
		delete theYBeamTags;
		theYBeamTags=0;
	}
	if(theSecXBeamTags!=0){
		delete theSecXBeamTags;
		theSecXBeamTags=0;
	}
	if(theSecYBeamTags!=0){
		delete theSecYBeamTags;
		theSecYBeamTags=0;
	}
	if(theColumnTags!=0){
		delete theColumnTags;

		theColumnTags=0;
	}
	if(theSlabTags!=0){
		delete theSlabTags;
		theSlabTags=0;
	}

}

int 
SIFCompartment::AddXBeam(int MemberTag){

	if(theXBeamTags->Size()==0){
		theXBeamTags->resize(1);
		(*theXBeamTags)(0)= MemberTag;
	}
	else{
		int XBeamTagSize = theXBeamTags->Size();
		for(int i=0; i<theXBeamTags->Size(); i++){
			if((*theXBeamTags)(i)== MemberTag)
				return 0;
		}
		theXBeamTags->resize(XBeamTagSize+1);
		(*theXBeamTags)(XBeamTagSize)= MemberTag;
	}
	return 0;
}

const ID&
SIFCompartment::getConnectedXBeams(void){
	
	return (*theXBeamTags);

}


int 
SIFCompartment::AddYBeam(int MemberTag){

	if(theYBeamTags->Size()==0){
		theYBeamTags->resize(1);
		(*theYBeamTags)(0)= MemberTag;
	}
	else {
		int YBeamTagSize = theYBeamTags->Size();
		for(int i=0; i<theYBeamTags->Size(); i++){
			if((*theYBeamTags)(i)== MemberTag)
				return 0;
		}
		theYBeamTags->resize(YBeamTagSize+1);
		(*theYBeamTags)(YBeamTagSize)= MemberTag;
	}
	return 0;
}

const ID&
SIFCompartment::getConnectedYBeams(void){

	return (*theYBeamTags);
}

int 
SIFCompartment::AddXBeamSec(int MemberTag){

	if(theSecXBeamTags->Size()==0){
		theSecXBeamTags->resize(1);
		(*theSecXBeamTags)(0)= MemberTag;
	}
	else{
		int SecXBeamTagSize = theSecXBeamTags->Size();
		for(int i=0; i<theSecXBeamTags->Size(); i++){
			if((*theSecXBeamTags)(i)== MemberTag)
				return 0;
		}
		theSecXBeamTags->resize(SecXBeamTagSize+1);
		(*theSecXBeamTags)(SecXBeamTagSize)= MemberTag;
	}
	return 0;
}

const ID&
SIFCompartment::getConnectedSecXBeams(void){
	
	return (*theSecXBeamTags);

}

/*
int 
SIFCompartment::AddSecYBeam(int MemberTag){

	if(theSecYBeamTags->Size()==0){
		theSecYBeamTags->resize(1);
		(*theSecYBeamTags)(0)= MemberTag;
	}
	else{
		int SecYBeamTagSize = theSecYBeamTags->Size();
		for(int i=0; i<theSecYBeamTags->Size(); i++){
			if((*theSecYBeamTags)(i)== MemberTag)
				return 0;
		}
		theYBeamTags->resize(SecYBeamTagSize+1);
		(*theSecYBeamTags)(SecYBeamTagSize)= MemberTag;
	}
	return 0;
}

const ID&
SIFCompartment::getConnectedSecYBeams(void){
	
	return (*theSecYBeamTags);

}
*/



int 
SIFCompartment::AddColumn(int MemberTag){

	if(theColumnTags->Size()==0){
		theColumnTags->resize(1);
		(*theColumnTags)(0)= MemberTag;
	}
	else {
		int ColumnTagSize = theColumnTags->Size();
		for(int i=0; i<theColumnTags->Size(); i++){
			if((*theColumnTags)(i)== MemberTag)
				return 0;
		}
		theColumnTags->resize(ColumnTagSize+1);
		(*theColumnTags)(ColumnTagSize)= MemberTag;
	}
	return 0;
}

const ID&
SIFCompartment::getConnectedColumns(void){

	return (*theColumnTags);
}


int 
SIFCompartment::AddSlab(int MemberTag){

	if(theSlabTags->Size()==0){
		theSlabTags->resize(1);
		(*theSlabTags)(0)= MemberTag;
	}
	else {
		int SlabTagSize = theSlabTags->Size();
		for(int i=0; i<theSlabTags->Size(); i++){
			if((*theSlabTags)(i)== MemberTag)
				return 0;
		}
		theSlabTags->resize( MemberTag+1);
		(*theSlabTags)(SlabTagSize)= MemberTag;
	}
}

const ID&
SIFCompartment::getConnectedSlabs(void){

	return (*theSlabTags);
}

/*
int 
SIFCompartment::AddFireAction(int theFireActionID){

	if(theFireActionID==0){
		opserr<<"WARNING::SIFCompartment "<<this->getTag()<<" received an empty definition of SIFfireAction"<<endln;
		return -1;
	}
	return 0;
}
*/

const ID& 
SIFCompartment::getCompartmentInfo(void){

return theCompInfo;

}


const Vector&
SIFCompartment::getCompartmentOrigin(){

	return theCompOrigin;

}

