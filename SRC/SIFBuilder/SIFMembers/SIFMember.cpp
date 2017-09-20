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
// This file constructs the class SIFMember which defines the 
// configuration of structural memebers.
// Created by Liming Jiang @UoE


#include <SIFMember.h>

SIFMember::SIFMember(int tag, int memberTypeTag, int typeTag,const ID& theMemInfo):TaggedObject(tag),
TypeTag(typeTag), MemberTypeTag(memberTypeTag), SIFSectionPtr(0),IntNodeTags(0),
ConnectedJoints(0),IntEleTags(0),AppliedFireActions(0),theMemberInfo(0), DamVec(0)
{
  theMemberInfo = ID();
  theMemberInfo=theMemInfo;
}


SIFMember::SIFMember( ):TaggedObject(0),
TypeTag(0), MemberTypeTag(0), SIFSectionPtr(0),ConnectedJoints(0),IntNodeTags(0),
IntEleTags(0),AppliedFireActions(0),theMemberInfo(0),DamVec(0)
{
  
}



SIFMember::~SIFMember()
{
	if(AppliedFireActions!=0)
		AppliedFireActions=0;
}


//General function for all SIFMembers to assign the SIFSection;
int
SIFMember::assignSection(SIFSection* theSection)
{
  if (theSection!=0) {
    SIFSectionPtr = theSection;
  } else {
    opserr<<"WARNING:: the SIFSection pointer failed to be assigned to SIFMember "
           <<this->getTag()<<endln;
    return -1;
  }
  
  return 0;
}


//General function for all SIFMembers to get the SIFSection Pointer;
SIFSection*
SIFMember::getSIFSectionPtr(void)
{
		
    return SIFSectionPtr;
}

//General function for all SIFMembers to get typeTag;
int
SIFMember::getTypeTag(void)
{
  return TypeTag;
  
}

int
SIFMember::getMemberTypeTag(void)
{
  return MemberTypeTag;
  //1:xBeam,2:YBeam,3:Column,10:Slab, 21 xSecBeam, 11: slab modelled as plane frame
}


//General function for all SIFMembers to get connected SIFjoints;
ID&
SIFMember::getConnectedJoints(void)
{
  return *ConnectedJoints;
}


int 
SIFMember::setIntNodeTags(const ID& nodeTag)
{
	if(nodeTag.Size()==0)
	{
		opserr<<"WARNNING::SIFMember "<<this->getTag()<<" failed to set NodeTag"<<endln;
		return -1;
	}
   IntNodeTags = nodeTag;
   return 0;
}

const ID&
SIFMember::getIntNodeTags()
{
	return IntNodeTags;

}

int 
SIFMember::setIntEleTags(const ID& eleTag)
{
	if(eleTag.Size()==0)
	{
		opserr<<"WARNNING::SIFMember "<<this->getTag()<<" failed to set EleTag"<<endln;
		return -1;
	}
   IntEleTags = eleTag;
   return 0;
}

const ID&
SIFMember::getIntEleTags()
{
	return IntEleTags;
}


int 
SIFMember::AddFireAction(int theFireActionID){

	if(theFireActionID==0){
		opserr<<"WARNING::SIFCompartment "<<this->getTag()<<" received an empty definition of SIFfireAction"<<endln;
		return -1;
	}
	if(AppliedFireActions==0){
		AppliedFireActions= ID(1);
		AppliedFireActions(0)= theFireActionID;
	}
	else{
		int NumFA= AppliedFireActions.Size();
		//check whether the fire action is already added;
		for(int i=0; i<NumFA; i++){
			if(AppliedFireActions(i)== theFireActionID){
				opserr<<"WARNING::SIFMember cannont add an existing fire action: "
						<<theFireActionID<<endln;
				return 0;
			}
		}
		//if it is a new one
		AppliedFireActions.resize(NumFA+1);
		AppliedFireActions(NumFA)= theFireActionID;
	}
	return 0;
}

int 
SIFMember::WipeFireAction(){
	if(AppliedFireActions!=0)
		AppliedFireActions=0;
	return 0;
}

const ID&
SIFMember::getAppliedFireActions(){
	return AppliedFireActions;
}

const ID& 
SIFMember::getMemberInfo(void){

return theMemberInfo;

}

double 
SIFMember::getPartialDamage(Vector& damageVec){
	double PartialDamage;
	//PartialDamage = UpperLoc - LowerLoc;
	if(DamVec!=0){
		PartialDamage = DamVec(1)- DamVec(0);
		damageVec  = DamVec;
		return PartialDamage;
	}
	else 
		return 0;
}

int  
SIFMember::setPartialDamage(const Vector& damageVec){
	DamVec = Vector();
	DamVec = damageVec;

	return 0;
}