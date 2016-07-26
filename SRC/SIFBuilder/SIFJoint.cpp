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
// This file constructs the class SIFJoint which holds the information for joints.
// Created by Liming Jiang @UoE, Praveen Kamath @UoE

#include <SIFJoint.h>

SIFJoint::SIFJoint(int tag, double crd1, double crd2, double crd3):TaggedObject(tag),
         Crd(0),NodeTag(0), theLoad(0),BCtypeTag(0)
{
		Crd= new Vector(3);
		(*Crd)(0) = crd1;
		(*Crd)(1) = crd2;
		(*Crd)(2) = crd3;

}

SIFJoint::~SIFJoint()
{
	//
}

const Vector &
SIFJoint::getCrds() const
{
    // return the vector of SIFJoints coordinates
    return *Crd;
}

int 
SIFJoint::setNodeTag(const ID& nodeTag)
{
	if(nodeTag.Size()==0)
	{
		opserr<<"WARNNING::SIFJoint failed to set NodeTag"<<endln;
		return -1;
	}
   NodeTag = nodeTag;
   return 0;
}

const ID&
SIFJoint::getNodeTag()
{
	return NodeTag;

}

void 
SIFJoint::setBCtype(int BC_Type)
{

   BCtypeTag = BC_Type;

}

int
SIFJoint::getBCtype()
{
	return BCtypeTag;

}

void 
SIFJoint::addLoad(const Vector& theLoadVec)
{
	if(theLoad==0){
		theLoad = theLoadVec;
	}

}
	
const Vector&
SIFJoint::getLoad()
{

	return theLoad;

}