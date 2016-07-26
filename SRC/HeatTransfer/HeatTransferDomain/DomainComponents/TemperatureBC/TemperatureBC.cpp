/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Fire & Heat Transfer modules developed by:                         **
**   Yaqiang Jiang (y.jiang@ed.ac.uk)                                 **
**   Asif Usmani (asif.usmani@ed.ac.uk)                               **
**                                                                    **
** ****************************************************************** */

//
// Written by Yaqiang Jiang (y.jiang@ed.ac.uk)
//
// Note: This class was adapted from SP_Constraint 
#include <TemperatureBC.h>
//#include <classTags.h>
//#include <Vector.h>
#include <HeatTransferDomain.h>
#include <HeatTransferNode.h>

static int numTMPBCs = 0;
static int nextTag = 0;


// constructor for object of type SP_Constraint
TemperatureBC::TemperatureBC(int node, int ndof, double value, bool ISconstant)
:HeatTransferDomainComponent(nextTag++), nodeTag(node), 
 dofNumber(ndof), valueR(value), valueC(value),
 isConstant(ISconstant), patternTag(-1)
{
  numTMPBCs++;

}

TemperatureBC::~TemperatureBC()
{
  numTMPBCs--;
  if (numTMPBCs == 0)
    nextTag = 0;
}


int
TemperatureBC::getNodeTag(void) const
{
    // return id of constrained node
    return nodeTag;
}


int
TemperatureBC::getDOF_Number(void) const
{
    //  return the number of the constrained DOF    
    return dofNumber;
}


double
TemperatureBC::getValue(void)
{
    // return the value of the constraint
    return valueC;
}

int
TemperatureBC::applyTemperatureBC(double factor)
{
    if (isConstant == false)
		valueC = factor * valueR;

	//// overwrite the initial value
	//HeatTransferDomain* theDomain = this->getDomain();
	//theNode = theDomain->getNode(nodeTag);

	//if (theNode == 0) {
	//	opserr << "TemperatureBC::applyTemperatureBC() - node " << nodeTag << " does not exist\n";
	//	return -1;
	//	}

	//int numNodeDOF = theNode->getNumberDOF();

	//if (dofNumber < 0 || numNodeDOF <= dofNumber) {
	//	opserr << "TemperatureBC::applyTemperatureBC() - dof number " << dofNumber++ << " at node " 
	//		<< nodeTag << " not valid\n";
	//	return -2;
	//	}

	//nodalResponses = new Vector(numNodeDOF);
	//if (nodalResponses == 0) {
	//	opserr << "TemperatureBC::applyTemperatureBC() - out of memory\n";
	//	return -2;
	//	}

	//*nodalResponses = theNode->getTrialTemperature();
	//(*nodalResponses)(dofNumber) = valueC;

	//theNode->setTrialTemperature(*nodalResponses);

    return 0;
}


void
TemperatureBC::setPatternTag(int tag)
{
    patternTag = tag;
}


int
TemperatureBC::getPatternTag(void) const
{
    return patternTag;
}


void 
TemperatureBC::setNodalValue(void)
{
    HeatTransferDomain* theDomain = this->getDomain();
	if (theDomain == 0) {
		opserr << "TemperatureBC::setNodalValue - failed as the TemperatureBC with tag "
			   << nextTag << "is not associated with a HeatTransferDomain\n";
		exit(-1);
		}
	HeatTransferNode* theNode = theDomain->getNode(nodeTag);
	theNode->setResponse(valueC, dofNumber);
}

