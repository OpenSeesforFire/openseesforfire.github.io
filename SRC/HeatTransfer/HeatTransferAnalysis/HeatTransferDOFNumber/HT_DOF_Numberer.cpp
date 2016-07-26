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
// Note: This class was adapted from DOF_Numberer
#include <HT_DOF_Numberer.h>
#include <HT_AnalysisModel.h>
#include <GraphNumberer.h>
#include <ID.h>
#include <HT_DOF_Group.h>
#include <HT_FE_Element.h>
#include <HT_FE_EleIter.h>
//#include <Channel.h>
//#include <Channel.h>
//#include <FEM_ObjectBroker.h>
//#include <Channel.h>
//#include <FEM_ObjectBroker.h>

#include <Graph.h>

#include <HeatTransferDomain.h>
#include <MP_TemperatureBC.h>
#include <HeatTransferNode.h>
#include <MP_TempBCIter.h>
#include <HT_DOF_GrpIter.h>
// Constructor


HT_DOF_Numberer::HT_DOF_Numberer(GraphNumberer& aGraphNumberer)
:MovableObject(999), theModel(0), 
 theGraphNumberer(&aGraphNumberer)
{

}    


HT_DOF_Numberer::~HT_DOF_Numberer() 
{
    if (theGraphNumberer != 0)
		delete theGraphNumberer;
}


void
HT_DOF_Numberer::setLinks(HT_AnalysisModel& the_model)
{
    theModel = &the_model;
}


int 
HT_DOF_Numberer::numberDOF(int lastDOF_Group) 
{
    // check we have a model and a numberer
    HeatTransferDomain* theDomain = 0;
    if (theModel != 0) 
		theDomain = theModel->getDomainPtr();

    if ((theModel == 0) || (theDomain == 0)) {
		opserr << "WARNING HT_DOF_Numberer::numberDOF - ";
		opserr << "Pointers are not set\n";
		return -1;
		}
    
    if ((theGraphNumberer == 0)) {
		opserr << "WARNING DOF_Numberer::numberDOF - ";
		opserr << "subclasses must provide own implementation\n";
		return -2;
		}    

    // check we cant do quick return
    
    if (theModel->getNumDOF_Groups() == 0)
		return 0;

    // we first number the dofs using the dof group graph

    const ID& orderedRefs = theGraphNumberer->
		number(theModel->getDOFGroupGraph(), lastDOF_Group);     

    // we now iterate through the DOFs first time setting -2 values

    int eqnNumber = 0;
    
    if (orderedRefs.Size() != theModel->getNumDOF_Groups()) {
		opserr << "WARNING HT_DOF_Numberer::numberDOF - ";
		opserr << "Incompatable Sizes\n";
		return -3;
		}
    int result = 0;
    
    int size = orderedRefs.Size();
    for (int i = 0; i < size; i++) {
		int dofTag = orderedRefs(i);
		HT_DOF_Group* dofPtr;	
		dofPtr = theModel->getDOF_GroupPtr(dofTag);
		if (dofPtr == 0) {
			opserr << "WARNING HT_DOF_Numberer::numberDOF - ";
			opserr << "HT_DOF_Group " << dofTag << "not in HT_AnalysisModel!\n";
			result = -4;
			} else {
				const ID& theID = dofPtr->getID();
				int idSize = theID.Size();
				for (int j=0; j < idSize; j++)
					if (theID(j) == -2) dofPtr->setID(j,eqnNumber++);
			}
		}

    // iterate throgh  the DOFs second time setting -3 values
    for (int k = 0; k < size; k++) {
		int dofTag = orderedRefs(k);
		HT_DOF_Group* dofPtr;	
		dofPtr = theModel->getDOF_GroupPtr(dofTag);
		if (dofPtr != 0) {
			const ID& theID = dofPtr->getID();
			int idSize = theID.Size();
			for (int j = 0; j < idSize; j++)
				if (theID(j) == -3) dofPtr->setID(j,eqnNumber++);
			}
		}

    // iterate through the DOFs one last time setting any -4 values
    HT_DOF_GrpIter &tDOFs = theModel->getDOFs();
	HT_DOF_Group *dofPtr;
    while ((dofPtr = tDOFs()) != 0) {
    	const ID &theID = dofPtr->getID();
    	int have4s = 0;
	for (int i=0; i<theID.Size(); i++)
	    if (theID(i) == -4) have4s = 1;

	if (have4s == 1) {
		int nodeID = dofPtr->getNodeTag();
		// loop through the MP_Constraints to see if any of the
		// DOFs are constrained, note constraint matrix must be diagonal
		// with 1's on the diagonal
		MP_TemperatureBCIter &theMPs = theDomain->getMPTemperatureBCs();
		MP_TemperatureBC *mpPtr;
		while ((mpPtr = theMPs()) != 0 ) {
			// note keep looping over all in case multiple constraints
			// are used to constrain a node -- can't assume intelli user
	    		if (mpPtr->getNodeConstrained() == nodeID) {
	    			int nodeRetained = mpPtr->getNodeRetained();
	    			HeatTransferNode *nodeRetainedPtr = theDomain->getNode(nodeRetained);
	    			HT_DOF_Group *retainedDOF = nodeRetainedPtr->getDOF_GroupPtr();
	    			const ID&retainedDOFIDs = retainedDOF->getID();
	    			const ID&constrainedDOFs = mpPtr->getConstrainedDOFs();
	    			const ID&retainedDOFs = mpPtr->getRetainedDOFs();
	    			for (int i=0; i<constrainedDOFs.Size(); i++) {
	    				int dofC = constrainedDOFs(i);
	    				int dofR = retainedDOFs(i);
	    				int dofID = retainedDOFIDs(dofR);
	    				dofPtr->setID(dofC, dofID);
	    			}
	    		}
		}		
	}	
    }

    eqnNumber--;
    int numEqn = eqnNumber+1;






    // iterate through the FE_Element getting them to set their IDs
    HT_FE_EleIter& theEle = theModel->getFEs();
    HT_FE_Element* elePtr;
    while ((elePtr = theEle()) != 0)
		elePtr->setID();

    // set the numOfEquation in the Model
    theModel->setNumEqn(numEqn);

    if (result == 0)
		return numEqn;

    return result;
}


HT_AnalysisModel*
HT_DOF_Numberer::getAnalysisModelPtr(void) const
{
    return theModel;
}


GraphNumberer*
HT_DOF_Numberer::getGraphNumbererPtr(void) const
{
    return theGraphNumberer;
}
