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
// Note: This class was adapted from PenaltyConstraintHandler   
#include <PenaltyBC_Handler.h>
#include <stdlib.h>

#include <HT_AnalysisModel.h>
#include <HeatTransferDomain.h>
#include <HT_FE_Element.h>
#include <HT_DOF_Group.h>
#include <HeatTransferNode.h>
#include <HeatTransferElement.h>
#include <HT_NodeIter.h>
#include <HT_ElementIter.h>
#include <TempBCIter.h>
#include <TemperatureBCIter.h>
#include <MP_TemperatureBCIter.h>
#include <MP_TempBCIter.h>
#include <TemperatureBC.h>
#include <MP_TemperatureBC.h>
#include <HeatTransferIntegrator.h>
#include <ID.h>
#include <Penalty_FE.h>
#include <HT_PenaltyMP_FE.h>



PenaltyBC_Handler::PenaltyBC_Handler(double penalty_num, double penalty_MPnum)
:alpha(penalty_num),alphaMP(penalty_MPnum)
{

}

PenaltyBC_Handler::~PenaltyBC_Handler()
{

}

int
PenaltyBC_Handler::handle()
{
    // first check links exist to a Domain and an AnalysisModel object
    HeatTransferDomain* theDomain = this->getDomainPtr();
    HT_AnalysisModel* theModel = this->getAnalysisModelPtr();
    HeatTransferIntegrator* theIntegrator = this->getIntegratorPtr();    
    
    if ((theDomain == 0) || (theModel == 0) || (theIntegrator == 0)) {
		opserr << "WARNING PenaltyBC_Handler::handle() - ";
		opserr << " setLinks() has not been called\n";
		return -1;
		}
    
    // get number ofelements and nodes in the domain 
    // and init the theFEs and theDOFs arrays
    int numTemps = 0;
    TemperatureBCIter& theTempBCs = theDomain->getAllTempBCs();
    TemperatureBC* tempPtr;
    while ((tempPtr = theTempBCs()) != 0)
		numTemps++;
    
    // initialse the DOF_Groups and add them to the AnalysisModel.
    //    : must of course set the initial IDs
    HT_NodeIter& theNod = theDomain->getNodes();
    HeatTransferNode* nodPtr;  
    HT_DOF_Group* dofPtr;
	MP_TemperatureBC* mpPtr;   
    
    int numDofGrp = 0;
    int count3 = 0;
    int countDOF =0;
    while ((nodPtr = theNod()) != 0) {
		if ((dofPtr = new HT_DOF_Group(numDofGrp++, nodPtr)) == 0) {
			opserr << "WARNING PenaltyBC_Handler::handle() ";
			opserr << "- ran out of memory";
			opserr << " creating HT_DOF_Group " << numDofGrp << endln;	
			return -4;    		
			}

		// initially set all the ID value to -2
		const ID& id = dofPtr->getID();
		for (int j = 0; j < id.Size(); j++) {
			dofPtr->setID(j,-2);
			countDOF++;
			}
		nodPtr->setDOF_GroupPtr(dofPtr);
		theModel->addHT_DOF_Group(dofPtr);
		}

    theModel->setNumEqn(countDOF);

    // create the HT_FE_Element for the HeatTransferElement and add to the HT_AnalysisModel
    HT_ElementIter& theEle = theDomain->getElements();
    HeatTransferElement* elePtr;

    int numFeEle = 0;
    HT_FE_Element* fePtr;
    while ((elePtr = theEle()) != 0) {
		// just a regular element .. create an FE_Element for it & add to AnalysisModel
		if ((fePtr = new HT_FE_Element(numFeEle++, elePtr)) == 0) {
			opserr << "WARNING PenaltyBC_Handler::handle() - ran out of memory";
			opserr << " creating HT_FE_Element " << elePtr->getTag() << endln; 
			return -5;
			}

		theModel->addHT_FE_Element(fePtr);
		}
    
    // create the Penalty_FE for the TemperatureBC and 
    // add to the HT_AnalysisModel
    TemperatureBCIter& theTemps = theDomain->getAllTempBCs();
    while ((tempPtr = theTemps()) != 0) {
		if ((fePtr = new Penalty_FE(numFeEle, *theDomain, *tempPtr, alpha)) == 0) {
			opserr << "WARNING PenaltyBC_Handler::handle()";
			opserr << " - ran out of memory";
			opserr << " creating Penalty_FE" << endln; 
			return -5;
			}		
		theModel->addHT_FE_Element(fePtr);
		numFeEle++;
		}	

	// create the PenaltyMP_FE for the MP_TemperatureBCs and 
    // add to the AnalysisModel    

    MP_TemperatureBCIter& theMPs = theDomain->getMPTemperatureBCs();
    while ((mpPtr = theMPs()) != 0) {
	if ((fePtr = new HT_PenaltyMP_FE(numFeEle, *theDomain, *mpPtr, alphaMP)) == 0) {
	    opserr << "WARNING PenaltyConstraintHandler::handle()";
	    opserr << " - ran out of memory";
	    opserr << " creating PenaltyMP_FE " << endln; 
	    return -5;
	}		
	theModel->addHT_FE_Element(fePtr);
	numFeEle++;
    }	        

    return 0;
}


void 
PenaltyBC_Handler::clearAll(void)
{
    // for the nodes reset the DOF_Group pointers to 0
    HeatTransferDomain* theDomain = this->getDomainPtr();
    if (theDomain == 0)
		return;

    HT_NodeIter& theNod = theDomain->getNodes();
    HeatTransferNode* nodPtr;
    while ((nodPtr = theNod()) != 0)
		nodPtr->setDOF_GroupPtr(0);
}    
