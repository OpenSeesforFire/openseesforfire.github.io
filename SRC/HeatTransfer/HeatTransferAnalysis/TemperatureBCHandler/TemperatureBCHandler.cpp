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
// Note: This class was adapted from ConstraintHandler
#include <TemperatureBCHandler.h>
#include <HeatTransferDomain.h>
#include <HT_AnalysisModel.h>
#include <HeatTransferIntegrator.h>
#include <HT_FE_EleIter.h>
#include <HT_FE_Element.h>

TemperatureBCHandler::TemperatureBCHandler()
:theDomainPtr(0),theAnalysisModelPtr(0),theIntegratorPtr(0)
{

}


TemperatureBCHandler::~TemperatureBCHandler()
{
    
}


int
TemperatureBCHandler::doneNumberingDOF(void)
{
    // iterate through the FE_Element getting them to set their IDs
    HT_FE_EleIter& theEle = theAnalysisModelPtr->getFEs();
    HT_FE_Element* elePtr;
    while ((elePtr = theEle()) != 0)
		elePtr->setID();
    return 0;
}


void 
TemperatureBCHandler::setLinks(HeatTransferDomain& theDomain, 
			          HT_AnalysisModel& theModel,
			          HeatTransferIntegrator& theIntegrator)
{
    theDomainPtr = &theDomain;
    theAnalysisModelPtr = &theModel;
    theIntegratorPtr = &theIntegrator;
}
	

int
TemperatureBCHandler::update(void)
{
    return 0;
}


int
TemperatureBCHandler::applyBCs(void)
{
    return 0;
}


HeatTransferDomain*
TemperatureBCHandler::getDomainPtr(void) const
{
    return theDomainPtr;
}

HT_AnalysisModel*
TemperatureBCHandler::getAnalysisModelPtr(void) const
{
    return theAnalysisModelPtr;
}

HeatTransferIntegrator*
TemperatureBCHandler::getIntegratorPtr(void) const
{
    return theIntegratorPtr;
}