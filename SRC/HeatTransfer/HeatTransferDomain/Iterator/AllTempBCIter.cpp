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
// Adapted by Yaqiang Jiang (y.jiang@ed.ac.uk)
//


#include <AllTempBCIter.h>

#include <HeatTransferDomain.h>
#include <BoundaryPattern.h>
#include <BoundaryPatternIter.h>
#include <TemperatureBC.h>
#include <TaggedObjectIter.h>
#include <TaggedObjectStorage.h>

// SingleDomAllSP_Iter(SingleDomAllain &theDomain):
//	constructor that takes the model, just the basic iter

AllTempBCIter::AllTempBCIter(HeatTransferDomain& domain)
:theDomain(&domain), doneDomainTemps(false)
{

}


AllTempBCIter::~AllTempBCIter()
{

}    


void
AllTempBCIter::reset(void)
{
    theDomainTemps = &(theDomain->getTemperatureBCs());
    theBoundaryPatterns = &(theDomain->getBoundaryPatterns());
    currentBoundaryPattern = (*theBoundaryPatterns)();
    if (currentBoundaryPattern != 0) {
		theBCPatternTemps = &(currentBoundaryPattern->getTemperatureBCs());
		}

    doneDomainTemps = false;
}


TemperatureBC*
AllTempBCIter::operator()(void)
{
    TemperatureBC* theRes = 0;

	// return TemperatureBCs added in HTDomain
    if (doneDomainTemps == false) {
		theRes = (*theDomainTemps)();
		if (theRes != 0)
			return theRes;
		else
			doneDomainTemps = true;
		}
    
	//  return TemperatureBCs added in BoundaryPattern
	while (currentBoundaryPattern != 0) {
		theRes = (*theBCPatternTemps)();
		if (theRes == 0) {
			currentBoundaryPattern = (*theBoundaryPatterns)();
			if (currentBoundaryPattern != 0)
				theBCPatternTemps = &(currentBoundaryPattern->getTemperatureBCs());
			} else
				return theRes;
		}

  return 0;
}