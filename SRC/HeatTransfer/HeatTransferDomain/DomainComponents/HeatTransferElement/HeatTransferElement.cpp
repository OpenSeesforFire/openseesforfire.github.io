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


#include <stdlib.h>
#include <algorithm>

#include <HeatTransferElement.h>
#include <HTElementResponse.h>
//#include <Renderer.h>
#include <Vector.h>
#include <Matrix.h>
#include <HeatTransferNode.h>
#include <HeatTransferDomain.h>
#include <Convection.h>
#include <Radiation.h>
#include <Information.h>

using std::list;

HeatTransferElement::HeatTransferElement(int tag) 
:HeatTransferDomainComponent(tag), theConvectionBCs(0), 
 theRadiationBCs(0), thePrescribedFluxBCs(0),
 convecFlag(false), radFlag(false),pFlag(false)
{
    // does nothing
}


HeatTransferElement::~HeatTransferElement() 
{
  //if (theConvectionBCs != 0)
	//  delete theConvectionBC;

  //if (theRadiationBCs != 0)
	//  delete theRadiationBC;
}


int
HeatTransferElement::commitState(void)
{
    return 0;
}


int
HeatTransferElement::update(void)
{
    return 0;
}


int
HeatTransferElement::revertToStart(void)
{
    return 0;
}


void 
HeatTransferElement::zeroFlux(void)
{
    // does nothing  
}


void 
HeatTransferElement::applyConvection(Convection* theConvection, double factor)
{
	if (theConvectionBCs.size() == 0) {
		convecFlag = true;
		theConvectionBCs.push_back(theConvection);
	} else {
		list<Convection*>::iterator result;
		result = find( theConvectionBCs.begin( ), theConvectionBCs.end( ), theConvection );
		if (result == theConvectionBCs.end())
			theConvectionBCs.push_back(theConvection);
	}	

}


void 
HeatTransferElement::applyRadiation(Radiation* theRadiation, double factor)
{
	if (theRadiationBCs.size() == 0) {
		radFlag = true;
		theRadiationBCs.push_back(theRadiation);
	} else {
		list<Radiation*>::iterator result;
		result = find( theRadiationBCs.begin( ), theRadiationBCs.end( ), theRadiation );
		if (result == theRadiationBCs.end())
			theRadiationBCs.push_back(theRadiation);
	}	
}


bool 
HeatTransferElement::hasConvection() const
{
    return convecFlag;		
}


bool 
HeatTransferElement::hasRadiation() const
{
    return radFlag;		
}

void 
HeatTransferElement::clearAllFluxBCs()
{
    if (theConvectionBCs.size() != 0) 
		theConvectionBCs.clear();

	if (theRadiationBCs.size() != 0) 
		theRadiationBCs.clear();

	if (thePrescribedFluxBCs.size() != 0) 
		thePrescribedFluxBCs.clear();
}



void
HeatTransferElement::removeConvection(Convection* theConvection)
{
	
	theConvectionBCs.remove(theConvection);
    // reset the flag if no convection bc exists
	if (theConvectionBCs.size() == 0)
		convecFlag = false;
}

void
HeatTransferElement::removeRadiation(Radiation* theRadiation)
{
    theRadiationBCs.remove(theRadiation);
    // reset the flag if no radiation bc exists
	if (theRadiationBCs.size() == 0)
		radFlag = false;
}

Response*
HeatTransferElement::setResponse(const char** argv, int argc, OPS_Stream& theHandler) {

	return 0;
}


int
HeatTransferElement::getResponse(int responseID, Information& eleInformation)
{
	return 0;
}