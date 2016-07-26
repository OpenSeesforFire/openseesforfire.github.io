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
// Note: This class was adapted from LoadPattern


#include <BoundaryPattern.h>
#include <stdlib.h>
#include <ID.h>
#include <TimeSeries.h>
#include <HeatFluxBC.h>
#include <TemperatureBC.h>
#include <MapOfTaggedObjects.h>
#include <HeatFluxBCIter.h>
#include <TempBCIter.h>
#include <TemperatureBCIter.h>


BoundaryPattern::BoundaryPattern(int tag)
:HeatTransferDomainComponent(tag),
 isConstant(1), factor(0.), theSeries(0), 
 theHeatFluxBCs(0), theTemperatureBCs(0),
 theFluxIter(0), theTempIter(0)
{
    theHeatFluxBCs = new MapOfTaggedObjects();
    theTemperatureBCs = new MapOfTaggedObjects();

    if (theHeatFluxBCs == 0 || theTemperatureBCs == 0) {
		opserr << " BoundaryPattern::BoundaryPattern() - ran out of memory\n";
		exit(-1);
		}    

    theFluxIter = new HeatFluxBCIter(theHeatFluxBCs);    
    theTempIter = new TempBCIter(theTemperatureBCs);
    
    if (theFluxIter == 0 || theTempIter == 0) {
		opserr << " BoundaryPattern::BoundaryPattern() - ran out of memory\n";
		exit(-1);
		}
}

    
// ~LoadPattern()
//	destructor

BoundaryPattern::~BoundaryPattern()
{
    if (theSeries != 0)
		delete theSeries;

    if (theHeatFluxBCs != 0)
		delete theHeatFluxBCs;

    if (theTemperatureBCs != 0)
		delete theTemperatureBCs;

    if (theFluxIter != 0)
		delete theFluxIter;

    if (theTempIter != 0)
		delete theTempIter;
}


void
BoundaryPattern::setTimeSeries(TimeSeries *theTimeSeries)
{
    // invoke the destructor on the old TimeSeries
    if (theSeries != 0)
		delete theSeries;

    // set the pointer to the new series object
    theSeries = theTimeSeries;
}


void
BoundaryPattern::setDomain(HeatTransferDomain* theDomain)
{
    if (theHeatFluxBCs != 0) {
		HeatFluxBC* fluxbc;
		HeatFluxBCIter& the_fluxbc_iter = this->getHeatFluxBCs();
		while ((fluxbc = the_fluxbc_iter()) != 0)
			fluxbc->setDomain(theDomain);    

		TemperatureBC* tempbc;
		TemperatureBCIter& the_tempbc_iter = this->getTemperatureBCs();
		while ((tempbc = the_tempbc_iter()) != 0)
			tempbc->setDomain(theDomain);
		}

    this->HeatTransferDomainComponent::setDomain(theDomain);
}


bool
BoundaryPattern::addHeatFluxBC(HeatFluxBC* fluxbc)
{
    HeatTransferDomain* theDomain = this->getDomain();
    
    bool result = theHeatFluxBCs->addComponent(fluxbc);
    if (result == true) {
		if (theDomain != 0)
			fluxbc->setDomain(theDomain);
		fluxbc->setPatternTag(this->getTag());
		} else
			opserr << "WARNING: BoundaryPattern::addHeatFluxBC() - FluxBC could not be added\n";	
    
    return result;
}

bool
BoundaryPattern::addTemperatureBC(TemperatureBC* tempbc)
{
    HeatTransferDomain* theDomain = this->getDomain();
    
    bool result = theTemperatureBCs->addComponent(tempbc);
    if (result == true) {
		if (theDomain != 0)
			tempbc->setDomain(theDomain);
		tempbc->setPatternTag(this->getTag());
		} else
			opserr << "WARNING: BoundaryPattern::addTemperatureBC() - load could not be added\n";
    return result;
}

    
HeatFluxBCIter&
BoundaryPattern::getHeatFluxBCs(void)  
{
    theFluxIter->reset();
    return *theFluxIter;    
}

TemperatureBCIter&
BoundaryPattern::getTemperatureBCs(void)  
{
    theTempIter->reset();
    return *theTempIter;    
}

void
BoundaryPattern::clearAll(void)
{
    theHeatFluxBCs->clearAll();
    theTemperatureBCs->clearAll();
}


HeatFluxBC*
BoundaryPattern::removeHeatFluxBC(int tag)
{
    TaggedObject* obj = theHeatFluxBCs->removeComponent(tag);
    if (obj == 0)
		return 0;
    HeatFluxBC* result = (HeatFluxBC *)obj;
    result->setDomain(0);
    return result;    
}    

TemperatureBC*
BoundaryPattern::removeTemperatureBC(int tag)
{
    TaggedObject* obj = theTemperatureBCs->removeComponent(tag);
    if (obj == 0)
		return 0;
    TemperatureBC* result = (TemperatureBC *)obj;
    result->setDomain(0);
    return result;    
}    


void
BoundaryPattern::applyBCs(double time)
{
  // first determine the load factor
    if (theSeries != 0 && isConstant != 0)
		factor = theSeries->getFactor(time);
    
    HeatFluxBC* fluxbc;
    HeatFluxBCIter& fluxIter = this->getHeatFluxBCs();
    while ((fluxbc = fluxIter()) != 0)
		fluxbc->applyFluxBC(factor);

    TemperatureBC* tempbc;
    TemperatureBCIter& tempIter = this->getTemperatureBCs();
    while ((tempbc = tempIter()) != 0)
		tempbc->applyTemperatureBC(factor);
}


void
BoundaryPattern::setFluxConstant(void) 
{
    isConstant = 0;
}


void
BoundaryPattern::unsetFluxConstant(void) 
{
    isConstant = 1;
}


double
BoundaryPattern::getFactor(void)
{
  if (theSeries != 0)
	  return factor;
  else
	  return 0.0;
}


int
BoundaryPattern::getNumHeatFluxBCs(void)const
{   
 return theHeatFluxBCs->getNumComponents();
}
