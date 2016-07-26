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
// Note: This class was adapted from EarthquakePattern

#include <FireImposedPattern.h>
#include <stdlib.h>
#include <ID.h>
#include <FireModel.h>
#include <HeatFluxBC.h>
//#include <TemperatureBC.h>
//#include <MapOfTaggedObjects.h>
#include <HeatFluxBCIter.h>
//#include <TempBCIter.h>
//#include <TemperatureBCIter.h>


FireImposedPattern::FireImposedPattern(int tag)
:BoundaryPattern(tag), the_firemodel(0)
{

}

    
// ~LoadPattern()
//	destructor

FireImposedPattern::~FireImposedPattern()
{
    
}


void
FireImposedPattern::setFireModel(FireModel* theFireModel)
{
    // invoke the destructor on the old TimeSeries
    if (the_firemodel != 0)
		delete the_firemodel;

    // set the pointer to the new series object
    the_firemodel = theFireModel;
}


void
FireImposedPattern::applyBCs(double time)
{
    HeatFluxBC* fluxbc;
    HeatFluxBCIter& fluxIter = this->getHeatFluxBCs();

	while ( (fluxbc = fluxIter()) != 0) {
		the_firemodel->applyFluxBC(fluxbc, time);
	}

}
