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

#include <HeatFluxBCIter.h>

#include <HeatFluxBC.h>
#include <TaggedObjectIter.h>
#include <TaggedObjectStorage.h>


HeatFluxBCIter::HeatFluxBCIter(TaggedObjectStorage* theStorage)
:myIter(theStorage->getComponents())
{

}


HeatFluxBCIter::~HeatFluxBCIter()
{

}    

void
HeatFluxBCIter::reset(void)
{
    myIter.reset();
}    


HeatFluxBC*
HeatFluxBCIter::operator()(void)
{
    TaggedObject* theComponent = myIter();
    if (theComponent == 0)
		return 0;
	else {
		HeatFluxBC* result = (HeatFluxBC *)theComponent;
		return result;
		}
}