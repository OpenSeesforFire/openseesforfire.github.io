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
// This file constructs the SIFCompartment iterator.
// Created by Liming Jiang @UoE

#include "SIFCompartmentIter.h"

#include <SIFCompartment.h>
#include <TaggedObjectIter.h>
#include <TaggedObjectStorage.h>


SIFCompartmentIter::SIFCompartmentIter(TaggedObjectStorage *theStorage)
  :myIter(theStorage->getComponents())
{
}


SIFCompartmentIter::~SIFCompartmentIter()
{
}    

void
SIFCompartmentIter::reset(void)
{
    myIter.reset();
}    


SIFCompartment*
SIFCompartmentIter::operator()(void)
{
    // check if we still have elements in the model
    // if not return 0, indicating we are done
    TaggedObject *theComponent = myIter();
    if (theComponent == 0)
	return 0;
    else {
	SIFCompartment* result = (SIFCompartment *)theComponent;
	return result;
    }
}

