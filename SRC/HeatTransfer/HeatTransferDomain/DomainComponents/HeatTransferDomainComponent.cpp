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
// Note: This class was adapted from DomainComponent 
#include <HeatTransferDomain.h>
#include <HeatTransferDomainComponent.h>

HeatTransferDomainComponent::HeatTransferDomainComponent(int tag)
  :TaggedObject(tag),theDomain(0)
{
    // does nothing else
}


HeatTransferDomainComponent::~HeatTransferDomainComponent()
{
    // does nothing
}


void
HeatTransferDomainComponent::setDomain(HeatTransferDomain* model)
{
    // sets the pointer 
    theDomain = model;
}


HeatTransferDomain *
HeatTransferDomainComponent::getDomain(void) const
{
    // returns the current pointer
    return theDomain;
}

