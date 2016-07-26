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

#include <HeatFluxBC.h>
#include <HeatTransferDomain.h>


HeatFluxBC::HeatFluxBC(int tag, int eTag, int fTag)
:HeatTransferDomainComponent(tag), elementTag(eTag),
 faceTag(fTag), patternTag(-1), theElement(0)
{
    // does nothing
}


HeatFluxBC::~HeatFluxBC()
{
    // does nothing
}


void
HeatFluxBC::setDomain(HeatTransferDomain* theDomain)
{
  this->HeatTransferDomainComponent::setDomain(theDomain);

  theElement = theDomain->getElement(elementTag);
  if (theElement == 0) {
	  opserr << "WARNING - HeatFluxBC::setDomain - no element with tag ";
	  opserr << elementTag << " exists in the domain\n";
	  }
}


int
HeatFluxBC::getFaceTag() const
{
	return faceTag;
}


void
HeatFluxBC::setPatternTag(int ptag)
{
    patternTag = ptag;
}


int
HeatFluxBC::getPatternTag(void) const
{
    return patternTag;
}


int
HeatFluxBC::getElementTag() const
{
	return elementTag;
}

