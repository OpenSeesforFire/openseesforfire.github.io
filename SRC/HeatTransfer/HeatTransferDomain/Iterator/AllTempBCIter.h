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

// This class is a modified version of SingleDomAllSP_Iter 
// for the heat transfer module
// Modified by: Yaqiang Jiang(y.jiang@ed.ac.uk)

// SingleDomAllSP_Iter
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A


#ifndef AllTempBCIter_h
#define AllTempBCIter_h

#include <TemperatureBCIter.h>

class TaggedObjectStorage;
class TaggedObjectIter;
class HeatTransferDomain;
class BoundaryPatternIter;
class BoundaryPattern;

class AllTempBCIter: public TemperatureBCIter
{
  public:
    AllTempBCIter(HeatTransferDomain& theDomain);
    virtual ~AllTempBCIter();
    
    virtual void reset(void);
    virtual TemperatureBC* operator()(void);
    
  private:
    HeatTransferDomain* theDomain;
    bool doneDomainTemps;
    TemperatureBCIter* theDomainTemps;
    BoundaryPatternIter* theBoundaryPatterns;
    BoundaryPattern* currentBoundaryPattern;
    TemperatureBCIter* theBCPatternTemps;
};

#endif
