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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:18 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/single/SingleDomMP_Iter.h,v $
                                                                        
                                                                        
// File: ~/domain/domain/single/SingleDomMP_Iter.h
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for 
// SingleDomMP_Iter. SingleDomMP_Iter is an iter for returning 
// the multiple point constraints  of an object of class SingleDomain. 
// SingleDomMP_Iters must be written for each subclass of SingleDomain.
// where the storage of the MP_Constraints changes.

#ifndef MP_TempBCIter_h
#define MP_TempBCIter_h

#include <MP_TemperatureBCIter.h>

class TaggedObjectStorage;
class TaggedObjectIter;

class MP_TempBCIter: public MP_TemperatureBCIter
{
  public:
    MP_TempBCIter(TaggedObjectStorage* theStorage);
    virtual ~MP_TempBCIter();
    
    virtual void reset(void);
    virtual MP_TemperatureBC* operator()(void);
    
  private:
    TaggedObjectIter &myIter;
};

#endif

