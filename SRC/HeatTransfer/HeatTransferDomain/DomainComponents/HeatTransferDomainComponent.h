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

#ifndef HeatTransferDomainComponent_h
#define HeatTransferDomainComponent_h

#include <TaggedObject.h>
#include <MovableObject.h>
#include <Vector.h>

class HeatTransferDomain;


class HeatTransferDomainComponent: public TaggedObject
{
  public:
    virtual ~HeatTransferDomainComponent();

    virtual void setDomain(HeatTransferDomain* theDomain);
    virtual HeatTransferDomain* getDomain(void) const;
	virtual void  Print(OPS_Stream&, int = 0) {return;};

  protected:
    HeatTransferDomainComponent(int tag);
    
  private:    
    HeatTransferDomain* theDomain; // a pointer to the enclosing Domain object
};

#endif
