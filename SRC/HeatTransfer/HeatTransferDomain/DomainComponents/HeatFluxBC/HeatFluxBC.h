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

#ifndef HeatFluxBC_h
#define HeatFluxBC_h

#include <ID.h>
#include <HeatTransferDomainComponent.h>

class HeatTransferElement;

class HeatFluxBC : public HeatTransferDomainComponent
{
    public:
	  HeatFluxBC(int tag, int eTag, int fTag);
	  virtual ~HeatFluxBC();
	  
	  virtual void setDomain(HeatTransferDomain* theDomain);
	  virtual void applyFluxBC(double factor) = 0;
	  int getFaceTag() const;
	  int getElementTag(void) const;

	  virtual void setPatternTag(int ptag);
      int getPatternTag(void) const;
	  virtual int getTypeTag() = 0;

	protected:
      HeatTransferElement* theElement;  // pointer to associated element, need to be determined private or protected??
	
    private:
	  int elementTag, faceTag, patternTag;
};

#endif
