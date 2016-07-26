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

#ifndef PrescribedSurfFlux_h
#define PrescribedSurfFlux_h

#include <HeatFluxBC.h>

class Vector;

class PrescribedSurfFlux : public HeatFluxBC
{
  public:
    PrescribedSurfFlux(int tag, int eleTag, int fTag, 
		               int numNPF, const Vector& data);

	PrescribedSurfFlux(int tag, int eleTag, int fTag, 
		               int numNPF,int fireType =0);
    ~PrescribedSurfFlux();

	virtual void applyFluxBC(double factor = 1.0);

	virtual const Vector& getData(void) const;
	virtual void setData(const Vector& theData);
	virtual int getFireType(void);
	int getTypeTag(){return type_tag;};

  protected:

  private:
	Vector* flux_data;
	int NPF;  // number of nodes per face
	int FireType;
	static const int type_tag = 3;
};

#endif