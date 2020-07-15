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

#ifndef Radiation_h
#define Radiation_h

#include <HeatFluxBC.h>
#include <Vector.h>


class Radiation : public HeatFluxBC
{
  public:
    Radiation(int tag, int eleTag, int fag, double epsilon, double sigma, double alpha, double qir);
    ~Radiation();

	virtual void applyFluxBC(double factor);

	virtual double getEmissionParameter(void) const;  // return product of emissivity and BZM_const
	virtual double getBLZMConstant(void) const;
	virtual double getIrradiation(void) const;
	virtual double getAbsorptivity(void) const;
	int getTypeTag(){return type_tag;};
	virtual void setIrradiation(double qir);

  protected:

  private:
	double emissivity, BZM_const, irradiation;
	double absorptivity;
	bool setIndex;
	static const int type_tag = 2;
};

#endif