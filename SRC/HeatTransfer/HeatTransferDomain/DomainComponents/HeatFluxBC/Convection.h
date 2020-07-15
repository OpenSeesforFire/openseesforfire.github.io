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

#ifndef Convection_h
#define Convection_h

#include <HeatFluxBC.h>
#include <Vector.h>


class Convection : public HeatFluxBC
{
  public:
    Convection(int tag, int eleTag, int fTag, double h, double T = 293.15, int fireType =0);
    ~Convection();

	void applyFluxBC(double factor);

	double getParameter(void) const;  // return convective heat transfer coefficient hc
	void setParameter(double);
	double getSurroundingTemp(void) const;
	void setSurroundingTemp(double T);
	int getTypeTag(){return type_tag;};
	int getFireType(void);

  protected:

  private:
	  int FireType;
	  double hc;
	  double Ta;
	  bool setIndex;
	  int type_tag = 1;
};

#endif