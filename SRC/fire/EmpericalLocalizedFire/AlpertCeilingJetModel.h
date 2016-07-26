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

#ifndef AlpertCeilingJetModel_h
#define AlpertCeilingJetModel_h
#include <FireModel.h>

class Vector;

class AlpertCeilingJetModel : public FireModel
{
    public:

	  //Note: Q in Alpert equation is in kW, instead of W in 
	  //Hasemi's model. crd1, crd2, crd3 are coordinates of the center
	  //Ta is ambient temperautre, Tf is cap tempeature in flamming region(1200 by default)
	  AlpertCeilingJetModel(int tag, double crd1, double crd2, double crd3, double Q,
		                    double H, double Ta, double Tf = 1200, 
							int centerLineTag = 3);

	  virtual ~AlpertCeilingJetModel();
	  
	  void applyFluxBC(HeatFluxBC* theFlux, double time);
	  double getGasTemperature(double xx1, double xx2, double xx3);
	  double getGasTemperature(HeatFluxBC* flux, double time);
	protected:

    private:
	  
	  double x1, x2, x3, q, h, T0, Tf; 
	  int centerLine;
};

#endif

