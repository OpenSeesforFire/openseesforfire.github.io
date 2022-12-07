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
// Added by Liming Jiang (liming.jiang@ed.ac.uk)
// Modified based on LocalisedFireEC1.cpp
//


#ifndef LocalizedFireSFPE_h
#define LocalizedFireSFPE_h
#include <FireModel.h>

class Vector;
class HeatTransferNode;

class LocalizedFireSFPE : public FireModel
{
    public:
	  // typeTag indicates the type of norminal fire given in EN1991-1-2
	  // default is 1 corresponding to the standard temperature-time curve,
	  // 2 is for external fire curve, 3 is for hydrocarbon curve. Default 
	  // value is 1.
	  LocalizedFireSFPE(int tag, double crd1, double crd2, double crd3, double D, double Q,
		               double HC, double HB, int centerLineTag = 3, double qsquare = 0);

	  //LocalizedFireSFPE(double crd1, double crd2, double crd3, const Vector& time,
		 //              const Vector& d, const Vector& Q, double H);


	  virtual ~LocalizedFireSFPE();
	  
	  void applyFluxBC(HeatFluxBC* theFlux, double time);
	  double getFlux(HeatTransferNode* the_node, double time,int FireType);
	protected:

    private:
	  double determineFireConvec(HeatTransferNode* the_node, double time, int FireType = 0);
	  double x1, x2, x3, d, ini_q, hc,hb;
	  int centerLine, FireType;
	  double Qsquare;
};

#endif

