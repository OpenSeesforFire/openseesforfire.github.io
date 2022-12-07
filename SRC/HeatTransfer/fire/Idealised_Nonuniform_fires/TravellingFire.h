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
// Added by Liming Jiang (liming.jiang@polyu.edu.hk)
// Based on LocalisedBurning fire model

#ifndef TravellingFire_h
#define TravellingFire_h
#include <FireModel.h>
#include <PathTimeSeriesThermal.h>
#include <Vector.h>

//class Vector;
class HeatTransferNode;

class TravellingFire : public FireModel
{
    public:
	  // typeTag indicates the type of norminal fire given in EN1991-1-2
	  // default is 1 corresponding to the standard temperature-time curve,
	  // 2 is for external fire curve, 3 is for hydrocarbon curve. Default 
	  // value is 1.
		TravellingFire(int tag, double D = 1, double Q = 1e6,
			double H = 3, int centerLineTag = 2, double smokeTemp = 293.15, PathTimeSeriesThermal* fireLocPath = 0);

	  //TravellingFire(double crd1, double crd2, double crd3, const Vector& time,
		 //              const Vector& d, const Vector& Q, double H);


	  virtual ~TravellingFire();
	  
	  void applyFluxBC(HeatFluxBC* theFlux, double time);
	  int setFirePars(double time,const Vector& firePars =0);
	  double getFirePars(int ParTag=1);
	  double getFireOut(double time, const Vector&);
	protected:

    private:
	  double getFlux(HeatTransferNode* the_node, double time);
	  PathTimeSeriesThermal* FireLocPath;
	  Vector fireLocs;
	  double  d, q, h;
	  double smokeT;
	  double maxq;
	  int centerLine;
};

#endif

