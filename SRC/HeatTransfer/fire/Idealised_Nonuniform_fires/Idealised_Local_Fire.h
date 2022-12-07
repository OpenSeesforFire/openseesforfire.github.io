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
// Added by Liming Jiang (y.jiang@ed.ac.uk)
//

#ifndef Idealised_Local_Fire_h
#define Idealised_Local_Fire_h
#include <FireModel.h>
#include <PathTimeSeriesThermal.h>

class Vector;
class HeatTransferNode;

class Idealised_Local_Fire : public FireModel
{
    public:
	  // typeTag indicates the type of norminal fire given in EN1991-1-2
	  // default is 1 corresponding to the standard temperature-time curve,
	  // 2 is for external fire curve, 3 is for hydrocarbon curve. Default 
	  // value is 1.
	  Idealised_Local_Fire(int tag, double crd1, double crd2, double crd3,
	  					double Q, double D1, double D2, int centerLineTag = 3);
	
	  Idealised_Local_Fire(int tag, double crd1, double crd2, double crd3, 
				  	  	double Q, double D1, double D2, double K1, double K2, int centerLineTag = 3);
  
	  Idealised_Local_Fire(int tag, double crd1, double crd2, double crd3,
                       double Q, double D1, double D2, double factor, int centerLineTag = 3);

	  Idealised_Local_Fire(int tag, PathTimeSeriesThermal* fireLocPath,
		  double Q, double D1, double D2, double factor, int centerLineTag = 3);
						

	  //Idealised_Local_Fire(double crd1, double crd2, double crd3, const Vector& time,
		 //              const Vector& d, const Vector& Q, double H);


	  virtual ~Idealised_Local_Fire();
	  
	  void applyFluxBC(HeatFluxBC* theFlux, double time);
	  double getFlux(HeatTransferNode* the_node, double time);
	protected:

    private:
	  PathTimeSeriesThermal* FireLocPath;
	  double x1, x2, x3, d1, d2,k1, k2, q, Factor;
	  //x1,x2,x2: Cordinates of fire origin
	  // d1,d2: length of constant zone, length of the whole zone
	  //k, reduction ratio of mid point for quadratic corelation
	  int centerLine;
  int CurveType; //1:Linear, 2:Quadratic, 3:exponential
};

#endif

