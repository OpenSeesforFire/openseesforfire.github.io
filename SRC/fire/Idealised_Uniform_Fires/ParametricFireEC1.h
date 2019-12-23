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

#ifndef ParametricFireEC1_h
#define ParametricFireEC1_h
#include <FireModel.h>


class ParametricFireEC1 : public FireModel
{
    public:
      ParametricFireEC1(int tag);
	  ParametricFireEC1(int tag, double I, double Av, double H, double At, double Af, double Qf, double tlim, double startTime = 0);

	  virtual ~ParametricFireEC1();
	  
	  void applyFluxBC(HeatFluxBC* theFlux, double time);
	  double getGasTemperature(double time);

	protected:

    private:

	  double ThI;  // thermal inertia of the compartment boundaries
	  double Avent;  // total  area of vertial openings on walls
	  double Hvent;  // weighted average of window heights on the walls
	  double Atotal;  // total area of the compartment(including walls)
	  double Afire;  // area of the floor with fire
	  double Qfire;  // total design fire related with Afire
	  double tlimit;  // time levels corresponds to different fire growth rate
	  double StartTime; //starttime of the fire action
};

#endif

