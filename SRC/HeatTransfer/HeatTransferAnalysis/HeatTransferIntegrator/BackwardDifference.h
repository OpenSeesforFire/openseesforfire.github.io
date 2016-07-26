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

#ifndef BackwardDifference_h
#define BackwardDifference_h


#include <HT_TransientIntegrator.h>
#include <Penalty_FE.h>

class HT_FE_Element;
class Vector;

class BackwardDifference: public HT_TransientIntegrator
{
  public:
	  BackwardDifference(bool vflag = true);
	  ~BackwardDifference();

	  // methods which define what the FE_Element and DOF_Groups add
	  // to the system of equation object.
	  int formEleTangent(HT_FE_Element* theEle);

	  int domainChanged(void);    
	  int newStep(double delta_t);    
	  int revertToLastStep(void);        
	  int update(const Vector& deltaT);
	  bool getUpdatingFlag(){return v_form;};
	  double getAlphaDeltat(){return alpha;};
   
  protected:
    
  private:
	  Vector *Tt, *Ttdot;  // response quantities at time t
	  Vector *T, *Tdot;   // response quantities at time t+deltaT
	  double alpha;
	  bool v_form; // use v-form updating algorithm, otherwise use d-form
};

#endif
