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

// This class is a modified version of PenaltyConstraintHandler 
// for the heat transfer module
// Modified by: Yaqiang Jiang(y.jiang@ed.ac.uk)

// PenaltyConstraintHandler
// Written: fmk 
// Created: May 1998.
// Revision: A
#ifndef PenaltyBC_Handler_h
#define PenaltyBC_Handler_h

#include <TemperatureBCHandler.h>

class HT_FE_Element;
class HT_DOF_Group;

class PenaltyBC_Handler: public TemperatureBCHandler
{
  public:
    PenaltyBC_Handler(double penalty_num, double penalty_MPnum);//edited by Liming Jiang FOR MP_TemperatureBC
    ~PenaltyBC_Handler();

    int handle();
    void clearAll(void);    
    
  protected:
    
  private:
    double alpha;
	double alphaMP;
};

#endif
