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
**   Yaqiang Jiang (y.jiang@ed.ac.uk)         
**   Liming JIang
**   Asif Usmani (asif.usmani@ed.ac.uk)                               **
**                                                                    **
** ****************************************************************** */

//
// Written by Liming Jiang (liming.jiang@polyu.edu.hk)
//

#ifndef StainlessSteelEC_h
#define StainlessSteelEC_h

#include <HeatTransferMaterial.h>

class StainlessSteelEC: public HeatTransferMaterial
{
    public:
		StainlessSteelEC(int tag);    
		virtual ~StainlessSteelEC();

		// method for this material to update itself according to its new parameters
		int setTrialTemperature(double T, int par = 0);
		const Matrix& getConductivity();
		double getRho();
		double getSpecificHeat();
		double getEnthalpy();
		// used for returning nodal enthalpy values
		double getEnthalpy(double NodalTemp);  
		HeatTransferMaterial* getCopy();

		int commitState();
		int revertToLastCommit();
		int revertToStart();
		void update();
		void  Print(OPS_Stream&, int = 0) {return;};

    protected:
    
    private:
		double rho, cp, enthalpy;
		double trial_temp;
		double ini_temp;  // keep a copy of initial temperature
		static double epsilon;
};


#endif



