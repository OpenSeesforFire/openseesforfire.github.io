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

#ifndef SimpleMaterial_h
#define SimpleMaterial_h


#include <HeatTransferMaterial.h>

class SimpleMaterial: public HeatTransferMaterial
{
    public:
		SimpleMaterial(int tag, double rho, double cp, double kc);     
		~SimpleMaterial();

		// method for this material to update itself according to its new parameters
		int setTrialTemperature(double T, int par = 0);  // seem of no use

		const Matrix& getConductivity();
		double getRho();
		double getSpecificHeat();
		double getEnthalpy();
		double getEnthalpy(double temp);
		HeatTransferMaterial* getCopy();

		void update();
		int commitState();
		int revertToLastCommit();
		int revertToStart();


    protected:
    
    private:
		double rho, cp;
		double trial_temp;
		//double kc;
		double ini_temp;  // keep a copy of initial temperature
};


#endif
