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
// Written by Liming Jiang (liming.jiang@polyu.edu.hk)
//

#ifndef TimberHTMaterial_h
#define TimberHTMaterial_h

#include <HeatTransferMaterial.h>
#include <HeatTransferDomain.h>
#include <Vector.h>

class TimberHTMaterial: public HeatTransferMaterial
{
    public:
		TimberHTMaterial(int tag,int typeTag, HeatTransferDomain* theDomain, Vector matPars);
		TimberHTMaterial(int tag, int typeTag, HeatTransferDomain* theDomain, Matrix thePars, Vector matPars);
		virtual ~TimberHTMaterial();

		// method for this material to update itself according to its new parameters
		int setTrialTemperature(double T, int par=0 );
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
		bool getIfHeatGen();
		double getHeatGen(double locy=0);

		const Vector& getPars();
    protected:
		
    private:
		
		int TypeTag;
		int trialphTag, PhaseTag;
		int TempTag;
		double rho, cp, enthalpy;
		double commit_rho, commit_cp, commit_k;
		double rho0, moist;
		double trial_temp;
		double ini_temp;  // keep a copy of initial temperature
		//static double charEndt;
		Vector MatPars;
		double pht1, pht2;
		double pht13, pht23;
		double HtComb, Qgen;
		double charTime;
		HeatTransferDomain* theHTDomain;
		const char* fileName;
		Matrix* thePars;
		double T1, T2, T3;
		double transt23;
		double commit_time;
		double current_time;
		double current_Qc;
		double commit_Qc;
		int determinePhase(double temp, double time);
};


#endif



