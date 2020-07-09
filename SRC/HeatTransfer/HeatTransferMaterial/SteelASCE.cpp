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

#include <SteelASCE.h>
#include <Matrix.h>
#include <OPS_Globals.h>
#include <cmath>

double SteelASCE::epsilon = 1e-5;

SteelASCE::SteelASCE(int tag)
:HeatTransferMaterial(tag), trial_temp(0.0), 
 ini_temp(0.0), rho(7850.0), cp(0.0), enthalpy(0.0)
{
    if ( k == 0){
		k = new Matrix(3,3);
		if (k == 0) {
			opserr << "SteelASCE::SteelASCE() - out of memory\n";
			exit(-1);
			}
		}
}


SteelASCE::~SteelASCE()
{
    if (k != 0)
		delete k;
}

int 
SteelASCE::setTrialTemperature(double temp, int par)
{
    trial_temp = temp - 273.15;
    return 0;
}


const Matrix& 
SteelASCE::getConductivity(void)
{
    if (trial_temp < 20.0) {
		// if ((20.0 - trial_temp) < epsilon) {
		//if (int(trial_temp + 0.5) == 20) {
			//double kc = 54.0 - 3.33 * 1e-2 * trial_temp;
			(*k)(0,0) = 47.56;
			(*k)(1,1) = 47.56;
			(*k)(2,2) = 47.56;
			//} else {
			//opserr << "CarbonSteelEN93::getConductivity() - trial temperature " 
			//	<< trial_temp << " out of bounds\n";
			//exit(-1);}
		} else if ((20.0 <= trial_temp) && (trial_temp <= 900.0)) {
			double kc = -0.022 * trial_temp + 48.0;
			(*k)(0,0) = kc;
			(*k)(1,1) = kc;
			(*k)(2,2) = kc;
		}else if ((900.0 < trial_temp) && (trial_temp <= 1300.0)) {
			(*k)(0,0) = 28.2;
			(*k)(1,1) = 28.2;
			(*k)(2,2) = 28.2;
			} else {
				opserr << "SteelASCE::getConductivity() - trial temperature "
					<< trial_temp << " out of bounds\n";
				exit(-1);
		 }

    return *k; // change
}


double  
SteelASCE::getRho(void)
{
    return rho;
}


double 
SteelASCE::getSpecificHeat(void)
{

    if (trial_temp < 20.0){
			  cp = 430.57;
		} else if ((20.0 <= trial_temp) && (trial_temp <= 650.0)) {
			cp = (0.004 * trial_temp + 3.3) * 1e6 / rho;
		} else if ((600.0 < trial_temp) && (trial_temp <= 725.0)) {
			cp = (0.068 * trial_temp - 38.3) * 1e6 / rho;
		} else if ((725.0 < trial_temp) && (trial_temp <= 800.0)) {
			cp = (-0.086 * trial_temp + 73.35) * 1e6 / rho;
		} else if ((800.0 < trial_temp) && (trial_temp <= 1300.0)) {
			cp = 4.55 * 1e6 / rho;
		} else {
			opserr << "SteelASCE::getSpecificHeat() - trial temperature " 
				   << trial_temp << " out of bounds, cap values used\n";
			cp = 579.62;
			}

    return cp;
}


double
SteelASCE::getEnthalpy()
{
 
    if ((0.0 <= trial_temp) && (trial_temp <= 650.0)) {
		static double p1 = 0.002 * 1e6;
		static double p2 = 3.3 * 1e6;
		enthalpy = p1 * pow(trial_temp,2) + p2 * trial_temp;
		} else if ((650.0 < trial_temp) && (trial_temp <= 725.0)) {
			static double p3 = 0.034 * 1e6;
			static double p4 = -38.3 * 1e6;
			static double p5 = 1.352 * 1e10;
			enthalpy = p3 * pow(trial_temp,2) + p4 * trial_temp + p5;
		} else if ((725.0 < trial_temp) && (trial_temp <= 800.0)) {
			static double p6 = -0.043 * 1e6;
			static double p7 = 73.35 * 1e6;
			static double p8 = 2.6953125 * 1e10;

			enthalpy = p6 * pow(trial_temp,2) + p7 * trial_temp - p8;
			} else if ((800.0 < trial_temp) && (trial_temp <= 1300.0)) {
				static double p9 = 4.55 * 1e6;
				static double p10 = 566.875 * 1e6;

				enthalpy = p9 * trial_temp + p10;
			} else {
				opserr << "SteelASCE::getEnthalpy() - nodal temperature " 
					<< trial_temp << " out of bounds\n";
				exit(-1);
				}

    return enthalpy;	
	
}


double
SteelASCE::getEnthalpy(double temp)
{
    double enthp;
	double nod_temp = temp - 273.15;
    
    if ((0.0 <= nod_temp) && (nod_temp <= 650.0)) {
		static double p1 = 0.002 * 1e6;
		static double p2 = 3.3 * 1e6;
		enthp = p1 * pow(nod_temp,2) + p2 * nod_temp;
		} else if ((650.0 < nod_temp) && (nod_temp <= 725.0)) {
			static double p3 = 0.034 * 1e6;
			static double p4 = -38.3 * 1e6;
			static double p5 = 1.352 * 1e10;
			enthp = p3 * pow(nod_temp,2) + p4 * nod_temp + p5;
		} else if ((725.0 < nod_temp) && (nod_temp <= 800.0)) {
			static double p6 = -0.043 * 1e6;
			static double p7 = 73.35 * 1e6;
			static double p8 = 2.6953125 * 1e10;

			enthp = p6 * pow(nod_temp,2) + p7 * nod_temp - p8;
			} else if ((800.0 < nod_temp) && (nod_temp <= 1300.0)) {
				static double p9 = 4.55 * 1e6;
				static double p10 = 566.875 * 1e6;

				enthp = p9 * nod_temp + p10;
			} else {
				opserr << "SteelASCE::getEnthalpy(double temp) - nodal temperature " 
					<< nod_temp << " out of bounds\n";
				exit(-1);
				}

    return enthp;	
}


HeatTransferMaterial*
SteelASCE::getCopy(void)
{
    SteelASCE* theCopy = new SteelASCE(this->getTag());
    theCopy->trial_temp = trial_temp;
    return theCopy;
}


void
SteelASCE::update()
{
    return; 
}


int 
SteelASCE::commitState(void)
{
    return 0;
}


int 
SteelASCE::revertToLastCommit(void)
{
    return 0;
}


int 
SteelASCE::revertToStart(void)
{
    return 0;
}


