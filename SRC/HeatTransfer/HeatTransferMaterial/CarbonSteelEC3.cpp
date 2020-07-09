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

#include <CarbonSteelEC3.h>
#include <Matrix.h>
#include <OPS_Globals.h>
#include <cmath>

double CarbonSteelEC3::epsilon = 1e-5;

CarbonSteelEC3::CarbonSteelEC3(int tag)
:HeatTransferMaterial(tag), trial_temp(0.0), 
 ini_temp(0.0), rho(7850.0), cp(0.0), enthalpy(0.0)
{
    if ( k == 0){
		k = new Matrix(3,3);
		if (k == 0) {
			opserr << "CarbonSteelEC3::CarbonSteelEN93() - out of memory\n";
			exit(-1);
			}
		}
}


CarbonSteelEC3::~CarbonSteelEC3()
{
    if (k != 0)
		delete k;
}

int 
CarbonSteelEC3::setTrialTemperature(double temp, int par)
{
    trial_temp = temp - 273.15;
    return 0;
}


const Matrix& 
CarbonSteelEC3::getConductivity(void)
{
    if (trial_temp < 20.0) {
		// if ((20.0 - trial_temp) < epsilon) {
		//if (int(trial_temp + 0.5) == 20) {
			//double kc = 54.0 - 3.33 * 1e-2 * trial_temp;
			(*k)(0,0) = 53.334000000000003;
			(*k)(1,1) = 53.334000000000003;
			(*k)(2,2) = 53.334000000000003;
			//} else {
			//opserr << "CarbonSteelEN93::getConductivity() - trial temperature " 
			//	<< trial_temp << " out of bounds\n";
			//exit(-1);}
		} else if ((20.0 <= trial_temp) && (trial_temp < 800.0)) {
			double kc = 54.0 - 3.33 * 1e-2 * trial_temp;
			(*k)(0,0) = kc;
			(*k)(1,1) = kc;
			(*k)(2,2) = kc;
		}else if ((800.0 <= trial_temp) && (trial_temp <= 1200.0)) {
			(*k)(0,0) = 27.3;
			(*k)(1,1) = 27.3;
			(*k)(2,2) = 27.3;
			} else {
				opserr << "CarbonSteelEC3::getConductivity() - trial temperature "
					<< trial_temp << " out of bounds\n";
				exit(-1);
		 }

    return *k; // change
}


double  
CarbonSteelEC3::getRho(void)
{
    return rho;
}


double 
CarbonSteelEC3::getSpecificHeat(void)
{
  //  if (((20.0 - trial_temp) <= epsilon) && (trial_temp < 600.0)) {
		//cp = 425.0 + 0.773 * trial_temp - 1.69 * 1e-3 * pow(trial_temp,2)
		//	 + 2.22 * 1e-6 * pow(trial_temp,3);
		//} else if (((600.0 - trial_temp) <= epsilon) && (trial_temp < 735.0)) {
		//	cp = 666.0 + 13002.0 / (738.0 - trial_temp);
		//} else if (((735.0 - trial_temp) <= epsilon) && (trial_temp < 900.0)) {
		//	cp = 545.0 + 17820.0 / (trial_temp - 731.0);
		//	} else if (((900.0 - trial_temp) <= epsilon) && (trial_temp <= 1200.0)) {
		//		cp = 650.0;
		//	} else {
		//		opserr << "CarbonSteelEN93::getSpecificHeat() - trial temperature" 
		//			   << trial_temp << " out of bounds\n";
		//		exit(-1);
		//		}

    if (trial_temp < 20.0){
		//if ((20.0 - trial_temp) < epsilon) {
		//if (int(trial_temp+0.5) == 20) {
			  cp = 439.7841332;
			//cp = 425.0 + 0.773 * trial_temp - 1.69 * 1e-3 * pow(trial_temp,2) 
			//	 + 2.22 * 1e-6 * pow(trial_temp,3);
			//} else {
			//	opserr << "CarbonSteelEN93::getSpecificHeat() - trial temperature " 
			//		   << trial_temp << " out of bounds\n";
			//	exit(-1);}
		} else if ((20.0 <= trial_temp) && (trial_temp < 600.0)) {
			cp = 425.0 + 0.773 * trial_temp - 1.69 * 1e-3 * pow(trial_temp,2) 
				 + 2.22 * 1e-6 * pow(trial_temp,3);
		} else if ((600.0 <= trial_temp) && (trial_temp < 735.0)) {
			cp = 666.0 + 13002.0 / (738.0 - trial_temp);
		} else if ((735.0 <= trial_temp) && (trial_temp < 900.0)) {
			cp = 545.0 + 17820.0 / (trial_temp - 731.0);
		} else if ((900.0 <= trial_temp) && (trial_temp <= 1200.0)) {
			cp = 650.0;
		} else {
			opserr << "CarbonSteelEC3::getSpecificHeat() - trial temperature " 
				   << trial_temp << " out of bounds, cap values used\n";
			cp = 650.0;
			}

    return cp;
}


double
CarbonSteelEC3::getEnthalpy()
{
    // The enthalpy data here is obtained by integrating the spesific heat over
    // the temperature range [20,1200] given in EN 1993-1-2
    if ((0.0 <= trial_temp) && (trial_temp < 600.0)) {
		static double p1 = 3336250.0;
		static double p2 = 3034.0250000000001;
		static double p3 = 4.4221666666666666;
		static double p4 = 0.0043567499999999995;
		enthalpy = p1 * trial_temp + p2 * pow(trial_temp,2) - p3 * pow(trial_temp,3)
			+ p4 * pow(trial_temp,4) - 67903929.746666670;
		} else if ((600.0 <= trial_temp) && (trial_temp < 735.0)) {
			static double p11 = 5228100.0;
			static double p12 = 102065700.0;
			enthalpy = p11 * trial_temp - p12 * log(738.0 - trial_temp) + 1585466.7064828873;
		} else if ((735.0 <= trial_temp) && (trial_temp < 900.0)) {
			static double p13 = 4278250.0;
			static double p14 = 139887000.0;
			enthalpy = p13 * trial_temp + p14 * log(trial_temp - 731.0) + 393670025.14099216;
			} else if ((900.0 <= trial_temp) && (trial_temp <= 1200.0)) {
				static double p15 = 5102500.0000000000;
				enthalpy = p15 * trial_temp + 369451166.67543614;
			} else {
				opserr << "CarbonSteelEC3::getEnthalpy() - nodal temperature " 
					<< trial_temp << " out of bounds\n";
				exit(-1);
				}

    return enthalpy;	
}


double
CarbonSteelEC3::getEnthalpy(double temp)
{
    double enthp;
	double nod_temp = temp - 273.15;
    
	// The temperature range is expanded to [0,1200] rather than the original one [20,1200] in Eurocode.
	// The reason is, for an analysis with initial temperature at 20, the solution could be lower than
	// 20 after initial iterations. Eventhough, the slope of H-T within the expanded temperature range is 
	// kept constant, the same as the heat capacity at T = 20;
	if ((0.0 <= nod_temp) && (nod_temp < 600.0)) {
		static double p1 = 3336250.0;
		static double p2 = 3034.0250000000001;
		static double p3 = 4.4221666666666666;
		static double p4 = 0.0043567499999999995;
		enthp = p1 * nod_temp + p2 * pow(nod_temp,2) - p3 * pow(nod_temp,3)
			+ p4 * pow(nod_temp,4) - 67903929.746666670;
		} else if ((600.0 <= nod_temp) && (nod_temp < 735.0)) {
			static double p11 = 5228100.0;
			static double p12 = 102065700.0;
			enthp = p11 * nod_temp - p12 * log(738.0 - nod_temp) + 1585466.7064828873;
		} else if ((735.0 <= nod_temp) && (nod_temp < 900.0)) {
			static double p13 = 4278250.0;
			static double p14 = 139887000.0;
			enthp = p13 * nod_temp + p14 * log(nod_temp - 731.0) + 393670025.14099216;
			} else if ((900.0 <= nod_temp) && (nod_temp <= 1300.0)) {
				static double p15 = 5102500.0000000000;
				enthp = p15 * nod_temp + 369451166.67543614;
			} else {
				opserr << "CarbonSteelEC3::getEnthalpy(double ) - nodal temperature " 
					<< nod_temp << " out of bounds\n";
				exit(-1);
				}

    return enthp;	
}


HeatTransferMaterial*
CarbonSteelEC3::getCopy(void)
{
    CarbonSteelEC3* theCopy = new CarbonSteelEC3(this->getTag());
    theCopy->trial_temp = trial_temp;
    return theCopy;
}


void
CarbonSteelEC3::update()
{
    return; 
}


int 
CarbonSteelEC3::commitState(void)
{
    return 0;
}


int 
CarbonSteelEC3::revertToLastCommit(void)
{
    return 0;
}


int 
CarbonSteelEC3::revertToStart(void)
{
    return 0;
}


