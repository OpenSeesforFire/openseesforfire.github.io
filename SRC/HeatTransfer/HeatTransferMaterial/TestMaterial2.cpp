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

#include <TestMaterial2.h>
#include <Matrix.h>
#include <OPS_Globals.h>
#include <cmath>


TestMaterial2::TestMaterial2(int tag)
:HeatTransferMaterial(tag), trial_temp(0.0), 
 ini_temp(0.0), rho(1.0), cp(0.0), enthalpy(0.0)
{
    if ( k == 0){
		k = new Matrix(3,3);
		if (k == 0) {
			opserr << "TestMaterial::CarbonSteelEN93() - out of memory\n";
			exit(-1);
			}
		}
}


TestMaterial2::~TestMaterial2()
{
    if (k != 0)
		delete k;
}

int 
TestMaterial2::setTrialTemperature(double temp)
{
    trial_temp = temp;
    return 0;
}


const Matrix& 
TestMaterial2::getConductivity(void)
{
    (*k)(0,0) = 53.0;
	(*k)(1,1) = 53.0;
	(*k)(2,2) = 53.0;

	return *k; // change
}


double  
TestMaterial2::getRho(void)
{

return 7850.0;
}


double 
TestMaterial2::getSpecificHeat(void)
{
if (trial_temp < 20.0){
		//if ((20.0 - trial_temp) < epsilon) {
		//if (int(trial_temp+0.5) == 20) {
			  cp = 439.80176;
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
			opserr << "CarbonSteelEN93::getSpecificHeat() - trial temperature " 
				   << trial_temp << " out of bounds, cap values used\n";
			cp = 3452443.816;
			}

    return cp;
}



double
TestMaterial2::getEnthalpy()
{
if (trial_temp < 20.0){
		//if (int(trial_temp+0.5) == 20) {
			enthalpy = 69048876.320000008;
			//} else {
				//opserr << "CarbonSteelEN93::getEnthalpy() - nodal temperature " 
					//<< trial_temp << " out of bounds\n";
				//exit(-1);}
		} else if ((20.0 <= trial_temp) && (trial_temp < 600.0)) {
			static double p1 = 3336250.0;
			static double p2 = 3034.0250000000001;
			static double p3 = 4.4221666666666666;
			static double p4 = 0.0043567499999999995;
			enthalpy = p1 * trial_temp + p2 * pow(trial_temp,2) - p3 * pow(trial_temp,3)
				+ p4 * pow(trial_temp,4) + 1144946.5733333379;
		} else if ((600.0 <= trial_temp) && (trial_temp < 735.0)) {
			static double p11 = 5228100.0;
			static double p12 = 102065700.0;
			enthalpy = p11 * trial_temp - p12 * log(738.0 - trial_temp) + 70634343.026482895;
			} else if ((735.0 <= trial_temp) && (trial_temp < 900.0)) {
				static double p13 = 4278250.0;
				static double p14 = 139887000.0;
				enthalpy = p13 * trial_temp + p14 * log(trial_temp - 731.0) + 462718901.46099216;
			} else if ((900.0 <= trial_temp) && (trial_temp <= 1200.0)) {
				static double p15 = 5102500.0000000000;
				enthalpy = p15 * trial_temp + 438500042.99543613;
				} else {
					opserr << "CarbonSteelEN93::getEnthalpy() - nodal temperature " 
						<< trial_temp << " out of bounds\n";
					exit(-1);
				}

    return enthalpy;	
}


double
TestMaterial2::getEnthalpy(double nod_temp)
{
    double enthp;
    if (nod_temp < 20.0){
		//if (int(nod_temp+0.5) == 20) {
			enthp = 69048876.320000008;
			//} else {
				//opserr << "CarbonSteelEN93::getEnthalpy(double ) - nodal temperature " 
				//	<< nod_temp << " out of bounds\n";
				//exit(-1);}
		} else if ((20.0 <= nod_temp) && (nod_temp < 600.0)) {
			static double p1 = 3336250.0;
			static double p2 = 3034.0250000000001;
			static double p3 = 4.4221666666666666;
			static double p4 = 0.0043567499999999995;
			enthp = p1 * nod_temp + p2 * pow(nod_temp,2) - p3 * pow(nod_temp,3)
				+ p4 * pow(nod_temp,4) + 1144946.5733333379;
		} else if ((600.0 <= nod_temp) && (nod_temp < 735.0)) {
			static double p11 = 5228100.0;
			static double p12 = 102065700.0;
			enthp = p11 * nod_temp - p12 * log(738.0 - nod_temp) + 70634343.026482895;
			} else if ((735.0 <= nod_temp) && (nod_temp < 900.0)) {
				static double p13 = 4278250.0;
				static double p14 = 139887000.0;
				enthp = p13 * nod_temp + p14 * log(nod_temp - 731.0) + 462718901.46099216;
			} else if ((900.0 <= nod_temp) && (nod_temp <= 1200.0)) {
				static double p15 = 5102500.0000000000;
				enthp = p15 * nod_temp + 438500042.99543613;
				} else {
					opserr << "CarbonSteelEN93::getEnthalpy(double ) - nodal temperature " 
						<< nod_temp << " out of bounds\n";
					exit(-1);
				}

    return enthp;	
}


HeatTransferMaterial*
TestMaterial2::getCopy(void)
{
    TestMaterial2* theCopy = new TestMaterial2(this->getTag());
    theCopy->trial_temp = trial_temp;
    return theCopy;
}


void
TestMaterial2::update()
{
    return; 
}


int 
TestMaterial2::commitState(void)
{
    return 0;
}


int 
TestMaterial2::revertToLastCommit(void)
{
    return 0;
}


int 
TestMaterial2::revertToStart(void)
{
    return 0;
}


