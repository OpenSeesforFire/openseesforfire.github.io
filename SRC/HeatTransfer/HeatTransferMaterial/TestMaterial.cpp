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

#include <TestMaterial.h>
#include <Matrix.h>
#include <OPS_Globals.h>
#include <cmath>


TestMaterial::TestMaterial(int tag)
:HeatTransferMaterial(tag), trial_temp(0.0), 
 ini_temp(0.0), rho(7850.0), cp(0.0), enthalpy(0.0)
{
    if ( k == 0){
		k = new Matrix(3,3);
		if (k == 0) {
			opserr << "TestMaterial::CarbonSteelEN93() - out of memory\n";
			exit(-1);
			}
		}
}


TestMaterial::~TestMaterial()
{
    if (k != 0)
		delete k;
}


int 
TestMaterial::setTrialTemperature(double temp)
{
    trial_temp = temp- 273.15;
    return 0;
}


const Matrix& 
TestMaterial::getConductivity(void)
{
    (*k)(0,0) = 1.08;
	(*k)(1,1) = 1.08;
	(*k)(2,2) = 1.08;
	return *k; // change
}


double  
TestMaterial::getRho(void)
{
    return 1.0;
}


double 
TestMaterial::getSpecificHeat(void)
{
    return 1.0;
}


double
TestMaterial::getEnthalpy()
{

  if ((-50.0 <= trial_temp) && (trial_temp < -5.1)) {
			enthalpy = trial_temp + 51.0;
		} else if ((-5.1 <= trial_temp) && (trial_temp < -0.1)) {
			//enthalpy = 8.021 * trial_temp + 121.9121;  // this is for 10degrees interval
			enthalpy = 15.042 * trial_temp + 122.6142; // this is for 5 degrees interval
			} else if ((-0.1 <= trial_temp) && (trial_temp < 15.0)) {
				enthalpy = trial_temp + 121.21;
			} else {
				opserr << "TestMaterial::getEnthalpy() - trial temperature " 
						<< trial_temp << " out of bounds\n";
				//	exit(-1);
				}

    return enthalpy;	
}


double
TestMaterial::getEnthalpy(double temp)
{
    double enthp;
	double nod_temp = temp - 273.15;

     if ((-50.0 <= nod_temp) && (nod_temp < -5.1)) {
			enthp = nod_temp + 51.0;
		} else if ((-5.1 <= nod_temp) && (nod_temp < -0.1)) {
			//enthp = 8.021 * nod_temp + 121.9121; // this is for 10degrees interval
			enthp = 15.042 * nod_temp + 122.6142; // this is for 5 degrees interval
			} else if ((-0.1 <= nod_temp) && (nod_temp < 15.0)) {
				enthp = nod_temp + 121.21;
			} else {
				opserr << "TestMaterial::getEnthalpy(double ) - nodal temperature " 
						<< nod_temp << " out of bounds\n";
					exit(-1);
				}

    return enthp;	
}


HeatTransferMaterial*
TestMaterial::getCopy(void)
{
    TestMaterial* theCopy = new TestMaterial(this->getTag());
    theCopy->trial_temp = trial_temp;
    return theCopy;
}


void
TestMaterial::update()
{
    return; 
}


int 
TestMaterial::commitState(void)
{
    return 0;
}


int 
TestMaterial::revertToLastCommit(void)
{
    return 0;
}


int 
TestMaterial::revertToStart(void)
{
    return 0;
}


