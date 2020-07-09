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

#include <StainlessSteelEC.h>
#include <Matrix.h>
#include <OPS_Globals.h>
#include <cmath>

double StainlessSteelEC::epsilon = 1e-5;

StainlessSteelEC::StainlessSteelEC(int tag)
:HeatTransferMaterial(tag), trial_temp(0.0), 
 ini_temp(0.0), rho(7850.0), cp(0.0), enthalpy(0.0)
{
    if ( k == 0){
		k = new Matrix(3,3);
		if (k == 0) {
			opserr << "StainlessSteelEC::CarbonSteelEN93() - out of memory\n";
			exit(-1);
			}
		}
}


StainlessSteelEC::~StainlessSteelEC()
{
    if (k != 0)
		delete k;
}

int 
StainlessSteelEC::setTrialTemperature(double temp, int par)
{
    trial_temp = temp - 273.15;
    return 0;
}


const Matrix& 
StainlessSteelEC::getConductivity(void)
{
	double lamda = 14.6 + 1.27*0.01*trial_temp;
	
	(*k)(0,0) = lamda;
	(*k)(1,1) = lamda;
	(*k)(2,2) = lamda;

    return *k; // change
}


double  
StainlessSteelEC::getRho(void)
{
    return rho;
}


double 
StainlessSteelEC::getSpecificHeat(void)
{

   cp = 450 + 0.280*trial_temp - 2.91*1e-4*trial_temp*trial_temp + 1.34*(1e-7)*trial_temp*trial_temp*trial_temp;

    return cp;
}


double
StainlessSteelEC::getEnthalpy()
{
	opserr << "StainlessSteelEC::getEnthalpy() should not be called" << endln;
    return -1;	
}


double
StainlessSteelEC::getEnthalpy(double temp)
{
	opserr << "StainlessSteelEC::getEnthalpy() should not be called" << endln;
	return -1;
}


HeatTransferMaterial*
StainlessSteelEC::getCopy(void)
{
    StainlessSteelEC* theCopy = new StainlessSteelEC(this->getTag());
    theCopy->trial_temp = trial_temp;
    return theCopy;
}


void
StainlessSteelEC::update()
{
    return; 
}


int 
StainlessSteelEC::commitState(void)
{
    return 0;
}


int 
StainlessSteelEC::revertToLastCommit(void)
{
    return 0;
}


int 
StainlessSteelEC::revertToStart(void)
{
    return 0;
}


