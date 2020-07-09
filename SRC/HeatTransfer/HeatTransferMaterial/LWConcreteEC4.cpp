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

#include <LWConcreteEC4.h>
#include <Matrix.h>
#include <OPS_Globals.h>
#include <cmath>


LWConcreteEC4::LWConcreteEC4(int tag, double moisture)
:HeatTransferMaterial(tag), trial_temp(0.0), ini_temp(0.0), 
 rho(1800.0), cp(840), enthalpy(0.0), moist(moisture)
{
    if ( k == 0){
		k = new Matrix(3,3);
		if (k == 0) {
			opserr << "CarbonSteelEC3::CarbonSteelEN93() - out of memory\n";
			exit(-1);
			}
		}
}


LWConcreteEC4::~LWConcreteEC4()
{
    if (k != 0)
		delete k;
}

int 
LWConcreteEC4::setTrialTemperature(double temp, int par)
{
    trial_temp = temp - 273.15;
    return 0;
}


const Matrix& 
LWConcreteEC4::getConductivity(void)
{
    if (trial_temp <= 800.0) {
		double kc = 1.0 - trial_temp / 1600.0;
		(*k)(0,0) = kc;
		(*k)(1,1) = kc;
		(*k)(2,2) = kc;
		} else {
			(*k)(0,0) = 0.5;
			(*k)(1,1) = 0.5;
			(*k)(2,2) = 0.5;
		 }

    return *k;
}


double  
LWConcreteEC4::getRho(void)
{
    return rho;
}


double 
LWConcreteEC4::getSpecificHeat(void)
{
    return cp;
}


double
LWConcreteEC4::getEnthalpy()
{
    opserr << "LWConcreteEC4::getEnthalpy() - no phase change implemented for this light weight concrete"
		  <<  " .\n";
	exit(-1);
}


double
LWConcreteEC4::getEnthalpy(double temp)
{
    opserr << "LWConcreteEC4::getEnthalpy(double ) - no phase change implemented for this light weight concrete"
		  <<  " .\n";
	exit(-1);

}


HeatTransferMaterial*
LWConcreteEC4::getCopy(void)
{
    LWConcreteEC4* theCopy = new LWConcreteEC4(this->getTag(), moist);
    theCopy->trial_temp = trial_temp;
    return theCopy;
}


void
LWConcreteEC4::update()
{
    return; 
}


int 
LWConcreteEC4::commitState(void)
{
    return 0;
}


int 
LWConcreteEC4::revertToLastCommit(void)
{
    return 0;
}


int 
LWConcreteEC4::revertToStart(void)
{
    return 0;
}


