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

#include <SimpleMaterial.h>
#include <Matrix.h>

#include <OPS_Globals.h>

SimpleMaterial::SimpleMaterial(int tag, double my_rho, double my_cp, double kc)
:HeatTransferMaterial(tag), trial_temp(0.0), 
 ini_temp(0.0), rho(my_rho), cp(my_cp)
{
    if ( k == 0){
		k = new Matrix(3,3);
		(*k)(0,0) = kc;
		(*k)(1,1) = kc;
		(*k)(2,2) = kc;
		}
}


SimpleMaterial::~SimpleMaterial()
{
    if (k != 0)
		delete k;
}

int 
SimpleMaterial::setTrialTemperature(double temp, int par)
{
    trial_temp = temp;
    return 0;
}


const Matrix& 
SimpleMaterial::getConductivity(void)
{
   
  //  if (trial_temp < 800.0) {
		//double thek = 1.0 - trial_temp / 1600.0;
		//(*k)(0,0) = thek;
		//(*k)(1,1) = thek;
		//(*k)(2,2) = thek;
		//} else {
		//	(*k)(0,0) = 0.5;
		//	(*k)(1,1) = 0.5;
		//	(*k)(2,2) = 0.5;
		//}

	return *k;
}


double  
SimpleMaterial::getRho(void)
{
    return rho;
}


double 
SimpleMaterial::getSpecificHeat(void)
{
    return cp;
}


double
SimpleMaterial::getEnthalpy(void)
{
    opserr << "SimpleMaterial::getEnthalpy() not implemented yet\n";
    exit(-1);
}

double 
SimpleMaterial::getEnthalpy(double temperature)
{
    opserr << "SimpleMaterial::getEnthalpy(double ) not implemented yet\n";
    exit(-1);
}

HeatTransferMaterial*
SimpleMaterial::getCopy(void)
{
    SimpleMaterial* theCopy = new SimpleMaterial(this->getTag(), rho, cp, (*k)(0,0));
    theCopy->trial_temp = trial_temp;
    return theCopy;
}


void
SimpleMaterial::update()
{
    return;
}


int 
SimpleMaterial::commitState(void)
{
    return 0;
}


int 
SimpleMaterial::revertToLastCommit(void)
{
    return 0;
}


int 
SimpleMaterial::revertToStart(void)
{
    return 0;
}


