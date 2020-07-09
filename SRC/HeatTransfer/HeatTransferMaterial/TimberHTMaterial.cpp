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
// Written by Liming Jiang (liming.jiang@ed.ac.uk)
//

#include <TimberHTMaterial.h>
#include <Matrix.h>
#include <OPS_Globals.h>
#include <cmath>

double TimberHTMaterial::epsilon = 1e-5;

TimberHTMaterial::TimberHTMaterial(int tag,int typeTag)
:HeatTransferMaterial(tag), trial_temp(0.0), 
 ini_temp(0.0), rho(7850.0), cp(0.0), enthalpy(0.0),TypeTag(typeTag), PhaseTag(0)
{
    if ( k == 0){
		k = new Matrix(3,3);
		if (k == 0) {
			opserr << "TimberHTMaterial::CarbonSteelEN93() - out of memory\n";
			exit(-1);
			}
		}
}


TimberHTMaterial::~TimberHTMaterial()
{
    if (k != 0)
		delete k;
}

int 
TimberHTMaterial::setTrialTemperature(double temp, int par)
{
    trial_temp = temp - 273.15;
    return 0;
}


const Matrix& 
TimberHTMaterial::getConductivity(void)
{
	double materialK = 0;
	if(PhaseTag ==0){
      materialK = 0.126;
	}
  else if(PhaseTag ==1){
    if (trial_temp <= 200.0)
      materialK = 0.121-0.000319*trial_temp;
    else if(trial_temp <= 700.0)
      materialK = 0.0468+0.000050*trial_temp;
    else
      materialK = 0.0468+0.000050*700;
  }
  else if(PhaseTag ==2){
    if (trial_temp <= 200.0)
      materialK = 0.207-0.000318*trial_temp;
    else if(trial_temp <= 400.0)
      materialK = 0.147-0.000035*trial_temp;
    else if(trial_temp <= 700.0)
      materialK = 0.0054+0.000321*trial_temp;
    else
      materialK = 0.0054+0.000321*700;
  }
  else
    opserr<<"TimberHTMaterial::unrecognised PhaseTag "<<PhaseTag;
  

	(*k)(0,0) = materialK;
	(*k)(1,1) = materialK;
	(*k)(2,2) = materialK;

    return *k; // change
}


double  
TimberHTMaterial::getRho(void)
{
  if (PhaseTag==0) {
    rho = 430.0;
  }
  else if (PhaseTag==1) {
    rho = 423.2;
  }
  else if (PhaseTag==2) {
    rho = 451.8;
  }
  return rho;
}


double 
TimberHTMaterial::getSpecificHeat(void)
{
  
  if(PhaseTag ==0){
      cp = 2300;
  }
  else if(PhaseTag ==1){
	if (trial_temp < 100.0)
      cp = 1627+ 22.3*trial_temp;
    else if(trial_temp < 400.0)
      cp = 4446 - 5.05*trial_temp;
    else if(trial_temp < 700.0)
      cp = -1336 + 9.37*trial_temp;
    else if(trial_temp < 1200.0)
      cp = -1336 + 9.37*700;
	else 
		opserr<<"SFRM Coating ,invalid temperature"<<trial_temp;
  }
  else if(PhaseTag ==2){
	if (trial_temp < 200.0)
      cp = 643+1.93*trial_temp;
    else if(trial_temp < 400.0)
      cp = 1241-0.924*trial_temp;
    else if(trial_temp < 600.0)
      cp = 195+ 1.71*trial_temp;
    else if(trial_temp < 700.0)
      cp = 1826-1.08*trial_temp;
    else if(trial_temp < 1200.0)
      cp = 1826-1.08*700;
	else 
		opserr<<"SFRM Coating ,invalid temperature"<<trial_temp;
  }
  else
    opserr<<"TimberHTMaterial::unrecognised PhaseTag "<<PhaseTag;
  // cp =170;
    return cp;
}


double
TimberHTMaterial::getEnthalpy()
{
    	return 0;	
}


double
TimberHTMaterial::getEnthalpy(double temp)
{
    
    return 0;	
}


HeatTransferMaterial*
TimberHTMaterial::getCopy(void)
{
    TimberHTMaterial* theCopy = new TimberHTMaterial(this->getTag(),PhaseTag);
    theCopy->trial_temp = trial_temp;
    return theCopy;
}


void
TimberHTMaterial::update()
{
    return; 
}


int 
TimberHTMaterial::commitState(void)
{
    return 0;
}


int 
TimberHTMaterial::revertToLastCommit(void)
{
    return 0;
}


int 
TimberHTMaterial::revertToStart(void)
{
    return 0;
}


