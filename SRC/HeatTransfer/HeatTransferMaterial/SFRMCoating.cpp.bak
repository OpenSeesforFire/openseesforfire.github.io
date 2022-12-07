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

#include <SFRMCoating.h>
#include <Matrix.h>
#include <OPS_Globals.h>
#include <cmath>

double SFRMCoating::epsilon = 1e-5;

SFRMCoating::SFRMCoating(int tag,int typeTag)
:HeatTransferMaterial(tag), trial_temp(0.0), 
 ini_temp(0.0), rho(7850.0), cp(0.0), enthalpy(0.0),TypeTag(typeTag)
{
    if ( k == 0){
		k = new Matrix(3,3);
		if (k == 0) {
			opserr << "SFRMCoating::CarbonSteelEN93() - out of memory\n";
			exit(-1);
			}
		}
}


SFRMCoating::~SFRMCoating()
{
    if (k != 0)
		delete k;
}

int 
SFRMCoating::setTrialTemperature(double temp, int par)
{
    trial_temp = temp - 273.15;
    return 0;
}


const Matrix& 
SFRMCoating::getConductivity(void)
{
	double materialK = 0;
	if(TypeTag ==1){
		if (trial_temp <= 300.0) 
            materialK = 0.0778-0.000054*trial_temp;
        else if(trial_temp <= 700.0)
            materialK = -0.08+0.000472*trial_temp;
        else
            materialK = -0.08+0.000472*700;
	}
  else if(TypeTag ==2){
    if (trial_temp <= 200.0)
      materialK = 0.121-0.000319*trial_temp;
    else if(trial_temp <= 700.0)
      materialK = 0.0468+0.000050*trial_temp;
    else
      materialK = 0.0468+0.000050*700;
  }
  else if(TypeTag ==3){
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
    opserr<<"SFRMCoating::unrecognised TypeTag "<<TypeTag;
  

	(*k)(0,0) = materialK;
	(*k)(1,1) = materialK;
	(*k)(2,2) = materialK;

    return *k; // change
}


double  
SFRMCoating::getRho(void)
{
  if (TypeTag==1) {
    rho = 298;
  }
  else if (TypeTag==2) {
    rho = 423.2;
  }
  else if (TypeTag==2) {
    rho = 451.8;
  }
  return rho;
}


double 
SFRMCoating::getSpecificHeat(void)
{
  
  if(TypeTag ==1){
    if (trial_temp <0)
          cp = 3236;
	else if (trial_temp < 100.0)
      cp = 3236 + 4.16*trial_temp;
    else if(trial_temp < 200.0)
        //cp = 3236 + 4.16 * trial_temp;
      cp = 3009+ 6.43 * trial_temp;
    else if(trial_temp < 400.0)
      cp = 6815 - 12.60*trial_temp;
    else if(trial_temp < 700.0)
      cp = 1645 - 0.36*trial_temp;
    else if(trial_temp < 1200.0)
      cp = 1645 - 0.36*700;
	else 
		opserr<<"SFRM Coating ,invalid temperature"<<trial_temp;
  }
  else if(TypeTag ==2){
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
  else if(TypeTag ==3){
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
    opserr<<"SFRMCoating::unrecognised TypeTag "<<TypeTag;
  // cp =170;
    return cp;
}


double
SFRMCoating::getEnthalpy()
{
    	return 0;	
}


double
SFRMCoating::getEnthalpy(double temp)
{
    
    return 0;	
}


HeatTransferMaterial*
SFRMCoating::getCopy(void)
{
    SFRMCoating* theCopy = new SFRMCoating(this->getTag(),TypeTag);
    theCopy->trial_temp = trial_temp;
    return theCopy;
}


void
SFRMCoating::update()
{
    return; 
}


int 
SFRMCoating::commitState(void)
{
    return 0;
}


int 
SFRMCoating::revertToLastCommit(void)
{
    return 0;
}


int 
SFRMCoating::revertToStart(void)
{
    return 0;
}


