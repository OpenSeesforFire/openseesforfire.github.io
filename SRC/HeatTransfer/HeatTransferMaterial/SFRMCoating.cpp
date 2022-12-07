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
      materialK = -0.08+0.000469*trial_temp;
    else
      materialK = -0.08+0.000469*700;
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
  else if (TypeTag == 4) {
        if (trial_temp <= 750.0)
            materialK = 0.11 + 0.00028 * trial_temp;
        else
            materialK = 0.32;
    }
  else if (TypeTag == 5) {
        if (trial_temp <= 50)
            materialK = 0.094;
        else if (trial_temp <= 100)
            materialK = (6.24e-4) * trial_temp + 0.0628;
        else if (trial_temp <= 200)
            materialK = (-3.33e-4) * trial_temp + 0.1585;
        else if (trial_temp <= 300)
            materialK = (2.95e-4) * trial_temp + 0.0329;
        else if (trial_temp <= 800)
            materialK = (1.362e-4) * trial_temp + 0.0805;
        else if (trial_temp <= 1000)
            materialK = (3.615e-4) * trial_temp - 0.0997;
        else if (trial_temp <= 1200)
            materialK = 0.2618;
        else
            materialK = 0.32;
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
  else if (TypeTag==3) {
    rho = 451.8;
  }
  else if (TypeTag == 4) {
      rho = 357.0;
  }
  else if (TypeTag == 5) {
      if (trial_temp <= 25)
          rho = 297.48;
      else if (trial_temp <= 1200)
          rho = 4.504e-10 * pow(trial_temp, 4) - 9.087e-07 * pow(trial_temp, 3) + 0.0008126 * pow(trial_temp, 2) - 0.3926 * trial_temp + 306.8;
      else if (trial_temp > 1200)
          rho = 369.54;
  }
  return rho;
}


double 
SFRMCoating::getSpecificHeat(void)
{
  if(TypeTag ==1){
	if (trial_temp < 200.0)
      cp = 3236 + 5.295*trial_temp;
    else if(trial_temp < 400.0)
      cp = 7089 - 13.97*trial_temp;
    else if(trial_temp < 700.0)
      cp = 1645 - 0.36*trial_temp;
    else if(trial_temp < 1200.0)
      cp = 1645 - 0.36*700;
	else 
		opserr<<"sfrm coating ,invalid temperature"<<trial_temp;
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
  else if (TypeTag == 4) {
      if (trial_temp <= 1200.0)
          cp = 1111;
      else
          opserr << "SFRM Coating ,invalid temperature" << trial_temp;
  }
  else if (TypeTag == 5) {
      if (trial_temp <= 25)
          cp = 841;
      else if (trial_temp <= 100)
          cp = 1.3653 * trial_temp + 806.8667;
      else if (trial_temp <= 140)
          cp = 191.14 * trial_temp - 18171;
      else if (trial_temp <= 200)
          cp = -123.0583 * trial_temp + 25817;
      else if (trial_temp <= 600)
          cp = 0.4882 * trial_temp + 1107.8;
      else if (trial_temp <= 1000)
          cp = 0.3 * trial_temp + 1220.8;
      else
          cp = 1520.8;
  
  }
  else
    opserr<<"SFRMCoating::unrecognised TypeTag "<<TypeTag;
    return cp;
}


double
SFRMCoating::getEnthalpy()
{
    if (TypeTag == 5) {
        if (trial_temp <= 25)
            enthalpy = 841 * trial_temp;
        else if (trial_temp <= 100)
            enthalpy = 102.4 * trial_temp + 18465;
        else if (trial_temp <= 500)
            enthalpy = 140.55 * trial_temp + 2.1688e6;
        else if (trial_temp <= 600)
            enthalpy = 380.6 * trial_temp + 2.0488e6;
        else if (trial_temp <= 1000)
            enthalpy = 236.3 * trial_temp + 2.1353e6;
        else if (trial_temp <= 1200)
            enthalpy = 2371600;
        
        return enthalpy;
    }
    else { return 0; }
    //if (trial_temp < 200.0)
    //    enthalpy = 3236 * trial_temp + 5.295 * 0.5 * pow(trial_temp, 2);
    //else if (trial_temp < 400.0)
    //    enthalpy = -385300 + 7089 * trial_temp - 13.97 * 0.5 * pow(trial_temp, 2);
    //else if (trial_temp < 700.0)
    //    enthalpy = 1088800 + 1645 * trial_temp - 0.36 * 0.5 * pow(trial_temp, 2);
    //else if (trial_temp < 1200.0)
    //    enthalpy = 88200 + (1645 - 0.36 * 700) * trial_temp;
    //else
    //    opserr << "SFRM Coating ,invalid temperature" << trial_temp << endln;

    //return enthalpy;
}


double
SFRMCoating::getEnthalpy(double temp)
{
    double temp_C = temp - 273.15;
    if (TypeTag == 5) {
        if (temp_C <= 25)
            enthalpy = 841 * temp_C;
        else if (temp_C <= 100)
            enthalpy = 102.4 * temp_C + 18465;
        else if (temp_C <= 500)
            enthalpy = 140.55 * temp_C + 2.1688e6;
        else if (temp_C <= 600)
            enthalpy = 380.6 * temp_C + 2.0488e6;
        else if (temp_C <= 1000)
            enthalpy = 236.3 * temp_C + 2.1353e6;
        else if (temp_C <= 1200)
            enthalpy = 2371600;

        return enthalpy;
    }
    else { return 0; }
    //double enthp = 0;
    //double nod_temp = temp - 273.15;

    //if (nod_temp < 200.0)
    //    enthp = 3236 * nod_temp + 4.295 * 0.5 * pow(nod_temp, 2);
    //else if (nod_temp < 400.0)
    //    enthp = -385300 + 7089 * nod_temp - 13.97 * 0.5 * pow(nod_temp, 2);
    //else if (nod_temp < 700.0)
    //    enthp = 1088800 + 1645 * nod_temp - 0.36 * 0.5 * pow(nod_temp, 2);
    //else if (nod_temp < 1200.0)
    //    enthp = 88200 + (1645 - 0.36 * 700) * nod_temp;
    //else
    //    opserr << "SFRM Coating ,invalid temperature" << nod_temp << endln;

    //return enthp;
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


