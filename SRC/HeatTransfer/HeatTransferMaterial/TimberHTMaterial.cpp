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

TimberHTMaterial::TimberHTMaterial(int tag,int typeTag, HeatTransferDomain* theDomain, Vector matPars)
:HeatTransferMaterial(tag), trial_temp(0.0), charTime(0.0),
 ini_temp(0.0), rho(0), cp(0.0), enthalpy(0.0),TypeTag(typeTag), PhaseTag(0), theHTDomain(theDomain)
{
    if ( k == 0){
		k = new Matrix(3,3);
		if (k == 0) {
			opserr << "TimberHTMaterial::CarbonSteelEN93() - out of memory\n";
			exit(-1);
			}
		}

    pht1 = 0;
    pht2 = 0;
    //phaseTag =0: Wet Wood
    // phaseTag =1: Dry Wood
    // phaseTag =2: Char 
    //PhaseTag =3: Ash
    MatPars = matPars;
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
    

    double time = theHTDomain->getCurrentTime();

    this->determinePhase(trial_temp, time);

    return 0;
}


const Matrix& 
TimberHTMaterial::getConductivity(void)
{
	double materialK = 0;
	if(PhaseTag ==0){
      materialK = 0.186;
      //wet wood
	}
  else if(PhaseTag ==1){
        materialK = 0.176;
    //dry wood
  }
  else if(PhaseTag ==2){
        materialK = 0.065;
    //char
  }
  else if (PhaseTag == 3) {
        materialK = 0.058;
        //ash
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

  if (PhaseTag == 0) {
      rho = 380.0;
      //wet wood
  }
  else if (PhaseTag == 1) {
      rho = 360.0;
      //dry wood
  }
  else if (PhaseTag == 2) {
      rho = 73.0;
      //char
  }
  else if (PhaseTag == 3) {
      rho = 5.7;
      //ash
  }
  else
      opserr << "TimberHTMaterial::unrecognised PhaseTag " << PhaseTag;

  return rho;
}


double 
TimberHTMaterial::getSpecificHeat(void)
{
  
  if(PhaseTag ==0){
      cp = 1764;
  }
  else if(PhaseTag ==1){
      cp = 1664;
  
      /*  Temperature based definition
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
      */

  }
  else if (PhaseTag == 2) {
      cp = 1219;
  }
  else if (PhaseTag == 3) {
      cp = 1244;
  }
  else
      opserr << "TimberHTMaterial::unrecognised PhaseTag " << PhaseTag;



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
    TimberHTMaterial* theCopy = new TimberHTMaterial(this->getTag(), TypeTag,theHTDomain, MatPars);
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


int
TimberHTMaterial::determinePhase(double temp, double time)
{
    double T1, T2, T3;
    double dt1, dt2, dt3;
    T1 = MatPars(0); dt1 = MatPars(1);
    T2 = MatPars(2); dt2 = MatPars(3);
    T3 = MatPars(4); dt3 = MatPars(5);

    if (temp < T1) {
        PhaseTag = 0;
        pht1 = 0;
        //Wet wood
    }
    else if (temp < T2)
    {
        if (pht1<1e-6 && PhaseTag<1) {
            pht1 = time;
        }
        else {
            pht2 = time;
        }

        if ((pht2 - pht1) > dt1) {
            PhaseTag = 1;
            pht1 = 0;
        }
        //dry wood
    }
    else if (temp < T3)
    {
        if (pht1 < 1e-6 && PhaseTag<2) {
            pht1 = time;
        }
        else {
            pht2 = time;
        }

        if ((pht2 - pht1) > dt2) {
            PhaseTag = 2;
            pht1 = 0;
            if (charTime < 1e-6)
                charTime = time;
        }

        //char
    }
    else {
        if (pht1 < 1e-6 && PhaseTag<3) {
            pht1 = time;
        }
        else {
            pht2 = time;
        }

        if ((pht2 - pht1) > dt3) {
            PhaseTag = 3;
            pht1 = 0;
        }
        //ash
    }

      
    return 0;
}


const Vector&
TimberHTMaterial::getPars() {
    static Vector pars(2);
    pars(0) = PhaseTag;
    pars(1) = charTime;

    return pars;

}