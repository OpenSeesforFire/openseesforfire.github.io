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
#include <fstream>
#include <iomanip>

using std::ios;
using std::ifstream;

//double TimberHTMaterial::epsilon = 1e-5;

TimberHTMaterial::TimberHTMaterial(int tag,int typeTag, HeatTransferDomain* theDomain, Vector matPars)
:HeatTransferMaterial(tag), trial_temp(0.0), charTime(0.0), thePars(0), HtComb(0.0), trialphTag(0), TempTag(0),
 ini_temp(0.0), rho(0), cp(0.0), enthalpy(0.0),TypeTag(typeTag), PhaseTag(0), theHTDomain(theDomain),
    commit_rho(0), commit_cp(0), commit_k(0), Qgen(0), commit_time(0), current_time(0), current_Qc(0), commit_Qc(0)
{
    if ( k == 0){
		k = new Matrix(3,3);
		if (k == 0) {
			opserr << "TimberHTMaterial::CarbonSteelEN93() - out of memory\n";
			exit(-1);
			}
		}
    rho0 = 0;
    moist = 0;
    pht1 = 0;
    pht2 = 0;
    pht13 = 0;
    pht23 = 0;
    //phaseTag =0: Wet Wood
    // phaseTag =1: Dry Wood
    // phaseTag =2: Char 
    //PhaseTag =3: Ash
    MatPars = matPars;
    if (typeTag == 0) {
        rho0 = MatPars(0);
        moist = MatPars(1);
    }

    T1 = 0; T2 = 0; T3 = 0;

}

TimberHTMaterial::TimberHTMaterial(int tag, int typeTag, HeatTransferDomain* theDomain, Matrix thepars, Vector matPars)
    :HeatTransferMaterial(tag), trial_temp(0.0), charTime(0.0), HtComb(0.0), trialphTag(0), TempTag(0),
    ini_temp(0.0), rho(0), cp(0.0), enthalpy(0.0), TypeTag(typeTag), PhaseTag(0), theHTDomain(theDomain),
    commit_rho(0), commit_cp(0), commit_k(0), Qgen(0), commit_time(0), current_time(0), current_Qc(0), commit_Qc(0)
{
    if (k == 0) {
        k = new Matrix(3, 3);
        if (k == 0) {
            opserr << "TimberHTMaterial::CarbonSteelEN93() - out of memory\n";
            exit(-1);
        }
    }

    thePars = new Matrix(thepars.noRows(), thepars.noCols());
    (*thePars) = thepars;
    MatPars = matPars;
    pht1 = 0;
    pht2 = 0;
    pht13 = 0;
    pht23 = 0;
    rho0 = 0;
    moist = 0;
    //phaseTag =0: Wet Wood
    // phaseTag =1: Dry Wood
    // phaseTag =2: Char 
    //PhaseTag =3: Ash

    T1 = (*thePars)(1, 0);
    T2 = (*thePars)(2, 0);
    T3 = (*thePars)(3, 0);

    if (MatPars.Size() == 1) {
        moist = MatPars(0);
        HtComb = 0;
        transt23 = 3600;
    }
    else if (MatPars.Size() == 2) {
        moist = MatPars(0);
        HtComb = MatPars(1);
        transt23 = 3600;
    }
    else if (MatPars.Size() == 3) {
        moist = MatPars(0);
        HtComb = MatPars(1);
        transt23 = MatPars(2);
        // MatPars(3) may be needed
    }
    else
        opserr << "Timber Material recieves incorrect material properties" << endln;
    transt23 = 0.05 * 19e6 * (*thePars)(0, 1) / HtComb;
#ifdef _DEBUG
   // opserr << "transt23 " << transt23 << endln;
#endif
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

    if (TypeTag!=0) {
        this->determinePhase(trial_temp, time);
    }
    

    return 0;
}


const Matrix& 
TimberHTMaterial::getConductivity(void)
{
    double materialK = 0;
    if (TypeTag == 0) {
        if (trial_temp <= 20)
            materialK = 0.12;
        else if (trial_temp <= 200)
            materialK = 0.12 + 0.03 * (trial_temp - 20) / 180;
        else if (trial_temp <= 350)
            materialK = 0.15 - 0.08 * (trial_temp - 200) / 150;
        else if (trial_temp <= 500)
            materialK = 0.07 + 0.02 * (trial_temp - 350) / 150;
        else if (trial_temp <= 800)
            materialK = 0.09 + 0.26 * (trial_temp - 500) / 300;
        else if (trial_temp <= 1200)
            materialK = 0.35 + 1.15 * (trial_temp - 800) / 400;
        else
            opserr << "TimberHTMaterial::getSpecificHeat recieves incorrect temperature: " << trial_temp << endln;
    }
    else {
        
        if (trialphTag == 0) {
            materialK = (*thePars)(0, 2) + (0.11 * (*thePars)(0, 2)) * (trial_temp - 20) / 75;
            //wet wood
        }
        else if (trialphTag == 10)
        {
            if (trial_temp <= 95)
                //materialK = (*thePars)(0, 2);
                materialK = commit_k;
            else if (trial_temp <= 125)
                materialK = 1.11 * (*thePars)(0, 2) + (0.03 * (*thePars)(0, 2)) * (trial_temp - 95) / 30;
            else
                materialK = 1.14 * (*thePars)(0, 2);
        }
        else if (trialphTag == 1) {
            if (trial_temp <= 125)
                materialK = commit_k;
            else if (trial_temp <= 200)
                materialK = 1.14 * (*thePars)(0, 2) + (0.11 * (*thePars)(0, 2)) * (trial_temp - 125) / 75;
            else if (trial_temp <= 300)
                materialK = 1.25 * (*thePars)(0, 2) - (0.44 * (*thePars)(0, 2)) * (trial_temp - 200) / 100;
                //+ ((*thePars)(2, 2) - (*thePars)(1, 2)) * (trial_temp - 125) / 175;
      
            //dry wood
        }
        else if (trialphTag == 2) {
            if (trial_temp <= 300)
                materialK = commit_k;
            else if (trial_temp <= 600)
                materialK = 0.81 * (*thePars)(0, 2) + ((*thePars)(3, 2) - 0.81 * (*thePars)(0, 2)) * (trial_temp - 300) / 300;
            else if (trial_temp <= 800)
                materialK = (*thePars)(3, 2) + (0.39 * (*thePars)(3, 2)) * (trial_temp - 600) / 200;
            else if (trial_temp <= 1200)
                materialK = 1.39 * (*thePars)(3, 2) + (1.50 - 1.39 * (*thePars)(3, 2)) * (trial_temp - 800) / 400;     
            //char
        }
        else if (trialphTag == 3) {
            if (trial_temp <= 600)
                materialK = commit_k;
            else {
                if (trial_temp <= 800)
                    materialK = (*thePars)(3, 2) + (0.39 * (*thePars)(3, 2)) * (trial_temp - 600) / 200;
                else if (trial_temp <= 1200)
                    materialK = 1.39 * (*thePars)(3, 2) + (1.50 - 1.39 * (*thePars)(3, 2)) * (trial_temp - 800) / 400;

                if (materialK < commit_k)
                    materialK = commit_k;    //lower temperature than previously commited, using the commited k
            }

            //ash
        }
        else
            opserr << "TimberHTMaterial::unrecognised trialphTag " << trialphTag;
    }
	
    if (materialK < 0)
        opserr << "incorrect conductivity" << endln;

	(*k)(0,0) = materialK;
	(*k)(1,1) = materialK;
	(*k)(2,2) = materialK;

    return *k; // change
}


double  
TimberHTMaterial::getRho(void)
{
    if (TypeTag == 0) {
        if (trial_temp <= 100)
            rho = rho0*(1 + moist);
        else if (trial_temp <= 200)
            rho = rho0;
        else if (trial_temp <= 250)
            rho = rho0*(1.00 - 0.07 * (trial_temp - 200) / 50);
        else if (trial_temp <= 300)
            rho = rho0 * (0.93 - 0.17 * (trial_temp - 250) / 50);
        else if (trial_temp <= 350)
            rho = rho0 * (0.76 - 0.24 * (trial_temp - 300) / 50);
        else if (trial_temp <= 400)
            rho = rho0 * (0.52 - 0.14 * (trial_temp - 350) / 50);
        else if (trial_temp <= 600)
            rho = rho0 * (0.38 - 0.1 * (trial_temp - 400) / 200);
        else if (trial_temp <= 800)
            rho = rho0 * (0.28 - 0.02 * (trial_temp - 600) / 200);
        else if (trial_temp <= 1200)
            rho = rho0 * (0.26 - 0.26 * (trial_temp - 800) / 400);
        else
            opserr << "TimberHTMaterial::getSpecificHeat recieves incorrect temperature: " << trial_temp << endln;
    }
    else {
        if (trialphTag == 0) {
            rho = (*thePars)(0, 1) - (0.0736 * (*thePars)(0, 1)) * (trial_temp - 20) / 75;
            //wet wood
        }
        else if (trialphTag == 10) {
            if (trial_temp <= 95)
                rho = commit_rho;
            else if (trial_temp <= 125)
                rho = 0.9264 * (*thePars)(0, 1) - (0.0294 * (*thePars)(0, 1)) * (trial_temp - 95) / 30;
            else
                rho = (*thePars)(0, 1);
        }
        else if (trialphTag == 1) {
            // if (trial_temp <= 200)
              //   rho = (*thePars)(1, 1);
            // else
            if (trial_temp <= 125)
                rho = commit_rho;
            else if (trial_temp <= 300)
                rho = 0.897 * (*thePars)(0, 1) - (0.215 * (*thePars)(0, 1)) * (trial_temp - 125) / 175;
            // rho = (*thePars)(1, 1) + ((*thePars)(2, 1) - (*thePars)(1, 1)) * (trial_temp - 125) / 175;

             //dry wood
        }
        else if (trialphTag == 2) {
            if (trial_temp <= 300)
                rho = commit_rho;
                //rho = 641.91 - 138.97 * (trial_temp - 125) / 175;
            else if (trial_temp <= 600)
                rho = 0.682 * (*thePars)(0, 1) - (0.431 * (*thePars)(0, 1)) * (trial_temp - 300) / 300;
            else if (trial_temp <= 1200)
                rho = 0.251 * (*thePars)(0, 1) - (0.251 * (*thePars)(0, 1)) * (trial_temp - 600) / 600;

            //char
        }
        else if (trialphTag == 3) {
            if (trial_temp <= 600)
                rho = commit_rho;
            else {
                if (trial_temp <= 1200)
                    rho = 0.251 * (*thePars)(0, 1) - (0.251 * (*thePars)(0, 1)) * (trial_temp - 600) / 600;
                if(rho> commit_rho)
                    rho = commit_rho;    //lower temperature than previously commited, using the commited rho
            }

            //ash
        }
        else
            opserr << "TimberHTMaterial::unrecognised trialphTag " << trialphTag;
    }

    if (rho < 0)
        opserr << "incorrect density" << endln;

  return rho;
}


double 
TimberHTMaterial::getSpecificHeat(void)
{
    if (TypeTag == 0) {
        if (trial_temp <= 20)
            cp = 1530.0;
        else if (trial_temp <= 99)
            cp = 1530.0 + 240.0 * (trial_temp - 20) / 79;
        else if (trial_temp <= 109)
            cp = 1770.0 + 11830 * (trial_temp - 99) / 10;
        else if (trial_temp <= 119)
            cp = 13600.0 - 100.0 * (trial_temp - 109) / 10;
           // cp = 1770.0 + 350.0 * (trial_temp - 100) / 20;
        else if (trial_temp <= 129)
            cp = 13500.0 - 11380.0 * (trial_temp - 119) / 10;
        else if (trial_temp <= 200)
            cp = 2120.0 - 120.0 * (trial_temp - 121) / 79;
        else if (trial_temp <= 250)
            cp = 2000.0 - 380.0 * (trial_temp - 200) / 50;
        else if (trial_temp <= 300)
            cp = 1620.0 - 910.0 * (trial_temp - 250) / 50;
        else if (trial_temp <= 350)
            cp = 710.0 + 140.0 * (trial_temp - 300) / 50;
        else if (trial_temp <= 400)
            cp = 850.0 + 150.0 * (trial_temp - 350) / 50;
        else if (trial_temp <= 600)
            cp = 1000.0 + 400.0 * (trial_temp - 400) / 200;
        else if (trial_temp <= 800)
            cp = 1400.0 + 250.0 * (trial_temp - 600) / 200;
        else if (trial_temp <= 1200)
            cp = 1650.0;
        else
            opserr << "TimberHTMaterial::getSpecificHeat recieves incorrect temperature: " << trial_temp << endln;
    }
    else {
        if (trialphTag == 0) {
            cp = (*thePars)(0, 3);
        }
        else if (trialphTag == 10) {
            double maxcp = 14208e3 / (*thePars)(0, 1);
            if (trial_temp <= 95)
                cp = (*thePars)(0, 3);
             else if (trial_temp <= 105)
                 cp = (*thePars)(0, 3) + (maxcp - (*thePars)(0, 3)) * (trial_temp - 95.0) / 10.0;
            else if (trial_temp <= 115)
                 cp = maxcp;
            else if (trial_temp <= 125)
                cp = maxcp- (maxcp - (*thePars)(1, 3)) * (trial_temp - 115.0) / 10.0;
            else
                cp = (*thePars)(1, 3);
        }
        else if (trialphTag == 1) {
            //if (trial_temp <= 200)
               // cp = (*thePars)(1, 3);
           // else
            cp = (*thePars)(1, 3);
            //cp = (*thePars)(1, 3) + ((*thePars)(2, 3) - (*thePars)(1, 3)) * (trial_temp - 125) / 175;

            //dry wood
        }
        else if (trialphTag == 2) {
            if (trial_temp <= 300)
                cp = commit_cp;               //if temperature is lower than 300, the commited cp is used
            else if (trial_temp <= 600)
                cp = (*thePars)(1, 3) - ((*thePars)(1, 3) - (*thePars)(3, 3)) * (trial_temp - 300) / 300;
            else  
                cp = (*thePars)(3, 3);
            

            //char
        }
        else if (trialphTag == 3) {
                cp = (*thePars)(3, 3);   
            //ash
        }
        else
            opserr << "TimberHTMaterial::unrecognised trialphTag " << trialphTag;
    }

    if (cp < 0)
        opserr << "incorrect specific heat" << endln;

    return cp;
}


double
TimberHTMaterial::getEnthalpy()
{
    
    // The temperature range is expanded other than the original one [20,1200] in Eurocode.
    // The reason is, for an analysis with initial temperature at 20, the solution could be lower than
    // 20 after initial iterations. Eventhough, the slope of H-T within the expanded temperature range is 
    // kept constant, the same as the heat capacity at T = 20;
    //if ((0.0 <= nod_temp) && (nod_temp <= 100.0)) {
    /*
      double maxcp = 13600;
    if (trial_temp <= 100.0) {
        double c1 = (*thePars)(0, 3)* (*thePars)(0, 1);
        double c2 = (*thePars)(0, 3) * (*thePars)(0, 1)*25;
        enthalpy = c1 * trial_temp - c2;
    }
    else if ((100.0 < trial_temp) && (trial_temp <= 115.0)) {
        double c = (*thePars)(0, 3) * (*thePars)(0, 1) * (100.0 - 25.0);
        double c11 = maxcp* (*thePars)(0, 1);
        double c12 = c- maxcp * (*thePars)(0, 1)*100;
        enthalpy = c11 * trial_temp + c12 ;
    }
    else if ((115.0 < trial_temp) && (trial_temp <= 125.0)) {
        double c = maxcp * (*thePars)(0, 1)*115.0+ (*thePars)(0, 3) * (*thePars)(0, 1) * (100.0 - 25.0) - maxcp * (*thePars)(0, 1) * 100.0;
        
        double c11 = (*thePars)(1, 1)*((*thePars)(1, 3)-maxcp)/20.0;
        double c12 = (*thePars)(1, 1)* maxcp-115.0*((*thePars)(1, 3)-maxcp)/10.0;
        double c13 = c11 * 115.0 * 115.0 + c12 * 115.0;
        enthalpy = c11 * trial_temp * trial_temp + c12 * trial_temp-c13+c;
    }
    else if (trial_temp <= 1200.0) {
        double c = maxcp * (*thePars)(0, 1) * 115.0 + (*thePars)(0, 3) * (*thePars)(0, 1) * (100.0 - 25.0) - maxcp * (*thePars)(0, 1) * 100.0;

        double c11 = (*thePars)(1, 1) * ((*thePars)(1, 3) - maxcp) / 20.0;
        double c12 = (*thePars)(1, 1) * maxcp - 115.0 * ((*thePars)(1, 3) - maxcp) / 10.0;
        double c13 = c11 * 115.0 * 115.0 + c12 * 115.0;
        double c20 = c11 * 125.0 * 125.0 + c12 * 125.0 - c13 + c;

        double c21 = (*thePars)(1, 3) * (*thePars)(1, 1);
        double c22 = c20 - (*thePars)(1, 3) * (*thePars)(1, 1) * 125.0;
        enthalpy = c21 * trial_temp + c22;
    }


    else
    */
  
        enthalpy = 0;
    


    
    return enthalpy;	
}


double
TimberHTMaterial::getEnthalpy(double temp)
{
    double enthp;

    /*
    double nod_temp = temp - 273.15;

    // The temperature range is expanded other than the original one [20,1200] in Eurocode.
    // The reason is, for an analysis with initial temperature at 20, the solution could be lower than
    // 20 after initial iterations. Eventhough, the slope of H-T within the expanded temperature range is 
    // kept constant, the same as the heat capacity at T = 20;
    //if ((0.0 <= nod_temp) && (nod_temp <= 100.0)) {
    double maxcp = 13600;
    if (nod_temp <= 100.0) {
        double c1 = (*thePars)(0, 3) * (*thePars)(0, 1);
        double c2 = (*thePars)(0, 3) * (*thePars)(0, 1) * 25;
        enthp = c1 * nod_temp - c2;
    }
    else if (nod_temp <= 115.0) {
        double c = (*thePars)(0, 3) * (*thePars)(0, 1) * (100.0 - 25.0);
        double c11 = maxcp * (*thePars)(0, 1);
        double c12 = c - maxcp * (*thePars)(0, 1) * 100;
        enthp = c11 * nod_temp + c12;
    }
    else if (nod_temp <= 125.0) {
        double c = maxcp * (*thePars)(0, 1) * 115.0 + (*thePars)(0, 3) * (*thePars)(0, 1) * (100.0 - 25.0) - maxcp * (*thePars)(0, 1) * 100.0;

        double c11 = (*thePars)(1, 1) * ((*thePars)(1, 3) - maxcp) / 20.0;
        double c12 = (*thePars)(1, 1) * maxcp - 115.0 * ((*thePars)(1, 3) - maxcp)* (*thePars)(1, 1) / 10.0;
        double c13 = c11 * 115.0 * 115.0 + c12 * 115.0;
        enthp = c11 * nod_temp * nod_temp + c12 * nod_temp - c13 + c;
    }
   
    else if (nod_temp <= 200.0) {
        double c = maxcp * (*thePars)(0, 1) * 115.0 + (*thePars)(0, 3) * (*thePars)(0, 1) * (100.0 - 25.0) - maxcp * (*thePars)(0, 1) * 100.0;

        double c11 = (*thePars)(1, 1) * ((*thePars)(1, 3) - maxcp) / 20.0;
        double c12 = (*thePars)(1, 1) * maxcp - 115.0 * ((*thePars)(1, 3) - maxcp) * (*thePars)(1, 1) / 10.0;
        double c13 = c11 * 115.0 * 115.0 + c12 * 115.0;
        double c20 = c11 * 125.0 * 125.0 + c12 * 125.0 - c13 + c;

        double c21 = (*thePars)(1, 3) * (*thePars)(1, 1);
        double c22 = c20 - (*thePars)(1, 3) * (*thePars)(1, 1) * 125.0;
        enthp = c21 * nod_temp + c22;
    }
  
    
    else
    
    */
    
        enthp = 0;

    return enthp;
}


HeatTransferMaterial*
TimberHTMaterial::getCopy(void)
{
    TimberHTMaterial* theCopy = 0;
    if(thePars!=0) {
        theCopy = new TimberHTMaterial(this->getTag(), TypeTag, theHTDomain, (*thePars), MatPars);
    }
    else
         theCopy = new TimberHTMaterial(this->getTag(), TypeTag,theHTDomain, MatPars);
    
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

    PhaseTag = trialphTag;
    commit_rho = rho;
    commit_cp = cp;
    commit_k = (*k)(0, 0);
    commit_time = current_time;
    commit_Qc = current_Qc;
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
    
    PhaseTag = 0;
    return 0;
}


int
TimberHTMaterial::determinePhase(double temp, double time)
{
    double trial_Qc = commit_Qc;
    double alp1 = 1.0;
    double alp2 = 1.5;
    double alpha = 0;
    double factor = 1.0;

    trialphTag = PhaseTag;
    
    //determine phase
    if (temp <95) {
        //<100oC
        if(trialphTag<1)
            trialphTag = 0;
       // if (trialphTag ==10)
         //   trialphTag = 0;
        if (trialphTag < 2)
            pht1 = 0;
        else if (trialphTag == 2) {
            //char layer back to <300oC
            pht1 = time - charTime;
        }
        
        //Wet wood
    }
    else if (temp < T2)
    {
        //<300oC
        if (trialphTag == 0 || trialphTag == 10) {

            // if (pht1 < 1e-6 && trialphTag == 0) {
             //    pht1 = time;
            // }
            // else {
            //     pht2 = time;
            // }

           //  if ((pht2 - pht1) > dt1) {
           //      trialphTag = 1;
             //    pht1 = 0;
            // } 
            if (temp < 125)
                trialphTag = 10; // evaporation
            else
                trialphTag = 1;  //dry wood
        }
        else if(trialphTag == 1)
            trialphTag = 1;
        pht1 = 0;
        //dry wood
        if(trialphTag < 2)
             TempTag = 0;
        else if (trialphTag == 2) {
            //char layer back to <300oC, minus the duration
            pht1 = time - charTime;
        }
    }
    else if (temp <=1200)
    {
        //<600oC
        if (trialphTag < 2|| trialphTag ==10) {
            //just over 300oC
            if (pht1 < 1e-6 ) {
                if (trialphTag == 1 || trialphTag == 10) {
                    pht1 = time;
                    pht2 = time;
                }
                //record the initial char time 
            }

           
            trialphTag = 2;// enter char stage
            TempTag = 0;
           //if (charTime < 1e-6)
               // charTime = time;
           
        }
        else if (trialphTag == 2)
        {
            //already char layer state
            pht2 = time;
            TempTag = 0;
            current_Qc = commit_Qc;
            //char to ash transition 2->3
             charTime = pht2 - pht1;

             if (trial_temp < 300)
                 alpha = 0;
             //else if (trial_temp < 400)
               //  alpha = (trial_temp - 300) / 100 * 0.1;   //combustion area at 400-600 range
             else if (trial_temp < 400)
                 alpha = ((trial_temp - 300) / 100) * alp1;   //combustion area at 400-600 range  alpha = (trial_temp - 400) / 200 * 0.43;
             else if (trial_temp < 600)
                 alpha = alp1; //linear transition from alph1 to alp2 if different
             else if (trial_temp < 800)
                 alpha = alp1 + ((trial_temp - 600) / 200) * (alp2 - alp1);
             else
                 alpha = alp2;

             //determine factor from abosorbed heat to generated heat
             if (trial_temp < 600)
                 factor = 0.05;
             else
                 factor = 0.01;

             current_Qc = current_Qc + alpha * HtComb / factor*(0.25);
            //if T>600 and timber has long time of combustion
             if (current_Qc > 19e6* (*thePars)(0, 1) && trial_temp > T3) {
                 trialphTag = 3;
                 pht13 = time; //reset the pht1 for transition to 3;
                 pht23 = time;
                 TempTag = 1;
             }
                
        }
        else if (trialphTag == 3) {
            pht23 = time;
            if(pht23-pht13>10)
                TempTag = 0;  //set a 10s short transition time for smooth change of thermal properties
        }
        
    }
    else {
        trialphTag = 3;  //over 1200oC
    }
      
    /*
    else {
        if (trialphTag < 3) {
            if ((pht2-pht1>1e-6)&&TempTag==0) {
                pht1 = time;
                pht2 = time;
                TempTag = 1;
            }

            if (TempTag == 1)
                pht2 = time;
            
            double dt3 = 30;
            if ((pht2 - pht1 > dt3) && TempTag==1) {
                trialphTag = 3;
                //pht1 = 0;
            }
        }
        //ash
    }
    */
    

      
    return 0;
}


bool
TimberHTMaterial::getIfHeatGen()
{
    if (trialphTag == 2|| trialphTag == 3)
        return true;
    else
        return false;
}


double
TimberHTMaterial::getHeatGen(double locy)
{
    double alp1 = 1.0;    //400-600
    double alp2 = 1.5;   //600-800
    double alp3 = 1.0;  // initial, locy<0.005
    double alp4 = 1.0;  // phTag =3, flaming

    //double Qgen = 0;
    double alpha = 0;
	double locRatio = 1; 

    if (trialphTag == 2)
    {
        if (TempTag == 0) {
            if (trial_temp < 300)
                alpha = 0;
            else if (trial_temp < 400)
                alpha = (trial_temp - 300) / 100 * alp1;   //combustion area at 400-600 range
            else if (trial_temp < 600)
                alpha = alp1; //linear transition from alph1 to alp2 if different
            else if (trial_temp < 800)
                alpha = alp1 + ((trial_temp - 600) / 200) * (alp2 - alp1);
            else
                alpha = alp2;
        }

        
//        if (locy < 0.015) {
//               if (transt23 - charTime < 200) {
//                    alpha = ((transt23 - charTime) / 200) * alp3;
//                }
//                else
//                    alpha = alp3;
//        }
      
       
    }
	else if(trialphTag ==3){
        //flames heating the ash layer
        if (TempTag == 1) {
            if (pht23 > pht13) {
                alpha = alp1 + ((pht23 - pht13) / 10) * (alp2 - alp1);  //transition to ash layer:alp2
            }
        }
        else {
            if (trial_temp < 400)
                alpha = 0;
            else if (trial_temp < 600)
                alpha = (1 + (trial_temp - 600) / 200) * alp2;
            else if (trial_temp < 800)
                alpha = (1 - (trial_temp - 600) / 200) * alp2;
            else
                alpha = 0;
        }

       
	}
    //Considering locy indicating depth effect
       //locy=0, exposed surface ,ratio =1; locy =100, deep layer, ratio=0. You can change 100 to large value
      // if (locRatio < 0)
        //   locRatio = 0;
      // else if (locRatio>1)
        //   locRatio = 1;  

     //  locRatio = 1;   currently not used
	//	alpha = locRatio * alpha;
        
        Qgen = alpha* HtComb;

        if (Qgen < -1e-5) {
            opserr << "incorrect Heat of generation" << endln;
            Qgen = 0;
        }
        

    return Qgen ;
}



const Vector&
TimberHTMaterial::getPars() {
    static Vector pars(2);
    pars(0) = PhaseTag;
    pars(1) = Qgen;
    ///pars(1) = commit_cp;

    return pars;

}



