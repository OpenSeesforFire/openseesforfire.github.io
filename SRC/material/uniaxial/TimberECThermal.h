/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
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
** ****************************************************************** */

// Added by Liming Jiang (UoE)
// Created: 06/13
//
// Description: This file contains the class definition for 
// TimberECThermal. TimberECThermal is modified from Concrete02Thermal
// TimberECThermal is dedicated to provide a concrete material which 
// strictly satisfy Eurocode regarding the temperature dependent properties.

// Concrete02 is written by FMK in the year of 2006 and based on Concr2.f




#ifndef TimberECThermal_h
#define TimberECThermal_h

#include <UniaxialMaterial.h>

class TimberECThermal : public UniaxialMaterial
{
public:
    TimberECThermal(int tag, double _fc, double _epsc0, double _fcu,
        double _Epos, double _Eneg, double _ft);

    TimberECThermal(void);

    virtual ~TimberECThermal();

    const char* getClassType(void) const { return "TimberECThermal"; };
    double getInitialTangent(void);
    UniaxialMaterial* getCopy(void);


    int setTrialStrain(double strain, double rate);     //JZ this function is no use, just for the definiation of pure virtual function.
    int setTrialStrain(double strain, double FiberTemperature, double strainRate); //***JZ

    double getStrain(void);
    double getStress(void);
    double getTangent(void);

    double getThermalElongation(void); //***JZ
    double getElongTangent(double, double&, double&, double);//***JZ //PK add to include max temp

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel,
        FEM_ObjectBroker& theBroker);

    void Print(OPS_Stream& s, int flag = 0);

    int getVariable(const char* variable, Information&);


protected:

private:

    double Temp;  // Timber temp
    double TempMax;  // layer state tag:: 0:wet wood, 1: dry wood, 2:char layer
    double TempMaxP; //committed fiber state
    double ThermalElongation; // eps(theata) = alpha * temperature
    double fcT;
    double epsc0T;
    double fcuT;
    double ftT;
    double EposT;
    double EnegT;


    // matpar : Concrete FIXED PROPERTIES
    double fc0;    // concrete compression strength           : mp(1)
    double epsc0; // strain at compression strength          : mp(2)
    double fcu0;   // stress at ultimate (crushing) strain    : mp(3)
//    double epscu; // ultimate (crushing) strain              : mp(4)       
//    double rat;   // ratio between unloading slope at epscu and original slope : mp(5)
    double Epos0;  // Modulus of Elasticity for compression   : mp(4)
    double Eneg0;  // Modulus of Elasticity for tension       : mp(5)
    double ft0;    // concrete tensile strength               : mp(6)
    double E;   // tension stiffening slope                : mp(7)



    // hstvP : Concerete HISTORY VARIABLES last committed step
//    double ecminP;  //  hstP(1)
//    double deptP;   //  hstP(2)
    double epsP;  //  = strain at previous converged step
    double sigP;  //  = stress at previous converged step
    double EposP;    //   Compression elasticity modulus at last converged step;
    double EnegP;

    double TempP; //

    // hstv : Concerete HISTORY VARIABLES  current step
//    double ecmin;
//    double dept;
    double sig;
//    double Epos;
//    double Eneg;
    double eps;


};


#endif

