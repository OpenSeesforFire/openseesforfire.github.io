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
#include <Vector.h>
#include <Matrix.h>

// Define material tag
#define MAT_TAG_TimberECThermal 1001

class TimberECThermal : public UniaxialMaterial
{
public:
    // Constructor that takes 12 parameters: E1, E2, E3, fc1, fc2, fc3, ft1, ft2, ft3, tempCoeff1, tempCoeff2, tempCoeff3
    TimberECThermal(int tag, double* params);

    // Default constructor
    TimberECThermal(void);

    // Destructor
    virtual ~TimberECThermal();

    // Get material type
    const char* getClassType(void) const { return "TimberECThermal"; };

    // Get a copy of the material
    UniaxialMaterial* getCopy(void);

    // Get initial tangent modulus
    double getInitialTangent(void);

    // Set trial strain (unused, only for pure virtual function definition)
    int setTrialStrain(double strain, double rate);

    // Set trial strain, considering temperature and strain rate
    int setTrialStrain(double strain, double FiberTemperature, double strainRate);

    // Get strain (uniaxial strain, averaged from three directional strains)
    double getStrain(void);

    // Get stress (uniaxial stress, averaged from three directional stresses)
    double getStress(void);

    // Get tangent modulus
    double getTangent(void);

    // Get thermal elongation
    double getThermalElongation(void);

    // Get tangent modulus of thermal elongation (including maximum temperature)
    double getElongTangent(double temperature, double& ET, double& Elong, double TempTmax);

    // Commit state
    int commitState(void);

    // Revert to last committed state
    int revertToLastCommit(void);

    // Revert to initial state
    int revertToStart(void);

    // Send self to channel
    int sendSelf(int commitTag, Channel& theChannel);

    // Receive self from channel
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

    // Print material information
    void Print(OPS_Stream& s, int flag = 0);

    // Get variable information
    int getVariable(const char* variable, Information& theInfo);

protected:

private:
    // Material parameters (three principal directions)
    double E1, E2, E3;           // Elastic modulus
    double fc1, fc2, fc3;        // Compressive strength
    double ft1, ft2, ft3;        // Tensile strength
    double tempCoeff1, tempCoeff2, tempCoeff3; // Temperature coefficients

    // Current stress and strain (three directions)
    double sig1, sig2, sig3;
    double eps1, eps2, eps3;

    // Positive and negative elastic moduli (three directions)
    double Epos1, Eneg1;
    double Epos2, Eneg2;
    double Epos3, Eneg3;

    // Temperature-related variables
    double Temp;         // Current temperature
    double TempMax;      // Maximum temperature reached
    double TempMaxP;     // Maximum temperature at commit
    double ThermalElongation; // Thermal elongation

    // Committed state variables (store the state from the previous step)
    double Csig1, Csig2, Csig3;
    double Ceps1, Ceps2, Ceps3;
    double CEpos1, CEneg1;
    double CEpos2, CEneg2;
    double CEpos3, CEneg3;
    double CTemp;
    double CTempMax;

    // Combined uniaxial stress and strain
    double sig; // Uniaxial stress
    double eps; // Uniaxial strain

    double E; // Combined tangent modulus

    // Helper function: compute stress
    void computeStress();

    // Helper function: update temperature effects
    void updateTemperatureEffects();

    // Set temperature
    void setTemperature(double temperature);
};

#endif













