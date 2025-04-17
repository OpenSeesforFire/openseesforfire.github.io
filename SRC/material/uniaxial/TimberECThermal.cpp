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



#include <stdlib.h>
#include <TimberECThermal.h>
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>
#include <elementAPI.h>
#include <math.h>

// Define Hill failure criterion coefficients, adjustable based on specific materials
#define HILL_A 1.0
#define HILL_B 1.0
#define HILL_C 1.0

void* OPS_TimberECThermal(void)
{
    // Pointer to the uniaxial material to be returned
    UniaxialMaterial* theMaterial = 0;

    int iData[1];
    double dData[12]; // Extended parameters to accommodate anisotropy and temperature effects
    int numData = 1;

    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid uniaxialMaterial TimberECThermal tag" << endln;
        return 0;
    }

    numData = OPS_GetNumRemainingInputArgs();

    if (numData != 12) { // Adjust the number of parameters as needed
        opserr << "Invalid #args, want: uniaxialMaterial TimberECThermal " << iData[0]
            << " E1 E2 E3 fc1 fc2 fc3 ft1 ft2 ft3 tempCoeff1 tempCoeff2 tempCoeff3\n";
        return 0;
    }

    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "Invalid #args, want: uniaxialMaterial TimberECThermal " << iData[0]
            << " E1 E2 E3 fc1 fc2 fc3 ft1 ft2 ft3 tempCoeff1 tempCoeff2 tempCoeff3\n";
        return 0;
    }

    // Parsing successful, allocate material
    theMaterial = new TimberECThermal(iData[0], dData);

    if (theMaterial == 0) {
        opserr << "WARNING could not create uniaxialMaterial of type TimberECThermal Material\n";
        return 0;
    }

    return theMaterial;
}

// Constructor
TimberECThermal::TimberECThermal(int tag, double* params) :
    UniaxialMaterial(tag, MAT_TAG_TimberECThermal),
    E1(params[0]), E2(params[1]), E3(params[2]),
    fc1(params[3]), fc2(params[4]), fc3(params[5]),
    ft1(params[6]), ft2(params[7]), ft3(params[8]),
    tempCoeff1(params[9]), tempCoeff2(params[10]), tempCoeff3(params[11]),
    sig1(0.0), sig2(0.0), sig3(0.0),
    eps1(0.0), eps2(0.0), eps3(0.0),
    Epos1(E1), Eneg1(E1),
    Epos2(E2), Eneg2(E2),
    Epos3(E3), Eneg3(E3),
    Temp(20.0), TempMax(20.0), TempMaxP(20.0),
    ThermalElongation(0.0),
    Csig1(0.0), Csig2(0.0), Csig3(0.0),
    Ceps1(0.0), Ceps2(0.0), Ceps3(0.0),
    CEpos1(Epos1), CEneg1(Eneg1),
    CEpos2(E2), CEneg2(Eneg2),
    CEpos3(E3), CEneg3(Eneg3),
    CTemp(20.0), CTempMax(20.0),
    sig(0.0), eps(0.0),
    E((E1 + E2 + E3) / 3.0)
{
    // Initialize other variables
}

TimberECThermal::TimberECThermal(void) :
    UniaxialMaterial(0, MAT_TAG_TimberECThermal),
    E1(0.0), E2(0.0), E3(0.0),
    fc1(0.0), fc2(0.0), fc3(0.0),
    ft1(0.0), ft2(0.0), ft3(0.0),
    tempCoeff1(0.0), tempCoeff2(0.0), tempCoeff3(0.0),
    sig1(0.0), sig2(0.0), sig3(0.0),
    eps1(0.0), eps2(0.0), eps3(0.0),
    Epos1(0.0), Eneg1(0.0),
    Epos2(0.0), Eneg2(0.0),
    Epos3(0.0), Eneg3(0.0),
    Temp(20.0), TempMax(20.0), TempMaxP(20.0),
    ThermalElongation(0.0),
    Csig1(0.0), Csig2(0.0), Csig3(0.0),
    Ceps1(0.0), Ceps2(0.0), Ceps3(0.0),
    CEpos1(0.0), CEneg1(0.0),
    CEpos2(0.0), CEneg2(0.0),
    CEpos3(0.0), CEneg3(0.0),
    CTemp(20.0), CTempMax(20.0),
    sig(0.0), eps(0.0),
    E(0.0)
{
    // Initialize other variables
}

// Destructor
TimberECThermal::~TimberECThermal(void)
{
    // No special operations needed
}

// Get a copy of the material
UniaxialMaterial* TimberECThermal::getCopy(void)
{
    double params[12] = { E1, E2, E3, fc1, fc2, fc3, ft1, ft2, ft3, tempCoeff1, tempCoeff2, tempCoeff3 };
    TimberECThermal* theCopy = new TimberECThermal(this->getTag(), params);

    // Copy current stress and strain state
    theCopy->sig1 = this->sig1;
    theCopy->sig2 = this->sig2;
    theCopy->sig3 = this->sig3;
    theCopy->eps1 = this->eps1;
    theCopy->eps2 = this->eps2;
    theCopy->eps3 = this->eps3;
    theCopy->Epos1 = this->Epos1;
    theCopy->Eneg1 = this->Eneg1;
    theCopy->Epos2 = this->Epos2;
    theCopy->Eneg2 = this->Eneg2;
    theCopy->Epos3 = this->Epos3;
    theCopy->Eneg3 = this->Eneg3;
    theCopy->Temp = this->Temp;
    theCopy->TempMax = this->TempMax;
    theCopy->sig = this->sig;
    theCopy->eps = this->eps;
    theCopy->E = this->E;

    return theCopy;
}

// Get initial tangent modulus
double TimberECThermal::getInitialTangent(void)
{
    return E;
}

// Set trial strain (unused, only for pure virtual function definition)
int TimberECThermal::setTrialStrain(double strain, double rate)
{
    opserr << "TimberECThermal::setTrialStrain(double strain, double rate) - should never be called\n";
    return -1;
}

// Set trial strain, considering temperature and strain rate
int TimberECThermal::setTrialStrain(double strain, double FiberTemperature, double strainRate)
{
    // Update temperature
    setTemperature(FiberTemperature);

    // Assign uniaxial strain to three principal directions, assuming uniform distribution
    eps1 = strain;
    eps2 = strain;
    eps3 = strain;

    // Compute stress and tangent modulus
    computeStress();

    return 0;
}

// Compute stress
void TimberECThermal::computeStress()
{
    // Tensile behavior: linear softening
    if (eps1 > 0) { // Tension
        if (eps1 <= ft1 / Epos1) {
            sig1 = Epos1 * eps1;
        }
        else {
            // Linear softening
            sig1 = ft1 - (eps1 - ft1 / Epos1) * (ft1 / (0.1 * (ft1 / Epos1)));
            Epos1 = 0.1 * Epos1; // Reduce stiffness
            if (sig1 < 0) sig1 = 0.0;
        }
    }
    else { // Compression
        if (fabs(sig1) < fc1) {
            sig1 = Eneg1 * eps1;
        }
        else {
            // Using Sandhass asymptotic failure model
            // Ideal elastoplastic: stress remains constant after exceeding limit
            sig1 = -fc1;
            Eneg1 = 1e-8; // Almost zero
        }
    }

    // Similar treatment for eps2 and eps3
    if (eps2 > 0) {
        if (eps2 <= ft2 / Epos2) {
            sig2 = Epos2 * eps2;
        }
        else {
            sig2 = ft2 - (eps2 - ft2 / Epos2) * (ft2 / (0.1 * (ft2 / Epos2)));
            Epos2 = 0.1 * Epos2;
            if (sig2 < 0) sig2 = 0.0;
        }
    }
    else {
        if (fabs(sig2) < fc2) {
            sig2 = Eneg2 * eps2;
        }
        else {
            sig2 = -fc2;
            Eneg2 = 1e-8;
        }
    }

    if (eps3 > 0) {
        if (eps3 <= ft3 / Epos3) {
            sig3 = Epos3 * eps3;
        }
        else {
            sig3 = ft3 - (eps3 - ft3 / Epos3) * (ft3 / (0.1 * (ft3 / Epos3)));
            Epos3 = 0.1 * Epos3;
            if (sig3 < 0) sig3 = 0.0;
        }
    }
    else {
        if (fabs(sig3) < fc3) {
            sig3 = Eneg3 * eps3;
        }
        else {
            sig3 = -fc3;
            Eneg3 = 1e-8;
        }
    }

    // Stress interaction weakening
    double interactionFactor = 1.0 - 0.1 * (fabs(sig1) + fabs(sig2) + fabs(sig3));
    if (interactionFactor < 0.0) interactionFactor = 0.0;
    sig1 *= interactionFactor;
    sig2 *= interactionFactor;
    sig3 *= interactionFactor;

    // Combined uniaxial stress (simple average)
    sig = (sig1 + sig2 + sig3) / 3.0;

    // Combined tangent modulus
    E = (Epos1 + Epos2 + Epos3) / 3.0;
}

// Update temperature effects
void TimberECThermal::updateTemperatureEffects()
{
    // This function has been replaced or supplemented by the getElongTangent function
    // Further implementation can be added as needed
}

// Set temperature
void TimberECThermal::setTemperature(double temperature)
{
    Temp = temperature;
    if (Temp > TempMax) {
        TempMax = Temp;
    }
    getElongTangent(Temp, E, ThermalElongation, TempMax);
    computeStress();
}

// Get strain
double TimberECThermal::getStrain(void)
{
    // Return uniaxial strain, averaged from three directional strains
    return eps;
}

// Get stress
double TimberECThermal::getStress(void)
{
    // Return uniaxial stress, averaged from three directional stresses
    return sig;
}

// Get tangent modulus
double TimberECThermal::getTangent(void)
{
    return E;
}

// Get thermal elongation
double TimberECThermal::getThermalElongation(void)
{
    return ThermalElongation;
}

// Get tangent modulus of thermal elongation (including maximum temperature)
double TimberECThermal::getElongTangent(double temperature, double& ET, double& Elong, double TempTmax)
{
    // Calculate reduction factors for elastic modulus and strength
    double factor_compression, factor_tension, factor_shear;

    if (temperature <= 20.0) {
        factor_compression = 1.0;
        factor_tension = 1.0;
        factor_shear = 1.0;
    }
    else if (temperature <= 100.0) {
        // Compression: linear decrease from 1 to 0.35
        // factor_compression = 1.0 - (1.0 - 0.35) * (temperature - 20.0) / 80.0;
        // Using tempCoeff1 = 0.008125 per °C
        factor_compression = 1.0 - tempCoeff1 * (temperature - 20.0); // tempCoeff1 = 0.008125

        // Tension: linear decrease from 1 to 0.50
        // factor_tension = 1.0 - (1.0 - 0.50) * (temperature - 20.0) / 80.0;
        // Using tempCoeff2 = 0.00625 per °C
        factor_tension = 1.0 - tempCoeff2 * (temperature - 20.0); // tempCoeff2 = 0.00625

        // Shear: linear decrease from 1 to 0.63
        // factor_shear = 1.0 - (1.0 - 0.63) * (temperature - 20.0) / 80.0;
        // Using tempCoeff3 = 0.004625 per °C
        factor_shear = 1.0 - tempCoeff3 * (temperature - 20.0); // tempCoeff3 = 0.004625
    }
    else if (temperature <= 300.0) {
        // Compression: linear decrease from 0.35 to 0
        factor_compression = 0.35 - 0.35 * (temperature - 100.0) / 200.0; // slope = -0.00175 per °C

        // Tension: linear decrease from 0.50 to 0
        factor_tension = 0.50 - 0.50 * (temperature - 100.0) / 200.0; // slope = -0.0025 per °C

        // Shear: linear decrease from 0.63 to 0
        factor_shear = 0.63 - 0.63 * (temperature - 100.0) / 200.0; // slope = -0.00315 per °C
    }
    else {
        // Above 300°C, factors are 0
        factor_compression = 0.0;
        factor_tension = 0.0;
        factor_shear = 0.0;
    }

    // Ensure factors are within [0,1] range
    factor_compression = (factor_compression < 0.0) ? 0.0 : ((factor_compression > 1.0) ? 1.0 : factor_compression);
    factor_tension = (factor_tension < 0.0) ? 0.0 : ((factor_tension > 1.0) ? 1.0 : factor_tension);
    factor_shear = (factor_shear < 0.0) ? 0.0 : ((factor_shear > 1.0) ? 1.0 : factor_shear);

    // Apply reduction factors to elastic modulus
    Eneg1 = E1 * factor_compression;
    Eneg2 = E2 * factor_compression;
    Eneg3 = E3 * factor_compression;

    Epos1 = E1 * factor_tension;
    Epos2 = E2 * factor_tension;
    Epos3 = E3 * factor_tension;

    // Apply reduction factors to strength
    fc1 = fc1 * factor_compression;
    fc2 = fc2 * factor_compression;
    fc3 = fc3 * factor_compression;

    ft1 = ft1 * factor_tension;
    ft2 = ft2 * factor_tension;
    ft3 = ft3 * factor_tension;

    // If there are shear strength variables, apply factor_shear here
    // For example:
    // fs1 = fs1 * factor_shear;
    // fs2 = fs2 * factor_shear;
    // fs3 = fs3 * factor_shear;

    // Calculate combined tangent modulus ET
    ET = (Epos1 + Epos2 + Epos3 + Eneg1 + Eneg2 + Eneg3) / 6.0;

    // Calculate thermal elongation, assuming a fixed thermal expansion coefficient alpha, e.g., 1e-5 /°C
    double alpha = 1e-5;
    Elong = alpha * temperature;

    // Update maximum temperature
    if (temperature > TempMax) {
        TempMax = temperature;
    }

    return 0.0;
}

// Commit state
int TimberECThermal::commitState()
{
    // Save current state as committed state
    Csig1 = sig1;
    Csig2 = sig2;
    Csig3 = sig3;
    Ceps1 = eps1;
    Ceps2 = eps2;
    Ceps3 = eps3;
    CEpos1 = Epos1;
    CEneg1 = Eneg1;
    CEpos2 = Epos2;
    CEneg2 = Eneg2;
    CEpos3 = Epos3;
    CEneg3 = Eneg3;
    CTemp = Temp;
    CTempMax = TempMax;
    return 0;
}

// Revert to last committed state
int TimberECThermal::revertToLastCommit()
{
    sig1 = Csig1;
    sig2 = Csig2;
    sig3 = Csig3;
    eps1 = Ceps1;
    eps2 = Ceps2;
    eps3 = Ceps3;
    Epos1 = CEpos1;
    Eneg1 = CEneg1;
    Epos2 = CEpos2;
    Eneg2 = CEneg2;
    Epos3 = CEpos3;
    Eneg3 = CEneg3;
    Temp = CTemp;
    TempMax = CTempMax;
    E = (Epos1 + Epos2 + Epos3) / 3.0;
    return 0;
}

// Revert to initial state
int TimberECThermal::revertToStart()
{
    sig1 = sig2 = sig3 = 0.0;
    eps1 = eps2 = eps3 = 0.0;
    Epos1 = E1;
    Eneg1 = E1;
    Epos2 = E2;
    Eneg2 = E2;
    Epos3 = E3;
    Eneg3 = E3;
    Temp = 20.0;
    TempMax = 20.0;
    E = (E1 + E2 + E3) / 3.0;
    return 0;
}

// Send self to channel
int TimberECThermal::sendSelf(int commitTag, Channel& theChannel)
{
    static Vector data(24);
    data(0) = E1;
    data(1) = E2;
    data(2) = E3;
    data(3) = fc1;
    data(4) = fc2;
    data(5) = fc3;
    data(6) = ft1;
    data(7) = ft2;
    data(8) = ft3;
    data(9) = tempCoeff1;
    data(10) = tempCoeff2;
    data(11) = tempCoeff3;

    data(12) = sig1;
    data(13) = sig2;
    data(14) = sig3;
    data(15) = eps1;
    data(16) = eps2;
    data(17) = eps3;
    data(18) = Epos1;
    data(19) = Eneg1;
    data(20) = Epos2;
    data(21) = Eneg2;
    data(22) = Epos3;
    data(23) = Eneg3;

    if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "TimberECThermal::sendSelf() - failed to sendSelf\n";
        return -1;
    }
    return 0;
}

// Receive self from channel
int TimberECThermal::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    static Vector data(24);
    if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "TimberECThermal::recvSelf() - failed to recvSelf\n";
        return -1;
    }

    E1 = data(0);
    E2 = data(1);
    E3 = data(2);
    fc1 = data(3);
    fc2 = data(4);
    fc3 = data(5);
    ft1 = data(6);
    ft2 = data(7);
    ft3 = data(8);
    tempCoeff1 = data(9);
    tempCoeff2 = data(10);
    tempCoeff3 = data(11);

    sig1 = data(12);
    sig2 = data(13);
    sig3 = data(14);
    eps1 = data(15);
    eps2 = data(16);
    eps3 = data(17);
    Epos1 = data(18);
    Eneg1 = data(19);
    Epos2 = data(20);
    Eneg2 = data(21);
    Epos3 = data(22);
    Eneg3 = data(23);

    E = (Epos1 + Epos2 + Epos3) / 3.0;

    return 0;
}

// Print material information
void TimberECThermal::Print(OPS_Stream& s, int flag)
{
    s << "TimberECThermal Material:\n";
    s << "E1: " << E1 << ", E2: " << E2 << ", E3: " << E3 << "\n";
    s << "fc1: " << fc1 << ", fc2: " << fc2 << ", fc3: " << fc3 << "\n";
    s << "ft1: " << ft1 << ", ft2: " << ft2 << ", ft3: " << ft3 << "\n";
    s << "Stress: [" << sig1 << ", " << sig2 << ", " << sig3 << "]\n";
    s << "Strain: [" << eps1 << ", " << eps2 << ", " << eps3 << "]\n";
    s << "Tangent Modulus: " << E << "\n";
    s << "Temperature: " << Temp << "°C\n";
    s << "Max Temperature: " << TempMax << "°C\n";
    s << "Thermal Elongation: " << ThermalElongation << "\n";
}

// Get variable information
int TimberECThermal::getVariable(const char* varName, Information& theInfo)
{
    if (strcmp(varName, "stress") == 0) {
        // Return stresses in three directions
        Vector stress(3);
        stress(0) = sig1;
        stress(1) = sig2;
        stress(2) = sig3;
        theInfo.setVector(stress);
        return 0;
    }
    else if (strcmp(varName, "strain") == 0) {
        // Return strains in three directions
        Vector strain(3);
        strain(0) = eps1;
        strain(1) = eps2;
        strain(2) = eps3;
        theInfo.setVector(strain);
        return 0;
    }
    else if (strcmp(varName, "tangent") == 0) {
        // Return tangent moduli in three directions
        Vector tangent(3);
        tangent(0) = Epos1;
        tangent(1) = Epos2;
        tangent(2) = Epos3;
        theInfo.setVector(tangent);
        return 0;
    }
    else if (strcmp(varName, "thermalElongation") == 0) {
        theInfo.theDouble = ThermalElongation;
        return 0;
    }
    else if (strcmp(varName, "tangentModulus") == 0) {
        theInfo.theDouble = E;
        return 0;
    }
    return -1;
}



