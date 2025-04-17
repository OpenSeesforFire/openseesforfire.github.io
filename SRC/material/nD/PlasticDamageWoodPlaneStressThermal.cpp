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

// Based on the work presented in:
// Lee, J., & Fenves, G. L. (2001). Return-mapping algorithm for plastic-damage models:
// 3-D and plane stress formulation. International Journal for Numerical Methods
// in Engineering, 50(2), 487-C506.
//
// Modified for SIF modelling by Liming Jiang [http://openseesforfire.github.io] 

#include <PlasticDamageWoodPlaneStressThermal.h>
#include <Channel.h>
#include <elementAPI.h>
#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>

static int numPlasticDamageWoodPlaneStressThermal = 0;

void* OPS_PlasticDamageWoodPlaneStressThermal()
{
    if (numPlasticDamageWoodPlaneStressThermal == 0) {
        numPlasticDamageWoodPlaneStressThermal++;
    }

    NDMaterial* theMaterial = 0;

    int numArgs = OPS_GetNumRemainingInputArgs();

    if (numArgs < 12) {
        opserr << "Want: nDMaterial PlasticDamageWoodPlaneStressThermal "
            << "$tag $EposL0 $EnegL0 $ftL0 $fcL0 "
            << "$EposT0 $EnegT0 $ftT0 $fcT0 $G0 $fs0 "
            << "$gf_L $gf_T\n";
        return 0;
    }

    int iData[1];
    double dData[12];

    int numData = 1;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "WARNING invalid integer tag for PlasticDamageWoodPlaneStressThermal\n";
        return 0;
    }

    numData = 12;
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING invalid data for PlasticDamageWoodPlaneStressThermal, tag: "
            << iData[0] << "\n";
        return 0;
    }

    theMaterial = new PlasticDamageWoodPlaneStressThermal(
        iData[0],
        dData[0], // EposL0
        dData[1], // EnegL0
        dData[2], // ftL0
        dData[3], // fcL0
        dData[4], // EposT0
        dData[5], // EnegT0
        dData[6], // ftT0
        dData[7], // fcT0
        dData[8], // G0
        dData[9], // fs0
        dData[10], // gf_L
        dData[11]  // gf_T
    );

    if (!theMaterial) {
        opserr << "WARNING: failed to create PlasticDamageWoodPlaneStressThermal material\n";
        return 0;
    }

    return theMaterial;
}

PlasticDamageWoodPlaneStressThermal::PlasticDamageWoodPlaneStressThermal(int tag,
    double EposL0_, double EnegL0_,
    double ftL0_, double fcL0_,
    double EposT0_, double EnegT0_,
    double ftT0_, double fcT0_,
    double G0_,
    double fs0_,
    double gf_L_, double gf_T_)
    : NDMaterial(tag, ND_TAG_PlasticDamageWood2D),
    EposL0(EposL0_), EnegL0(EnegL0_),
    ftL0(ftL0_), fcL0(fcL0_),
    EposT0(EposT0_), EnegT0(EnegT0_),
    ftT0(ftT0_), fcT0(fcT0_),
    G0(G0_),
    fs0(fs0_),
    EposL(EposL0_), EnegL(EnegL0_),
    ftL(ftL0_), fcL(fcL0_),
    EposT(EposT0_), EnegT(EnegT0_),
    ftT(ftT0_), fcT(fcT0_),
    G(G0_),
    fs(fs0_),

    D_L(0.0), D_T(0.0),
    gf_L(gf_L_), gf_T(gf_T_),
    D_LCommit(0.0), D_TCommit(0.0),

    eps(3), sig(3),
    epsCommit(3), sigCommit(3),
    C(3, 3), Ce0(3, 3),
    TempAndElong(2),
    TempMax(0.0),
    TempMaxP(0.0)
{
    eps.Zero();
    sig.Zero();

    epsCommit.Zero();
    sigCommit.Zero();

    initializeTangent();
}

PlasticDamageWoodPlaneStressThermal::PlasticDamageWoodPlaneStressThermal()
    : NDMaterial(0, ND_TAG_PlasticDamageWood2D),
    EposL0(0.0), EnegL0(0.0),
    ftL0(0.0), fcL0(0.0),
    EposT0(0.0), EnegT0(0.0),
    ftT0(0.0), fcT0(0.0),
    G0(0.0),
    fs0(0.0),
    EposL(0.0), EnegL(0.0),
    ftL(0.0), fcL(0.0),
    EposT(0.0), EnegT(0.0),
    ftT(0.0), fcT(0.0),
    G(0.0),
    fs(0.0),
    D_L(0.0), D_T(0.0),
    gf_L(0.0), gf_T(0.0),
    D_LCommit(0.0), D_TCommit(0.0),
    eps(3), sig(3),
    epsCommit(3), sigCommit(3),
    C(3, 3), Ce0(3, 3),
    TempAndElong(2),
    TempMax(0.0),
    TempMaxP(0.0)
{
}

PlasticDamageWoodPlaneStressThermal::~PlasticDamageWoodPlaneStressThermal()
{
}

void PlasticDamageWoodPlaneStressThermal::initializeTangent()
{
    Ce0.Zero();
    Ce0(0, 0) = EposL0;
    Ce0(1, 1) = EposT0;
    Ce0(2, 2) = G0;

    C = Ce0;
}

int
PlasticDamageWoodPlaneStressThermal::setTrialStrain(const Vector& strain)
{
    if (strain.Size() != 3) {
        opserr << "PlasticDamageWoodPlaneStressThermal::setTrialStrain - strain vector size must be 3\n";
        return -1;
    }

    eps = strain;

    if (TempMaxP > TempMax) {
        TempMax = TempMaxP;
    }

    double Temp = TempMax;
    updateMaterialProperties(Temp);

    double epsL = eps(0);
    double sigmaL = computeStressL(epsL);

    double epsT = eps(1);
    double sigmaT = computeStressT(epsT);

    double gammaLT = eps(2);
    double tauLT = computeShear(gammaLT);

    sig(0) = sigmaL;
    sig(1) = sigmaT;
    sig(2) = tauLT;

    updateDamage();

    updateTangent();

    return 0;
}

int
PlasticDamageWoodPlaneStressThermal::setTrialStrain(const Vector& strain, const Vector& rate)
{
    return this->setTrialStrain(strain);
}

int
PlasticDamageWoodPlaneStressThermal::setTrialStrainIncr(const Vector& dStrain)
{
    if (dStrain.Size() != 3) {
        opserr << "PlasticDamageWoodPlaneStressThermal::setTrialStrainIncr - strain increment vector size must be 3\n";
        return -1;
    }

    eps += dStrain;
    return this->setTrialStrain(eps);
}

int
PlasticDamageWoodPlaneStressThermal::setTrialStrainIncr(const Vector& dStrain, const Vector& rate)
{
    return this->setTrialStrainIncr(dStrain);
}

// L
double
PlasticDamageWoodPlaneStressThermal::computeStressL(double epsL)
{
    double sigma = 0.0;
    double epsc0 = -1.17633 * (fcL / EnegL);

    if (epsL < 0) {
        if (epsL >= -0.85 * (fcL / EnegL)) {
            sigma = EnegL * epsL;
        }
        else if (epsL >= -0.925 * (fcL / EnegL)) {
            sigma = ((-0.15 * fcL) / (0.85 * fcL / EnegL + 0.925 * epsc0))
                * (epsL + (0.85 * fcL / EnegL)) - 0.85 * fcL;
        }
        else if (epsL >= -1.075 * epsc0) {
            sigma = -fcL;
        }
        else if (epsL >= -1.700 * epsc0) {
            sigma = (0.15 * fcL / (0.625 * epsc0)) * (epsL - 1.075 * epsc0) - fcL;
        }
        else {
            sigma = -0.85 * fcL;
        }
    }
    else {
        if (epsL <= (ftL / EposL)) {
            sigma = EposL * epsL;
        }
        else if (epsL <= 26 * (ftL / EposL)) {
            sigma = -0.036 * EposL * (epsL - 26 * (ftL / EposL)) + 0.1 * ftL;
        }
        else {
            sigma = 0.1 * ftL;
        }
    }

    // L damage
    sigma *= (1.0 - D_L);

    return sigma;
}

// T
double
PlasticDamageWoodPlaneStressThermal::computeStressT(double epsT)
{
    double sigma = 0.0;
    double epsc0 = -1.17633 * (fcT / EnegT);

    if (epsT < 0) {
        if (epsT >= -0.85 * (fcT / EnegT)) {
            sigma = EnegT * epsT;
        }
        else if (epsT >= -0.925 * (fcT / EnegT)) {
            sigma = ((-0.15 * fcT) / (0.85 * fcT / EnegT + 0.925 * epsc0))
                * (epsT + (0.85 * fcT / EnegT)) - 0.85 * fcT;
        }
        else if (epsT >= -1.075 * epsc0) {
            sigma = -fcT;
        }
        else if (epsT >= -1.700 * epsc0) {
            sigma = (0.15 * fcT / (0.625 * epsc0)) * (epsT - 1.075 * epsc0) - fcT;
        }
        else {
            sigma = -0.85 * fcT;
        }
    }
    else {
        if (epsT <= (ftT / EposT)) {
            sigma = EposT * epsT;
        }
        else if (epsT <= 26 * (ftT / EposT)) {
            sigma = -0.036 * EposT * (epsT - 26 * (ftT / EposT)) + 0.1 * ftT;
        }
        else {
            sigma = 0.1 * ftT;
        }
    }

    // T damage
    sigma *= (1.0 - D_T);

    return sigma;
}

// LT
double
PlasticDamageWoodPlaneStressThermal::computeShear(double gammaLT)
{
    double tau = 0.0;
    double gammaCrit = fs / G;

    if (gammaLT >= 0.0) {
        if (gammaLT <= gammaCrit) {
            tau = G * gammaLT;
        }
        else if (gammaLT <= 26 * gammaCrit) {
            double tmp = gammaLT - 26 * gammaCrit;
            tau = -0.036 * G * tmp + 0.1 * fs;
        }
        else {
            tau = 0.1 * fs;
        }
    }
    else {
        double gammaAbs = fabs(gammaLT);
        double tauPos = 0.0;
        if (gammaAbs <= gammaCrit) {
            tauPos = G * gammaAbs;
        }
        else if (gammaAbs <= 26 * gammaCrit) {
            double tmp = gammaAbs - 26 * gammaCrit;
            tauPos = -0.036 * G * tmp + 0.1 * fs;
        }
        else {
            tauPos = 0.1 * fs;
        }

        tau = -tauPos;
    }

    // LT no damage
    return tau;
}

void PlasticDamageWoodPlaneStressThermal::updateTangent()
{
    double E_L = (eps(0) < 0) ? EnegL : EposL;
    double E_T = (eps(1) < 0) ? EnegT : EposT;

    C.Zero();
    C(0, 0) = (1.0 - D_L) * E_L;
    C(1, 1) = (1.0 - D_T) * E_T;
    C(2, 2) = G;
}

void PlasticDamageWoodPlaneStressThermal::updateMaterialProperties(double Temp)
{
    if (Temp > TempMax) {
        TempMax = Temp;
    }

    if (Temp <= 20.0) {
        ftL = ftL0;
        EposL = EposL0;
        fcL = fcL0;
        EnegL = EnegL0;
    }
    else if (Temp <= 100.0) {
        ftL = (1.0 - 0.35 * (Temp - 20.0) / 80.0) * ftL0;
        EposL = (1.0 - 0.50 * (Temp - 20.0) / 80.0) * EposL0;
        fcL = fcL0 * (1.0 - ((Temp - 20.0) / 80.0) * 0.75);
        EnegL = EnegL0 * (1.0 - 0.65 * (Temp - 20.0) / 80.0);
    }
    else if (Temp <= 300.0) {
        ftL = (0.65 - 0.65 * (Temp - 100.0) / 200.0) * ftL0;
        EposL = (0.50 - 0.50 * (Temp - 100.0) / 200.0) * EposL0;
        fcL = fcL0 * (0.25 - ((Temp - 100.0) / 200.0) * 0.25);
        EnegL = EnegL0 * (0.35 - 0.35 * (Temp - 100.0) / 200.0);
    }
    else {
        ftL = 1.0e-5;
        EposL = 1.0e-8;
        fcL = 1.0e-5;
        EnegL = 1.0e-8;
    }

    if (Temp <= 20.0) {
        ftT = ftT0;
        EposT = EposT0;
        fcT = fcT0;
        EnegT = EnegT0;
    }
    else if (Temp <= 100.0) {
        ftT = (1.0 - 0.35 * (Temp - 20.0) / 80.0) * ftT0;
        EposT = (1.0 - 0.50 * (Temp - 20.0) / 80.0) * EposT0;
        fcT = fcT0 * (1.0 - ((Temp - 20.0) / 80.0) * 0.75);
        EnegT = EnegT0 * (1.0 - 0.65 * (Temp - 20.0) / 80.0);
    }
    else if (Temp <= 300.0) {
        ftT = (0.65 - 0.65 * (Temp - 100.0) / 200.0) * ftT0;
        EposT = (0.50 - 0.50 * (Temp - 100.0) / 200.0) * EposT0;
        fcT = fcT0 * (0.25 - ((Temp - 100.0) / 200.0) * 0.25);
        EnegT = EnegT0 * (0.35 - 0.35 * (Temp - 100.0) / 200.0);
    }
    else {
        ftT = 1.0e-5;
        EposT = 1.0e-8;
        fcT = 1.0e-5;
        EnegT = 1.0e-8;
    }

    if (Temp <= 20.0) {
        fs = fs0;
    }
    else if (Temp <= 100.0) {
        double ratio = 1.0 - ((Temp - 20.0) / 80.0) * (1.0 - 0.63);
        fs = fs0 * ratio;
    }
    else if (Temp <= 300.0) {
        double ratio100 = 0.63;
        double ratio300 = 1.0e-5 / fs0;
        double alpha = (Temp - 100.0) / 200.0;
        double ratio = ratio100 + (ratio300 - ratio100) * alpha;
        fs = fs0 * ratio;
    }
    else {
        fs = 1.0e-5;
    }

    initializeTangent();
}

void PlasticDamageWoodPlaneStressThermal::updateDamage()
{
    double E_L = (eps(0) < 0) ? EnegL : EposL;
    double E_T = (eps(1) < 0) ? EnegT : EposT;
    double W_L = 0.5 * E_L * eps(0) * eps(0);
    double W_T = 0.5 * E_T * eps(1) * eps(1);

    if (W_L < gf_L) {
            D_L = std::min(W_L / gf_L, 1.0);
        }
    else {
            D_L = std::min(1.0, W_L / gf_L);
        }

    if (W_T < gf_T) {
        D_T = std::min(W_T / gf_T, 1.0);
    }
    else {
        D_T = std::min(1.0, W_T / gf_T);
    }

}

const Matrix&
PlasticDamageWoodPlaneStressThermal::getTangent(void)
{
    return C;
}

const Matrix&
PlasticDamageWoodPlaneStressThermal::getInitialTangent(void)
{
    return Ce0;
}

const Vector&
PlasticDamageWoodPlaneStressThermal::getStress(void)
{
    return sig;
}

const Vector&
PlasticDamageWoodPlaneStressThermal::getStrain(void)
{
    return eps;
}

int
PlasticDamageWoodPlaneStressThermal::commitState(void)
{
    epsCommit = eps;
    sigCommit = sig;
    TempMaxP = TempMax;

    D_LCommit = D_L;
    D_TCommit = D_T;

    return 0;
}

int
PlasticDamageWoodPlaneStressThermal::revertToLastCommit(void)
{
    eps = epsCommit;
    sig = sigCommit;
    TempMax = TempMaxP;

    D_L = D_LCommit;
    D_T = D_TCommit;

    updateTangent();
    return 0;
}

int
PlasticDamageWoodPlaneStressThermal::revertToStart(void)
{
    eps.Zero();
    sig.Zero();

    epsCommit.Zero();
    sigCommit.Zero();

    TempMax = 0.0;
    TempMaxP = 0.0;

    D_L = 0.0;
    D_T = 0.0;

    D_LCommit = 0.0;
    D_TCommit = 0.0;

    C = Ce0;
    return 0;
}

NDMaterial*
PlasticDamageWoodPlaneStressThermal::getCopy(void)
{
    PlasticDamageWoodPlaneStressThermal* theCopy =
        new PlasticDamageWoodPlaneStressThermal(this->getTag(),
            EposL0, EnegL0,
            ftL0, fcL0,
            EposT0, EnegT0,
            ftT0, fcT0,
            G0,
            fs0,
            gf_L, gf_T
        );

    theCopy->EposL = EposL;
    theCopy->EnegL = EnegL;
    theCopy->ftL = ftL;
    theCopy->fcL = fcL;
    theCopy->EposT = EposT;
    theCopy->EnegT = EnegT;
    theCopy->ftT = ftT;
    theCopy->fcT = fcT;
    theCopy->G = G;

    theCopy->fs = fs;

    theCopy->D_L = D_L;
    theCopy->D_T = D_T;

    theCopy->D_LCommit = D_LCommit;
    theCopy->D_TCommit = D_TCommit;

    theCopy->eps = eps;
    theCopy->sig = sig;

    theCopy->epsCommit = epsCommit;
    theCopy->sigCommit = sigCommit;

    theCopy->C = C;
    theCopy->Ce0 = Ce0;

    theCopy->TempAndElong = TempAndElong;
    theCopy->TempMax = TempMax;
    theCopy->TempMaxP = TempMaxP;

    return theCopy;
}

NDMaterial*
PlasticDamageWoodPlaneStressThermal::getCopy(const char* type)
{
    if ((strcmp(type, "PlaneStress") == 0) || (strcmp(type, "PlaneStress2D") == 0)) {
        return this->getCopy();
    }
    return 0;
}

const char*
PlasticDamageWoodPlaneStressThermal::getType(void) const
{
    return "PlaneStress2D";
}

int
PlasticDamageWoodPlaneStressThermal::getOrder(void) const
{
    return 3;
}

int
PlasticDamageWoodPlaneStressThermal::sendSelf(int commitTag, Channel& theChannel)
{
    static Vector data(22);
    data(0) = EposL0;
    data(1) = EnegL0;
    data(2) = ftL0;
    data(3) = fcL0;
    data(4) = EposT0;
    data(5) = EnegT0;
    data(6) = ftT0;
    data(7) = fcT0;
    data(8) = G0;
    data(9) = fs0;

    data(10) = gf_L;
    data(11) = gf_T;

    data(12) = eps(0);
    data(13) = eps(1);
    data(14) = eps(2);

    data(15) = sig(0);
    data(16) = sig(1);
    data(17) = sig(2);

    data(18) = TempMax;
    data(19) = TempMaxP;

    data(20) = D_L;
    data(21) = D_T;


    return theChannel.sendVector(this->getDbTag(), commitTag, data);
}

int
PlasticDamageWoodPlaneStressThermal::recvSelf(int commitTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    static Vector data(22);
    int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "PlasticDamageWoodPlaneStressThermal::recvSelf - failed to receive Vector\n";
        return res;
    }

    EposL0 = data(0);
    EnegL0 = data(1);
    ftL0 = data(2);
    fcL0 = data(3);
    EposT0 = data(4);
    EnegT0 = data(5);
    ftT0 = data(6);
    fcT0 = data(7);
    G0 = data(8);
    fs0 = data(9);

    gf_L = data(10);
    gf_T = data(11);

    EposL = EposL0;
    EnegL = EnegL0;
    ftL = ftL0;
    fcL = fcL0;
    EposT = EposT0;
    EnegT = EnegT0;
    ftT = ftT0;
    fcT = fcT0;
    G = G0;
    fs = fs0;

    eps(0) = data(12);
    eps(1) = data(13);
    eps(2) = data(14);

    sig(0) = data(15);
    sig(1) = data(16);
    sig(2) = data(17);

    TempMax = data(18);
    TempMaxP = data(19);

    D_L = data(20);
    D_T = data(21);

    updateTangent();

    return 0;
}

void
PlasticDamageWoodPlaneStressThermal::Print(OPS_Stream& s, int flag)
{
    s << "PlasticDamageWoodPlaneStressThermal, tag: " << this->getTag() << "\n";
    s << "L direction material parameters (basic):\n";
    s << "  EposL0: " << EposL0 << ", EnegL0: " << EnegL0 << "\n";
    s << "  ftL0: " << ftL0 << ", fcL0: " << fcL0 << "\n";
    s << "T direction material parameters (basic):\n";
    s << "  EposT0: " << EposT0 << ", EnegT0: " << EnegT0 << "\n";
    s << "  ftT0: " << ftT0 << ", fcT0: " << fcT0 << "\n";
    s << "Shear modulus G0: " << G0 << "\n";
    s << "Shear strength fs0: " << fs0 << "\n";
    s << "Fracture Energies:\n";
    s << "  gf_L: " << gf_L << ", gf_T: " << gf_T << "\n";
    s << "Current material parameters (affected by temperature):\n";
    s << "  EposL: " << EposL << ", EnegL: " << EnegL << "\n";
    s << "  ftL: " << ftL << ", fcL: " << fcL << "\n";
    s << "  EposT: " << EposT << ", EnegT: " << EnegT << "\n";
    s << "  ftT: " << ftT << ", fcT: " << fcT << "\n";
    s << "  G: " << G << "\n";
    s << "  fs: " << fs << "\n";
    s << "Damage Variables:\n";
    s << "  D_L: " << D_L << ", D_T: " << D_T << "\n";
    s << "Stress: " << sig;
    s << "Strain: " << eps;
    s << "Temperature dependence: TempMax = " << TempMax
        << ", TempMaxP = " << TempMaxP << "\n";
}

const Vector&
PlasticDamageWoodPlaneStressThermal::getTempAndElong(void)
{
    return TempAndElong;
}

double
PlasticDamageWoodPlaneStressThermal::setThermalTangentAndElongation(double& tempT, double& ET, double& Elong)
{
    if (tempT > TempMax) {
        TempMaxP = tempT;
    }

    TempAndElong(0) = tempT;
    TempAndElong(1) = Elong;

    updateMaterialProperties(tempT);

    ET = 0.5 * (EposL + EposT);
    Elong = 1.0e-5 * tempT;

    return 0;
}
