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

// $Revision: 1.6 $
// $Date: 2006-08-04 18:18:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicPlaneStress2D.h,v $

#ifndef PlasticDamageWoodPlaneStressThermal_h
#define PlasticDamageWoodPlaneStressThermal_h

#ifndef ND_TAG_PlasticDamageWood2D
#define ND_TAG_PlasticDamageWood2D 9999 
#endif

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

class PlasticDamageWoodPlaneStressThermal : public NDMaterial
{
public:

    PlasticDamageWoodPlaneStressThermal(int tag,
        double EposL0, double EnegL0,
        double ftL0, double fcL0,
        double EposT0, double EnegT0,
        double ftT0, double fcT0,
        double G0,
        double fs0,
        double gf_L, double gf_T);

    PlasticDamageWoodPlaneStressThermal();

    ~PlasticDamageWoodPlaneStressThermal();

    int setTrialStrain(const Vector& strain);
    int setTrialStrain(const Vector& strain, const Vector& rate);
    int setTrialStrainIncr(const Vector& dStrain);
    int setTrialStrainIncr(const Vector& dStrain, const Vector& rate);

    const Matrix& getTangent(void);
    const Matrix& getInitialTangent(void);

    const Vector& getStress(void);
    const Vector& getStrain(void);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    NDMaterial* getCopy(void);
    NDMaterial* getCopy(const char* type);

    const char* getType(void) const;
    int getOrder(void) const;

    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

    void Print(OPS_Stream& s, int flag = 0);

    const Vector& getTempAndElong(void);
    double setThermalTangentAndElongation(double& tempT, double& ET, double& Elong);

private:

    // L
    double EposL0;
    double EnegL0;
    double ftL0;
    double fcL0;

    // T
    double EposT0;
    double EnegT0;
    double ftT0;
    double fcT0;

    double G0;

    double fs0;

    double EposL;
    double EnegL;
    double ftL;
    double fcL;

    double EposT;
    double EnegT;
    double ftT;
    double fcT;

    double G;

    double fs;

    double D_L;
    double D_T;

    double gf_L;
    double gf_T;

    double D_LCommit;
    double D_TCommit;

    Vector eps;
    Vector sig;

    Vector epsCommit;
    Vector sigCommit;

    Matrix C;
    Matrix Ce0;

    Vector TempAndElong;
    double TempMax;
    double TempMaxP;

    double computeStressL(double epsL);
    double computeStressT(double epsT);
    double computeShear(double gammaLT);

    void updateTangent();
    void initializeTangent();
    void updateMaterialProperties(double Temp);

    void updateDamage();
};

#endif
