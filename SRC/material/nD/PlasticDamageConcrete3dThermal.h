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

#ifndef PlasticDamageConcrete3dThermal_h
#define PlasticDamageConcrete3dThermal_h

// Written: Thanh Do
// Created: 07/16
//
// Description: 
//
// What: "@(#) ElasticIsotropicThreeDimesnional.h, revA"

#include <ElasticIsotropicMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class PlasticDamageConcrete3dThermal : public NDMaterial
{
  public:
  PlasticDamageConcrete3dThermal(int tag, 
			  double E, 
			  double nu, 
			  double ft,
			  double fc,  
			  double At = 0.8, 
			  double Ac = 1.2, 
			  double Dt = 0.5,
			  double Dc = 0.4);
    PlasticDamageConcrete3dThermal();
    ~PlasticDamageConcrete3dThermal();

    const char *getClassType(void) const {return "PlasticDamageConcrete3dThermal";};

    int setTrialStrain (const Vector &v);
    int setTrialStrain (const Vector &v, const Vector &r);
    int setTrialStrainIncr (const Vector &v);
    int setTrialStrainIncr (const Vector &v, const Vector &r);
    const Matrix &getTangent (void);
    const Matrix &getInitialTangent (void);
    
    const Vector &getStress (void);
    const Vector &getStrain (void);
    
    int commitState (void);
    int revertToLastCommit (void);
    int revertToStart (void);

	void StrsDecA(const Vector &sigT, Vector &sigPr, Matrix &Peigvec, double &sigMax);
	void StrsInvar(const Vector &sigT, double &I1, double &J2);

	//Temperature
	double setThermalTangentAndElongation(double &TempT, double &, double &);//L. Jiang [SIF]
	const Vector& getTempAndElong(void);  ///L. Jiang [SIF]

    NDMaterial*getCopy(const char *type);
    NDMaterial *getCopy (void);
    const char *getType (void) const;
    int getOrder (void) const;

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    void Print(OPS_Stream &s, int flag =0);       

  protected:

  private:
    // parameters
    double E;     // elastic modulus
    double nu;    // Poisson ratio 
    double ft;    // tensile yield strength
    double fc;    // compressive yield strength

	double E0;     // elastic modulus
	double ft0;    // tensile yield strength
	double fc0;    // compressive yield strength
    
    double At;    // damage parameter
    double Ac;    // damage parameter

    // current state variables
    double kt;    // positive damage variable
    double kc;    // negative damage variable
	double dt;    // tensile damage parameter
	double dc;    // compressive damage parameter
    double Dbart;    // positive damage variable
    double Dbarc;    // negative damage variable

    Vector eps;   // strain
    Vector sig;   // stress
    Vector sige;  // effective stress
    Vector eps_p; // plastic strain
    Vector sigeP; // effective stress

    // committed state variables
    double ktCommit;
    double kcCommit;
    double dtCommit; 
    double dcCommit; 



    Vector epsCommit;
    Vector sigCommit;
    Vector sigeCommit;  
    Vector eps_pCommit; 
    Vector sigePCommit; 

    // tangent matrices
	Matrix Ce0;
    Matrix Ce; 
    Matrix C; 
    Matrix Ccommit; 

	Vector TempAndElong;
};

#endif
