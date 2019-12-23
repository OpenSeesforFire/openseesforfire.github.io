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


#ifndef J2PlaneStressThermal_h
#define J2PlaneStressThermal_h


// Written for SIF modelling by Liming Jiang [http://openseesforfire.github.io] 

#include <ElasticIsotropicMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>


//  double sig[3];
//  double km[9];
//  struct2_T Pres;
//  double eps[3];
//  struct3_T Past;
//  double Deps[3];
//} struct1_T;

class J2PlaneStressThermal : public NDMaterial
{
public:
	J2PlaneStressThermal(int tag,
		double e,
		double nu,
		double fy_0,
		double fy_infty,
		double d0,
		double H0
		);
	J2PlaneStressThermal();
	~J2PlaneStressThermal();

	const char *getClassType(void) const { return "J2PlaneStressThermal"; };

	int setTrialStrain(const Vector &v);
	int setTrialStrain(const Vector &v, const Vector &r);
	int setTrialStrainIncr(const Vector &v);
	int setTrialStrainIncr(const Vector &v, const Vector &r);
	const Matrix &getTangent(void);
	const Matrix &getInitialTangent(void);

	const Vector &getStress(void);
	const Vector &getStrain(void);

	//Temperature
	double setThermalTangentAndElongation(double &TempT, double &, double &);//L. Jiang [SIF]
	void StrsPr(const Vector&sig, Matrix &Peig, Vector &sigPr);
	const Vector& getTempAndElong(void);  ///L. Jiang [SIF]

	double q(double xi);
	//double qprime(double xi);

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	NDMaterial*getCopy(const char *type);
	NDMaterial *getCopy(void);
	const char *getType(void) const;
	int getOrder(void) const;

	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel,
		FEM_ObjectBroker &theBroker);
	void Print(OPS_Stream &s, int flag = 0);

protected:

private:
	// parameters
	double E;     // elastic modulus
	double nu;    // Poisson ratio 
	double fy;    // tensile yield strength
	double fy_inf;    // compressive yield strength
	double ECommit;    //Commited Modulus;

	double fyt;
	double E0;     // elastic modulus
	double fy0;    // tensile yield strength
	double fy0_inf;    // compressive yield strength
	double d;
	double H;

	double xi_p;

	Vector eps;   // strain
	Vector sig;   // stress
	Vector sige;  // effective stress
	Vector eps_p; // plastic strain
	Vector sigeP; // effective stress

	double Temp;
	double TempT;

	double kxi;
	double kxi_Commit;
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
	int Cchange;
	int Tchange;

	Vector TempAndElong;
};

#endif
