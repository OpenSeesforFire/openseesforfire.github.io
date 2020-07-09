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


#ifndef PlasticDamageConcretePlaneStressThermal_h
#define PlasticDamageConcretePlaneStressThermal_h

 // Based on PlasticDamageConcretePlaneStress
 
// Modified for SIF modelling by Liming Jiang [http://openseesforfire.github.io] 

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

class PlasticDamageConcretePlaneStressThermal : public NDMaterial
{
  public:
  PlasticDamageConcretePlaneStressThermal(int tag, 
	  double E,
	  double nu,
	  double ft,
	  double fc,
	  double Gt ,
	  double Gc ,
	  double At ,
	  double Ac,
	  double Dt,
	  double Dc
	 );
  PlasticDamageConcretePlaneStressThermal();
  ~PlasticDamageConcretePlaneStressThermal();
  
  const char *getClassType(void) const {return "PlasticDamageConcrete3d";};
  
  int setTrialStrain (const Vector &v);
  int setTrialStrain (const Vector &v, const Vector &r);
  int setTrialStrainIncr (const Vector &v);
  int setTrialStrainIncr (const Vector &v, const Vector &r);
  const Matrix &getTangent (void);
  const Matrix &getInitialTangent (void);
  
  const Vector &getStress (void);
  const Vector &getStrain (void);
  
  //Temperature
  double setThermalTangentAndElongation(double &TempT, double &, double &);//L. Jiang [SIF]
  void StrsPr(const Vector&sig, Matrix &Peig, Vector &sigPr);
  const Vector& getTempAndElong( void);  ///L. Jiang [SIF]
  
  int commitState (void);
  int revertToLastCommit (void);
  int revertToStart (void);
  
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
	 double ECommit;    //Commited Modulus;

	 double E0;     // elastic modulus
	 double ft0;    // tensile yield strength
	 double fc0;    // compressive yield strength

	 double At;    // damage parameter
	 double Ac;    // damage parameter
	 double gt;    // Tensile strain energy for full damage
	 double gc;    // Compressivestrain energy  for full damage
	 double gt0;
	 double gc0;


				   // current state variables
	 double kt;    // positive damage variable
	 double kc;    // negative damage variable
	 double dt;    // tensile damage parameter
	 double dc;    // compressive damage parameter
	 double Dbart;    // positive damage variable
	 double Dbarc;    // negative damage variable

	 double xi_p;  

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

	 double Temp;
	 double TempT;

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
	 double epsLitsp;

	 Vector TempAndElong;
};

#endif
