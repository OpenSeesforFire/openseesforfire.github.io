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
                                                                        
#include <ElasticIsotropic3DThermal.h>           
#include <Channel.h>

Vector ElasticIsotropic3DThermal::sigma(6);
Matrix ElasticIsotropic3DThermal::D(6,6);

ElasticIsotropic3DThermal::ElasticIsotropic3DThermal
(int tag, double e, double nu, double rho, double alpha, bool Soft) :
 ElasticIsotropicMaterial (tag, ND_TAG_ElasticIsotropic3DThermal, e, nu, rho, alpha),
 epsilon(6), Cepsilon(6),E0T(e),Alpha(alpha),Temp(0),ThermalElong(0),softening(true)
{
  E=E0T;
  epsilon.Zero();
  Cepsilon.Zero();
}

ElasticIsotropic3DThermal::ElasticIsotropic3DThermal():
 ElasticIsotropicMaterial (0, ND_TAG_ElasticIsotropic3DThermal,0.0, 0.0, 0.0),
	 Temp(0),ThermalElong(0) , epsilon(6), Cepsilon(6),softening(true)
{
  epsilon.Zero();
  Cepsilon.Zero();
}

ElasticIsotropic3DThermal::~ElasticIsotropic3DThermal ()
{

}

int
ElasticIsotropic3DThermal::setTrialStrain (const Vector &strain)
{
  epsilon = strain;
  return 0;
}

int
ElasticIsotropic3DThermal::setTrialStrain (const Vector &strain, const Vector &rate)
{
  epsilon = strain;
  return 0;
}

int
ElasticIsotropic3DThermal::setTrialStrainIncr (const Vector &strain)
{
  epsilon += strain;
  return 0;
}

int
ElasticIsotropic3DThermal::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  epsilon += strain;
  return 0;
}

const Matrix&
ElasticIsotropic3DThermal::getTangent (void)
{
  double mu2 = E/(1.0+v);
  double lam = v*mu2/(1.0-2.0*v);
  double mu = 0.50*mu2;
  mu2 += lam;

  D(0,0) = D(1,1) = D(2,2) = mu2;
  D(0,1) = D(1,0) = D(0,2) = D(2,0) = D(1,2) = D(2,1) = lam;
  D(3,3) = mu;
  D(4,4) = mu;
  D(5,5) = mu;

  return D;
}

const Matrix&
ElasticIsotropic3DThermal::getInitialTangent (void)
{
  //  return this->getTangent();
  double mu2 = E/(1.0+v);
  double lam = v*mu2/(1.0-2.0*v);
  double mu = 0.50*mu2;
  mu2 += lam;

  D(0,0) = D(1,1) = D(2,2) = mu2;
  D(0,1) = D(1,0) = D(0,2) = D(2,0) = D(1,2) = D(2,1) = lam;
  D(3,3) = mu;
  D(4,4) = mu;
  D(5,5) = mu;

  return D;
}

const Vector&
ElasticIsotropic3DThermal::getStress (void)
{
  double mu2 = E/(1.0+v);
  double lam = v*mu2/(1.0-2.0*v);
  double mu = 0.50*mu2;
  mu2 += lam;

  double eps0 = epsilon(0);
  double eps1 = epsilon(1);
  double eps2 = epsilon(2);

  D(0,0) = D(1,1) = D(2,2) = mu2;
  D(0,1) = D(1,0) = D(0,2) = D(2,0) = D(1,2) = D(2,1) = lam;
  D(3,3) = mu;
  D(4,4) = mu;
  D(5,5) = mu;

  //sigma = D*epsilon;
  sigma(0) = mu2*eps0 + lam*(eps1+eps2);
  sigma(1) = mu2*eps1 + lam*(eps0+eps2);
  sigma(2) = mu2*eps2 + lam*(eps0+eps1);
  sigma(3) = mu*epsilon(3);
  sigma(4) = mu*epsilon(4);
  sigma(5) = mu*epsilon(5);
	
  return sigma;
}
double 
ElasticIsotropic3DThermal::setThermalTangentAndElongation(double &TempT, double&ET, double&Elong)
{

if(softening){
	double Temp = TempT;
  double epsc0 = -0.0025;
 double fc;
  // compression strength, at elevated temperature
  //   strain at compression strength, at elevated temperature
  //   ultimate (crushing) strain, at elevated temperature
  if (Temp <= 0) {
    fc = 1;
    epsc0 = -0.0025;
    //fcu = fcuT;
    //epscu = -0.02;
    //Ets = EtsT;  jz what is there the statement?
  }
  else if (Temp <= 80) {
    fc = 1;
    epsc0 = -(0.0025 + (0.004-0.0025)*(Temp - 0)/(80 - 0));
  }
  else if (Temp <= 180) {
    fc = (1 - (Temp - 80)*0.05/100);
    epsc0 = -(0.0040 + (0.0055-0.0040)*(Temp - 80)/100);
  }
  else if (Temp <= 280) {
    fc = (0.95 - (Temp - 180)*0.1/100);
    epsc0 = -(0.0055 + (0.0070-0.0055)*(Temp - 180)/100);
   
  }
  else if (Temp <= 380) {
    fc = (0.85 - (Temp - 280)*0.1/100);
    epsc0 = -(0.0070 + (0.0100-0.0070)*(Temp - 280)/100);
  }
  else if (Temp <= 480) {
    fc = (0.75 - (Temp - 380)*0.15/100);
    epsc0 = -(0.0100 + (0.0150-0.0100)*(Temp - 380)/100);

  }
  else if (Temp <= 580) {
    fc = (0.60 - (Temp - 480)*0.15/100);
    epsc0 = -(0.0150 + (0.0250-0.0150)*(Temp - 480)/100);
  }
  else if (Temp <= 680) {
    fc = (0.45 - (Temp - 580)*0.15/100);
    epsc0 = -0.0250;
  }
  else if (Temp <= 780) {
    fc = (0.30 - (Temp - 680)*0.15/100);
    epsc0 = -0.0250;

  }
  else if (Temp <= 880) {
    fc = (0.15 - (Temp - 780)*0.07/100);
    epsc0 = -0.0250;
  }
  else if (Temp <= 980) {
    fc =(0.08 - (Temp - 880)*0.04/100);
    epsc0 = -0.0250;
  }
  else if (Temp <= 1080) {
    fc = (0.04 - (Temp - 980)*0.03/100);
    epsc0 = -0.0250;
  }
  else  {
    opserr << "the temperature is invalid\n";
    
  }
  
  double ThermalElongation;
  if (Temp <= 1) {
		  ThermalElongation = Temp  * 9.213e-6;
  }
  else if (Temp <= 680) {
    ThermalElongation = -1.8e-4 + 9e-6 *(Temp+20) + 2.3e-11 *(Temp+20)*(Temp+20)*(Temp+20);
  }
  else if (Temp <= 1180) {
    ThermalElongation = 14e-3;
  }
  else {
    opserr << "the temperature is invalid\n";
  }
  
  E = -fc*0.0025/epsc0*E0T;
  Elong = ThermalElongation;

}
else{
  ET =E0T;
  Elong = Alpha*TempT;
}
   Temp = TempT;
   ThermalElong = Elong ;
  return 0;
}

const Vector&
ElasticIsotropic3DThermal::getStrain (void)
{
  return epsilon;
}

int
ElasticIsotropic3DThermal::commitState (void)
{
  Cepsilon=epsilon;
  return 0;
}

int
ElasticIsotropic3DThermal::revertToLastCommit (void)
{
  epsilon=Cepsilon;
  return 0;
}

int
ElasticIsotropic3DThermal::revertToStart (void)
{
  epsilon.Zero();
  Cepsilon.Zero();
  return 0;
}

NDMaterial*
ElasticIsotropic3DThermal::getCopy (void)
{
  ElasticIsotropic3DThermal *theCopy =
    new ElasticIsotropic3DThermal (this->getTag(), E, v, rho,Alpha);
  
  theCopy->epsilon = epsilon;
  theCopy->Cepsilon = Cepsilon;
  
  return theCopy;
}

const char*
ElasticIsotropic3DThermal::getType (void) const
{
  return "ThreeDimensionalThermal";
}

int
ElasticIsotropic3DThermal::getOrder (void) const
{
  return 6;
}

int 
ElasticIsotropic3DThermal::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(10);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = v;
  data(3) = rho;
  data(4) = Cepsilon(0);
  data(5) = Cepsilon(1);
  data(6) = Cepsilon(2);
  data(7) = Cepsilon(3);
  data(8) = Cepsilon(4);
  data(9) = Cepsilon(5);
  
  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ElasticIsotropic3DThermal::sendSelf -- could not send Vector\n";
    return res;
  }
  
  return res;
}

int 
ElasticIsotropic3DThermal::recvSelf(int commitTag, Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
  static Vector data(10);
  
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ElasticIsotropic3DThermal::sendSelf -- could not send Vector\n";
    return res;
  }

  this->setTag((int)data(0));
  E = data(1);
  v = data(2);
  rho = data(3);
  epsilon(0)=data(4);
  epsilon(1)=data(5);
  epsilon(2)=data(6);
  epsilon(3)=data(7);
  epsilon(4)=data(8);
  epsilon(5)=data(9);

  Cepsilon = epsilon;

  return res;
}

//send back TempAndElong(Liming,UoE)
const Vector& 
ElasticIsotropic3DThermal::getTempAndElong( void)
{
	static Vector TempElong = Vector(2);
	TempElong(0) = Temp;
	TempElong(1) = ThermalElong;
  return TempElong;
}
