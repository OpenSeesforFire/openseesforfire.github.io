//===============================================================================
//# COPYRIGHT (C): Woody's license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//# PROJECT:           Object Oriented Finite Element Program
//# PURPOSE:           Elastic Isotropic Material implementation:
//# CLASS:             ElasticIsotropic3DThermalSteel
//#
//# VERSION:           0.61803398874989 (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhaohui Yang, Boris Jeremic (jeremic@ucdavis.edu)
//# PROGRAMMER(S):     Zhaohui Yang, Boris Jeremic
//#
//#
//# DATE:              10Oct2000
//# UPDATE HISTORY:    22Nov2002 small fixes 
//#                    Aug2006   Z.Cheng
//#
//===============================================================================
                                                                        
#include <ElasticIsotropic3DThermalSteel.h>


Matrix ElasticIsotropic3DThermalSteel::D(6,6);	  // global for ElasticIsotropic3DThermalSteel only
Vector ElasticIsotropic3DThermalSteel::sigma(6);	 // global for ElasticIsotropic3DThermalSteel onyl

Tensor ElasticIsotropic3DThermalSteel::Dt(4, def_dim_4, 0.0);
stresstensor ElasticIsotropic3DThermalSteel::Stress;

ElasticIsotropic3DThermalSteel::ElasticIsotropic3DThermalSteel
(int tag, double e, double nu, double rho):
 ElasticIsotropicMaterial (tag, ND_TAG_ElasticIsotropic3DThermalSteel, e, nu, rho),
 epsilon(6),E0T(e) 
{
	E = E0T;
	ThermalElongation = 0.0;
	// Set up the elastic constant matrix for 3D elastic isotropic 
	D.Zero();
}

ElasticIsotropic3DThermalSteel::ElasticIsotropic3DThermalSteel():
 ElasticIsotropicMaterial (0, ND_TAG_ElasticIsotropic3DThermalSteel, 0.0, 0.0, 0.0),
 epsilon(6)
{

}

ElasticIsotropic3DThermalSteel::~ElasticIsotropic3DThermalSteel ()
{

}

int
ElasticIsotropic3DThermalSteel::setTrialStrain (const Vector &v)
{
	epsilon = v;

	return 0;
}

int
ElasticIsotropic3DThermalSteel::setTrialStrain (const Vector &v, const Vector &r)
{
	epsilon = v;

	return 0;
}

int
ElasticIsotropic3DThermalSteel::setTrialStrainIncr (const Vector &v)
{
	epsilon += v;

	return 0;
}

int
ElasticIsotropic3DThermalSteel::setTrialStrainIncr (const Vector &v, const Vector &r)
{
	epsilon += v;

	return 0;
}

const Matrix&
ElasticIsotropic3DThermalSteel::getTangent (void)
{
   double mu2 = E/(1.0+v);
   double lam = v*mu2/(1.0-2.0*v);
   double mu  = 0.50*mu2;

   mu2 += lam;

   D(0,0) = D(1,1) = D(2,2) = mu2;
   D(0,1) = D(1,0) = lam;
   D(0,2) = D(2,0) = lam;
   D(1,2) = D(2,1) = lam;
   D(3,3) = mu;
   D(4,4) = mu;
   D(5,5) = mu;

   return D;
}

const Vector&
ElasticIsotropic3DThermalSteel::getStress (void)
{
       double mu2 = E/(1.0+v);
       double lam = v*mu2/(1.0-2.0*v);
       double mu  = 0.50*mu2;

       mu2 += lam;
       double eps0 = epsilon(0);
       double eps1 = epsilon(1);

       double eps2 = epsilon(2);

       sigma(0) = mu2*eps0 + lam*(eps1+eps2);
       sigma(1) = mu2*eps1 + lam*(eps2+eps0);
       sigma(2) = mu2*eps2 + lam*(eps0+eps1);

       //J.Jiang modify to subtract elongation from mechanical strain
	   //sigma(0) = mu2*eps0 + lam*(eps1+eps2) - E/(1.0+v)/(1.0-2.0*v)*ThermalElongation;
       //sigma(1) = mu2*eps1 + lam*(eps2+eps0) - E/(1.0+v)/(1.0-2.0*v)*ThermalElongation;
       //sigma(2) = mu2*eps2 + lam*(eps0+eps1) - 2*v*E/(1.0+v)/(1.0-2.0*v)*ThermalElongation;

       sigma(3) = mu*epsilon(3);
       sigma(4) = mu*epsilon(4);
       sigma(5) = mu*epsilon(5);

	return sigma;
}

double 
ElasticIsotropic3DThermalSteel::setThermalTangentAndElongation(double &TempT, double&ET, double&Elong)
{

	// EN 1992 pt 1-2-1. Class N hot rolled  reinforcing steel at elevated temperatures
  if (TempT <= 1000) {
		E = E0T;
  }
  else if (TempT <= 200) {
      E = E0T*(1 - (TempT - 100)*0.1/100);
  }
  else if (TempT <= 300) {
      E = E0T*(0.9 - (TempT - 200)*0.1/100);
  }
  else if (TempT <= 400) {
		E = E0T*(0.8 - (TempT - 300)*0.1/100);
  }
  else if (TempT <= 500) {
	  E = E0T*(0.7 - (TempT - 400)*0.1/100);
  }
  else if (TempT <= 600) {
	  E = E0T*(0.6 - (TempT - 500)*0.29/100);
  }
  else if (TempT <= 700) {
	  E = E0T*(0.31 - (TempT - 600)*0.18/100);
  }
  else if (TempT <= 800) {
	  E = E0T*(0.13 - (TempT - 700)*0.04/100);
  }
  else if (TempT <= 900) {
	  E = E0T*(0.09 - (TempT - 800)*0.02/100);
  }
  else if (TempT <= 1000) {
	  E = E0T*(0.0675 - (TempT - 900)*(0.00675 - 0.0045)/100);
  }
  else if (TempT <= 1100) {
	  E = E0T*(0.045 - (TempT - 1000)*(0.0045 - 0.00225)/100);
  }
  else if (TempT <= 1200) {
	  E = E0T*(0.0225 - (TempT - 1100)*0.0225/100);
  }
  else  {
      opserr << "the temperature is invalid\n"; 
  } 

  
  //steel
  if (TempT <= 20) {
	  ThermalElongation = 0.0;
  }
  else if (TempT <= 750) {
      ThermalElongation = -2.416e-4 + 1.2e-5 *TempT + 0.4e-8 *TempT*TempT;
	  //ThermalElongation = 12e-6*TempT;
  }
  else if (TempT <= 860) {
      ThermalElongation = 11e-3;
	  //ThermalElongation = 12e-6*TempT;
  }
  else if (TempT <= 1200) {
      ThermalElongation = -6.2e-3 + 2e-5*TempT;
	  //ThermalElongation = 12e-6*TempT;
  }
  else {
	  opserr << "the temperature is invalid\n";
  }

  //ET = E;  
  ET =E;
  Elong = ThermalElongation;

  return 0;
}

const Vector&
ElasticIsotropic3DThermalSteel::getStrain (void)
{
  return epsilon;
}

int
ElasticIsotropic3DThermalSteel::setTrialStrain (const Tensor &v)
{
    Strain = v;
    return 0;
}

int
ElasticIsotropic3DThermalSteel::setTrialStrain (const Tensor &v, const Tensor &r)
{
    Strain = v;
    return 0;
}

int
ElasticIsotropic3DThermalSteel::setTrialStrainIncr (const Tensor &v)
{
    //opserr << " before set Tri St Incr " << Strain;
    //opserr << " Strain Incr " << v << endlnn;
    Strain = Strain + v;
    //opserr << " after setTrialStrainIncr  " << Strain << endlnn;
    return 0;
}

int
ElasticIsotropic3DThermalSteel::setTrialStrainIncr (const Tensor &v, const Tensor &r)
{
    Strain = Strain + v;
    return 0;
}

const Tensor&
ElasticIsotropic3DThermalSteel::getTangentTensor (void)
{
  setInitElasticStiffness();  
  return Dt;
}

const stresstensor&
ElasticIsotropic3DThermalSteel::getStressTensor (void)
{
  setInitElasticStiffness();
  Stress = Dt("ijkl") * Strain("kl");
  return Stress;
}

const straintensor&
ElasticIsotropic3DThermalSteel::getStrainTensor (void)
{
    return Strain;
}

int
ElasticIsotropic3DThermalSteel::commitState (void)
{
    return 0;
}

int
ElasticIsotropic3DThermalSteel::revertToLastCommit (void)
{
	return 0;
}

int
ElasticIsotropic3DThermalSteel::revertToStart (void)
{
	return 0;
}

NDMaterial*
ElasticIsotropic3DThermalSteel::getCopy (void)
{
	ElasticIsotropic3DThermalSteel *theCopy =
		new ElasticIsotropic3DThermalSteel (this->getTag(), E, v, rho);
	theCopy->epsilon = this->epsilon;
	theCopy->Strain = this->Strain;

	return theCopy;
}

const char*
ElasticIsotropic3DThermalSteel::getType (void) const
{
	return "ThreeDimensional";
}

int
ElasticIsotropic3DThermalSteel::getOrder (void) const
{
	return 6;
}

void
ElasticIsotropic3DThermalSteel::Print(OPS_Stream &s, int flag)
{
	s << "ElasticIsotropic3DThermalSteel" << endln;
	s << "\ttag: " << this->getTag() << endln;
	s << "\tE: " << E << endln;
	s << "\tv: " << v << endln;
	s << "\trho: " << rho << endln;
}


//================================================================================
void ElasticIsotropic3DThermalSteel::setInitElasticStiffness(void)
{        				       
    // Kronecker delta tensor
    tensor I2("I", 2, def_dim_2);

    tensor I_ijkl = I2("ij")*I2("kl");

    //I_ijkl.null_indices();
    tensor I_ikjl = I_ijkl.transpose0110();
    tensor I_iljk = I_ijkl.transpose0111();
    tensor I4s = (I_ikjl+I_iljk)*0.5;
    
    // Building elasticity tensor
    Dt = I_ijkl*( E*v / ( (1.0+v)*(1.0 - 2.0*v) ) ) + I4s*( E / (1.0 + v) );

    return;

}


