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
** ****************************************************************** */
                                                                        
// $Revision: 1.12 $
// $Date: 2008-10-20 22:23:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/J2PlasticityThermal.cpp,v $

// Written: Ed "C++" Love
//
// J2 isotropic hardening material class
// 
//  Elastic Model
//  sigma = K*trace(epsilion_elastic) + (2*G)*dev(epsilon_elastic)
//
//  Yield Function
//  phi(sigma,q) = || dev(sigma) ||  - sqrt(2/3)*q(xi) 
//
//  Saturation Isotropic Hardening with linear term
//  q(xi) = simga_infty + (sigma_y - sigma_infty)*exp(-delta*xi) + H*xi 
//
//  Flow Rules
//  \dot{epsilon_p} =  gamma * d_phi/d_sigma
//  \dot{xi}        = -gamma * d_phi/d_q 
//
//  Linear Viscosity 
//  gamma = phi / eta  ( if phi > 0 ) 
//
//  Backward Euler Integration Routine 
//  Yield condition enforced at time n+1 
//
//  set eta := 0 for rate independent case

//Modified for SIF modelling by Liming Jiang [http://openseesforfire.github.io]
//

#include <J2PlasticityThermal.h>
#include <J2PlaneStress.h>
#include <J2PlaneStrain.h>
#include <J2AxiSymm.h>
#include <J2PlateFiber.h>

#include <J2ThreeDimensional.h> 
#include <J2ThreeDimensionalThermal.h> 
#include <string.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

//parameters
const double J2PlasticityThermal :: one3   = 1.0 / 3.0 ;
const double J2PlasticityThermal :: two3   = 2.0 / 3.0 ;
const double J2PlasticityThermal :: four3  = 4.0 / 3.0 ;
const double J2PlasticityThermal :: root23 = sqrt( 2.0 / 3.0 ) ;

double J2PlasticityThermal::initialTangent[3][3][3][3] ;   //material tangent
double J2PlasticityThermal::IIdev[3][3][3][3] ; //rank 4 deviatoric 
double J2PlasticityThermal::IbunI[3][3][3][3] ; //rank 4 I bun I 


void*
OPS_J2PlasticityThermal(void)
{
    NDMaterial* theMaterial = 0;

    int numArgs = OPS_GetNumRemainingInputArgs();

    if (numArgs < 5) {
        opserr << "Want: nDMaterial J2PlasticityThermal $tag $typeTag $E $nu $fy $fyinf \n";
        return 0;
    }

    int iData[2];
    double dData[10];
    dData[4] = 0;
    dData[5] = 0;



    int numData = 2;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "WARNING invalid integer tag: nDMaterial EasticIsotropic \n";
        return 0;
    }

    numData = numArgs - 2;;
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING invalid data: nDMaterial EasticIsotropic : " << iData[0] << "\n";
        return 0;
    }

    double K = dData[0] / 2 / (1 - dData[1]);
    double G = dData[0] / 2 / (1 + dData[1]);


    theMaterial = new J2PlasticityThermal(iData[0], 0,iData[1],
        K, G, dData[2], dData[3], dData[4], dData[5]);

    return theMaterial;
}


//zero internal variables
void J2PlasticityThermal :: zero ( ) 
{
  xi_n = 0.0 ;
  xi_nplus1 = 0.0 ;
  
  epsilon_p_n.Zero( ) ;
  epsilon_p_nplus1.Zero( ) ;

  stress.Zero();
  strain.Zero();
}


//null constructor
J2PlasticityThermal ::  J2PlasticityThermal( ) : 
NDMaterial( ),
epsilon_p_n(3,3),
epsilon_p_nplus1(3,3),
stress(3,3),
strain(3,3),
TempAndElong(2)
{ 
  bulk        = 0.0 ;
  shear       = 0.0 ;
  sigma_y     = 0.0 ;
  bulk_0 = 0.0;
  shear_0 = 0.0;
  sigma_0 = 0.0;

  sigma_infty = 0.0 ;
  sigma_infty0 = 0.0;
  delta       = 0.0 ;
  Hard        = 0.0 ;
  eta         = 0.0 ;
  rho = 0.0;

  this->zero( ) ;     // or (*this).zero( ) 

  int i, j, k, l ;

  //zero rank4 IIdev and IbunI 
  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++ )  {
      for ( k = 0; k < 3; k++ ) {
	for ( l = 0; l < 3; l++)  { 

	  IbunI[i][j][k][l] = 0.0 ;

	  IIdev[i][j][k][l] = 0.0 ;

	} // end for l
      } // end for k
    } // end for j
  } // end for i


  //form rank4 IbunI 

  IbunI [0][0] [0][0] = 1.0 ;
  IbunI [0][0] [1][1] = 1.0 ;
  IbunI [0][0] [2][2] = 1.0 ;
  IbunI [1][1] [0][0] = 1.0 ;
  IbunI [1][1] [1][1] = 1.0 ;
  IbunI [1][1] [2][2] = 1.0 ;
  IbunI [2][2] [0][0] = 1.0 ;
  IbunI [2][2] [1][1] = 1.0 ;
  IbunI [2][2] [2][2] = 1.0 ;

  //form rank4 IIdev

  IIdev [0][0] [0][0] =  two3 ; // 0.666667 
  IIdev [0][0] [1][1] = -one3 ; //-0.333333 
  IIdev [0][0] [2][2] = -one3 ; //-0.333333 
  IIdev [0][1] [0][1] = 0.5 ;
  IIdev [0][1] [1][0] = 0.5 ;
  IIdev [0][2] [0][2] = 0.5 ;
  IIdev [0][2] [2][0] = 0.5 ;
  IIdev [1][0] [0][1] = 0.5 ;
  IIdev [1][0] [1][0] = 0.5 ;
  IIdev [1][1] [0][0] = -one3 ; //-0.333333 
  IIdev [1][1] [1][1] =  two3 ; // 0.666667 
  IIdev [1][1] [2][2] = -one3 ; //-0.333333 
  IIdev [1][2] [1][2] = 0.5 ;
  IIdev [1][2] [2][1] = 0.5 ;
  IIdev [2][0] [0][2] = 0.5 ;
  IIdev [2][0] [2][0] = 0.5 ;
  IIdev [2][1] [1][2] = 0.5 ;
  IIdev [2][1] [2][1] = 0.5 ;
  IIdev [2][2] [0][0] = -one3 ; //-0.333333 
  IIdev [2][2] [1][1] = -one3 ; //-0.333333 
  IIdev [2][2] [2][2] =  two3 ; // 0.666667 

  ThermalElongation = 0.0;
  plastic_integrator();
}


//full constructor
J2PlasticityThermal :: J2PlasticityThermal(int    tag,
			     int classTag,
                 int typeTag,
			     double K,
			     double G,
			     double yield0,
			     double yield_infty,
			     double d,
			     double H,
			     double viscosity,
			     double r) 
: 
  NDMaterial(tag, classTag),
  epsilon_p_n(3,3),
  epsilon_p_nplus1(3,3),
  stress(3,3),
  strain(3,3),
	TempAndElong(2)
{
  bulk        = K ;
  shear       = G ;
  sigma_y     = yield0 ;
  bulk_0 = K;
  shear_0 = G;
  sigma_0 = yield0;
  TypeTag = typeTag;

  sigma_infty = yield_infty ;
  sigma_infty0 = yield_infty;
  delta       = d ;
  Hard        = H ;
  eta         = viscosity ;
  rho = r;

  this->zero( ) ;

  int i, j, k, l ;

  //zero rank4 IIdev and IbunI 
  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++ )  {
      for ( k = 0; k < 3; k++ ) {
	for ( l = 0; l < 3; l++)  { 

	  IbunI[i][j][k][l] = 0.0 ;

	  IIdev[i][j][k][l] = 0.0 ;

	} // end for l
      } // end for k
    } // end for j
  } // end for i


  //form rank4 IbunI 

  IbunI [0][0] [0][0] = 1.0 ;
  IbunI [0][0] [1][1] = 1.0 ;
  IbunI [0][0] [2][2] = 1.0 ;
  IbunI [1][1] [0][0] = 1.0 ;
  IbunI [1][1] [1][1] = 1.0 ;
  IbunI [1][1] [2][2] = 1.0 ;
  IbunI [2][2] [0][0] = 1.0 ;
  IbunI [2][2] [1][1] = 1.0 ;
  IbunI [2][2] [2][2] = 1.0 ;

  //form rank4 IIdev

  IIdev [0][0] [0][0] =  two3 ; // 0.666667 
  IIdev [0][0] [1][1] = -one3 ; //-0.333333 
  IIdev [0][0] [2][2] = -one3 ; //-0.333333 
  IIdev [0][1] [0][1] = 0.5 ;
  IIdev [0][1] [1][0] = 0.5 ;
  IIdev [0][2] [0][2] = 0.5 ;
  IIdev [0][2] [2][0] = 0.5 ;
  IIdev [1][0] [0][1] = 0.5 ;
  IIdev [1][0] [1][0] = 0.5 ;
  IIdev [1][1] [0][0] = -one3 ; //-0.333333 
  IIdev [1][1] [1][1] =  two3 ; // 0.666667 
  IIdev [1][1] [2][2] = -one3 ; //-0.333333 
  IIdev [1][2] [1][2] = 0.5 ;
  IIdev [1][2] [2][1] = 0.5 ;
  IIdev [2][0] [0][2] = 0.5 ;
  IIdev [2][0] [2][0] = 0.5 ;
  IIdev [2][1] [1][2] = 0.5 ;
  IIdev [2][1] [2][1] = 0.5 ;
  IIdev [2][2] [0][0] = -one3 ; //-0.333333 
  IIdev [2][2] [1][1] = -one3 ; //-0.333333 
  IIdev [2][2] [2][2] =  two3 ; // 0.666667 

  ThermalElongation = 0;
  plastic_integrator();
}


//elastic constructor
J2PlasticityThermal :: 
J2PlasticityThermal(   int    tag, 
                int  classTag,
                int typeTag,
                double K, 
                double G ) :
NDMaterial(tag, classTag),
epsilon_p_n(3,3),
epsilon_p_nplus1(3,3),
stress(3,3),
strain(3,3),
TempAndElong(2)
{
  bulk        = K ;
  shear       = G ; 
  bulk_0 = K;
  shear_0 = G;
  sigma_y     = 1.0e16*shear ;
  sigma_0 = 1.0e16*shear;
  sigma_infty = sigma_y ;
  sigma_infty0 = sigma_y;
  delta       = 0.0 ;
  Hard        = 0.0 ;
  eta         = 0.0 ;

  this->zero( ) ;

  int i, j, k, l ;

  //zero rank4 IIdev and IbunI 
  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++ )  {
      for ( k = 0; k < 3; k++ ) {
	for ( l = 0; l < 3; l++)  { 

	  IbunI[i][j][k][l] = 0.0 ;

	  IIdev[i][j][k][l] = 0.0 ;

	} // end for l
      } // end for k
    } // end for j
  } // end for i


  //form rank4 IbunI 

  IbunI [0][0] [0][0] = 1.0 ;
  IbunI [0][0] [1][1] = 1.0 ;
  IbunI [0][0] [2][2] = 1.0 ;
  IbunI [1][1] [0][0] = 1.0 ;
  IbunI [1][1] [1][1] = 1.0 ;
  IbunI [1][1] [2][2] = 1.0 ;
  IbunI [2][2] [0][0] = 1.0 ;
  IbunI [2][2] [1][1] = 1.0 ;
  IbunI [2][2] [2][2] = 1.0 ;

  //form rank4 IIdev

  IIdev [0][0] [0][0] =  two3 ; // 0.666667 
  IIdev [0][0] [1][1] = -one3 ; //-0.333333 
  IIdev [0][0] [2][2] = -one3 ; //-0.333333 
  IIdev [0][1] [0][1] = 0.5 ;
  IIdev [0][1] [1][0] = 0.5 ;
  IIdev [0][2] [0][2] = 0.5 ;
  IIdev [0][2] [2][0] = 0.5 ;
  IIdev [1][0] [0][1] = 0.5 ;
  IIdev [1][0] [1][0] = 0.5 ;
  IIdev [1][1] [0][0] = -one3 ; //-0.333333 
  IIdev [1][1] [1][1] =  two3 ; // 0.666667 
  IIdev [1][1] [2][2] = -one3 ; //-0.333333 
  IIdev [1][2] [1][2] = 0.5 ;
  IIdev [1][2] [2][1] = 0.5 ;
  IIdev [2][0] [0][2] = 0.5 ;
  IIdev [2][0] [2][0] = 0.5 ;
  IIdev [2][1] [1][2] = 0.5 ;
  IIdev [2][1] [2][1] = 0.5 ;
  IIdev [2][2] [0][0] = -one3 ; //-0.333333 
  IIdev [2][2] [1][1] = -one3 ; //-0.333333 
  IIdev [2][2] [2][2] =  two3 ; // 0.666667 

  ThermalElongation = 0.0;
}


//destructor
J2PlasticityThermal :: ~J2PlasticityThermal( ) 
{  } 



NDMaterial*
J2PlasticityThermal :: getCopy (const char *type)
{
    if (strcmp(type,"PlaneStress2D") == 0 || strcmp(type,"PlaneStress") == 0)
    {
	J2PlaneStress  *clone ;
	clone = new J2PlaneStress(this->getTag(), bulk, shear, sigma_y,
				  sigma_infty, delta, Hard, eta, rho) ;
	return clone ;
    }
    else if ((strcmp(type,"ThreeDimensional") == 0) ||
	     (strcmp(type,"3D") == 0))
    {
	J2ThreeDimensional  *clone ;
	clone = new J2ThreeDimensional(this->getTag(), bulk, shear, sigma_y,
				       sigma_infty, delta, Hard, eta, rho) ;
	return clone ;	
    }
	else if ((strcmp(type, "ThreeDimensionalThermal") == 0) ||
		(strcmp(type, "3DThermal") == 0))
	{
		J2ThreeDimensionalThermal  *clone;
		clone = new J2ThreeDimensionalThermal(this->getTag(), TypeTag, bulk, shear, sigma_y,
			sigma_infty, delta, Hard, eta, rho);
		return clone;
	}
    // Handle other cases
    else
    {
      return NDMaterial::getCopy(type);
    }
}

//print out material data
void J2PlasticityThermal :: Print( OPS_Stream &s, int flag )
{
  s << endln ;
  s << "J2-Plasticity : " ; 
  s << this->getType( ) << endln ;
  s << "Bulk Modulus =   " << bulk        << endln ;
  s << "Shear Modulus =  " << shear       << endln ;
  s << "sigma_y =        " << sigma_y     << endln ;
  s << "Sigma_infty =    " << sigma_infty << endln ;
  s << "Delta =          " << delta       << endln ;
  s << "H =              " << Hard        << endln ;
  s << "Eta =            " << eta         << endln ;
  s << "Rho =            " << rho         << endln ;
  s << endln ;
}


//--------------------Plasticity-------------------------------------

//plasticity integration routine
void J2PlasticityThermal :: plastic_integrator( )
{
  const double tolerance = (1.0e-6) ;

  const double dt = ops_Dt ; //time step

  static Matrix dev_strain(3,3) ; //deviatoric strain

  static Matrix dev_stress(3,3) ; //deviatoric stress
 
  static Matrix normal(3,3) ;     //normal to yield surface

  double NbunN ; //normal bun normal   

  double norm_tau = 0.0 ;   //norm of deviatoric stress 
  double inv_norm_tau = 0.0 ;

  double phi = 0.0 ; //trial value of yield function

  double trace = 0.0 ; //trace of strain

  double gamma = 0.0 ; //consistency parameter

  double resid = 1.0 ; 
  double tang  = 0.0 ;
  
  double theta = 0.0 ; 
  double theta_inv = 0.0 ;

  double c1 = 0.0 ; 
  double c2 = 0.0 ;
  double c3 = 0.0 ;

  int i,j,k,l;
  int ii, jj ; 

  int iteration_counter ;
  const int max_iterations = 50 ;

  //compute the deviatoric strains

  trace = strain(0,0) + strain(1,1) + strain(2,2) ;
 
  dev_strain = strain ;
  for ( i = 0; i < 3; i++ )
    dev_strain(i,i) -= ( one3*trace ) ;
   
  //compute the trial deviatoric stresses

  //   dev_stress = (2.0*shear) * ( dev_strain - epsilon_p_n ) ;
  dev_stress = dev_strain;
  dev_stress -= epsilon_p_n;
  dev_stress *= 2.0 * shear;

  //compute norm of deviatoric stress

  norm_tau = 0.0 ;
  for ( i = 0; i < 3; i++ ){
    for ( j = 0; j < 3; j++ ) 
      norm_tau += dev_stress(i,j)*dev_stress(i,j) ;
  } //end for i 
   
  norm_tau = sqrt( norm_tau ) ;

  if ( norm_tau > tolerance ) {
    inv_norm_tau = 1.0 / norm_tau ;
    normal =  inv_norm_tau * dev_stress ;
  }
  else {
    normal.Zero( ) ;
    inv_norm_tau = 0.0 ;
  } //end if 

  //compute trial value of yield function

  phi = norm_tau -  root23 * q(xi_n) ;

  // check if phi > 0 
  
  if ( phi > 0.0 ) { //plastic

     //solve for gamma 
     gamma = 0.0 ;
     resid = 1.0 ;
     iteration_counter = 0 ;
     while ( fabs(resid) > tolerance ) {

        resid = norm_tau 
              - (2.0*shear) * gamma 
              - root23 * q( xi_n + root23*gamma ) 
              - (eta/dt) * gamma ;

        tang =  - (2.0*shear)  
                - two3 * qprime( xi_n + root23*gamma )
                - (eta/dt) ;

        gamma -= ( resid / tang ) ;

	iteration_counter++ ;

	if ( iteration_counter > max_iterations ) {
	    //opserr << "More than " << max_iterations ;
 	    //opserr << " iterations in constituive subroutine J2-plasticity \n" ;
	    break ;
	} //end if 
	
     } //end while resid

     gamma *= (1.0 - 1e-08) ;

     //update plastic internal variables

     epsilon_p_nplus1 = epsilon_p_n + gamma*normal ;

     xi_nplus1 = xi_n + root23*gamma ;

     //recompute deviatoric stresses 

     dev_stress = (2.0*shear) * ( dev_strain - epsilon_p_nplus1 ) ;

     //compute the terms for plastic part of tangent

     theta =  (2.0*shear)  
           +  two3 * qprime( xi_nplus1 )
           +  (eta/dt) ;

     theta_inv = 1.0/theta ;

  }
  else { //elastic 

    //update history variables -- they remain unchanged

    epsilon_p_nplus1 = epsilon_p_n ;

    xi_nplus1 = xi_n ;

    //no extra tangent terms to compute 
    
    gamma = 0.0 ; 
    theta = 0.0 ;
    theta_inv = 0.0 ;

  } //end if phi > 0


  //add on bulk part of stress

  stress = dev_stress ;
  for ( i = 0; i < 3; i++ )
     stress(i,i) += bulk*trace ;

  //compute the tangent

  c1 = -4.0 * shear * shear ;
  c2 = c1 * theta_inv ;
  c3 = c1 * gamma * inv_norm_tau ;

  for ( ii = 0; ii < 6; ii++ ) {
    for ( jj = 0; jj < 6; jj++ )  {

          index_map( ii, i, j ) ;
          index_map( jj, k, l ) ;

          NbunN  = normal(i,j)*normal(k,l) ; 

          //elastic terms
          tangent[i][j][k][l]  = bulk * IbunI[i][j][k][l] ;

          tangent[i][j][k][l] += (2.0*shear) * IIdev[i][j][k][l] ;

          //plastic terms 
          tangent[i][j][k][l] += c2 * NbunN ;

	  tangent[i][j][k][l] += c3 * (  IIdev[i][j][k][l] - NbunN ) ;

          //minor symmetries 
          tangent [j][i][k][l] = tangent[i][j][k][l] ;
          tangent [i][j][l][k] = tangent[i][j][k][l] ;
          tangent [j][i][l][k] = tangent[i][j][k][l] ;

    } // end for jj
  } // end for ii

  return ;
} 




double
J2PlasticityThermal::setThermalTangentAndElongation(double &tempT, double&ET, double&Elong)
{
    double E0 = 2e11;
    double TempT = tempT;
    //double E00; //Initial tangent 
    double E = E0;
    double fy;
    double fy0 = sigma_y;

    // EN 1992&1993
    //typeTag:3   EC3 Structural Steel
    //typeTag:21  EC2 Reinforcing Steel EC2 NHotRolled
    //typeTag:22  EC2 Reinforcing Steel EC2 NCold formed
    //typeTag:23  EC2 Reinforcing Steel EC2 X

    double FyRfactors[12];
    double FpRfactors[12];
    double E0Rfactors[12];

    if (TypeTag == 0 || TypeTag == 3) {
        double FyRfEC3[12] = { 1.0, 1.0 ,1.0, 1.0 ,0.78, 0.47, 0.23, 0.11, 0.06, 0.04 ,0.02, 0.0 };
        double FpRfEC3[12] = { 1.0, 0.807 ,0.613, 0.420 ,0.36, 0.18, 0.075, 0.050, 0.0375, 0.025 ,0.0125, 0.0 };
        double E0RfEC3[12] = { 1.0, 0.9, 0.8 ,0.7, 0.6 ,0.31, 0.13, 0.09, 0.0675, 0.045, 0.0225 , 0.0 };
        for (int i = 0; i < 12; i++) {
            FyRfactors[i] = FyRfEC3[i];
            FpRfactors[i] = FpRfEC3[i];
            E0Rfactors[i] = E0RfEC3[i];
        }

    }
    else if (TypeTag == 21) {
        double FyRfEC21[12] = { 1.0, 1.0 ,1.0, 1.0 ,0.78, 0.47, 0.23, 0.11, 0.06, 0.04 ,0.02, 0.0 };
        double FpRfEC21[12] = { 1.0, 0.81 ,0.61, 0.42 ,0.36, 0.18, 0.07, 0.05, 0.04, 0.02 ,0.01, 0.0 };
        double E0RfEC21[12] = { 1.0, 0.9, 0.8 ,0.7, 0.6 ,0.31, 0.13, 0.09, 0.07, 0.04, 0.02 , 0.0 };
        for (int i = 0; i < 12; i++) {
            FyRfactors[i] = FyRfEC21[i];
            FpRfactors[i] = FpRfEC21[i];
            E0Rfactors[i] = E0RfEC21[i];
        }

    }
    else if (TypeTag == 22) {
        double FyRfEC22[12] = { 1.0, 1.0 ,1.0, 0.94 ,0.67, 0.40, 0.12, 0.11, 0.08, 0.05 ,0.03, 0.0 };
        double FpRfEC22[12] = { 0.96 ,0.92, 0.81 ,0.63, 0.44, 0.26, 0.08, 0.06, 0.05 ,0.03, 0.02, 0.0 };
        double E0RfEC22[12] = { 1.0, 0.87, 0.72 ,0.56, 0.40 ,0.24, 0.08, 0.06, 0.05, 0.03, 0.02 , 0.0 };
        for (int i = 0; i < 12; i++) {
            FyRfactors[i] = FyRfEC22[i];
            FpRfactors[i] = FpRfEC22[i];
            E0Rfactors[i] = E0RfEC22[i];
        }
    }
    else if (TypeTag == 23) {
        double FyRfEC23[12] = { 1.0, 1.0 ,1.0, 0.90 ,0.70, 0.47, 0.23, 0.11, 0.06, 0.04 ,0.02, 0.0 };
        double FpRfEC23[12] = { 1.00 ,0.87, 0.74 ,0.70, 0.51, 0.18, 0.07, 0.05, 0.04 ,0.02, 0.01, 0.0 };
        double E0RfEC23[12] = { 1.0, 0.95, 0.90 ,0.75, 0.60 ,0.31, 0.13, 0.09, 0.07, 0.04, 0.02 , 0.0 };
        for (int i = 0; i < 12; i++) {
            FyRfactors[i] = FyRfEC23[i];
            FpRfactors[i] = FpRfEC23[i];
            E0Rfactors[i] = E0RfEC23[i];
        }
    }
    else
        opserr << "WARNING SteelECThermal received an invalid typeTag: " << TypeTag << endln;

    //Now Updating modulus, strengths
    for (int i = 0; i < 13; i++) {
        if (TempT <= 80 + 100 * i)
        {
            if (i == 0) {
                fy = fy0 * (1.0 - TempT * (1.0 - FyRfactors[0]) / 80);
                E = E0 * (1.0 - TempT * (1.0 - E0Rfactors[0]) / 80);
            }
            else if (i == 12) {
                opserr << "Warning:The temperature " << TempT << " for SteelECthermal is out of range\n";
                return -1;
            }
            else {
                fy = fy0 * (FyRfactors[i - 1] - (TempT + 20 - 100 * i) * (FyRfactors[i - 1] - FyRfactors[i]) / 100);
                E = E0 * (E0Rfactors[i - 1] - (TempT + 20 - 100 * i) * (E0Rfactors[i - 1] - E0Rfactors[i]) / 100);
            }
            break;
        }
    }
#ifdef _BDEBUG
    //opserr<<", TempT:"<<TempT<< " fy: "<< fy<< " fp: "<< fp <<" E0T:  "<< E0<<endln;
#endif
    //E = E0;
    // caculation of thermal elongation
    double ThermalElongation = 0;
    if (TempT <= 1) {
        ThermalElongation = TempT * 1.2164e-5;
    }
    else if (TempT <= 730) {
        ThermalElongation = -2.416e-4 + 1.2e-5 * (TempT + 20) + 0.4e-8 * (TempT + 20) * (TempT + 20);
    }
    else if (TempT <= 840) {
        ThermalElongation = 11e-3;
    }
    else if (TempT <= 1180) {
        ThermalElongation = -6.2e-3 + 2e-5 * (TempT + 20);
    }
    else {
        opserr << " SteelEC Temperature " << TempT << " is invalid\n";
        return -1;
    }


    sigma_y = fy / fy0 * sigma_0;
    sigma_infty = fy / fy0 * sigma_infty0;
    bulk = bulk_0 * E / E0;
    shear = shear_0 * E / E0;
    //H = fy / fy0 * 200e6;
    //HT = H * fy / fy0;
  //  sigma_infty = fy_inf - (fy_inf - fy) * exp(-d * kxi_Commit) + HT * kxi_Commit;


	
	//ThermalElongation = 12e-6*(tempT);
	TempAndElong(0) = TempT - 20;
	TempAndElong(1) = ThermalElongation;
	//bulk = bulk_0;
	//shear = shear_0;
	//sigma_y = sigma_0;

	//ET = 2E11;
	//ET = E;  
	//ET = 3.84e10;
	Elong = ThermalElongation;
	this->plastic_integrator();
	return 0;
}




// set up for initial elastic
void J2PlasticityThermal :: doInitialTangent( )
{
  int ii,jj,i,j,k,l;

  //compute the deviatoric strains
  for ( ii = 0; ii < 6; ii++ ) {
    for ( jj = 0; jj < 6; jj++ )  {

          index_map( ii, i, j ) ;
          index_map( jj, k, l ) ;

          //elastic terms
          initialTangent[i][j][k][l]  = bulk * IbunI[i][j][k][l] ;
          initialTangent[i][j][k][l] += (2.0*shear) * IIdev[i][j][k][l] ;

          //minor symmetries 
          //minor symmetries 
          initialTangent [j][i][k][l] = initialTangent[i][j][k][l] ;
          initialTangent [i][j][l][k] = initialTangent[i][j][k][l] ;
          initialTangent [j][i][l][k] = initialTangent[i][j][k][l] ;

    } // end for jj
  } // end for ii

  return ;
} 



//hardening function
double J2PlasticityThermal :: q( double xi ) 
{
//  q(xi) = simga_infty + (sigma_y - sigma_infty)*exp(-delta*xi) + H*xi 

 return    sigma_infty
         + (sigma_y - sigma_infty)*exp(-delta*xi)
         + Hard*xi ;
}


//hardening function derivative
double J2PlasticityThermal :: qprime( double xi )
{
  return  (sigma_y - sigma_infty) * (-delta) * exp(-delta*xi)
         + Hard ;
}


//matrix_index ---> tensor indices i,j
void J2PlasticityThermal :: index_map( int matrix_index, int &i, int &j )
{
  switch ( matrix_index+1 ) { //add 1 for standard tensor indices

    case 1 :
      i = 1 ; 
      j = 1 ;
      break ;
 
    case 2 :
      i = 2 ;
      j = 2 ; 
      break ;

    case 3 :
      i = 3 ;
      j = 3 ;
      break ;

    case 4 :
      i = 1 ;
      j = 2 ;
      break ;

    case 5 :
      i = 2 ;
      j = 3 ;
      break ;

    case 6 :
      i = 3 ;
      j = 1 ;
      break ;


    default :
      i = 1 ;
      j = 1 ;
      break ;

  } //end switch

i-- ; //subtract 1 for C-indexing
j-- ;

return ; 
}


NDMaterial*
J2PlasticityThermal::getCopy (void)
{
  opserr << "J2PlasticityThermal::getCopy -- subclass responsibility\n"; 
  exit(-1);
  return 0;
}

const char*
J2PlasticityThermal::getType (void) const
{
    opserr << "J2PlasticityThermal::getType -- subclass responsibility\n";
    exit(-1);
    return 0;
}

int
J2PlasticityThermal::getOrder (void) const
{
    opserr << "J2PlasticityThermal::getOrder -- subclass responsibility\n";
    exit(-1);
    return 0;
}


int 
J2PlasticityThermal::commitState( ) 
{
  epsilon_p_n = epsilon_p_nplus1 ;
  xi_n        = xi_nplus1 ;

  return 0;
}

int 
J2PlasticityThermal::revertToLastCommit( ) 
{
  return 0;
}


int 
J2PlasticityThermal::revertToStart( ) {

  // added: C.McGann, U.Washington for InitialStateAnalysis
  if (ops_InitialStateAnalysis) {
	// do nothing, keep state variables from last step
  } else {
	// normal call for revertToStart (not initialStateAnalysis)
    this->zero( ) ;
  }

  return 0;
}

int
J2PlasticityThermal::sendSelf(int commitTag, Channel &theChannel)
{
  // we place all the data needed to define material and it's state
  // int a vector object
  static Vector data(10+9);
  int cnt = 0;
  data(cnt++) = this->getTag();
  data(cnt++) = bulk;
  data(cnt++) = shear;
  data(cnt++) = sigma_y;
  data(cnt++) = sigma_infty;
  data(cnt++) = delta;
  data(cnt++) = Hard;
  data(cnt++) = eta;
  data(cnt++) = rho;

  data(cnt++) = xi_n;

  for (int i=0; i<3; i++) 
    for (int j=0; j<3; j++) 
      data(cnt++) = epsilon_p_n(i,j);


  // send the vector object to the channel
  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "J2PlasticityThermal::sendSelf - failed to send vector to channel\n";
    return -1;
  }

  return 0;
}

int
J2PlasticityThermal::recvSelf (int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  // recv the vector object from the channel which defines material param and state
  static Vector data(10+9);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "J2PlasticityThermal::recvSelf - failed to recv vector from channel\n";
    return -1;
  }

  // set the material parameters and state variables
  int cnt = 0;
  this->setTag(data(cnt++));
  bulk = data(cnt++);
  shear = data(cnt++);
  sigma_y = data(cnt++);
  sigma_infty = data(cnt++);
  delta = data(cnt++);
  Hard = data(cnt++);
  eta = data(cnt++);
  rho = data(cnt++);

  xi_n = data(cnt++);

  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++) 
      epsilon_p_n(i,j) = data(cnt++);

  epsilon_p_nplus1 = epsilon_p_n;
  xi_nplus1        = xi_n;

  return 0;
}
