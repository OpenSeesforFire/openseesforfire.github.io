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

//Written by Liming Jiang
//Based on Lee&Fenves 1998 & 2001, and Omid et al 2013.
// Damaged Plasticity material model for concrete (speically modified for elevated temperature behaviour)
// Created: 07/16
                                                                        
#include <PlasticDamageConcrete3dThermal.h>           
#include <Channel.h>
#include <cmath>
#include <elementAPI.h>

static Vector Iv6(6); 
static Matrix Ivp(6,6); 
static Matrix Idp(6,6); 
static Matrix I(6,6);
static Matrix Id(6,6); 

static int matTag;

void *
OPS_PlasticDamageConcrete3dThermal(void)
{
  NDMaterial *theMaterial = 0;
  
  int numArgs = OPS_GetNumRemainingInputArgs();
  
  if (numArgs < 5 || numArgs > 13) {
    opserr << "Want: nDMaterial PlasticDamageConcrete3dThermal $tag $E $nu $ft $fc <$At $Ac $Dbart $Dbarc $gt $gc>\n";
    return 0;	
  }
  
  int iData[1];
  double dData[8];
  dData[4] = 0.6;
  dData[5] = 1.2;
  dData[6] = 0.4;
  dData[7] = 0.5;
  
  int numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: nDMaterial EasticIsotropic \n";
    return 0;
  }
  
  numData = numArgs - 1;;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid data: nDMaterial EasticIsotropic : " << iData[0] <<"\n";
    return 0;
  }  
  
  theMaterial = new PlasticDamageConcrete3dThermal(iData[0], 
					    dData[0], dData[1], dData[2],dData[3], 
					    dData[4], dData[5], dData[6],dData[7]);
  matTag = iData[0];
  return theMaterial;
}


PlasticDamageConcrete3dThermal::PlasticDamageConcrete3dThermal(int tag,
						 double _e, 
						 double _nu, 
						 double _ft,
						 double _fc,  
						 double _At, 
						 double _Ac,
						 double _Dt,
						 double _Dc)
  :NDMaterial(tag,ND_TAG_PlasticDamageConcrete3dThermal),
 E(_e), nu(_nu), ft(_ft), fc(_fc), At(_At), Ac(_Ac),Dbarc(_Dc),Dbart(_Dt),
 eps(6), sig(6), sige(6), eps_p(6), sigeP(6), TempAndElong(2),
 epsCommit(6), sigCommit(6), sigeCommit(6), eps_pCommit(6), sigePCommit(6),
 Ce(6,6), C(6,6), Ccommit(6,6),ft0(_ft),fc0(_fc),E0(_e),Ce0(6,6)
{
  eps.Zero();
  sig.Zero();
  sige.Zero();
  eps_p.Zero();
  sigeP.Zero();
  Ce.Zero();

  // additional material parameters
  double G   = E/2/(1+nu);        //shear modulus
  double  K   = E/3/(1-2*nu);     // bulk  modulus

  Iv6.Zero(); Iv6(0)=1.;Iv6(1)=1.;Iv6(2)=1.;

  Ivp.Zero(); 
  for (int i=0; i<3; i++) 
    for (int j=0; j<3; j++)
      Ivp(i,j)=1.;

  Idp.Zero(); 
  I.Zero(); 
  Id.Zero(); 
  for (int i=0; i<6; i++) {
    Idp(i,i) = 1.;
    if (i<3) {
      I(i,i) = 1.0;
      Id(i,i) = 1.0;
    } else {
      I(i,i) = 0.5;
      Id(i,i) = 0.5;
    }
  }
  for (int i=0; i<3; i++) 
    for (int j=0; j<3; j++) {
      Id(i,j)=Idp(i,j)-1/3.;
      Idp(i,j) = Id(i,j);
    }

  Ce0.addMatrix(0.0, Ivp, K);
  Ce0.addMatrix(1.0,  Id, 2.*G);
  Ce = Ce0;
  C = Ce;

  kt = 0;
  kc = 0;
  ktCommit = 0;
  kcCommit = 0;

  dc = 0;
  dt = 0;
  dcCommit = 0;
  dtCommit = 0;
  
 // double f2c = 1.16*fc;
  //double k = sqrt(2.0)*(f2c - fc)/(2.*f2c - fc);

  this->commitState();
}

PlasticDamageConcrete3dThermal::PlasticDamageConcrete3dThermal()
  :NDMaterial (0, ND_TAG_PlasticDamageConcrete3dThermal),
   eps(6), sig(6), sige(6), eps_p(6), sigeP(6), TempAndElong(2),
   epsCommit(6), sigCommit(6), sigeCommit(6), eps_pCommit(6), sigePCommit(6),
   Ce(6,6), C(6,6), Ccommit(6,6), At(0), Ac(0), Dbarc(0), Dbart(0)
{

}

PlasticDamageConcrete3dThermal::~PlasticDamageConcrete3dThermal ()
{

}

/*
%% Stress invariant function: octahedral normal and shear stresses
function [sigoct, tauoct] = StrsInvar (sig)

% normal stress
sigoct = (sig(1) + sig(2) + sig(3))/3;

% shear stress
J2 = ((sig(1) - sig(2))^2 + (sig(1) - sig(3))^2 + (sig(2) - sig(3))^2)/6 +...
     (sig(4))^2 + (sig(5))^2 + (sig(6))^2;
tauoct = (2/3*J2)^0.5;
end

%% Stress decomposition function: algebraic approach
function [sigpos, signeg, Qpos, Qneg] = StrsDecA (sig)

% positive and negative stress tensors
sigpos = (sig + abs(sig))/2;
signeg = sig - sigpos;

% projection tensors
Qpos = diag((1+sign(sig))/2);
Qneg = eye(6) - Qpos;
end
*/

void
PlasticDamageConcrete3dThermal::StrsInvar(const Vector &sigT, double &I1, double &J2) {
	//6 stresses-> mean hydrostatic stress: octahedral stress:(s1+s2+s3)/3
  I1 = sigT(0) + sigT(1) + sigT(2);
  if (sigT.Size() == 6)
	  J2 = (pow((sigT(0) - sigT(1)), 2) + pow((sigT(0) - sigT(2)), 2) + pow((sigT(1) - sigT(2)), 2)) / 6. +
	  pow(sigT(3), 2) + pow(sigT(4), 2) + pow(sigT(5), 2);
  else if (sigT.Size() == 3)
	  J2 = (pow((sigT(0) - sigT(1)), 2) + pow((sigT(0) - sigT(2)), 2) + pow((sigT(1) - sigT(2)), 2)) / 6.;
  //the shear stress on the octhedral plane is tauoct = sqrt(2/3*J2)
  
}


void PlasticDamageConcrete3dThermal::StrsDecA(const Vector &sigT, Vector &sigPr, Matrix &Peigvec, double &sigMax) {

  Peigvec.Zero();
  Vector TsigT(6);
  TsigT = sigT;
  int     rot, its, i, j , k ;

  double  g, h, aij, sm, thresh, t, c, s, tau ;

  Matrix  v(3,3) ;
  Vector  d(3) ;
  Vector  a(3) ;
  Vector  b(3) ; 
  Vector  z(3) ;

double tol = 1.0e-08 ;

  // set dataPtr 
 // double *dataPtr = data;

  // set v = M 
  v(0,0) = TsigT(0); v(1,1) = TsigT(1); v(2,2) = TsigT(2);
  v(0,1) = TsigT(3); v(0,2) = TsigT(4); v(1,2) = TsigT(5);
  v(1, 0) = TsigT(3); v(2, 0) = TsigT(4); v(2, 1) = TsigT(5);
  //opserr << v;
  //.... move array into one-d arrays
  a(0) = v(0,1) ;
  a(1) = v(1,2) ;
  a(2) = v(2,0) ;


  for ( i = 0; i < 3; i++ ) {
    d(i) = v(i,i) ;
    b(i) = v(i,i) ;
    z(i) = 0.0 ;

    for ( j = 0; j < 3; j++ ) 
      v(i,j) = 0.0 ;

    v(i,i) = 1.0 ;

  } //end for i

  rot = 0 ;
  its = 0 ;

  sm = fabs(a(0)) + fabs(a(1)) + fabs(a(2)) ;

  while ( sm > tol ) {
    //.... set convergence test and threshold
    if ( its < 3 ) 
      thresh = 0.011*sm ;
    else
      thresh = 0.0 ;
      
    //.... perform sweeps for rotations
    for ( i = 0; i < 3; i++ ) {

      j = (i+1)%3;
      k = (j+1)%3;

      aij  = a(i) ;

      g    = 100.0 * fabs(aij) ;

      if ( fabs(d(i)) + g != fabs(d(i))  ||
	   fabs(d(j)) + g != fabs(d(j))     ) {

	if ( fabs(aij) > thresh ) {

	  a(i) = 0.0 ; 
	  h    = d(j) - d(i) ; 

	  if( fabs(h)+g == fabs(h) )
	    t = aij / h ;
	  else {
	    //t = 2.0 * sign(h/aij) / ( fabs(h/aij) + sqrt(4.0+(h*h/aij/aij)));
	    double hDIVaij = h/aij;
	    if (hDIVaij > 0.0) 
	      t = 2.0 / (  hDIVaij + sqrt(4.0+(hDIVaij * hDIVaij)));
	    else
	      t = - 2.0 / (-hDIVaij + sqrt(4.0+(hDIVaij * hDIVaij)));
	  }

	  //.... set rotation parameters

	  c    = 1.0 / sqrt(1.0 + t*t) ;
	  s    = t*c ;
	  tau  = s / (1.0 + c) ;

	  //.... rotate diagonal terms

	  h    = t * aij ;
	  z(i) = z(i) - h ;
	  z(j) = z(j) + h ;
	  d(i) = d(i) - h ;
	  d(j) = d(j) + h ;

	  //.... rotate off-diagonal terms

	  h    = a(j) ;
	  g    = a[k] ;
	  a(j) = h + s*(g - h*tau) ;
	  a(k) = g - s*(h + g*tau) ;

	  //.... rotate eigenvectors

	  for ( k = 0; k < 3; k++ ) {
	    g      = v(k,i) ;
	    h      = v(k,j) ;
	    v(k,i) = g - s*(h + g*tau) ;
	    v(k,j) = h + s*(g - h*tau) ;
	  } // end for k

	  rot = rot + 1 ;

	} // end if fabs > thresh 
      } //else
      else 
	a(i) = 0.0 ;

    }  // end for i

    //.... update the diagonal terms
    for ( i = 0; i < 3; i++ ) {
      b(i) = b(i) + z(i) ;
      d(i) = b(i) ;
      z(i) = 0.0 ;
    } // end for i

    its += 1 ;

    sm = fabs(a(0)) + fabs(a(1)) + fabs(a(2)) ;

  } //end while sm
  //static Vector  dd(3) ;
  static Vector  dd(3);
  if (d(0)>d(1))
  {
	  if (d(0)>d(2))
	  {
		  dd(0) = 1;
		  if (d(1)>d(2))
		  {
			  dd(1) = 2;
			  dd(2) = 3;
		  }
		  else
		  {
			  dd(1) = 3;
			  dd(2) = 2;
		  }
	  }
	  else
	  {
		  dd(0) = 3;
		  dd(1) = 1;
		  dd(2) = 2;
	  }
  }
  else
  {
	  if (d(1)>d(2))
	  {
		  dd(0) = 2;
		  if (d(0)>d(2))
		  {
			  dd(1) = 1;
			  dd(2) = 3;
		  }
		  else
		  {
			  dd(1) = 3;
			  dd(2) = 1;
		  }
	  }
	  else
	  {
		  dd(0) = 3;
		  dd(1) = 2;
		  dd(2) = 1;
	  }
  }
  
  
	sigMax = d(dd(0)-1);
	//Peigvec = v;
	for (int i = 0;i < 3;i++) {
		sigPr(i) = d(dd(i) - 1);
		for (int j = 0;j < 3;j++)
			Peigvec(i, j) = v(i, dd(j) - 1);
	}
  
	
	/*
	Matrix PeigT(3,3);
	PeigT.addMatrixTranspose(0, Peigvec, 1);
	opserr << "v" << v;
	opserr <<"sig:"<< sigT << endln << sigPr << endln;
	opserr <<"Pig:"<< Peigvec << Peigvec*PeigT << endln;
	opserr << "ok" << endln;
	*/
	

 // return 0;
}

int
PlasticDamageConcrete3dThermal::setTrialStrain(Vector const&v1, Vector const&v2){
  return this->setTrialStrain(v1);
}

int
PlasticDamageConcrete3dThermal::setTrialStrain(const Vector &strain)
{
	//  opserr << "PlasticDamageConcrete3dThermal::setTrialStrain: " << strain << endln;

	// Update material parameters
	double G = E / 2 / (1 + nu);        //shear modulus
	double K = E / 3 / (1 - 2 * nu);     // bulk  modulus


	// bunch of Vectors and Matrices used in the method

	static Vector Deps(6);// increment of strain;
	static Vector Deps_p(6);// increment of plastic strain;
	static Vector Eeps_tr(6);

	static Vector sige_tr(6);

	static Matrix Cbar(6, 6);

	Vector sigmaPr1(3);
	Vector delta_epsPr(3);



	double ftbar = ft;
	double fcbar = fc;

	double f2c = 1.16*fc;
	double alpha = (1.16 - 1) / (2 * 1.16 - 1);
	double alphap = 0.2;
	double gamma = 3;  //gamma for concrete

	double dtotal; // overall consideration of tensile and compressive damage


	double tol = 1.0e-5;

	kt = ktCommit;
	kc = kcCommit;

	double gt = 5000;
	double gc = 50000;
	/*
	Matrix PeigT(3, 3); Matrix Product(3, 3);
	PeigT.addMatrixTranspose(0, Peig, 1);

	opserr << Peig << PeigT;
	Product = Peig*PeigT;
	opserr << Product;
	opserr << SigPr;
	*/
	eps_p = eps_pCommit;
	sigeP = sigePCommit;

	// current strain
	eps = strain;

	// incremental strain
	Eeps_tr = eps - eps_p;
	Deps = eps - epsCommit;

	// trial effective stress
	sige_tr = sigeP + Ce*Deps;
	//sige_tr = Ce*Eeps_tr;

#ifdef _DEBUG
	//if (this->getTag() == 103) {
		//opserr << this->getTag() << "--kt--: " << kt << "--kc--: " << kc << endln;
		//opserr << "---eps:" << eps << "--epsp:" << eps_p << endln;
	//}

#endif

	//Spectral transforamtion;
	Vector Sigtr(6); Vector SigPr(3); Matrix Peig(3, 3);
	double sigMax;
	double fyield; double sigpmp; double sigpmn;
	double beta;  double lamda;
	double I1, J2;

	StrsDecA(sige_tr, SigPr, Peig, sigMax);
	StrsInvar(sige_tr, I1, J2);

	if (sigMax > 0) {
		sigpmp = sigMax;
		sigpmn = 0;
	}
	else {
		sigpmp = 0;
		sigpmn = -sigMax;
	}


	//determine beta (considering updated strength)
	beta = fcbar / ftbar*(1 - alpha) - (1 + alpha);
	//Yield surface
	fyield = (alpha*I1 + sqrt(3 * J2) + beta*sigpmp) / (1 - alpha) - fcbar;
	/*
	if (this->getTag() == 113) {
	opserr << this->getTag() << "-Before-----kt--: " << kt << "--kc--: " << kc << endln;
			opserr << "---eps:" << eps << "--epsp:" << eps_p;
			opserr << "---sigTrial:" << sige_tr <<" Thermal elong"<<TempAndElong(1)<< endln<<endln;
	}
	*/
	if (fyield > tol) {

		//plastic
		
		

		double eta;


		if (sigMax > 0)
			eta = beta;
		else
			eta = 0;
		//determine lamda
		lamda = alpha*I1 + sqrt(3 * J2) + eta*sigMax - (1 - alpha)*fcbar;
		lamda = lamda / (9 * K*alphap*alpha + sqrt(6)*G + eta*(2 * G*sigMax / (sqrt(2 * J2)) - 2 * G*I1 / 3 / (sqrt(2 * J2)) + 3 * K*alphap));

		sigmaPr1 = SigPr - lamda*(2 * G*SigPr / sqrt(2 * J2) - 2 * G*I1 / 3 / sqrt(2 * J2) + 3 * K*alphap);
		/*
		StrsInvar(sigmaPr1, I1, J2);
		sigMax = sigmaPr1(0);
		if (sigMax > 0) {
			sigpmp = sigMax;
			sigpmn = 0;
		}
		else {
			sigpmp = 0;
			sigpmn = -sigMax;
		}

		fyield = (alpha*I1 + sqrt(3 * J2) + beta*sigpmp - gamma*sigpmn)/(1-alpha) - fcbar;
		*/
#ifdef _DEBUG
		//opserr << sigmaPr1;
		if (sigmaPr1(0) != sigmaPr1(0))
			opserr << "error" << endln;
#endif
		//now determine damage variables
		double Qrest, Qresc;
		double ktn1;
		double kcn1;
		double deltaktn1, deltakcn1;
		//double gt, gc;  //failure energy

		//calculate increment of eps_principle
		delta_epsPr = lamda*(SigPr / sqrt(2 * J2) + alphap - I1 / 3 / sqrt(2 * J2));

		double hsigmaSum = 0; double sigmaSum = 0;
		for (int i = 0;i < 3;i++) {
			if (sigmaPr1(i) > 0)
				hsigmaSum = hsigmaSum + sigmaPr1(i);

			sigmaSum = sigmaSum + fabs(sigmaPr1(i));
		}
		double gamma_sigma = hsigmaSum / sigmaSum;

		double ft_r, fc_r;
		double dtbt, dcbc;
		double fit, fic;

		fit = 1 + At*(2 + At)*kt;
		fic = 1 + Ac*(2 + Ac)*kc;
		if (fit < 0)
			opserr << "fit " << fit << " kt " << kt;
		if (fic < 0)
			opserr << "fic " << fic << " kc " << kc;

		dcbc = log(1 - Dbarc) / log((1 + Ac) / 2 / Ac);
		dtbt = log(1 - Dbart) / (log((1 + At) - sqrt(1 + At*At)) - log(2 * At));

		dt = 1 - pow((1 / At*(1 + At - sqrt(fit))), dtbt);
		dc = 1 - pow((1 / Ac*(1 + Ac - sqrt(fic))), dcbc);

		ft_r = (1 - dt)*ftbar;
		fc_r = (1 - dc)*fcbar;

		//calculate incremental k (tension& compression)
		if (delta_epsPr(0) > 0)
			deltaktn1 = gamma_sigma*ft_r / gt*delta_epsPr(0);
		else
			deltaktn1 = 0;
		if (delta_epsPr(2) < 0)
			deltakcn1 = -(1 - gamma_sigma)*fc_r / gc*delta_epsPr(2); //fc_r is positive
		else
			deltakcn1 = 0;

		ktn1 = kt + deltaktn1;
		kcn1 = kc + deltakcn1;

		double deltakt, deltakc;
		double dQdkc, dQdkt;
		Qrest = deltaktn1;
		Qresc = deltakcn1;

		double Qnorm = sqrt(Qrest*Qrest + Qresc*Qresc);

		if (ktn1 < 0)
			ktn1 = 0;
		else if (ktn1 >= 1)
			ktn1 = 1;
		if (kcn1 < 0)
			kcn1 = 0;
		else if (kc >= 1)
			kcn1 = 1;

		if (kt >= 1 || kc >= 1)
			Qnorm = 0;
		//damage completed, jumping to stress process
		int  count=0;
		while (Qnorm > 1e-4) {

			fit = 1 + At*(2 + At)*ktn1;
			fic = 1 + Ac*(2 + Ac)*kcn1;
			double ddt_dkt_u = dtbt *pow((1 + At - sqrt(fit))/At, (dtbt - 1));

			double ddt_dkt = ddt_dkt_u*(1 / 2.0 / sqrt(fit)*At * (2 + At));

			double ddc_dkc_u = dcbc / pow((1 + Ac - sqrt(fic)) / Ac, (dcbc - 1));
			double ddc_dkc= ddc_dkc_u*(1 / 2.0 / sqrt(fic)*Ac * (2 + Ac));
			dQdkt = -1 + gamma_sigma / gt*delta_epsPr(0)*ftbar*(-ddt_dkt);
			dQdkc = -1- (1 - gamma_sigma) / gc*delta_epsPr(2)*fcbar*(-ddc_dkc);
			
			if (fabs(Qrest) > tol& dQdkt != 0)
				deltakt = -Qrest / dQdkt;
			else
				deltakt = 0;

			if (fabs(Qresc) > tol& dQdkc != 0)
				deltakc = -Qresc / dQdkc;
			else
				deltakc = 0;

			ktn1 = ktn1 + deltakt;
			kcn1 = kcn1 + deltakc;

			 //double dht_dkt = (-alpha)*dDc_dkc/

			 //double check
			if (ktn1 < 0)
				ktn1 = 0;
			else if (ktn1 > 1)
				ktn1 = 1;
			if (kcn1 < 0)
				kcn1 = 0;
			else if (kc > 1)
				kcn1 = 1;

			fit = 1 + At*(2 + At)*ktn1;
			fic = 1 + Ac*(2 + Ac)*kcn1;
           
			if (fit < 0)
				opserr << "fit " << fit << " kt " << ktn1;
			if (fic < 0)
				opserr << "fic " << fic << " kc " << kcn1;

			dcbc = log(1 - Dbarc) / log((1 + Ac) / 2 / Ac);
			dtbt = log(1 - Dbart) / (log((1 + At) - sqrt(1 + At*At)) - log(2 * At));

			dt = 1 - pow((1 / At*(1 + At - sqrt(fit))), dtbt);
			dc = 1 - pow((1 / Ac*(1 + Ac - sqrt(fic))), dcbc);

			if (dc != dc|| dt != dt)
				opserr << "invalid" << endln;

			ft_r = (1 - dt)*ftbar;
			fc_r = (1 - dc)*fcbar;

			//calculate incremental k (tension& compression)
			if (delta_epsPr(0) > 0)
				deltaktn1 = gamma_sigma*ft_r / gt*delta_epsPr(0);
			else
				deltaktn1 = 0;
			if (delta_epsPr(2) < 0)
				deltakcn1 = -(1 - gamma_sigma)*fc_r / gc*delta_epsPr(2);
			else
				deltakcn1 = 0;

			
			Qrest = -ktn1 + kt + deltaktn1;
			Qresc = -kcn1 + kc + deltakcn1;
			if (ktn1 != ktn1)
				opserr << ktn1;

			//opserr << "ktn1: " << ktn1;
			Qnorm = sqrt(Qrest*Qrest + Qresc*Qresc);
			count++;
			if (count > 10) {
				//opserr << "count: " << count;
				break;
			}
				
		}
		if (ktn1 > kt)
			kt = ktn1;
		
		if(kcn1>kc)
			kc = kcn1;

		if (kt < 0)
			kt = 0;
		else if (kt > 1)
			kt = 1;
		if (kc < 0)
			kc = 0;
		else if (kc > 1)
			kc = 1;

		//obtain effective stress
		Matrix PeigT(3, 3);
		Matrix tempSigPr(3, 3);
		Matrix tempEpspPr(3, 3);
		Matrix Sign1(3, 3);
		Matrix Depspn1(3, 3);

		tempSigPr.Zero();
		tempEpspPr.Zero();
		tempSigPr(0, 0) = sigmaPr1(0); tempSigPr(1, 1) = sigmaPr1(1); tempSigPr(2, 2) = sigmaPr1(2);
		tempEpspPr(0, 0) = delta_epsPr(0); tempEpspPr(1, 1) = delta_epsPr(1); tempEpspPr(2, 2) = delta_epsPr(2);

		PeigT.addMatrixTranspose(0, Peig, 1);


		Sign1 = Peig*tempSigPr*PeigT;
		Depspn1 = Peig*tempEpspPr*PeigT;
		//#ifdef _DEBUG
		//if (this->getTag() == 113) {
		//	opserr << "lamda"<<lamda<<" fyield:"<<fyield<<"  sigmaPr1:" << sigmaPr1 << endln;
		//	opserr << "epsp:" << delta_epsPr << endln;

		//}
		//#endif
		sige(0) = Sign1(0, 0);sige(1) = Sign1(1, 1);sige(2) = Sign1(2, 2);
		sige(3) = Sign1(0, 1);sige(4) = Sign1(0, 2);sige(5) = Sign1(1, 2);

		Deps_p(0) = Depspn1(0, 0);Deps_p(1) = Depspn1(1, 1);Deps_p(2) = Depspn1(2, 2);
		Deps_p(3) = Depspn1(0, 1);Deps_p(4) = Depspn1(0, 2);Deps_p(5) = Depspn1(1, 2);

		//update plastic strain;
		eps_p = eps_pCommit + Deps_p;

		//obtain incremental plastic strain

		//opserr << "Deps_p" << Deps_p;
		Cbar = Ce;
		if (sige(0) != sige(0))
			opserr << "invalid" << endln;
	}
	else {
		//elastic
		
		dt = dtCommit;
		dc = dcCommit;

		sige = sige_tr;
		Cbar = Ce;
	}
 
	dtotal = 1-(1 - dc)*(1 - dt);   //not adding cyclic effect yet

	if (dtotal == 1)
		dtotal = 0.999; //set up residual stiffness

	sig = (1 - dtotal)*sige;

	C = (1 - dtotal)* Ce;
	/*
	if (this->getTag() == 113) {
	opserr << "-After-kt--: " << kt << "--kc--: " << kc << endln;
			opserr << "--epsp:" << eps_p << endln;
			opserr << "---sig:" << sig << endln;
	}
	*/		
#ifdef _DEBUG
	if (this->getTag() == 113) {
	opserr << this->getTag() << "--kt--: " << kt << "--kc--: " << kc << endln;
		opserr << "---eps:" << eps << "--epsp:" << eps_p << endln;
		opserr << "---sig:" << sig << endln<<" Thermal elong"<<TempAndElong(1)<<endln;
	}
#endif


     return 0;
}

int
PlasticDamageConcrete3dThermal::setTrialStrainIncr (const Vector &strain)
{
  eps += strain;
  this->setTrialStrain(eps);
  return 0;
}

int
PlasticDamageConcrete3dThermal::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  eps += strain;
  this->setTrialStrain(eps);
  return 0;
}

const Matrix&
PlasticDamageConcrete3dThermal::getTangent (void)
{
  return C;
}

const Matrix&
PlasticDamageConcrete3dThermal::getInitialTangent (void)
{
  return Ce;
}

const Vector&
PlasticDamageConcrete3dThermal::getStress (void)
{
	//if (sig(0) > 1e10)
		//opserr << "Stress too big" << sig << endln;
	//opserr << sig;
  return sig;
}

const Vector&
PlasticDamageConcrete3dThermal::getStrain (void)
{
  return eps;
}

int
PlasticDamageConcrete3dThermal::commitState (void)
{
  ktCommit = kt;
  kcCommit = kc;

  epsCommit = eps;
  sigCommit = sig;
  sigeCommit = sige;
  eps_pCommit = eps_p;
  sigePCommit = sige;
  if(this->getTag()==113)
  	opserr<<endln<<endln;

  return 0;
}

int
PlasticDamageConcrete3dThermal::revertToLastCommit (void)
{
  C = Ccommit;
  kt = ktCommit;
  kc = kcCommit;

  eps = epsCommit;
  sig = sigCommit;
  sige = sigeCommit;
  eps_p = eps_pCommit;
  sigeP = sigePCommit;
  return 0;
}

int
PlasticDamageConcrete3dThermal::revertToStart (void)
{
  eps.Zero();
  sig.Zero();
  sige.Zero();
  eps_p.Zero();
  sigeP.Zero();
  Ce.Zero();

  return 0;
}

NDMaterial*
PlasticDamageConcrete3dThermal::getCopy(const char *type)
{
  if (strcmp(type,"ThreeDimensional") == 0 || strcmp(type,"3D") == 0) {
	  //for debugging
	  matTag++;
	  
    PlasticDamageConcrete3dThermal *theCopy =
      new PlasticDamageConcrete3dThermal (matTag, E0, nu, ft0, fc0, At, Ac, Dbart, Dbarc);
  
    return theCopy;  
  } else {
    return 0;
  }
}

NDMaterial*
PlasticDamageConcrete3dThermal::getCopy (void)
{
  PlasticDamageConcrete3dThermal *theCopy =
    new PlasticDamageConcrete3dThermal (this->getTag(), E0, nu, ft0, fc0, At, Ac, Dbart, Dbarc);
  
  return theCopy;
}

const char*
PlasticDamageConcrete3dThermal::getType (void) const
{
  return "ThreeDimensional";
}

int
PlasticDamageConcrete3dThermal::getOrder (void) const
{
  return 6;
}


int 
PlasticDamageConcrete3dThermal::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(10);

  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "PlasticDamageConcrete3dThermal::sendSelf -- could not send Vector\n";
    return res;
  }
  
  return res;
}

int 
PlasticDamageConcrete3dThermal::recvSelf(int commitTag, Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
  static Vector data(10);
  
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "PlasticDamageConcrete3dThermal::sendSelf -- could not send Vector\n";
    return res;
  }

  return res;
}

void 
PlasticDamageConcrete3dThermal::Print(OPS_Stream &s, int flag) {
  opserr << "PlasticDamageConcrete3dThermal: " << this->getTag();
  opserr << "strain: " << eps;
  opserr << "strain: " << sig;
  opserr << "tangent: " << C;
}       



//send back TempAndElong(Liming,UoE)
const Vector&
PlasticDamageConcrete3dThermal::getTempAndElong(void)
{
	return TempAndElong;
}

//Set TemperatureAndElongation

double
PlasticDamageConcrete3dThermal::setThermalTangentAndElongation(double &tempT, double&ET, double&Elong)
{
	
	double Temp = tempT;
	double epsc0 = 0.0025;
	if (Temp <= 80) {
		ft = ft0;
	}
	else if (Temp <= 580) {
		ft = (1.0 - 1.0*(Temp - 80) / 500)*ft0;
		//Ets = (1.0 - 1.0*(Temp -80)/500)*fc0 * 1.5 / epsc0T;
		//Ets = (1.0 - 1.0*(Temp -80)/500)*EtsT;
		if (ft<0.01*ft0)
			ft = 0.01*ft0;
	}
	else {
		ft = 0.01*ft0;
		//Ets = 1.0e-10;
		//ft = 0;
		//Ets = 0;
	}

	// compression strength, at elevated temperature
	//   strain at compression strength, at elevated temperature
	//   ultimate (crushing) strain, at elevated temperature
	if (Temp <= 0) {
		fc = fc0;
		epsc0 = 0.0025;
		//fcu = fcuT;
		//epscu = -0.02;
		//Ets = EtsT;  jz what is there the statement?
	}
	else if (Temp <= 80) {
		fc = fc0;
		epsc0 = (0.0025 + (0.004 - 0.0025)*(Temp - 0) / (80 - 0));
	}
	else if (Temp <= 180) {
		fc = fc0*(1 - (Temp - 80)*0.05 / 100);
		epsc0 = (0.0040 + (0.0055 - 0.0040)*(Temp - 80) / 100);
	}
	else if (Temp <= 280) {
		fc = fc0*(0.95 - (Temp - 180)*0.1 / 100);
		epsc0 = (0.0055 + (0.0070 - 0.0055)*(Temp - 180) / 100);

	}
	else if (Temp <= 380) {
		fc = fc0*(0.85 - (Temp - 280)*0.1 / 100);
		epsc0 = (0.0070 + (0.0100 - 0.0070)*(Temp - 280) / 100);
	}
	else if (Temp <= 480) {
		fc = fc0*(0.75 - (Temp - 380)*0.15 / 100);
		epsc0 = -(0.0100 + (0.0150 - 0.0100)*(Temp - 380) / 100);

	}
	else if (Temp <= 580) {
		fc = fc0*(0.60 - (Temp - 480)*0.15 / 100);
		epsc0 = (0.0150 + (0.0250 - 0.0150)*(Temp - 480) / 100);
	}
	else if (Temp <= 680) {
		fc = fc0*(0.45 - (Temp - 580)*0.15 / 100);
		epsc0 = 0.0250;
	}
	else if (Temp <= 780) {
		fc = fc0*(0.30 - (Temp - 680)*0.15 / 100);
		epsc0 = 0.0250;

	}
	else if (Temp <= 880) {
		fc = fc0*(0.15 - (Temp - 780)*0.07 / 100);
		epsc0 = 0.0250;
	}
	else if (Temp <= 980) {
		fc = fc0*(0.08 - (Temp - 880)*0.04 / 100);
		epsc0 = 0.0250;
	}
	else if (Temp <= 1080) {
		fc = fc0*(0.04 - (Temp - 980)*0.03 / 100);
		epsc0 = 0.0250;
	}
	else {
		opserr << "the temperature is invalid\n";

	}

	double ThermalElongation = 0;
	if (Temp <= 1) {
		ThermalElongation = Temp  * 9.213e-6;
	}
	else if (Temp <= 680) {
		ThermalElongation = -1.8e-4 + 9e-6 *(Temp + 20) + 2.3e-11 *(Temp + 20)*(Temp + 20)*(Temp + 20);
	}
	else if (Temp <= 1180) {
		ThermalElongation = 14e-3;
	}
	else {
		opserr << "the temperature is invalid\n";
	}

	E = fc*0.0025/epsc0/fc0*E0;
	//Ce = -fc*0.0025 / epsc0 / fc0*Ce0;
	//Es = -fc*0.0025/epsc0/fc0*Es0;
	Elong = ThermalElongation;
	// this->setTempInitials();
	// E = E0;
	//Es = Es0;
	TempAndElong(0) = Temp ;
	TempAndElong(1) = ThermalElongation;
	ET = E;
	
	return 0;
}


