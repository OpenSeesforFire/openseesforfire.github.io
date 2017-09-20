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

// Based on PlasticDamageConcretePlaneStress

// Modified for SIF modelling by Liming Jiang [http://openseesforfire.github.io] 

#include <J2PlaneStressThermal.h>           
#include <Channel.h>
#include <cmath>
#include <elementAPI.h>

//#ifdef _DEBUG

#include <fstream>
#include <string>
#include <sstream>
using std::ios;
//#endif

//#define _DEBUG_PDC_PlaneStress 1

static Vector Iv6(6);
static Matrix Ivp(6, 6);
static Matrix Idp(6, 6);
static Matrix I(6, 6);
static Matrix Id(6, 6);

static int matTag;

//static const signed char b_A[3] = { -1, 1, 0 };
//static const signed char c_a[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
//static const signed char iv1[3] = { 0, 0, 1 };

void *
OPS_J2PlaneStressThermal(void)
{
	NDMaterial *theMaterial = 0;

	int numArgs = OPS_GetNumRemainingInputArgs();

	if (numArgs < 5 || numArgs > 11) {
		opserr << "Want: nDMaterial J2PlaneStressThermal $tag $E $nu $ft $fc <$beta $Ap $An $Bn>\n";
		return 0;
	}

	int iData[1];
	double dData[10];
	dData[4] = 0;
	dData[5] = 0;
	dData[6] = 0.0;


	int numData = 1;
	if (OPS_GetInt(&numData, iData) != 0) {
		opserr << "WARNING invalid integer tag: nDMaterial EasticIsotropic \n";
		return 0;
	}

	numData = numArgs - 1;;
	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "WARNING invalid data: nDMaterial EasticIsotropic : " << iData[0] << "\n";
		return 0;
	}

	theMaterial = new J2PlaneStressThermal(iData[0],
		dData[0], dData[1], dData[2], dData[3],dData[4], dData[5]);
	matTag = iData[0];
	return theMaterial;
}

J2PlaneStressThermal::J2PlaneStressThermal(int tag,
	double e,
	double nu,
	double fy_0,
	double fy_infty,
	double d0,
	double H0
)
	:NDMaterial(tag, 100),
	E(e), nu(nu), fy(fy_0), fy_inf(fy_infty),
	eps(3), sig(3), sige(3), eps_p(3), sigeP(3), TempAndElong(2),
	epsCommit(3), sigCommit(3), sigeCommit(3), eps_pCommit(3),d(d0),H(H0),
	Ce(3, 3), C(3, 3), Ce0(3, 3), Ccommit(3, 3), fy0(fy_0), fy0_inf(fy_infty), E0(e)
{
	eps.Zero();
	sig.Zero();
	sige.Zero();
	eps_p.Zero();
	sigeP.Zero();
	Ce.Zero();


						   //  some useful constants
						   //  shear modulus
	double G = E / 2.0 / (1.0 + nu);

	// initial tangent
	Ce0(0, 0) = E / (1.0 - nu * nu);
	Ce0(0, 1) = nu * E / (1.0 - nu * nu);
	Ce0(0, 2) = 0.0;
	Ce0(1, 0) = nu * E / (1.0 - nu * nu);
	Ce0(1, 1) = E / (1.0 - nu * nu);
	Ce0(1, 2) = 0.0;
	Ce0(2, 0) = 0.0;
	Ce0(2, 1) = 0.0;
	Ce0(2, 2) = G;

	Ce = Ce0;
	C = Ce0;

	xi_p = 0;

	TempAndElong(0) = 0;
	TempAndElong(1) = 0;

	Cchange = -1;
	Tchange = 0;
	Temp = 0;
	TempT = 0;

	this->commitState();
}

J2PlaneStressThermal::J2PlaneStressThermal()
	:NDMaterial(0, 0),
	eps(3), sig(3), sige(3), eps_p(3), sigeP(3), TempAndElong(2),
	epsCommit(3), sigCommit(3), sigeCommit(3), eps_pCommit(3), sigePCommit(3),
	Ce(3, 3), C(3, 3), Ccommit(3, 3)
{

}

J2PlaneStressThermal::~J2PlaneStressThermal()
{

}



void
J2PlaneStressThermal::StrsPr(const Vector&sig, Matrix &Peig, Vector &SigPr)
{
	double sigpr1, sigpr2, fi;
	//  rotation angle: principal effective stress directions
	if (fabs(sig(0) - sig(1)) > 1.0E-14) {
		fi = 0.5 * atan(2.0 * sig(2) / (sig(0) - sig(1)));
	}
	else if (fabs(sig(2)) > 1.0E-14) {
		fi = 0.78539816339744828;  // fi = pi/4.0
	}
	else {
		fi = 0.0;
	}

	if (fabs(fi) < 1.0E-14) {
		fi = 0.0;
	}

	//  principal effective stresses, sigp_effective
	double c = cos(fi);
	double s = sin(fi);

	sigpr1 = (sig(0) * (c * c) + sig(1) * (s * s)) + 2.0 *
		sig(2) * c * s;
	sigpr2 = (sig(0) * (s * s) + sig(1) * (c * c)) - 2.0 *
		sig(2) * c * s;

	if (sigpr1 < sigpr2) {
		SigPr(0) = sigpr2;
		SigPr(1) = sigpr1;
		Peig(0, 0) = sin(fi);
		Peig(0, 1) = cos(fi);
		Peig(1, 0) = -cos(fi);
		Peig(1, 1) = sin(fi);

	}
	else {
		SigPr(0) = sigpr1;
		SigPr(1) = sigpr2;
		Peig(0, 0) = cos(fi);
		Peig(0, 1) = sin(fi);
		Peig(1, 0) = sin(fi);
		Peig(1, 1) = -cos(fi);

	}


	//------end  of principal stressess


	//strsMises = sqrt(sigpr1*sigPr1 + sigpr2*sigpr2 - sigpr1*sigpr2);
}

int
J2PlaneStressThermal::setTrialStrain(Vector const&v1, Vector const&v2) {
	return this->setTrialStrain(v1);
}

int
J2PlaneStressThermal::setTrialStrain(const Vector &strain)
{
	double tol = 1.0e-6;
	double xi = 0;

	double Fres = 0;


	static Vector Deps(3);// increment of strain;
	static Vector Deps_p(3);// increment of plastic strain;
	static Vector Eeps_tr(3);

	static Vector sige_tr(3);

	static Matrix Cbar(3, 3);

	Vector sigPr(2);
	Vector delta_epsPr(2);

	//  tolerance for function sign evaluation

	eps_p = eps_pCommit;
	sigeP = sigePCommit;

	// current strain
	eps = strain;

	// incremental strain
	Eeps_tr = eps - eps_p;
	Deps = eps - epsCommit;


	Ce = E / E0*Ce0;

	// trial effective stress
	sige_tr = Ce*Eeps_tr;
	//sige_tr = sigeP+ Ce*Deps;

#ifdef _sDEBUG
	if (this->getTag() == 109) {
		opserr << this->getTag() << "--kt--: " << kt << "--kc--: " << kc << endln;
		opserr << "---eps:" << eps << "--epsp:" << eps_p << endln;
	}

#endif

	//Spectral transforamtion;
	Vector Sigtr(3); Vector SigtrPr(2); Matrix Peig(2, 2);
	double fyield; double sigpmp;
	double beta;  double lamda;
	double I1tr, sigmaVM;

	//int findex =0;  //index of yielding criteria: 0: elastic; 1: Prinicipal sigma; 2: Drucker Prager

	this->StrsPr(sige_tr, Peig, SigtrPr);


	//determine beta (considering updated strength)

	//Yield surface
	I1tr = SigtrPr(0) + SigtrPr(1);

	sigmaVM = sqrt(SigtrPr(0)*SigtrPr(0) + SigtrPr(1)*SigtrPr(1) - SigtrPr(0)*SigtrPr(1));

	fyield = sigmaVM - fy;


	if (fyield > tol) {
		// opserr << (alpha*I1tr + sigmaVM + beta*sigpmp) / (1 - alpha);
		//plastic


		double nuwave = (1 - 2 * nu) / (1 - nu);
		double Kwave = E / 2 / (1 - nu);
		double G = E / 2.0 / (1.0 + nu);

		double Snorm, I1;
		double sigma1 = SigtrPr(0);
		double sigma2 = SigtrPr(1);

		double zeta1;
		double zeta2;

		double deltaxi;
		double Xres;
		double xia;
		double xib;
		double xi_a;
		double xi_b;

		double ft_r, fc_r;

		double Qnorm;
		double deltakt, deltakc;
		double dQdkc, dQdkt;

		double Qrest, Qresc;
		double ktn1;
		double kcn1;
		double deltaktn1, deltakcn1;


		zeta1 = 0;
		zeta2 =  - 2.0 * nuwave* fy;

		I1 = sigma1 + sigma2;
		Snorm = sqrt(2.0 / 3.0 * (sigma1*sigma1 + sigma2*sigma2 - sigma1*sigma2));
		
		xia = (2.0*nuwave - 3.0)*fy / zeta2;
		

		xib = 1.0 - 3.0 * sqrt(6.0)*G / ( 2.0 * sqrt(6.0)*nuwave*G);

		if (xia < xib) {
			if (xia < 0)
				xi_a = 0;
			else
				xi_a = xia;

			if (xib < 1)
				xi_b = xib;
			else
				xi_b = 1;
		}
		else {
			if (xib < 0)
				xi_a = 0;
			else
				xi_a = xib;

			if (xia < 1)
				xi_b = xia;
			else
				xi_b = 1;
		}

		//iteration to determine lamda
		//-----------------------------------------------------------------------
		//variables
		double newbeta;
		for (int i = 0; i < 40; i++) {

			//residual of X = xi*(s+2*G*lamda)-s;

			if (xi < tol) {
				xi = (xi_a + xi_b) / 2;
			}


			double  lamdap1 =  zeta2*xi + (2.0 * nuwave - 3.0)*fy;
			double  lamdap2 =  - sqrt(6.0)*G*(3.0 - 2.0 * nuwave + 2.0 * nuwave*xi);
			lamda = (1.0 - xi) / xi*lamdap1 / lamdap2;
#ifdef _SDEBUG
			if (lamda < 0) {
				opserr << this->getTag() << " meets an nagetive lamda: " << lamda << "sigPr " << SigtrPr << endln;
				return(-1);
			}
#endif

			double deltasig;

			deltasig = xi*(nuwave*(1 - xi)*I1tr ) / (3.0 - 2.0 * nuwave*(1 - xi));

			sigma1 = SigtrPr(0)*xi + deltasig;
			sigma2 = SigtrPr(1)*xi + deltasig;


			I1 = sigma1 + sigma2;
			Snorm = sqrt(2.0 / 3.0 * (sigma1*sigma1 + sigma2*sigma2 - sigma1*sigma2));
			Xres = xi*(Snorm + 2.0 * G*lamda) - Snorm;


			//returning to different zone


			Fres = sqrt(3.0 / 2.0)*Snorm  - fy;


			if ((fabs(Xres) / Snorm < tol) && (fabs(Fres) / fy < tol) && lamda > 0) {
					break;
			}


			//if (i > 15) {

			//xi = xi_a ;
			//}

			double dlamda1_dxi =  zeta2;
			double dlamda2_dxi = - sqrt(6.0)*G * 2.0 * nuwave;

			double dlamdadxi = (-1.0 / xi / xi)*lamdap1 / lamdap2 + (1.0 / xi - 1)*(dlamda1_dxi*lamdap2 - dlamda2_dxi*lamdap1) / lamdap2 / lamdap2;

			double dsigdxi = (nuwave*(1.0 - 2.0*xi)*I1tr )*(3.0 - 2.0 * nuwave*(1.0 - xi)) - 2.0*nuwave*(nuwave*(1.0 - xi)*xi*I1tr);
			dsigdxi = dsigdxi / (3.0 - 2.0 * nuwave*(1.0 - xi)) / (3.0 - 2.0 * nuwave*(1 - xi));

			double dsigma1dxi = SigtrPr(0) + dsigdxi;
			double dsigma2dxi = SigtrPr(1) + dsigdxi;
			double dsnormdxi = ((2.0 * sigma1 - sigma2)*dsigma1dxi + (2.0 * sigma2 - sigma1)*dsigma2dxi) / 3.0 / Snorm;

			double  dXdxi = Snorm + 2.0 * G*lamda + (xi - 1.0)*dsnormdxi + 2.0*G*xi*dlamdadxi;

			deltaxi = -Xres / dXdxi;
			xi = xi + deltaxi;


#ifdef _SDEBUG
			opserr << "sigma1: " << sigma1 << "  sigma2: " << sigma2 << endln;
			opserr << "  Fres: " << Fres << " sigma1Res: " << sigma1Res << endln;
#endif
			//if (this->getTag()==177) {

			//opserr<< "Material: " <<this->getTag()<<" ," << SigtrPr(0)<< ", "<< SigtrPr(1)<< "..,sigPr:"<<sigma1<<sigma2<< ", ft: "<<ftbar<<", fc: "<<fcbar<<endln;
			//}

		}
		//end of determining lamda, principal stresses;


		//now determine damage variables
		//-----------------------------------------------------------------------

		//double gt, gc;  //failure energy
		sigPr(0) = sigma1; sigPr(1) = sigma2;

		//calculate increment of eps_principle
		if (lamda>0)
			delta_epsPr = lamda*(sigPr / Snorm - (sigPr(0) + sigPr(1)) / 3 / Snorm);

		

		//obtain effective stress
		Matrix PeigT(2, 2);

		Matrix tempEpspPr(2, 2);

		Matrix Depspn1(2, 2);


		tempEpspPr.Zero();

		tempEpspPr(0, 0) = delta_epsPr(0); tempEpspPr(1, 1) = delta_epsPr(1);


		PeigT.addMatrixTranspose(0, Peig, 1);


		Depspn1 = Peig*tempEpspPr*PeigT;


		Deps_p(0) = Depspn1(0, 0); Deps_p(1) = Depspn1(1, 1); Deps_p(2) = Depspn1(0, 1);

#ifdef _SDEBUG
		if (this->getTag() == 62) {
			opserr << "lamda" << lamda << "  sigmaPr1:" << sigPr << endln;
			opserr << "epsp:" << delta_epsPr << endln;
			opserr << "sige:" << sige << endln;

		}
#endif
		//update plastic strain;
		eps_p = eps_pCommit + Deps_p;

		//  double s0 = 0.01;
		//double s_sigma = s0 + (1 - s0)* gamma_sigma;
		//dt =  s_sigma*dt;

		//obtain incremental plastic strain

#ifdef _SDEBUG
		if (this->getTag() == 282) {
			opserr << this->getTag() << "--kt--: " << kt << "--kc--: " << kc << endln;
			opserr << "---eps:" << eps << "--epsp:" << eps_pCommit << " Eeps_tr" << Eeps_tr << endln;
			opserr << "  Sig: " << sig << " DeltaPeps " << delta_epsPr << endln;
		}


#endif

		//opserr << "Deps_p" << Deps_p;
		Cbar = Ce;

	}
	else {
		//elastic

		sige = sige_tr;
		Cbar = Ce;

		sigPr = SigtrPr;

	}

	// dtotal =1-(1 - dc)*(1 - dt);   //not adding cyclic effect yet
	Matrix PeigT(2, 2);
	Matrix tempSigPr(2, 2);
	Matrix Sign1(2, 2);


	tempSigPr.Zero();


	tempSigPr(0, 0) = sigPr(0); tempSigPr(1, 1) = sigPr(1);

	PeigT.addMatrixTranspose(0, Peig, 1);

	Sign1 = Peig*tempSigPr*PeigT;

	sige(0) = Sign1(0, 0); sige(1) = Sign1(1, 1); sige(2) = Sign1(0, 1);

	sig = sige;

	//adding remaini
	//if (dtotal > 0.99)
	//dtotal = 0.99;

	C = Cbar;

#ifdef _SDEBUG
	if (this->getTag() == 109) {
		opserr << this->getTag() << "--kt--: " << kt << "--kc--: " << kc << endln;
		opserr << "---eps:" << eps << "--epsp:" << eps_p << endln;
		opserr << "  Sig: " << sig << " DeltaPeps " << delta_epsPr << endln;
	}
#endif

	
	///////////////////////////////////////////////output for debug////////////////////////////// 
	//#endif

	return 0;
}

/*
//hardening function
	double J2PlaneStressThermal::q(double xi)
	{
		//  q(xi) = simga_infty + (sigma_0 - sigma_infty)*exp(-delta*xi) + H*xi 

		return    sigma_infty
			+ (sigma_0 - sigma_infty)*exp(-delta*xi)
			+ Hard*xi;
	}


//hardening function derivative
	double J2PlaneStressThermal::qprime(double xi)
	{
		return  (sigma_0 - sigma_infty) * (-delta) * exp(-delta*xi)
			+ Hard;
	}
	*/
int
J2PlaneStressThermal::setTrialStrainIncr(const Vector &dStrain)
{
	eps += dStrain;
	this->setTrialStrain(eps);
	return 0;
}

int
J2PlaneStressThermal::setTrialStrainIncr(const Vector &dStrain, const Vector &rate)
{
	eps += dStrain;
	this->setTrialStrain(eps);
	return 0;
}

const Matrix&
J2PlaneStressThermal::getTangent(void)
{
	return C;
}

const Matrix&
J2PlaneStressThermal::getInitialTangent(void)
{
	return Ce0;
}

const Vector&
J2PlaneStressThermal::getStress(void)
{
	return sig;
}

const Vector&
J2PlaneStressThermal::getStrain(void)
{
	return eps;
}

int
J2PlaneStressThermal::commitState(void)
{

	ECommit = E;
	epsCommit = eps;
	sigCommit = sig;
	sigeCommit = sige;
	eps_pCommit = eps_p;
	sigePCommit = sige;
	if (Cchange >0)
		Cchange = 0;
	return 0;
}

int
J2PlaneStressThermal::revertToLastCommit(void)
{
	C = Ccommit;
	E = ECommit;

	eps = epsCommit;
	sig = sigCommit;
	sige = sigeCommit;
	eps_p = eps_pCommit;
	sigeP = sigePCommit;
	return 0;
}

int
J2PlaneStressThermal::revertToStart(void)
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
J2PlaneStressThermal::getCopy(const char *type)
{
	if (strcmp(type, "PlaneStress") == 0 || strcmp(type, "PlaneStress2D") == 0) {
		//for debugging
		matTag++;
		J2PlaneStressThermal *theCopy =
			new J2PlaneStressThermal(matTag, E0, nu, fy, fy,d,H);
		return theCopy;
	}
	else {
		return 0;
	}
}

NDMaterial*
J2PlaneStressThermal::getCopy(void)
{
	J2PlaneStressThermal *theCopy =
		new J2PlaneStressThermal(this->getTag(), E0, nu, fy0, fy0,d,H);
	return theCopy;
}

const char*
J2PlaneStressThermal::getType(void) const
{
	return "PlaneStress2d";
}

int
J2PlaneStressThermal::getOrder(void) const
{
	return 6;
}


int
J2PlaneStressThermal::sendSelf(int commitTag, Channel &theChannel)
{
	static Vector data(10);

	int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "J2PlaneStressThermal::sendSelf -- could not send Vector\n";
		return res;
	}

	return res;
}

int
J2PlaneStressThermal::recvSelf(int commitTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	static Vector data(10);

	int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "J2PlaneStressThermal::sendSelf -- could not send Vector\n";
		return res;
	}

	return res;
}

void
J2PlaneStressThermal::Print(OPS_Stream &s, int flag) {
	//opserr << "J2PlaneStressThermal: " << this->getTag();
	//opserr << "strain: " << eps;
	//opserr << "strain: " << sig;
	//opserr << "tangent: " << this->getTangent();
}

//send back TempAndElong(Liming,UoE)
const Vector&
J2PlaneStressThermal::getTempAndElong(void)
{
	return TempAndElong;
}

//Set TemperatureAndElongation

double
J2PlaneStressThermal::setThermalTangentAndElongation(double &tempT, double&ET, double&Elong)
{

	double TempT = tempT + 20;
	double E00; //Initial tangent 
	ET = 2E11;
	E00 = 2E11;

	// EN 1992 pt 1-2-1. Class N hot rolled  reinforcing steel at elevated temperatures

	if (TempT <= 100) {

	}
	else if (TempT <= 200) {

		E = E0*(1 - (TempT - 100)*0.1 / 100);
		
		fy = fy0;

		H = 0.01*E / 2.8;

	}
	else if (TempT <= 300) {

		E = E0*(0.9 - (TempT - 200)*0.1 / 100);

		fy = fy0;
		
		H = 0.01*E / 2.8;

	}
	else if (TempT <= 400) {
		E = E0*(0.8 - (TempT - 300)*0.1 / 100);
		
		fy = fy0;

		H = 0.01*E / 2.8;
	}
	else if (TempT <= 500) {
		E = E0*(0.7 - (TempT - 400)*0.1 / 100);


		fy = fy0*(1 - (TempT - 400)*0.22 / 100);

		H = 0.01*E / 2.8;
	}
	else if (TempT <= 600) {
		E = E0*(0.6 - (TempT - 500)*0.29 / 100);

		fy = fy0*(0.78 - (TempT - 500)*0.31 / 100);

		H = 0.01*E / 2.8;
	}
	else if (TempT <= 700) {
		E = E0*(0.31 - (TempT - 600)*0.18 / 100);


		fy = fy0*(0.47 - (TempT - 600)*0.24 / 100);

		H = 0.01*E / 2.8;
	}
	else if (TempT <= 800) {
		E = E0*(0.13 - (TempT - 700)*0.04 / 100);


		fy = fy0*(0.23 - (TempT - 700)*0.12 / 100);

		H = 0.01*E / 2.8;
	}
	else if (TempT <= 900) {
		E = E0*(0.09 - (TempT - 800)*0.02 / 100);


		fy = fy0*(0.11 - (TempT - 800)*0.05 / 100);

		H = 0.01*E / 2.8;
	}
	else if (TempT <= 1000) {
		E = E0*(0.0675 - (TempT - 900)*(0.00675 - 0.0045) / 100);


		fy = fy0*(0.06 - (TempT - 900)*0.02 / 100);

		H = 0.01*E / 2.8;
	}

	else {
		opserr << "the temperature is invalid\n";
	}
	double ThermalElongation = 0;
	// Calculate thermal elongation 
	if (TempT <= 20) {
		ThermalElongation = 0.0;
	}
	else if (TempT <= 750) {
		ThermalElongation = -2.416e-4 + 1.2e-5 *TempT + 0.4e-8 *TempT*TempT;

	}
	else if (TempT <= 860) {
		ThermalElongation = 11e-3;

	}
	else if (TempT <= 1200) {
		ThermalElongation = -6.2e-3 + 2e-5*TempT;

	}
	else {
		opserr << "the temperature is invalid\n";
	}

	//ThermalElongation = 0;
	TempAndElong(0) = Temp;
	TempAndElong(1) = ThermalElongation;
	//bulk = bulk_0;
	//shear = shear_0;
	//sigma_y = sigma_0;

	//ET = 2E11;
	//ET = E;  
	//ET = 3.84e10;
	Elong = ThermalElongation;
	//this->plastic_integrator();
	Tchange++;
	return 0;
}
