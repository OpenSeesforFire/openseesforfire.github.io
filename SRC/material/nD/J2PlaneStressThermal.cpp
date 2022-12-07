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

	if (numArgs < 5) {
		opserr << "Want: nDMaterial J2PlaneStressThermal $tag $typeTag $E $nu $fy $fyinf \n";
		return 0;
	}

	int iData[2];
	double dData[10];
	dData[4] = 0;
	dData[5] = 0;
	


	int numData = 2;
	if (OPS_GetInt(&numData, iData) != 0) {
		opserr << "WARNING invalid integer tag: nDMaterial J2PlaneStressThermal \n";
		return 0;
	}

	numData = numArgs - 2;;
	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "WARNING invalid data: nDMaterial J2PlaneStressThermal : " << iData[0] << "\n";
		return 0;
	}

	theMaterial = new J2PlaneStressThermal(iData[0], iData[1],
		dData[0], dData[1], dData[2], dData[3],dData[4], dData[5]);
	matTag = iData[0];
	return theMaterial;
}

J2PlaneStressThermal::J2PlaneStressThermal(int tag, int typeTag,
	double e,
	double nu,
	double fy_0,
	double fy_infty,
	double d0,
	double H0
)
	:NDMaterial(tag, ND_TAG_J2PlaneStressThermal),
	TypeTag(typeTag),E(e), nu(nu), fy(fy_0), fy_inf(fy_infty),
	eps(3), sig(3), sige(3), eps_p(3), sigeP(3), TempAndElong(4),
	epsCommit(3), sigCommit(3), sigeCommit(3), eps_pCommit(3),d(d0),H(H0),HT(H0),
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
	TempAndElong(2) = 0;
	TempAndElong(3) = 0;

	Cchange = -1;
	Tchange = 0;
	Temp = 0;
	TempT = 0;
	fyt = fy_0;

	kxi = 0;
	kxi_Commit = 0;

	this->commitState();
}

J2PlaneStressThermal::J2PlaneStressThermal()
	:NDMaterial(0, 0),
	eps(3), sig(3), sige(3), eps_p(3), sigeP(3), TempAndElong(4),
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
	double tol = 1.0e-8;
	double xi = 0;
	double root23 = sqrt(2.0 / 3.0);

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


	//Spectral transforamtion;
	Vector Sigtr(3); Vector SigtrPr(2); Matrix Peig(2, 2);
	double fyield =0; double sigpmp;
	double beta;  double lamda =0.0;
	double I1tr, sigmaVM;

	//int findex =0;  //index of yielding criteria: 0: elastic; 1: Prinicipal sigma; 2: Drucker Prager

	this->StrsPr(sige_tr, Peig, SigtrPr);
	//opserr <<  "  sigmaPr1:" << SigtrPr << endln;

	//determine beta (considering updated strength)

	//Yield surface
	I1tr = SigtrPr(0) + SigtrPr(1);

	sigmaVM = sqrt(SigtrPr(0)*SigtrPr(0) + SigtrPr(1)*SigtrPr(1) - SigtrPr(0)*SigtrPr(1));
	fyt = fy_inf + (fy - fy_inf) * exp(-d * kxi) + HT * kxi;

	fyield = sigmaVM - fyt;

	//sigmaV = sqrt(3*J2)
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
	



		double Qnorm;

	
		I1 = sigma1 + sigma2;
		Snorm = sqrt(2.0 / 3.0 * (sigma1*sigma1 + sigma2*sigma2 - sigma1*sigma2));
		double phi = Snorm - root23 * fyt;
		double resid = 1.0;
		double tang = 0;
		double iteration_counter = 0;
		double kxi_trial = 0;
		double f_trial = 0;
		double f_trial_p = 0;
		lamda = (Snorm - root23 * fyt) / (2.0 * G);
		while (fabs(resid) > 1e-8) {

			kxi_trial = kxi + root23 * lamda;
			f_trial = fy_inf + (fy - fy_inf) * exp(-d * kxi_trial) + HT * kxi_trial;
			f_trial_p= -d * (fy - fy_inf) * exp(-d * kxi) + HT;
			resid = Snorm
				- (2.0 * G) * lamda
				- root23 * f_trial;
				//- (eta / dt) * gamma;


			tang = -(2.0 * G)
				- 2.0 / 3.0 * f_trial_p;
				//- (eta / dt);

			lamda -= (resid / tang);

			iteration_counter++;

			if (iteration_counter > 100) {
				//opserr << "More than 100" ;
				//opserr << " iterations in constituive subroutine J2-plasticity \n";
				break;
			} //end if 

		} //end while resid


		//iteration to determine lamda
		//-----------------------------------------------------------------------
		//variables
		
		kxi = kxi_trial;
		fyt = fy_inf+ (fy - fy_inf)*exp(-d*kxi)+ HT*kxi;
		Snorm = root23* fyt;
		//now determine damage variables
		//-----------------------------------------------------------------------
		double devS1 = (sigma1 - I1tr / 3.0) / (1 + 2 * G * lamda / Snorm);
		double devS2 = (sigma2 - I1tr / 3.0) / (1 + 2 * G * lamda / Snorm);
		
		sigma1 = 2 * devS1 + devS2;
		sigma2 = devS1 + 2*devS2;
	
		
		sigPr(0) = sigma1; sigPr(1) = sigma2;

		//calculate increment of eps_principle
		if (lamda >0)
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
	double J2PlaneStressThermal::q(double kxi)
	{
		//  q(xi) = simga_infty + (sigma_0 - sigma_infty)*exp(-delta*xi) + H*xi 

		return    
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
	kxi_Commit = kxi;
	if (Cchange >0)
		Cchange = 0;
	return 0;
}

int
J2PlaneStressThermal::revertToLastCommit(void)
{
	C = Ccommit;
	E = ECommit;
	kxi = kxi_Commit;
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
			new J2PlaneStressThermal(matTag, TypeTag, E0, nu, fy0, fy0_inf,d,H);
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
		new J2PlaneStressThermal(this->getTag(), TypeTag, E0, nu, fy0, fy0_inf,d,H);
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
	//TempAndElong(0) = kxi;
	//TempAndElong(1) = fyt;
	return TempAndElong;
}

//Set TemperatureAndElongation

double
J2PlaneStressThermal::setThermalTangentAndElongation(double &tempT, double&ET, double&Elong)
{
	
	double TempT = tempT;
	double E00; //Initial tangent 
	ET = E0;
	
	if (TempT < -100 || TempT>1200)
		opserr << "WARNING J2PlaneStressThermal received an invalid Temperature: " << TempT << endln;
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
		for (int i = 0; i<12; i++) {
			FyRfactors[i] = FyRfEC3[i];
			FpRfactors[i] = FpRfEC3[i];
			E0Rfactors[i] = E0RfEC3[i];
		}

	}
	else if (TypeTag == 21) {
		double FyRfEC21[12] = { 1.0, 1.0 ,1.0, 1.0 ,0.78, 0.47, 0.23, 0.11, 0.06, 0.04 ,0.02, 0.0 };
		double FpRfEC21[12] = { 1.0, 0.81 ,0.61, 0.42 ,0.36, 0.18, 0.07, 0.05, 0.04, 0.02 ,0.01, 0.0 };
		double E0RfEC21[12] = { 1.0, 0.9, 0.8 ,0.7, 0.6 ,0.31, 0.13, 0.09, 0.07, 0.04, 0.02 , 0.0 };
		for (int i = 0; i<12; i++) {
			FyRfactors[i] = FyRfEC21[i];
			FpRfactors[i] = FpRfEC21[i];
			E0Rfactors[i] = E0RfEC21[i];
		}

	}
	else if (TypeTag == 22) {
		double FyRfEC22[12] = { 1.0, 1.0 ,1.0, 0.94 ,0.67, 0.40, 0.12, 0.11, 0.08, 0.05 ,0.03, 0.0 };
		double FpRfEC22[12] = { 0.96 ,0.92, 0.81 ,0.63, 0.44, 0.26, 0.08, 0.06, 0.05 ,0.03, 0.02, 0.0 };
		double E0RfEC22[12] = { 1.0, 0.87, 0.72 ,0.56, 0.40 ,0.24, 0.08, 0.06, 0.05, 0.03, 0.02 , 0.0 };
		for (int i = 0; i<12; i++) {
			FyRfactors[i] = FyRfEC22[i];
			FpRfactors[i] = FpRfEC22[i];
			E0Rfactors[i] = E0RfEC22[i];
		}
	}
	else if (TypeTag == 23) {
		double FyRfEC23[12] = { 1.0, 1.0 ,1.0, 0.90 ,0.70, 0.47, 0.23, 0.11, 0.06, 0.04 ,0.02, 0.0 };
		double FpRfEC23[12] = { 1.00 ,0.87, 0.74 ,0.70, 0.51, 0.18, 0.07, 0.05, 0.04 ,0.02, 0.01, 0.0 };
		double E0RfEC23[12] = { 1.0, 0.95, 0.90 ,0.75, 0.60 ,0.31, 0.13, 0.09, 0.07, 0.04, 0.02 , 0.0 };
		for (int i = 0; i<12; i++) {
			FyRfactors[i] = FyRfEC23[i];
			FpRfactors[i] = FpRfEC23[i];
			E0Rfactors[i] = E0RfEC23[i];
		}
	}
	else
		opserr << "WARNING J2PlaneStressThermal received an invalid typeTag: " << TypeTag << endln;

	//Now Updating modulus, strengths
	for (int i = 0; i<13; i++) {
		if (TempT <= 80 + 100 * i)
		{
			if (i == 0) {
				fy = fy0*(1.0 - TempT*(1.0 - FyRfactors[0]) / 80);
				E = E0*(1.0 - TempT*(1.0 - E0Rfactors[0]) / 80);
			}
			else if (i == 12) {
				opserr << "Warning:The temperature " << TempT << " for J2PlaneStressThermal is out of range\n";
				return -1;
			}
			else {
				fy = fy0*(FyRfactors[i - 1] - (TempT + 20 - 100 * i)*(FyRfactors[i - 1] - FyRfactors[i]) / 100);
				E = E0*(E0Rfactors[i - 1] - (TempT + 20 - 100 * i)*(E0Rfactors[i - 1] - E0Rfactors[i]) / 100);
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
		ThermalElongation = -2.416e-4 + 1.2e-5 *(TempT + 20) + 0.4e-8 *(TempT + 20)*(TempT + 20);
	}
	else if (TempT <= 840) {
		ThermalElongation = 11e-3;
	}
	else if (TempT <= 1180) {
		ThermalElongation = -6.2e-3 + 2e-5*(TempT + 20);
	}
	else {
		opserr << " J2PlaneStressThermal Temperature " << TempT << " is invalid\n";
		return -1;
	}

	
	fy_inf = fy / fy0*fy0_inf;
	//H = fy / fy0 * 200e6;
	HT = H * fy / fy0;
	fyt = fy_inf -(fy_inf-fy)*exp(-d*kxi_Commit) + HT*kxi_Commit;
	
	

	//ThermalElongation = ThermalElongation;
	//(TempT -20)*1.0e-5;
	//TempAndElong(0) = tempT;
	//TempAndElong(1) = fy;
	//bulk = bulk_0;
	//shear = shear_0;
	//sigma_y = sigma_0;

	//ET = 2E11;
	//ET = E;  
	//ET = 3.84e10;
	//ThermalElongation = 0;
	Elong = ThermalElongation;
	//this->plastic_integrator();
    TempAndElong(0) = TempT;
	TempAndElong(1) = fyt;
	//TempAndElong(1) = ThermalElongation;
	Tchange++;
	return 0;
}
