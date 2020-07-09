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

//Based on the work presented in Lee, J., & Fenves, G. L. (2001). Return-mapping algorithm for plastic-damage models: 3-D and plane stress formulation. 
//International Journal for Numerical Methods in Engineering, 50(2), 487-C506. 
// Modified for SIF modelling by Liming Jiang [http://openseesforfire.github.io] 

#include <PlasticDamageConcretePlaneStressThermal.h>           
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
OPS_PlasticDamageConcretePlaneStressThermal(void)
{
	NDMaterial *theMaterial = 0;

	int numArgs = OPS_GetNumRemainingInputArgs();

	if (numArgs < 5 || numArgs > 11) {
		opserr << "Want: nDMaterial PlasticDamageConcretePlaneStressThermal $tag $E $nu $ft $fc <$gt $gc $An $Bn>\n";
		return 0;
	}

	int iData[1];
	double dData[10];
	dData[4] = 0;
	dData[5] = 0;
	dData[6] = 0.4;
	dData[7] = 4.0;
	dData[8] = 0.4;
	dData[9] = 0.1;
	

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

	if (dData[4] < 1E-6)
		dData[4] = dData[2] * dData[2] / dData[0] * 4.0;

	if (dData[5] < 1E-6)
		dData[5] = dData[3] * dData[3] / dData[0] * 6.0;

	theMaterial = new PlasticDamageConcretePlaneStressThermal(iData[0],
		dData[0], dData[1], dData[2], dData[3],
		dData[4], dData[5], dData[6], dData[7], dData[8], dData[9]);
	matTag = iData[0];
	return theMaterial;
}

PlasticDamageConcretePlaneStressThermal::PlasticDamageConcretePlaneStressThermal(int tag,
	double _e,
	double _nu,
	double _ft,
	double _fc,
	double Gt,
	double Gc,
	double _At,
	double _Ac,
	double _Dt,
	double _Dc
	)
	:NDMaterial(tag, 100),
	E(_e), nu(_nu), ft(_ft), At(_At), Ac(_Ac), Dbarc(_Dc), Dbart(_Dt),gt0(Gt),gc0(Gc),
	eps(3), sig(3), sige(3), eps_p(3), sigeP(3), TempAndElong(2),
	epsCommit(3), sigCommit(3), sigeCommit(3), eps_pCommit(3),
	Ce(3, 3), C(3, 3), Ce0(3, 3), Ccommit(3, 3), ft0(_ft), E0(_e)
{
	eps.Zero();
	sig.Zero();
	sige.Zero();
	eps_p.Zero();
	sigeP.Zero();
	Ce.Zero();

	fc = _fc / 1.562491022;  //for maximum stress
	fc0 = _fc / 1.562491022;  //for maximum stress

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


	kt = 0;
	kc = 0;
	ktCommit = 0;
	kcCommit = 0;

	gt = gt0;
	gc = gc0;
	dc = 0;
	dt = 0;
	dcCommit = 0;
	dtCommit = 0;

	TempAndElong(0) = 0;
	TempAndElong(1) = 0;

	Cchange = -1;
	Tchange = 0;
	Temp = 0;
	TempT = 0;
	epsLitsp = 0;

	this->commitState();
}

PlasticDamageConcretePlaneStressThermal::PlasticDamageConcretePlaneStressThermal()
	:NDMaterial(0, 0),
	eps(3), sig(3), sige(3), eps_p(3), sigeP(3), TempAndElong(2),
	epsCommit(3), sigCommit(3), sigeCommit(3), eps_pCommit(3), sigePCommit(3),
	Ce(3, 3), C(3, 3), Ccommit(3, 3), At(0), Ac(0), Dbarc(0), Dbart(0), epsLitsp(0)
{

}

PlasticDamageConcretePlaneStressThermal::~PlasticDamageConcretePlaneStressThermal()
{

}



void
PlasticDamageConcretePlaneStressThermal::StrsPr(const Vector&sig, Matrix &Peig, Vector &SigPr)
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
PlasticDamageConcretePlaneStressThermal::setTrialStrain(Vector const&v1, Vector const&v2) {
	return this->setTrialStrain(v1);
}

int
PlasticDamageConcretePlaneStressThermal::setTrialStrain(const Vector &strain)
{
	//double tempT= 580;
	//double ET;
	//double Elong;
	//this->setThermalTangentAndElongation(tempT, ET, Elong);


	double tol = 1.0e-6;
	double xi = 0;
	double f2c = 1.16 * fc;
	double Fres=0;
	double Qnorm = 0;
	

	static Vector Deps(3);// increment of strain;
	static Vector Deps_p(3);// increment of plastic strain;
	static Vector Eeps_tr(3);

	static Vector sige_tr(3);

	static Matrix Cbar(3, 3);

	Vector sigPr(2);
	Vector delta_epsPr(2);

	//  tolerance for function sign evaluation

//	if (this->getTag() == 109 && TempAndElong(0) > 286)
		//opserr << endln;

	double ftbar; double fcbar;  //effective uniaxial strengths
	double dtbt, dcbc;    //constants
	double fit, fic;

	//f2c = 1.16*fc;
	double alpha = (1.16 - 1) / (2 * 1.16 - 1);
	double alphap = 0.3;


	double dtotal; // overall consideration of tensile and compressive damage


	kt = ktCommit;
	kc = kcCommit;
	dt = dtCommit;
	dc = dcCommit;

	double gamma_sigma;

	//effective ft & fc;
	dcbc = log(1 - Dbarc) / log((1 + Ac) / 2 / Ac);
	dtbt = log(1 - Dbart) / (log((1 + At) - sqrt(1 + At*At)) - log(2 * At));


	fit = 1 + At*(2 + At)*kt;
	fic = 1 + Ac*(2 + Ac)*kc;

	//fcbar = fc; 
	ftbar = ft;
	//ftbar = ft*pow((1.0 / At*(1 + At - sqrt(fit))), (1.0 - dtbt))*sqrt(fit);
	fcbar = fc*pow((1.0 / Ac*(1 + Ac - sqrt(fic))), (1.0 - dcbc))*sqrt(fic);


	eps_p = eps_pCommit;
	sigeP = sigePCommit;

	// current strain
	eps = strain;

	// incremental strain
	Eeps_tr = eps - eps_p;
	Deps = eps - epsCommit;

	//if (Tchange==1|| Tchange == 2) {
		//Ce = ECommit / E0*Ce0;
		//Eeps_tr = epsCommit - eps_pCommit;
		//Tchange++;
	//}
	//else
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

	if (SigtrPr(0) > 0) {
		sigpmp = SigtrPr(0);
		beta = fcbar / ftbar*(1 - alpha) - (1 + alpha);
	}
	else {
		sigpmp = 0;
		beta = 0;
	}
	//determine beta (considering updated strength)

	//Yield surface
	I1tr = SigtrPr(0) + SigtrPr(1);

	sigmaVM = sqrt(SigtrPr(0)*SigtrPr(0) + SigtrPr(1)*SigtrPr(1) - SigtrPr(0)*SigtrPr(1));

	fyield = (alpha*I1tr + sigmaVM + beta*sigpmp) / (1 - alpha) - fcbar;


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

		double deltakt, deltakc;
		double dQdkc, dQdkt;

		double Qrest, Qresc;
		double ktn1;
		double kcn1;
		double deltaktn1, deltakcn1;


		ktn1 = kt;
		kcn1 = kc;


		
			zeta1 = beta*nuwave*(SigtrPr(0) - SigtrPr(1));
			zeta2 = (3.0 * alpha + beta*nuwave)*I1tr + (3.0 - 2.0 * nuwave)*beta*SigtrPr(0) - 2.0 * nuwave*(1.0 - alpha)* fcbar;

			I1 = sigma1 + sigma2;
			Snorm = sqrt(2.0 / 3.0 * (sigma1*sigma1 + sigma2*sigma2 - sigma1*sigma2));
			if (zeta1 < 1e-8)
				xia = (2.0*nuwave - 3.0)*(1.0 - alpha)*fcbar / zeta2;
			else
				xia = (sqrt(zeta2*zeta2 - 4.0 * zeta1*(2.0 * nuwave - 3.0)*(1.0 - alpha)*fcbar) - zeta2) / 2.0 / zeta1;

			xib = 1.0 - 3.0 * sqrt(6.0)*G / (6.0 * (2.0 * alpha + beta)*Kwave*alphap + 2.0 * sqrt(6.0)*nuwave*G);

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



				if (i == 20) {
					//try different yielding surface
					if (beta < tol) {
						sigpmp = SigtrPr(0);
						beta = fcbar / ftbar*(1 - alpha) - (1 + alpha);
					}
					else
					{
						sigpmp = 0;
						beta = 0;
					}
					zeta1 = beta*nuwave*(SigtrPr(0) - SigtrPr(1));
					zeta2 = (3.0 * alpha + beta*nuwave)*I1tr + (3.0 - 2.0 * nuwave)*beta*SigtrPr(0) - 2.0 * nuwave*(1.0 - alpha)* fcbar;
					if (zeta1 < tol)
						xia = (2.0*nuwave - 3.0)*(1.0 - alpha)*fcbar / zeta2;
					else
						xia = (sqrt(zeta2*zeta2 - 4.0 * zeta1*(2.0 * nuwave - 3.0)*(1.0 - alpha)*fcbar) - zeta2) / 2.0 / zeta1;

					xib = 1.0 - 3.0 * sqrt(6.0)*G / (6.0 * (2.0 * alpha + beta)*Kwave*alphap + 2.0 * sqrt(6.0)*nuwave*G);

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

					xi = (xi_a + xi_b) / 2;

				}

				double  lamdap1 = zeta1*xi*xi + zeta2*xi + (2.0 * nuwave - 3.0)*(1.0 - alpha)*fcbar;
				double  lamdap2 = 6.0 * (2.0 * alpha + beta)*(1.0 - xi)*Kwave*alphap - sqrt(6.0)*G*(3.0 - 2.0 * nuwave + 2.0 * nuwave*xi);
				lamda = (1.0 - xi) / xi*lamdap1 / lamdap2;
#ifdef _SDEBUG
				if (lamda < 0) {
					opserr << this->getTag() << " meets an nagetive lamda: " << lamda << "sigPr " << SigtrPr << endln;
					return(-1);
				}
#endif

				double deltasig;

				deltasig = xi*(nuwave*(1 - xi)*I1tr - 6.0 * lamda*alphap*Kwave) / (3.0 - 2.0 * nuwave*(1 - xi));

				sigma1 = SigtrPr(0)*xi + deltasig;
				sigma2 = SigtrPr(1)*xi + deltasig;


				I1 = sigma1 + sigma2;
				Snorm = sqrt(2.0 / 3.0 * (sigma1*sigma1 + sigma2*sigma2 - sigma1*sigma2));
				Xres = xi*(Snorm + 2.0 * G*lamda) - Snorm;


				//returning to different zone


				Fres = (alpha*I1 + sqrt(3.0 / 2.0)*Snorm + beta*sigma1) / (1 - alpha) - fcbar;


				if ((fabs(Xres) / Snorm < tol) && (fabs(Fres) / fcbar < tol) && lamda > 0) {
					if (sigma1 > 0 && beta > 0)
						break;
					else if (sigma1 < 0 && beta == 0)
						break;
				}
					

				//if (i > 15) {

					//xi = xi_a ;
				//}

				double dlamda1_dxi = zeta1 * 2.0 * xi + zeta2;
				double dlamda2_dxi = -6.0 * (2.0 * alpha + beta)*Kwave*alphap - sqrt(6.0)*G * 2.0 * nuwave;

				double dlamdadxi = (-1.0 / xi / xi)*lamdap1 / lamdap2 + (1.0 / xi - 1)*(dlamda1_dxi*lamdap2 - dlamda2_dxi*lamdap1) / lamdap2 / lamdap2;

				double dsigdxi = (nuwave*(1.0 - 2.0*xi)*I1tr - 6.0 * lamda*alphap*Kwave - 6.0*alphap*Kwave*xi*dlamdadxi)*(3.0 - 2.0 * nuwave*(1.0 - xi)) - 2.0*nuwave*(nuwave*(1.0 - xi)*xi*I1tr - 6.0*lamda*alphap*Kwave*xi);
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
			//-//---//--end of determining lamda, principal stresses------///--//--//


			//now determine damage variables
			//-----------------------------------------------------------------------

			//double gt, gc;  //failure energy
			sigPr(0) = sigma1; sigPr(1) = sigma2;

			//calculate increment of eps_principle
			if(lamda>0)
				delta_epsPr = lamda*(sigPr / Snorm + alphap - (sigPr(0) + sigPr(1)) / 3 / Snorm);

			double hsigmaSum = 0; double sigmaSum = 0;

			for (int i = 0; i < 2; i++) {
				if (sigPr(i) > 0)
					hsigmaSum = hsigmaSum + sigPr(i);

				sigmaSum = sigmaSum + fabs(sigPr(i));
			}

			gamma_sigma = hsigmaSum / sigmaSum;


			for (int k = 0; k < 50; k++) {
			
				if (ktn1 < 0)
					ktn1 = 0;
				else if (ktn1 > 1)
					ktn1 = 1;
				if (kcn1 < 0)
					kcn1 = 0;
				else if (kcn1 > 1)
					kcn1 = 1;

				fit = 1 + At*(2 + At)*ktn1;
				fic = 1 + Ac*(2 + Ac)*kcn1;

				if (fit < 0)
					opserr << "fit " << fit << " kt " << ktn1;
				if (fic < 0)
					opserr << "fic " << fic << " kc " << kcn1;



				ft_r = ft/At*((1+At)*sqrt(fit)-fit);
				fc_r = fc / Ac*((1 + Ac)*sqrt(fic) - fic);

				//calculate incremental k (tension& compression)
				if (delta_epsPr(0) > 0)
					deltaktn1 = gamma_sigma*ft_r / gt*delta_epsPr(0);
				else
					deltaktn1 = 0;
				if (delta_epsPr(1) < 0)
					deltakcn1 = -(1 - gamma_sigma)*fc_r / gc*delta_epsPr(1);
				else
					deltakcn1 = 0;
				
				//double ddt_dkt_u = dtbt *pow((1 + At - sqrt(fit)) / At, (dtbt - 1));

				//double ddt_dkt = ddt_dkt_u*(1 / 2.0 / sqrt(fit)*At * (2 + At));

				double dfit_dkt = At*( 2+At);
				double dfic_dkc = Ac*(2 + Ac);
				dQdkt = -1 + gamma_sigma / gt*delta_epsPr(0)*ft/At*((1+At)/2/sqrt(fit)-1)*dfit_dkt;
				dQdkc = -1 - (1 - gamma_sigma) / gc*delta_epsPr(1)*fc / Ac*((1 + Ac) / 2 / sqrt(fic) - 1)*dfic_dkc;

				Qrest = -ktn1 + ktCommit + deltaktn1;
				Qresc = -kcn1 + kcCommit + deltakcn1;

				if (fabs(Qrest) > tol&& dQdkt != 0)
					deltakt = -Qrest / dQdkt;
				else
					deltakt = 0;

				if (fabs(Qresc) > tol&& dQdkc != 0)
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

				Qrest = -ktn1 + ktCommit + deltaktn1;
				Qresc = -kcn1 + kcCommit + deltakcn1;
				if (ktn1 != ktn1)
					opserr << ktn1;

				Qnorm = sqrt(Qrest*Qrest + Qresc*Qresc);

				if (Qnorm < 1e-10) {
					//opserr << "count: " << count;
					break;
				}
				//if(i>40)
				 // opserr<<this->getTag() << " plastic iteration  "<<i;

				//ftbar = ft*pow((1.0 / At*(1 + At - sqrt(fit))), (1.0 - dtbt))*sqrt(fit);
				fcbar = fc*pow((1.0 / Ac*(1 + Ac - sqrt(fic))), (1.0 - dcbc))*sqrt(fic);
		
			}
			//end of local iteration for i

			//determine the damage parameters----------//
			dt = 1.0 - 1.0/At*((1 + At)*sqrt(fit) - fit);
			dc = 1.0 - 1.0 / Ac*((1 + Ac)*sqrt(fic) - fic);

			if (dc != dc || dt != dt) {
				opserr << "invalid" << dt << "," << dc << endln;
				return -1;
			}
			if (sigPr(0) > 0) {
				sigpmp = sigPr(0);
				beta = fcbar / ftbar*(1 - alpha) - (1 + alpha);
			}
			else {
				sigpmp = 0;
				beta = 0;
			}

		//	if (TempAndElong(0) > 391 &&this->getTag()==2589)
				//opserr << "Temp " << TempAndElong(0)<< endln;

			Fres = (alpha*I1 + sqrt(3.0 / 2.0)*Snorm + beta*sigpmp) / (1 - alpha) - fcbar;

		if (ktn1 > kt)
			kt = ktn1;
		if (kcn1>kc)
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


		dt = dtCommit;
		dc = dcCommit;

		sige = sige_tr;
		Cbar = Ce;

		sigPr = SigtrPr;
		
	}

	// dtotal =1-(1 - dc)*(1 - dt);   //not adding cyclic effect yet
	Matrix PeigT(2, 2);
	Matrix tempSigPr(2, 2);
	Matrix Sign1(2, 2);


	tempSigPr.Zero();
	//sigPr(0) = (1 - dt) * sigPr(0);
	//sigPr(1) = (1 - dc) * sigPr(1);

	//dt = 0; dc = 0;

	
	if (sigPr(0) > 0 && sigPr(1) > 0) {
		sigPr(0) = (1 - dt) * sigPr(0);
		sigPr(1) = (1 - dt) * sigPr(1);
	}
	 else if (sigPr(0) > 0 && sigPr(1) < 0) {
		 sigPr(0) = (1 - dt)*sigPr(0);
		 sigPr(1) = (1 - dc)*sigPr(1);
	 }
	 else if (sigPr(0) < 0 && sigPr(1) < 0) {
		sigPr(0) = (1 - dc) * sigPr(0);
		sigPr(1) = (1 - dc) * sigPr(1);
	}



	 tempSigPr(0, 0) = sigPr(0); tempSigPr(1, 1) = sigPr(1);

	 PeigT.addMatrixTranspose(0, Peig, 1);

	 Sign1 = Peig*tempSigPr*PeigT;

	 sige(0) = Sign1(0, 0); sige(1) = Sign1(1, 1); sige(2) = Sign1(0, 1);
		
	 sig = sige;

	//adding remaini
	//if (dtotal > 0.99)
	//dtotal = 0.99;

	C =  Cbar;

#ifdef _SDEBUG
	if (this->getTag() == 109) {
		opserr << this->getTag() << "--kt--: " << kt << "--kc--: " << kc << endln;
		opserr << "---eps:" << eps << "--epsp:" << eps_p << endln;
		opserr << "  Sig: " << sig << " DeltaPeps " << delta_epsPr << endln;
	}
#endif

	///////////////////////////////////////////////output for debug//////////////////////////////
	//#ifdef _DEBUG
	//if (this->getTag() == 8512) {
		//opserr << sig(0) << " " << sig(1) << " " << sig(2) << endln;
	//}

#ifdef _DEBUG
  if (this->getTag() == 2589|| this->getTag() == 2599) {
	  
	std::string recorder = "ShellData/Concrete";
	std::string recorderSuffix = ".out";
	int Rtag = this->getTag();
	std::stringstream f;
	f << recorder << Rtag << recorderSuffix;
	std::string filename;
	filename = f.str();
	//static const char* fname = filename.c_str();
	const char* fname = filename.c_str();
	
	std::ofstream outpt;
	if(Tchange==0){
		outpt.open(fname, ios::out);
		Tchange ++;
	} 
	else 
		outpt.open(fname, ios::out | ios::app);
   
//outpt.open(fname);
	
	if (Tchange == 1) {
		outpt << endln << TempAndElong(0) << endln;
		Tchange++;
	}
	
 // outpt<< "kt: "<< kt << ", kc: " << kc<< "...eps:" << eps(0)<<", "<<eps(1)<<", "<<eps(2) ;
	//outpt << "...epsp:" << eps_p(0) << ", " << eps_p(1) << ", " << eps_p(2) << "..Fres " << Fres << endln;
//	outpt << "..SigPr: " << sigPr(0) << "," << sigPr(1) << "...sig: " << sig(0) << ", " << sig(1) << ", " << sig(2) << "...fc " << fcbar << "...ftbar " << ftbar << "..C00 " << C(0, 0) << "..C22" << C(2, 2) << endln;

	//opserr << fyield<<" "<< kt << "  " << kc << " Fres " << Fres << " sigma " << sig(0)<< endln;
	outpt<< fyield << " " << kt << "  " << kc<< " " << eps(0)<<" "<<eps(1)<<"  "<<eps(2) ;
	outpt << "  " << eps_p(0) << "  " << eps_p(1) << "  " << eps_p(2) << "   " << Fres;
  outpt << "   "<<sigPr(0)<<"  "<<sigPr(1)<<"   " << sig(0)<<"  " <<sig(1)<<"   "<<sig(2)<< "    "<<fcbar<<"   "<<ftbar<<"   "<< delta_epsPr(0)<<"   "<< delta_epsPr(1) <<"  "<< Qnorm <<endln;
	
  }
  ///////////////////////////////////////////////output for debug////////////////////////////// 
 #endif
  
   if (sig(0) != sig(0))
	   opserr << "invalid sig" <<sig(0)<< endln;
   
   // TempAndElong(0) = this->getTag();
   //TempAndElong(0)= kt;
   //TempAndElong(1) = kc;

	if (sig(0) != sig(0))
		opserr << "invalid sig" << sig(0) << endln;

	// TempAndElong(0) = kc;
	TempAndElong(0) = sigPr(0);
	TempAndElong(1) = sigPr(1);


	//-================================
	return 0;
}


int
PlasticDamageConcretePlaneStressThermal::setTrialStrainIncr(const Vector &dStrain)
{
	eps += dStrain;
	this->setTrialStrain(eps);
	return 0;
}

int
PlasticDamageConcretePlaneStressThermal::setTrialStrainIncr(const Vector &dStrain, const Vector &rate)
{
	eps += dStrain;
	this->setTrialStrain(eps);
	return 0;
}

const Matrix&
PlasticDamageConcretePlaneStressThermal::getTangent(void)
{
	return C;
}

const Matrix&
PlasticDamageConcretePlaneStressThermal::getInitialTangent(void)
{
	return Ce0;
}

const Vector&
PlasticDamageConcretePlaneStressThermal::getStress(void)
{
	return sig;
}

const Vector&
PlasticDamageConcretePlaneStressThermal::getStrain(void)
{
	return eps;
}

int
PlasticDamageConcretePlaneStressThermal::commitState(void)
{
	ktCommit = kt;
	kcCommit = kc;
	ECommit = E;
	dtCommit = dt;
	dcCommit = dc;

	epsCommit = eps;
	sigCommit = sig;
	sigeCommit = sige;
	eps_pCommit = eps_p;
	sigePCommit = sige;
	if(Cchange >0 )
		Cchange = 0;
	return 0;
}

int
PlasticDamageConcretePlaneStressThermal::revertToLastCommit(void)
{
	C = Ccommit;
	kt = ktCommit;
	kc = kcCommit;
	E = ECommit ;
	dt = dtCommit;
	dc = dcCommit;

	eps = epsCommit;
	sig = sigCommit;
	sige = sigeCommit;
	eps_p = eps_pCommit;
	sigeP = sigePCommit;
	return 0;
}

int
PlasticDamageConcretePlaneStressThermal::revertToStart(void)
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
PlasticDamageConcretePlaneStressThermal::getCopy(const char *type)
{
	if (strcmp(type, "PlaneStress") == 0 || strcmp(type, "PlaneStress2D") == 0) {
		//for debugging
		matTag++;
		PlasticDamageConcretePlaneStressThermal *theCopy =
			new PlasticDamageConcretePlaneStressThermal(matTag, E0, nu, ft0, fc0*1.562491022,gt,gc, At, Ac, Dbart, Dbarc);
		return theCopy;
	}
	else {
		return 0;
	}
}

NDMaterial*
PlasticDamageConcretePlaneStressThermal::getCopy(void)
{
	PlasticDamageConcretePlaneStressThermal *theCopy =
		new PlasticDamageConcretePlaneStressThermal(this->getTag(), E0, nu, ft0, fc0*1.562491022, gt,gc,At, Ac, Dbart, Dbarc);
	return theCopy;
}

const char*
PlasticDamageConcretePlaneStressThermal::getType(void) const
{
	return "PlaneStress2d";
}

int
PlasticDamageConcretePlaneStressThermal::getOrder(void) const
{
	return 6;
}


int
PlasticDamageConcretePlaneStressThermal::sendSelf(int commitTag, Channel &theChannel)
{
	static Vector data(10);

	int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "PlasticDamageConcretePlaneStressThermal::sendSelf -- could not send Vector\n";
		return res;
	}

	return res;
}

int
PlasticDamageConcretePlaneStressThermal::recvSelf(int commitTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	static Vector data(10);

	int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "PlasticDamageConcretePlaneStressThermal::sendSelf -- could not send Vector\n";
		return res;
	}

	return res;
}

void
PlasticDamageConcretePlaneStressThermal::Print(OPS_Stream &s, int flag) {
	//opserr << "PlasticDamageConcretePlaneStressThermal: " << this->getTag();
	//opserr << "strain: " << eps;
	//opserr << "strain: " << sig;
	//opserr << "tangent: " << this->getTangent();
}

//send back TempAndElong(Liming,UoE)
const Vector&
PlasticDamageConcretePlaneStressThermal::getTempAndElong(void)
{
	return TempAndElong;
}

//Set TemperatureAndElongation

double
PlasticDamageConcretePlaneStressThermal::setThermalTangentAndElongation(double &tempT, double&ET, double&Elong)
{
	bool Lits = true;
	double Eps_lits = 0;
	Temp = tempT;
	double epsc0 = -0.0025;
	double epscu = -0.02;
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
		epsc0 = -0.0025;
		//fcu = fcuT;
		epscu = -0.02;
		//Ets = EtsT;  jz what is there the statement?
	}
	else if (Temp <= 80) {
		fc = fc0;
		epsc0 = -(0.0025 + (0.004 - 0.0025)*(Temp - 0) / (80 - 0));
		epscu = -(0.02 + (0.025 - 0.02)*(Temp - 0) / (80 - 0));
	}
	else if (Temp <= 180) {
		fc = fc0*(1 - (Temp - 80)*0.05 / 100);
		epsc0 = -(0.0040 + (0.0055 - 0.0040)*(Temp - 80) / 100);
		epscu = -(0.0225 + (0.025 - 0.0225)*(Temp - 80) / 100);
	}
	else if (Temp <= 280) {
		fc = fc0*(0.95 - (Temp - 180)*0.1 / 100);
		epsc0 = -(0.0055 + (0.0070 - 0.0055)*(Temp - 180) / 100);
		epscu = -(0.025 + (0.0275 - 0.025)*(Temp - 180) / 100);

	}
	else if (Temp <= 380) {
		fc = fc0*(0.85 - (Temp - 280)*0.1 / 100);
		epsc0 = -(0.0070 + (0.0100 - 0.0070)*(Temp - 280) / 100);
		epscu = -(0.0275 + (0.03 - 0.0275)*(Temp - 280) / 100);
	}
	else if (Temp <= 480) {
		fc = fc0*(0.75 - (Temp - 380)*0.15 / 100);
		epsc0 = -(0.0100 + (0.0150 - 0.0100)*(Temp - 380) / 100);
		epscu = -(0.03 + (0.0325 - 0.03)*(Temp - 380) / 100);
	}
	else if (Temp <= 580) {
		fc = fc0*(0.60 - (Temp - 480)*0.15 / 100);
		epsc0 = -(0.0150 + (0.0250 - 0.0150)*(Temp - 480) / 100);
		epscu = -(0.0325 + (0.035 - 0.0325)*(Temp - 480) / 100);

	}
	else if (Temp <= 680) {
		fc = fc0*(0.45 - (Temp - 580)*0.15 / 100);
		epsc0 = -0.0250;
		epscu = -(0.035 + (0.0375 - 0.035)*(Temp - 580) / 100);
	}
	else if (Temp <= 780) {
		fc = fc0*(0.30 - (Temp - 680)*0.15 / 100);
		epsc0 = -0.0250;
		epscu = -(0.0375 + (0.04 - 0.0375)*(Temp - 680) / 100);

	}
	else if (Temp <= 880) {
		fc = fc0*(0.15 - (Temp - 780)*0.07 / 100);
		epsc0 = -0.0250;
		epscu = -(0.04 + (0.0425 - 0.04)*(Temp - 780) / 100);
	}
	else if (Temp <= 980) {
		fc = fc0*(0.08 - (Temp - 880)*0.04 / 100);
		epsc0 = -0.0250;
		epscu = -(0.0425 + (0.045 - 0.0425)*(Temp - 880) / 100);
	}
	else if (Temp <= 1080) {
		fc = fc0*(0.04 - (Temp - 980)*0.03 / 100);
		epsc0 = -0.0250;
		epscu = -(0.045 + (0.0475 - 0.045)*(Temp - 980) / 100);
	}
	else {
		opserr << "Material temperature " << Temp << " is invalid\n";

	}
	
	double ThermalElongation = 0;
	int matType =1; //Siliceous aggregates
//  int matType =2;//Calcareous aggregates
	if (matType == 1) {
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
			opserr << "Material temperature " << Temp << " is invalid\n";
		}
	}
	else if(matType ==2){
		if (Temp <= 1) {
			ThermalElongation = Temp  * 6e-6;
		}
		else if (Temp <= 680) {
			ThermalElongation = -1.2e-4 + 6e-6 *(Temp + 20) + 1.4e-11 *(Temp + 20)*(Temp + 20)*(Temp + 20);
		}
		else if (Temp <= 1180) {
			ThermalElongation = 12e-3;
		}
		else {
			opserr << "Material temperature "<<Temp<<" is invalid\n";
		}
	
	}
	
	//ThermalElongation = Temp*1.0e-5;
	double stif = -fc*0.0025 / epsc0 / fc0;

	E = stif*E0;
//if (E < 0)
	//	E = 1e-4*E0;

		//E = -fc*0.0025 / epsc0 / fc0*E0;
	//Elong = (1 - dt)* ThermalElongation;
	//Es = -fc*0.0025/epsc0/fc0*Es0;

	if (Lits) {
		if (sigCommit(0) < 0 && sigCommit(1) < 0) {
			double Sig_factor = -(sigCommit(0) + sigCommit(1)) / 2 / fc0 / 1.562491022;
			Eps_lits = (4.12e-5 * Temp - 1.72e-7 * Temp * Temp + 3.3e-10 * Temp * Temp * Temp) * Sig_factor;
			if (Eps_lits > epsLitsp)
				epsLitsp = Eps_lits;

		}
		if (epsLitsp > 0)
			ThermalElongation = ThermalElongation - epsLitsp;
	}


   //else
	//	ThermalElongation = 0;


	Elong = ThermalElongation;
	// this->setTempInitials();
	// E = E0;
	//Es = Es0;
	//fc = fc0;
	//ft = ft0;
	gt = gt0*ft / ft0*E0 / E;
	gc = gc0*fc*epscu / fc0/(-0.02);
	
//TempAndElong(0) = kt;
//TempAndElong(1) = dt;
	//TempAndElong(0) = Temp;
	//TempAndElong(1) = ThermalElongation;
	//TempAndElong(1) = E;
	ET = E;
	Tchange = 1;

	return 0;
}
