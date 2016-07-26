#include <QuadFour.h>
#include <HeatTransferNode.h>
#include <HeatTransferMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
//#include <Renderer.h>
#include <HeatTransferDomain.h>
#include <cmath>
#include <Radiation.h>
#include <Convection.h>
#include <PrescribedSurfFlux.h>

double QuadFour::matrixData[16];
Matrix QuadFour::K(matrixData, 4, 4);
Vector QuadFour::Q(4);
double QuadFour::shp[3][4];
double QuadFour::pts[4][2];
double QuadFour::wts[4];
double QuadFour::shp2[2];
double QuadFour::pts2[2];
double QuadFour::wts2[2];
int QuadFour::npface[4][2];

using std::list;

QuadFour::QuadFour(int tag, int nd1, int nd2, int nd3, int nd4, 
				   HeatTransferMaterial& m,
				   bool phaseChange)
:HeatTransferElement(tag), theMaterial(0), connectedExternalNodes(4), Qp(4),
 phaseTransformation(phaseChange)
{
	pts[0][0] = -0.5773502691896258;
	pts[0][1] = -0.5773502691896258;
	pts[1][0] =  0.5773502691896258;
	pts[1][1] = -0.5773502691896258;
	pts[2][0] =  0.5773502691896258;
	pts[2][1] =  0.5773502691896258;
	pts[3][0] = -0.5773502691896258;
	pts[3][1] =  0.5773502691896258;

	wts[0] = 1.0;
	wts[1] = 1.0;
	wts[2] = 1.0;
	wts[3] = 1.0;

	pts2[0] = -0.5773502691896258;
	pts2[1] =  0.5773502691896258;
	
	wts2[0] = 1.0;
	wts2[1] = 1.0;

	npface[0][0] = 0;
    npface[0][1] = 1;
    npface[1][0] = 1;
    npface[1][1] = 2;
    npface[2][0] = 2;
    npface[2][1] = 3;
    npface[3][0] = 3;
    npface[3][1] = 0;

    // Allocate arrays of pointers to NDMaterials
    theMaterial = new HeatTransferMaterial*[4];
    
	if (theMaterial == 0) {
		opserr << "QuadFour::QuadFour - failed allocate material model pointer\n";
		exit(-1);
		}

	int i;
	for (i = 0; i < 4; i++) {

		// Get copies of the material model for each integration point
		theMaterial[i] = m.getCopy();

		// Check allocation
		if (theMaterial[i] == 0) {
			opserr << "QuadFour::QuadFour -- failed to get a copy of material model\n";
			exit(-1);
			}
		}

    // Set connected external node IDs
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;
    connectedExternalNodes(2) = nd3;
    connectedExternalNodes(3) = nd4;
    
    for (i=0; i<4; i++)
      theNodes[i] = 0;
}

QuadFour::~QuadFour()
{    
    for (int i = 0; i < 4; i++) {
		if (theMaterial[i])
			delete theMaterial[i];
		}

  // Delete the array of pointers to NDMaterial pointer arrays
	if (theMaterial)
		delete[] theMaterial;
}


int
QuadFour::getNumExternalNodes() const
{
    return 4;
}


const ID&
QuadFour::getExternalNodes()
{
    return connectedExternalNodes;
}


HeatTransferNode**
QuadFour::getNodePtrs(void) 
{
    return theNodes;
}


int
QuadFour::getNumDOF()
{
    return 4;
}


void
QuadFour::setDomain(HeatTransferDomain* theDomain)
{
	// Check HeatTransferDomain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
		theNodes[0] = 0;
		theNodes[1] = 0;
		theNodes[2] = 0;
		theNodes[3] = 0;
		return;
		}

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    int Nd3 = connectedExternalNodes(2);
    int Nd4 = connectedExternalNodes(3);

    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);
    theNodes[2] = theDomain->getNode(Nd3);
    theNodes[3] = theDomain->getNode(Nd4);

    if (theNodes[0] == 0 || theNodes[1] == 0 || theNodes[2] == 0 || theNodes[3] == 0) {
		return;
		}

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    int dofNd3 = theNodes[2]->getNumberDOF();
    int dofNd4 = theNodes[3]->getNumberDOF();
    
    if (dofNd1 != 1 || dofNd2 != 1 || dofNd3 != 1 || dofNd4 != 1) {
		return;
		}

    this->HeatTransferDomainComponent::setDomain(theDomain);
}


int
QuadFour::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
	if ((retVal = this->HeatTransferElement::commitState()) != 0) {
		opserr << "QuadFour::commitState () - failed in base class";
		}    

	// Loop over the integration points and commit the material states
	for (int i = 0; i < 4; i++)
		retVal += theMaterial[i]->commitState();

    return retVal;
}


int
QuadFour::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < 4; i++)
		retVal += theMaterial[i]->revertToLastCommit();

    return retVal;
}


int
QuadFour::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < 4; i++)
		retVal += theMaterial[i]->revertToStart();

    return retVal;
}


int
QuadFour::update()
{
	const Vector& Temp1 = theNodes[0]->getTrialTemperature();
	const Vector& Temp2 = theNodes[1]->getTrialTemperature();
	const Vector& Temp3 = theNodes[2]->getTrialTemperature();
	const Vector& Temp4 = theNodes[3]->getTrialTemperature();
	
	static double T[4];

	T[0] = Temp1(0);
	T[1] = Temp2(0);
	T[2] = Temp3(0);
	T[3] = Temp4(0);

	double T_hat;

	int ret = 0;

	// Loop over the integration points
	for (int i = 0; i < 4; i++) {

		// Determine Jacobian for this integration point
		this->shapeFunction(pts[i][0], pts[i][1]);

		// Interpolate nodal temperatures
		// T_hat = NT;
		T_hat = 0;
		for (int beta = 0; beta < 4; beta++) {
			T_hat += shp[2][beta] * T[beta];
			}

		// Set the temperature in material model
		ret += theMaterial[i]->setTrialTemperature(T_hat);
		}

	return ret;
}


const Matrix&
QuadFour::getConductionTangent()
{
	K.Zero();

	double dvol;

	// opserr << "the element tag is " << this->getTag() << " \n";
	// Loop over the integration points
	for (int i = 0; i < 4; i++) {
		// Determine Jacobian for this integration point
		dvol = this->shapeFunction(pts[i][0], pts[i][1]);
		dvol *= wts[i];
		// Get the material information--heat conductivity
		const Matrix& D = theMaterial[i]->getConductivity();

		// Perform numerical integration
		// K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;

		double D00 = D(0,0); double D11 = D(1,1);  // Just deal with OrthotropicMaterial 

		// Form elemental tangent  
		for (int m = 0; m < 4; m++) {
			for (int n = 0; n < 4; n++) {
				K(m,n) += (shp[0][m] * D00 * shp[0][n] + shp[1][m] * D11 * shp[1][n]) * dvol;
				}
			}
		}   
	
	//opserr << this->getTag() << "conguction tang:" << endln;
	//opserr << K << endln;
    return K;
}


const Matrix&
QuadFour::getCapacityTangent()
{
	K.Zero();
	if (phaseTransformation == true){
		static double enth[4];
		static double entht[4];
		static double rc[4]; // rc = rho times specific heat;
		const Vector& Temp1 = theNodes[0]->getTrialTemperature();
		const Vector& Temp2 = theNodes[1]->getTrialTemperature();
		const Vector& Temp3 = theNodes[2]->getTrialTemperature();
		const Vector& Temp4 = theNodes[3]->getTrialTemperature();

		const Vector& Tmp1 = theNodes[0]->getTemperature();
		const Vector& Tmp2 = theNodes[1]->getTemperature();
		const Vector& Tmp3 = theNodes[2]->getTemperature();
		const Vector& Tmp4 = theNodes[3]->getTemperature();

		static double T[4];
		static double Tt[4];

		T[0] = Temp1(0);
		T[1] = Temp2(0);
		T[2] = Temp3(0);
		T[3] = Temp4(0);

		Tt[0] = Tmp1(0);
		Tt[1] = Tmp2(0);
		Tt[2] = Tmp3(0);
		Tt[3] = Tmp4(0);

		// getting the nodal enthalpy values
		//enth[0] = theMaterial[0]->getEnthalpy(T[0]);
		//enth[1] = theMaterial[0]->getEnthalpy(T[1]);
		//enth[2] = theMaterial[0]->getEnthalpy(T[2]);
		//enth[3] = theMaterial[0]->getEnthalpy(T[3]);

		entht[0] = theMaterial[0]->getEnthalpy(Tt[0]);
		entht[1] = theMaterial[0]->getEnthalpy(Tt[1]);
		entht[2] = theMaterial[0]->getEnthalpy(Tt[2]);
		entht[3] = theMaterial[0]->getEnthalpy(Tt[3]);

		double dHdx, dHdy, dTdx, dTdy, Ht, Tcmt, Ttrial;	
		double rcdvol;
 
		// Loop over the integration points
		for (int i = 0; i < 4; i++) {
			// Determine Jacobian for this integration point
			rcdvol = this->shapeFunction(pts[i][0], pts[i][1]);
			dHdx = 0;
			dHdy = 0;
			dTdx = 0;
			dTdy = 0;
			Ht = 0;
			Tcmt = 0;
			Ttrial = 0;
			for (int beta = 0; beta < 4; beta++) {
				dHdx += shp[0][beta] * enth[beta];
				dHdy += shp[1][beta] * enth[beta];
				dTdx += shp[0][beta] * T[beta];
				dTdy += shp[1][beta] * T[beta];
                // get committed enthalpy values at integration point,
				// needed by Morgan's method
				Ht += shp[2][beta] * entht[beta];
				Tcmt += shp[2][beta] * Tt[beta];
				Ttrial += shp[2][beta] * T[beta];
				}
			double dT = Ttrial - Tcmt;
			if (dT == 0) {
				double rhoi = theMaterial[i]->getRho(); 
				rc[i] = rhoi * (theMaterial[i]->getSpecificHeat()); 
				} else {
					// using Morgan's approximation
					double H = theMaterial[i]->getEnthalpy();
					rc[i] = (H - Ht)/dT;
				}


			//if ((dTdx == 0) && (dTdy ==0)) {
			//	double rhoi = theMaterial[i]->getRho(); 
			//	rc[i] = rhoi * (theMaterial[i]->getSpecificHeat()); 
			//	} else {
			//		double dT = dTdx * dTdx + dTdy * dTdy;
			//		// using Lemmon's approximation
			//		//double dH = dHdx * dHdx + dHdy * dHdy;
			//		//rc[i] = sqrt(dH/dT);
			//		
			//		// using DelGuidice's approximation
			//		double dH = dHdx * dTdx + dHdy * dTdy;
			//		rc[i] = dH/dT;
			//	}

			rcdvol *= (rc[i] * wts[i]);

			for (int m = 0; m < 4; m++) {
				for (int n = 0; n < 4; n++) {
					K(m,n) += shp[2][m] * rcdvol * shp[2][n];
					}
				}
			}
		} else {
			static double rhoi[4];
			static double cpi[4];
			//opserr << this->getTag() << " capacity tangent\n";

			for (int i = 0; i < 4; i++) {
				// enable to account for temperature dependent density
				rhoi[i] = theMaterial[i]->getRho(); 
				cpi[i] = theMaterial[i]->getSpecificHeat(); // get specific heat
				}

			double rcdvol;

			// Loop over the integration points
			for (int i = 0; i < 4; i++) {

				// Determine Jacobian for this integration point
				rcdvol = this->shapeFunction(pts[i][0], pts[i][1]);
				rcdvol *= (rhoi[i] * cpi[i] * wts[i]);

				for (int m = 0; m < 4; m++) {
					for (int n = 0; n < 4; n++) {
						K(m,n) += shp[2][m] * rcdvol * shp[2][n];
						}
					}
				}
		}

	return K;
}


const Matrix&
QuadFour::getRadiationTangent()
{
    K.Zero();
	list<Radiation*>::iterator rad_iter = theRadiationBCs.begin();

	while (rad_iter != theRadiationBCs.end()) {
		if (*rad_iter != 0){
			Radiation* theRadiationBC = *rad_iter;
			int tag = theRadiationBC->getFaceTag();
			double eps_sgm = theRadiationBC->getParameter();  // get the product of epsilon and sigma
			double ds;

			for (int i = 0; i < 2; i++) {
				ds = this->shapeFunction(pts2[i],tag);
				ds *= wts2[i];

				double T_hat = this->getInterpolatedTemp(tag);
				double hr = eps_sgm * 4 * pow(T_hat, 3);

				for (int m = 0; m < 2; m++) {
					int mm = npface[tag-1][m];
					for (int n = 0; n < 2; n++) {
						int nn = npface[tag-1][n];
						K(mm, nn) += shp2[m] * hr * shp2[n] * ds;
						}
					}
				}
			}
		rad_iter++;
		}

	//opserr << this->getTag() << " radiation K:" << endln;
	//opserr << K << endln;
	return K;
}


const Matrix& 
QuadFour::getConvectionTangent()
{
    K.Zero();
	list<Convection*>::iterator convec_iter = theConvectionBCs.begin();

	while (convec_iter != theConvectionBCs.end()) {
		if (*convec_iter != 0) {
			Convection* theConvectionBC = *convec_iter;
			int tag = theConvectionBC->getFaceTag();
			double h = theConvectionBC->getParameter();  // get the covective heat transfer coefficient

			double ds;

			for (int i = 0; i < 2; i++){
				ds = this->shapeFunction(pts2[i], tag);
				ds *= wts2[i];

				for (int m = 0; m < 2; m++) {
					int mm = npface[tag-1][m];
					for (int n = 0; n < 2; n++) {
						int nn = npface[tag-1][n];
						K(mm, nn) += shp2[m] * h * shp2[n] * ds;
						}
					}
				}
			}
		convec_iter++;
		}
	//opserr << this->getTag() << " convection tang:" << endln;
	//opserr << K << endln;

	return K;
}


void
QuadFour::zeroFlux()
{
   Qp.Zero();
   return;
}


int 
QuadFour::addPrecribedSurfFlux(PrescribedSurfFlux* theFlux, double factor)
{
    int tag = theFlux->getFaceTag();
	double data = theFlux->getData();
	data *= -factor;

	double ds = this->shapeFunction(0, tag);
	ds *= 2;

	int fnode1 = npface[tag-1][0];
	int fnode2 = npface[tag-1][1];

	data *= 0.5 * ds;

	// For constant flux along an element face, half of the flux is asigned to each node 
	// of the element. See P619 Burnett's book 
	Qp(fnode1) = data;
	Qp(fnode2) = data;
	
	return 0;
}


const Vector&
QuadFour::get_Q_Transient()
{
    Q.Zero();

    this->getCapacityTangent();

	const Vector& Tdot1 = theNodes[0]->getTrialTdot();
	const Vector& Tdot2 = theNodes[1]->getTrialTdot();
	const Vector& Tdot3 = theNodes[2]->getTrialTdot();
	const Vector& Tdot4 = theNodes[3]->getTrialTdot();

	static double Td[4];

	Td[0] = Tdot1(0);
	Td[1] = Tdot2(0);
	Td[2] = Tdot3(0);
	Td[3] = Tdot4(0);
    
	for (int m = 0; m < 4; m++) {
		for (int n = 0; n < 4; n++) {
			Q(m) += (-K(m, n) * Td[n]);
			}
		}

	return Q;	
}


const Vector&
QuadFour::get_Q_Conduction()
{
	Q.Zero();	
	
	this->getConductionTangent();

	const Vector& Temp1 = theNodes[0]->getTrialTemperature();
	const Vector& Temp2 = theNodes[1]->getTrialTemperature();
	const Vector& Temp3 = theNodes[2]->getTrialTemperature();
	const Vector& Temp4 = theNodes[3]->getTrialTemperature();

	static double Ttrial[4];

	Ttrial[0] = Temp1(0);
	Ttrial[1] = Temp2(0);
	Ttrial[2] = Temp3(0);
	Ttrial[3] = Temp4(0);

	for (int m = 0; m < 4; m++) {
		for (int n = 0; n < 4; n++) {
			Q(m) += (-K(m, n) * Ttrial[n]);
			}
		}

	Q.addVector(1.0, Qp, 1.0);  // addPrecribedSurfFlux, for efficiency reason, may
                                // need to test the existence of prescribed flux
	return Q;	
}


const Vector&
QuadFour::get_Q_Radiation()
{
    Q.Zero();
	list<Radiation*>::iterator rad_iter = theRadiationBCs.begin();

	while (rad_iter != theRadiationBCs.end()) {
		if (*rad_iter != 0){
			Radiation* theRadiationBC = *rad_iter;
			int tag = theRadiationBC->getFaceTag();
			double eps_sgm = theRadiationBC->getParameter();  // get the product of epsilon and sigma

			double alpha = theRadiationBC->getAbsorptivity();
			double g = theRadiationBC->getIrradiation();  // get G

			double ds;

			for (int i = 0; i < 2; i++) {
				ds = this->shapeFunction(pts2[i], tag);
				ds *= wts2[i];

				double T_hat = this->getInterpolatedTemp(tag);
				double q_ir = alpha * g - eps_sgm * pow(T_hat, 4); 

				for (int m = 0; m < 2; m++) {
					int mm = npface[tag-1][m];
					Q(mm) += shp2[m] * q_ir * ds;
					}
				}
			}
		rad_iter++;
		}

	//opserr << this->getTag() << " radiation Q:" << endln;
	//opserr << Q << endln;
	return Q;
}


const Vector&
QuadFour::get_Q_Convection()
{
    Q.Zero();
	list<Convection*>::iterator convec_iter = theConvectionBCs.begin();

	while (convec_iter != theConvectionBCs.end()) {
		if (*convec_iter != 0) {
			Convection* theConvectionBC = *convec_iter;	
			int tag = theConvectionBC->getFaceTag();
			double h = theConvectionBC->getParameter();
			double Ta = theConvectionBC->getAmbientTemp();

			double ds;

			for (int i = 0; i < 2; i++) {
				ds = this->shapeFunction(pts2[i], tag);
				ds *= wts2[i];

				double T_hat = this->getInterpolatedTemp(tag);
				double qc = h * (Ta - T_hat);

				for (int m = 0; m < 2; m++) {
					int mm = npface[tag-1][m];
					Q(mm) += shp2[m] * qc * ds;
					}
				}
			}
		convec_iter++;
		}
    //opserr << Q <<endln;

	//opserr << this->getTag() << " convection Q:" << endln;
	//opserr << Q << endln;
	return Q;
}


double 
QuadFour::shapeFunction(double xi, double eta)
{
	const Vector& nd1Crds = theNodes[0]->getCrds();
	const Vector& nd2Crds = theNodes[1]->getCrds();
	const Vector& nd3Crds = theNodes[2]->getCrds();
	const Vector& nd4Crds = theNodes[3]->getCrds();

	double oneMinuseta = 1.0 - eta;
	double onePluseta = 1.0 + eta;
	double oneMinusxi = 1.0 - xi;
	double onePlusxi = 1.0 + xi;

	shp[2][0] = 0.25 * oneMinusxi * oneMinuseta;	// N_1
	shp[2][1] = 0.25 * onePlusxi * oneMinuseta;		// N_2
	shp[2][2] = 0.25 * onePlusxi * onePluseta;		// N_3
	shp[2][3] = 0.25 * oneMinusxi * onePluseta;		// N_4

	double J[2][2];

	J[0][0] = 0.25 * (-nd1Crds(0) * oneMinuseta + nd2Crds(0) * oneMinuseta +
				nd3Crds(0) * (onePluseta) - nd4Crds(0) * (onePluseta));
	//J[0][1] = 0.25 * (-nd1Crds(0) * oneMinusxi - nd2Crds(0) * onePlusxi +
	//			nd3Crds(0) * onePlusxi + nd4Crds(0) * oneMinusxi);
	//J[1][0] = 0.25 * (-nd1Crds(1) * oneMinuseta + nd2Crds(1) * oneMinuseta +
	//			nd3Crds(1) * onePluseta - nd4Crds(1) * onePluseta);

	J[1][0] = 0.25 * (-nd1Crds(0) * oneMinusxi - nd2Crds(0) * onePlusxi +
		nd3Crds(0) * onePlusxi + nd4Crds(0) * oneMinusxi);
	J[0][1] = 0.25 * (-nd1Crds(1) * oneMinuseta + nd2Crds(1) * oneMinuseta +
		nd3Crds(1) * onePluseta - nd4Crds(1) * onePluseta);

	J[1][1] = 0.25 * (-nd1Crds(1) * oneMinusxi - nd2Crds(1) * onePlusxi +
				nd3Crds(1) * onePlusxi + nd4Crds(1) * oneMinusxi);

	double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
	
	if (detJ < 0) {
		opserr << "QuadFour::shapeFunction(xi,eta) -- detJ < 0 ";
		exit(-1);
		}

	double oneOverdetJ = 1.0/detJ;
	double L[2][2];

	// L = inv(J)
	L[0][0] =  J[1][1] * oneOverdetJ;
	L[1][0] = -J[0][1] * oneOverdetJ;
	L[0][1] = -J[1][0] * oneOverdetJ;
	L[1][1] =  J[0][0] * oneOverdetJ;

    double L00 = 0.25 * L[0][0];
    double L10 = 0.25 * L[1][0];
    double L01 = 0.25 * L[0][1];
    double L11 = 0.25 * L[1][1];
	
	double L00oneMinuseta = L00 * oneMinuseta;
	double L00onePluseta  = L00 * onePluseta;
	double L01oneMinusxi  = L01 * oneMinusxi;
	double L01onePlusxi   = L01 * onePlusxi;

	double L10oneMinuseta = L10 * oneMinuseta;
	double L10onePluseta  = L10 * onePluseta;
	double L11oneMinusxi  = L11 * oneMinusxi;
	double L11onePlusxi   = L11 * onePlusxi;

	// See Cook, Malkus, Plesha p. 169 for the derivation of these terms
    shp[0][0] = -L00oneMinuseta - L01oneMinusxi;	// N_1,1
    shp[0][1] =  L00oneMinuseta - L01onePlusxi;		// N_2,1
    shp[0][2] =  L00onePluseta  + L01onePlusxi;		// N_3,1
    shp[0][3] = -L00onePluseta  + L01oneMinusxi;	// N_4,1
	
    shp[1][0] = -L10oneMinuseta - L11oneMinusxi;	// N_1,2
    shp[1][1] =  L10oneMinuseta - L11onePlusxi;		// N_2,2
    shp[1][2] =  L10onePluseta  + L11onePlusxi;		// N_3,2
    shp[1][3] = -L10onePluseta  + L11oneMinusxi;	// N_4,2

    return detJ;
}


double 
QuadFour::shapeFunction(double xi, int faceTag)
{
	if (faceTag < 1 || faceTag > 4) {
		opserr << "QuadFour::shapeFunction -- improper faceTag: " << faceTag << "for QuadFour\n";
		exit(-1);
		}

	double oneMinusxi = 1.0 - xi;
	double onePlusxi = 1.0 + xi;

	shp2[0] = 0.5 * oneMinusxi;	// N_1
	shp2[1] = 0.5 * onePlusxi;	// N_2

	int node1 = npface[faceTag-1][0];
	int node2 = npface[faceTag-1][1];

	const Vector& nd1Crds = theNodes[node1]->getCrds();
	const Vector& nd2Crds = theNodes[node2]->getCrds();

	double J00, J01;

	J00 = 0.5 * (-nd1Crds(0) + nd2Crds(0));  // dx/dxi
	J01 = 0.5 * (-nd1Crds(1) + nd2Crds(1));  // dy/dxi

	double JGamma = sqrt(pow(J00, 2) + pow(J01, 2));

	if (JGamma < 0){
		opserr << "QuadFour::shapeFunction(xi,eta,faceTag) -- JGamma < 0 ";
		exit(-1);
		}
  
    return JGamma;  // return the bounday Jacobian
}


double 
QuadFour::getInterpolatedTemp(int tag)
{
	int node1 = npface[tag-1][0];
	int node2 = npface[tag-1][1];

	const Vector& Temp1 = theNodes[node1]->getTrialTemperature();
	const Vector& Temp2 = theNodes[node2]->getTrialTemperature();
	static double T[2];
	T[0] = Temp1(0);
	T[1] = Temp2(0);

	double T_hat;
	T_hat = shp2[0] * T[0] + shp2[1] * T[1];

	return T_hat;
}


