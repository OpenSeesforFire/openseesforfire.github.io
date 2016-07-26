/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Fire & Heat Transfer modules developed by:                         **
**   Yaqiang Jiang (y.jiang@ed.ac.uk)                                 **
**   Asif Usmani (asif.usmani@ed.ac.uk)                               **
**                                                                    **
** ****************************************************************** */

//
// Written by Yaqiang Jiang (y.jiang@ed.ac.uk)
//

#ifndef BrickEight_h
#define BrickEight_h

#include <HeatTransferElement.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class HeatTransferNode;
class HeatTransferMaterial;

class BrickEight : public HeatTransferElement
{
  public:
    BrickEight(int tag, int nd1, int nd2, int nd3, int nd4,
               int nd5, int nd6, int nd7, int nd8,
		       HeatTransferMaterial& m, 
			   bool phaseTransformation = false);
    ~BrickEight();

    const char* getClassType(void) const {return "BrickEight";};

    int getNumExternalNodes(void) const;
    const ID& getExternalNodes(void);
    HeatTransferNode** getNodePtrs(void);
    int getNumDOF(void);
	const ID& getNodesOnFace(int faceTag);
	int getNumNodesperFace(void) {return 4;};

    void setDomain(HeatTransferDomain* theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);

    // public methods to obtain tangent matrices    
	const Matrix& getConductionTangent(void);  // K_k
	const Matrix& getCapacityTangent(void);  // K_c   
	const Matrix& getRadiationTangent(void);  // K_qr
	const Matrix& getConvectionTangent(void);  // K_qc

    void zeroFlux();
    int addPrecribedSurfFlux(PrescribedSurfFlux* theFlux, double factor);

	// public methods for calculating heat vectors
    const Vector& get_Q_Transient();  //computeQcp(void); getTransient 
    const Vector& get_Q_Conduction(); //computeQk(void);  getConductionFlux  
    const Vector& get_Q_Radiation();  //computeQqr(void); getRadiation  
	const Vector& get_Q_Convection();  //computeQqc(void); getConvection  

	//const Vector& getEleResidual(void);            

  protected:
    
  private:
    // private attributes - a copy for each object of the class
    HeatTransferMaterial** theMaterial;  // pointer to material objects
    
    ID connectedExternalNodes;  // Tags of quad nodes

    HeatTransferNode* theNodes[8];

    static double matrixData[64];  // array data for matrix
    static Matrix K;  // Element stiffness matrix
    static Vector Q;  // Flux vector
	Vector Qp;  // Stores the PrecribedSurfFlux

	//quadrature data
    static const double root3 ;
    static const double one_over_root3 ;    

    static double shp[4][8];  // Stores shape functions and derivatives (overwritten)
    static double pts[8][3];  // Stores quadrature points
    static const double wts[];  // Stores quadrature weights

	static double shp2[4];  // Shape function for surface quadrature
	static double pts2[4][2];
	static const double wts2[];

	static const int npface[6][4];
	bool phaseTransformation;

    // private member functions - only objects of this class can call these
    double shapeFunction(double xi, double eta, double zeta);
	double shapeFunction(double xi, double eta, int faceTag);
	double getInterpolatedTemp(int tag);  // return interpolated temperature at quadrature point at boundary
};

#endif