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

// This code is modified from the class QuadFour (Written by Yaqiang Jiang)
// Created by Liming Jiang (liming.jiang@ed.ac.uk), 2014
////////////////////////////////////////////////////////

#ifndef LineTwo_h
#define LineTwo_h

#include <HeatTransferElement.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class HeatTransferNode;
class HeatTransferMaterial;

class LineTwo : public HeatTransferElement
{
  public:
    LineTwo(int tag, int nd1, int nd2,
		     HeatTransferMaterial& m, 
			 bool phaseTransformation = false);
    ~LineTwo();

    const char* getClassType(void) const {return "LineTwo";};

    int getNumExternalNodes(void) const;
    const ID& getExternalNodes(void);
    HeatTransferNode** getNodePtrs(void);
    int getNumDOF(void);
	const ID& getNodesOnFace(int faceTag);
	int getNumNodesperFace(void) {return 1;};

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

	// public methods for calculating flux vectors
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

    HeatTransferNode* theNodes[2];

    static double matrixData[4];  // array data for matrix
    static Matrix K;  // Element stiffness matrix
    static Vector Q;  // Flux vector
	Vector Qp;  // Stores the PrecribedSurfFlux

    static double shp[2];  // Stores shape functions 
    static double pts[2];  // Stores quadrature points
    static double wts[2];  // Stores quadrature weights
	static double shpd[2]; //shape funtion derivatives (overwritten)

	//static double shp2[2];  // Shape function and derivatives for line quadrature
	//static double pts2[2];
	//static double wts2[2];
	//static int npface[2];
	bool phaseTransformation;

    // private member functions - only objects of this class can call these
    double shapeFunction(double xi);
	//double shapeFunction(double xi, int faceTag);
	
};

#endif