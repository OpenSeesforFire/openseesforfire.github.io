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

// $Revision: 1.15 $
// $Date: 2009-08-07 20:01:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/IGAQuad/IGAQuad.h,v $

// Written: MHS
// Created: Feb 2000
// Revised: Dec 2000 for efficiency
//
// Description: This file contains the class definition for IGAQuad.

#ifndef IGAQuad_h
#define IGAQuad_h

#ifndef _bool_h
#include "bool.h"
#endif

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
//typedef struct {
//    double J2;
//    Vector WXg;
//    Vector Xi;
//    Matrix* N;
//}resultDersBasisFunsAtGPS;
//typedef struct {
//    Vector R0;
//    Matrix R1;
//}resultRationalize;
class Node;
class NDMaterial;
class Response;

class IGAQuad : public Element
{
public:
    IGAQuad(int tag, int numCPs1, int ex, int nex, int ey, int ney, int* obfs1, int ndof1,
        int* CPIds, Vector KnotVect_x1, Vector KnotVect_y1, int* numMults,
        NDMaterial& m, const char* type,
        double t, double pressure = 0.0,
        double rho = 0.0,
        double b1 = 0.0, double b2 = 0.0);
    IGAQuad();
    ~IGAQuad();

    const char* getClassType(void) const { return "IGAQuad"; };

    int getNumExternalNodes(void) const;
    const ID& getExternalNodes(void);
    Node** getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain* theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix& getTangentStiff(void);
    const Matrix& getInitialStiff(void);
    const Matrix& getMass(void);

    void zeroLoad();
    int addLoad(ElementalLoad* theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector& accel);

    const Vector& getResistingForce(void);
    const Vector& getResistingForceIncInertia(void);

    // public methods for element output
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker
        & theBroker);

    int displaySelf(Renderer&, int mode, float fact, const char** displayModes = 0, int numModes = 0);
    void Print(OPS_Stream& s, int flag = 0);

    Response* setResponse(const char** argv, int argc,
        OPS_Stream& s);

    int getResponse(int responseID, Information& eleInformation);

    int setParameter(const char** argv, int argc, Parameter& param);
    int updateParameter(int parameterID, Information& info);

    // RWB; PyLiq1 & TzLiq1 need to see the excess pore pressure and initial stresses.
    friend class PyLiq1;
    friend class TzLiq1;

protected:

private:
    // private attributes - a copy for each object of the class

    NDMaterial** theMaterial; // pointer to the ND material objects

    ID connectedExternalNodes; // Tags of quad nodes

    Node* theNodes[9];
    Matrix SpanIdxList; // the list of knot span index in x/y directions 
    Matrix refinedKnotVect;// the refined knot vector in x/y directions
    int numCPs;
    int ndof;
    int obfs[2];
    Vector KnotVect_x;
    Vector KnotVect_y;
    ID eleIdInfo;// stores ex, nex, ey, ney

    static double matrixData[324];  // array data for matrix
    static Matrix K;		// Element stiffness, damping, and mass Matrix
    static Vector P;		// Element resisting force vector
    Vector Q;		        // Applied nodal loads
    double b[2];		// Body forces

    Matrix* N0n; // stores B-Spline Basis Function and 1st order 

    double appliedB[2]; // Body forces applied with load pattern, C.McGann, U.Washington
    int applyLoad;      // flag for body force in load, C.McGann, U.Washington

    Vector pressureLoad;	// Pressure load at nodes

    double thickness;	        // Element thickness
    double pressure;	        // Normal surface traction (pressure) over entire element
                     // Note: positive for outward normal
    double rho;
    static double shp[3][9];	// Stores shape functions and derivatives (overwritten)
    static double pts[9][2];	// Stores quadrature points
    static double wts[9];		// Stores quadrature weights

    //int numKnots;             //number of knots inside IGA element (to be refined later)
    Matrix* N0nx; // stores shape function and derivatives in both direction
    Matrix* N0ny;

    // private member functions - only objects of this class can call these
    Vector shapeFunction(int qx, int qy, ID eleIdInfo, Vector KnotVect_x, Vector KnotVect_y);
    void setPressureLoadAtNodes(void);
    int findSpan(int order, int ncp, int eleId, int nele, Vector KnotVect);
    void DerBasisFuns(double Idx, Vector Pts, int obf, int n, Vector KnotVect, Matrix** N0n);
    void calcDersBasisFunsAtGPs(int obf, int ncp, Vector KnotVect, int d, int NGPs, int Idx, double* J2, Vector* Xg, Matrix** N0n);
    void Rationalize(Vector WeightsCP, Vector N0, Matrix N1, Vector* R0, Matrix* R1);
    void GaussRule(int NGPs1, Vector* Xg, Vector* WXg);
    Matrix* Ki;
    Vector dotProduct(Vector VecA, Vector VecB);
    Vector dotDivide(Vector VecA, Vector VecB);
};

#endif
