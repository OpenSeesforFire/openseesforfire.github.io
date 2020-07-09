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
                                                                        
// $Revision: 1.5 $
// $Date: 2007-07-27 19:23:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/BeamColumnJoint3dThermal.h,v $
                                                                        
// Written: NM (nmitra@calpoly.edu)
// Created: April 2002
// Last Revised: January 2007
//
// Description: This file contains the class implementation for beam-column joint.
// This element is a 4 noded 12 dof (3 dof at each node) finite area super-element.
// The element takes in 13 different material types in order to simulate
// the inelastic action observed in a reinforced beam column joint.
// Details about this element can be found in 
// Mitra and Lowes, J. Struct. Eng. ASCE, Volume 133, Number 1 (January 2007), pp. 105-120
// Updates: Several concerning Joint formulation (presently a revised formulation for joints)
                                                                        
#ifndef BeamColumnJoint3dThermal_h
#define BeamColumnJoint3dThermal_h

#include <Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>
#include <FileStream.h>
#include <OPS_Stream.h>

class Node;
class Channel;
class FEM_ObjectBroker;
class Response;
class Renderer;
class UniaxialMaterial;

class BeamColumnJoint3dThermal : public Element
{
 public:
  // default constructor
  BeamColumnJoint3dThermal(); 
  
  // defined constructor
  BeamColumnJoint3dThermal(int tag,int Nd1, int Nd2, 
		    UniaxialMaterial& theMat1, UniaxialMaterial& theMat2,
		    UniaxialMaterial& theMat3);
  
  BeamColumnJoint3dThermal(int tag,int Nd1, int Nd2,
		    UniaxialMaterial& theMat1, UniaxialMaterial& theMat2,
		    UniaxialMaterial& theMat3, double Hgtfac, double Wdtfac);
  
  // default destructor
  ~BeamColumnJoint3dThermal();




  ////////////// public methods to obtain information about dof & connectivity    
  bool	isSubdomain(void) { return false; };

  // get number of external nodes
  int getNumExternalNodes(void) const;

  // return connected external nodes
  const ID& getExternalNodes(void);
  Node** getNodePtrs(void);

  // return number of DOFs
  int getNumDOF(void);

  // set domain performs check on dof and associativity with node
  void setDomain(Domain* theDomain);

  //////////////////////////// public methods to set the state of the element    

  // commit state
  int commitState(void);

  // revert to last commit
  int revertToLastCommit(void);

  // revert to start
  int revertToStart(void);

  // determine current strain and set strain in material
  int update(void);

  //////////////////////// public methods to obtain stiffness, mass, damping and 
  ////////////////////////////////////// residual information    

  // returns converged tangent stiffness matrix
  const Matrix& getTangentStiff(void);
  const Matrix& getInitialStiff(void);

  // not required for this element formulation
  const Matrix& getDamp(void);
  const Matrix& getMass(void);

  // not required for this element formulation
  void zeroLoad(void);
  int addLoad(ElementalLoad* theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector& accel);

  // get converged residual
  const Vector& getResistingForce(void);

  // get converged residual with inertia terms
  const Vector& getResistingForceIncInertia(void);

  // public methods for element output for parallel and database processing
  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel,
      FEM_ObjectBroker& theBroker);

  // display element graphically
  int displaySelf(Renderer&, int mode, float fact, const char** displayModes = 0, int numModes = 0);

  // print out element data
  void Print(OPS_Stream& s, int flag = 0);

  // implemented to print into file
  const char* getClassType(void) const { return "BeamColumnJoint2d"; };

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& eleInformation);

  int setParameter(char** argv, int argc, Information& info);
  int updateParameter(int parameterID, Information& info);

  void updateDir(const Vector& x, const Vector& y);

protected:

private:
    int elemType;
    // private methods
    void   setUp(int Nd1, int Nd2, const Vector& x, const Vector& y);
   // void   checkDirection(ID& dir) const;

    void   setTran1d(int e, int n);
    double computeCurrentStrain1d(int mat, const Vector& diff) const;

    // material info
    UniaxialMaterial** MaterialPtr;  // pointer to the 13 different materials

    // node info
    ID  connectedExternalNodes;   // contains the tags of the end nodes
    Node* nodePtr[2];             // pointers to four nodes
    int dimension;                      // = 1, 2, or 3 dimensions
    int nDOF;	                        // number of dof for ZeroLength
    Matrix transformation;		// transformation matrix for orientation


    ID* dir1d;     	   // array of directions 0-5 for 1d materials
    Matrix* t1d; 	   // hold the transformation matrix

    // vector pointers to initial disp and vel if present
    Vector* d0;
    Vector* v0;
    Vector deform;


    Matrix* theMatrix;  // pointer to objects matrix (a class wide Matrix)
    Vector* theVector;  // pointer to objects vector (a class wide Vector)
    Vector* theLoad;    // pointer to the load vector 

    int mInitialize;  // tag to fix bug in recvSelf/setDomain when using database command

};

#endif
