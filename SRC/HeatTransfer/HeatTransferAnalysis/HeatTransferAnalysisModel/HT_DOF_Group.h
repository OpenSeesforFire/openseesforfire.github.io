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

// This class is a modified version of DOF_Group 
// for the heat transfer module
// Modified by: Yaqiang Jiang(y.jiang@ed.ac.uk)

// DOF_Group
// Written: fmk 
// Created: 11/96
// Revision: A

#ifndef HT_DOF_Group_h
#define HT_DOF_Group_h

#include <ID.h>
#include <TaggedObject.h>

class HeatTransferNode;
class Vector;
class Matrix;
class HT_TransientIntegrator;  // need to be confirmed after the two classes have been developed
class HeatTransferIntegrator;

class HT_DOF_Group: public TaggedObject
{
  public:
    HT_DOF_Group(int tag, HeatTransferNode* myNode);
    //HT_DOF_Group(int tag, int ndof);  // ignore this constructor at the moment 
    virtual ~HT_DOF_Group();    

    virtual void setID(int dof, int value);
    virtual void setID(const ID& values);
    virtual const ID& getID(void) const;
    virtual int doneID(void);    

    virtual int getNodeTag(void) const;
    virtual int getNumDOF(void) const;    
    virtual int getNumFreeDOF(void) const;
    virtual int getNumConstrainedDOF(void) const;

    // methods to obtain committed responses from the nodes
    virtual const Vector& getCommittedTemp(void);
    virtual const Vector& getCommittedTdot(void);
    
    // methods to update the trial response at the nodes
    virtual void setNodeTemp(const Vector& T);
    virtual void setNodeTot(const Vector& Tdot);

    virtual void incrNodeTemp(const Vector& T);
    virtual void incrNodeTdot(const Vector& Tdot);

    virtual void resetNodePtr(void);
	
	virtual void  Print(OPS_Stream&, int = 0) {return;};
  
   protected:  
    HeatTransferNode* myNode;
	Vector* storage; 
    
  private:
    // private variables - a copy for each object of the class        
    ID 	myID;
    int numDOF;
	static Vector errVect;
    static Vector** theVectors; // may be changed and removed from constructor, since this one and  
	                            // 'storage' are only used by HT_DOF_Group::setNodeTemp method
    static int numDOFs;  // number of objects    
};

#endif
