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

// This class is a modified version of Node 
// for the heat transfer module
// Modified by: Yaqiang Jiang(y.jiang@ed.ac.uk)

// Node
// Written: fmk 
// Created: 11/96
// Revision: A

#ifndef HeatTransferNode_h
#define HeatTransferNode_h

#include <HeatTransferDomainComponent.h>
#include <TemperatureBC.h>

class HeatTransferElement;
class Vector;
class Matrix;
//class Channel;
//class Renderer;
class HT_DOF_Group;
//class TemperatureBC;

class HeatTransferNode : public HeatTransferDomainComponent
{
  public:
    // constructors
    HeatTransferNode(int tag, int ndof, double Crd1);
    HeatTransferNode(int tag, int ndof, double Crd1, double Crd2);
    HeatTransferNode(int tag, int ndof, double Crd1, double Crd2, double Crd3);      

    // destructor
    virtual ~HeatTransferNode();

    // public methods dealing with the DOF at the node
    virtual int  getNumberDOF(void) const;    
    virtual void setDOF_GroupPtr(HT_DOF_Group* theDOF_Grp);
    virtual HT_DOF_Group* getDOF_GroupPtr(void);

    // public methods for obtaining the nodal coordinates
    virtual const Vector& getCrds(void) const;

    // public methods for obtaining committed and trial 
    // response quantities of the node
    virtual const Vector& getTemperature(void);
    virtual const Vector& getTdot(void);  // Return temperature derivative wrt time
    virtual const Vector& getIncrTemperature(void);
    virtual const Vector& getIncrDeltaTemperature(void);

    virtual const Vector& getTrialTemperature(void);    
    virtual const Vector& getTrialTdot(void);       

    // public methods for updating the trial response quantities 
    virtual int setTrialTemperature(const Vector& );    
    virtual int setTrialTdot(const Vector& );    
	virtual int setResponse(double value, int dof_number);

    virtual int incrTrialTemperature(const Vector& );    
    virtual int incrTrialTdot(const Vector& );        

    // public methods dealing with the committed state of the node
    virtual int commitState();
    virtual int revertToLastCommit();    
    virtual int revertToStart();     
	//virtual void  Print(OPS_Stream&, int = 0) {return;};
    
	// need to overrite the commitTemp
	//friend int TemperatureBC::applyTemperatureBC(double );    

  protected:

  private:
    // priavte methods used to create the Vector objects 
    // for the committed and trial response quantaties.
    int createTemp(void);
    int createTdot(void);

    // private data associated with each node object
    int numberDOF;                    // number of dof at Node
    HT_DOF_Group* theDOF_GroupPtr;  // pointer to associated HeatTransferDOF_Group
    Vector* Crd;                      // original nodal coords
    Vector* commitTemp;  // commited quantities
	Vector* commitTdot; 
    Vector* trialTemp;  // trial quantities
	Vector* trialTdot;
    Vector* incrTemp;  // the difference between current trial and last committed temperature vector  
    Vector* incrDeltaTemp;  // the difference between current and last trial temperature vector
    
    double *T, *Tdot; // double arrays holding the temperature and its derivative 
};

#endif
