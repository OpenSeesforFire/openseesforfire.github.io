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

// This class is a modified version of FE_Element 
// for the heat transfer module
// Modified by: Yaqiang Jiang(y.jiang@ed.ac.uk)

// FE_Element
// Written: fmk 
// Created: 11/96
// Revision: A

#ifndef HT_FE_Element_h
#define HT_FE_Element_h

#include <ID.h>
#include <Matrix.h>
#include <Vector.h>
#include <TaggedObject.h>

class HeatTransferIntegrator;
class HeatTransferElement;
class HeatTransferAnalysis;
class HT_AnalysisModel;

class HT_FE_Element: public TaggedObject
{
  public:
    HT_FE_Element(int tag, HeatTransferElement* theElement);
    HT_FE_Element(int tag, int numDOF_Group, int ndof);
    virtual ~HT_FE_Element();    

    // public methods for setting/obtaining mapping information
    virtual const ID& getDOFtags(void) const;
    virtual const ID& getID(void) const;
    void setModel(HT_AnalysisModel& theModel);
    virtual int setID(void);
    
    // methods to form and obtain the tangent and residual
    virtual const Matrix& getTangent(HeatTransferIntegrator* theIntegrator);
    virtual const Vector& getResidual(HeatTransferIntegrator* theIntegrator);

    // methods to allow integrator to build tangent
    virtual void  zeroTangent(void);
    virtual void  addMkAndMqToTang(double fact = 1.0);
	virtual void  addMcToTang(double fact = 1.0);
    
    // methods to allow integrator to build residual    
    virtual void  zeroResidual(void);    
    virtual void  addQkAndQqToResidual(double fact = 1.0);
    virtual void  addQcToResidual(double fact = 1.0);  
	virtual void  addQtotalToResidual(double fact = 1.0);  
	virtual void  Print(OPS_Stream&, int = 0) {return;};
    
  protected:  
    // protected variables - a copy for each object of the class        
    ID myDOF_Groups;
    ID myID;

  private:
    // private variables - a copy for each object of the class    
    int numDOF;
    HT_AnalysisModel* theModel;
    HeatTransferElement* myEle;
    Vector* theResidual;
    Matrix* theTangent;
    HeatTransferIntegrator* theIntegrator; // need for Subdomain
    
    // static variables - single copy for all objects of the class	
    static Matrix errMatrix;
    static Vector errVector;
    static Matrix **theMatrices; // array of pointers to class wide matrices
    static Vector **theVectors;  // array of pointers to class widde vectors
    static int numFEs;           // number of objects
};

#endif
