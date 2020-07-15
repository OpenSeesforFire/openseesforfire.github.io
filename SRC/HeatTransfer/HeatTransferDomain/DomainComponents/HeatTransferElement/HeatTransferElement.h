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

// This class is a modified version of Element 
// for the heat transfer module
// Modified by: Yaqiang Jiang(y.jiang@ed.ac.uk)

// Element
// Written: fmk 
// Created: 11/96
// Revision: A

#ifndef HeatTransferElement_h
#define HeatTransferElement_h
#include <list>
#include <ID.h>
#include <HeatTransferDomainComponent.h>
#include <Convection.h>
#include <Radiation.h>

class Matrix;
class Vector;
class Information;
class Response;
class HeatTransferNode;
class PrescribedSurfFlux;

class HeatTransferElement : public HeatTransferDomainComponent
{
  public:
    HeatTransferElement(int tag);    
    virtual ~HeatTransferElement();

    // methods dealing with nodes and number of external dof
    virtual int getNumExternalNodes(void) const = 0;
    virtual const ID& getExternalNodes(void)  = 0;	
    virtual HeatTransferNode** getNodePtrs(void)  = 0;	
    virtual int getNumDOF(void) = 0;
	virtual const ID& getNodesOnFace(int faceTag) = 0;
	virtual int getNumNodesperFace(void)=0;
	// methods dealing with committed state and update
	virtual int commitState(void); // called when a converged solution has been obtained for a time step   
	virtual int revertToLastCommit(void) = 0; // called when the solutionAlgorithm has failed to converge 
	                                          // to a solution        
	virtual int revertToStart(void); // called when model is reset to initial conditions               
	virtual int update(void); // called when a new trial step has been set at the nodes

	// methods to return the current linearized tangent matrices
	virtual const Matrix& getConductionTangent(void) = 0; // K_k
	virtual const Matrix& getCapacityTangent(void) = 0; // K_cp   
	virtual const Matrix& getRadiationTangent(void) = 0; // K_qr
	// the convection heat transfer coefficient hc is also a function of the solid surface temperature Ts
	// see p571 in Incropera's book(6ed)
	virtual const Matrix& getConvectionTangent(void) = 0; // K_qc

    // methods for applying boundaryFluxes
    virtual void zeroFlux(void);	
    virtual int addPrecribedSurfFlux(PrescribedSurfFlux* theFlux, double factor) = 0;
	virtual void applyConvection(Convection* theConvection, double factor);
	virtual void applyRadiation(Radiation* theRadiation, double factor);
	virtual bool hasConvection() const; // used by HeatTransferFE_Element to decide
	                                    // if ConvectionBC exists
	virtual bool hasRadiation() const; 
    
	virtual void clearAllFluxBCs();
	virtual void removeConvection(Convection* theConvection); // remove a specific convection bc
	virtual void removeRadiation(Radiation* theRadiation);
	//void removeConvection(int tag);
	//void removeRadiation(int tag);

    // methods for calculating flux vectors
	virtual const Vector& get_Q_Transient() = 0;  // computeQcp(void); getTransient 
    virtual const Vector& get_Q_Conduction() = 0; // computeQk(void);  getConductionFlux  
    virtual const Vector& get_Q_Radiation() = 0;  // computeQqr(void); getRadiation  
	virtual const Vector& get_Q_Convection() = 0;  // computeQqc(void); getConvection  
    
	// method for obtaining residual at element level
    //virtual const Vector &getEleResidual(void) = 0;
	virtual Response* setResponse(const char** argv, int argc,OPS_Stream& theHandler);
	virtual int getResponse(int responseID, Information& eleInformation);

  protected:
	bool convecFlag, radFlag, pFlag; // flags to indicate BCs' existence
	std::list<Convection*> theConvectionBCs;
	std::list<Radiation*> theRadiationBCs;
	std::list<PrescribedSurfFlux*> thePrescribedFluxBCs;
  private:

};


#endif
