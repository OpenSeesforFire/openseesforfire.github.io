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

// This class is a modified version of IncrementalIntegrator 
// for the heat transfer module
// Modified by: Yaqiang Jiang(y.jiang@ed.ac.uk)

// IncrementalIntegrator
// Written: fmk 
// Created: Tue Sept 17 15:54:47: 1996
// Revision: A

#ifndef HeatTransferIntegrator_h
#define HeatTransferIntegrator_h

#include <HeatTransferIntegrator.h>

class LinearSOE;
class HT_AnalysisModel;
class HT_ConvergenceTest;
class HT_FE_Element;
class HT_DOF_Group;
class Vector;

class HeatTransferIntegrator
{
    public:
		HeatTransferIntegrator();
		virtual ~HeatTransferIntegrator();

		virtual void setLinks(HT_AnalysisModel& theModel,
			                  LinearSOE& theSOE,
			                  HT_ConvergenceTest* theTest);

		virtual int domainChanged(void);

		// methods to set up the system of equations
		virtual int formTangent();    
		virtual int formUnbalance(void);

		// pure virtual methods to define the FE_ELe and DOF_Group contributions
		virtual int formEleTangent(HT_FE_Element *theEle) = 0;
		virtual int formEleResidual(HT_FE_Element *theEle) = 0;

		// methods to update the domain
		virtual int newStep(double delta_t);
		virtual int update(const Vector& deltaT) = 0;
		virtual int commit(void);
		virtual int revertToLastStep(void);
		virtual int initialize(void);
     
    protected:
		LinearSOE* getLinearSOE(void) const;
		HT_AnalysisModel* getModel(void) const;
		HT_ConvergenceTest* getConvergenceTest(void) const;

		virtual int formElementResidual(void);            
		int statusFlag;

    private:
		LinearSOE* the_soe;
		HT_AnalysisModel* the_model;
		HT_ConvergenceTest* the_test;
};

#endif
