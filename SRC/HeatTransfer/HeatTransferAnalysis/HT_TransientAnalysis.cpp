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
// Note: This class was adapted from DirectIntegrationAnalysis
#include <HT_TransientAnalysis.h>
#include <HT_SolutionAlgorithm.h>
#include <HT_AnalysisModel.h>
#include <LinearSOE.h>
#include <HT_DOF_Numberer.h>
#include <TemperatureBCHandler.h>
#include <HT_ConvergenceTest.h>
#include <HT_TransientIntegrator.h>
#include <HeatTransferDomain.h>
#include <HT_FE_Element.h>
#include <HT_DOF_Group.h>
#include <HT_FE_EleIter.h>
#include <HT_NodeIter.h>
#include <HT_DOF_GrpIter.h>
#include <Matrix.h>
#include <ID.h>
#include <Graph.h>
#include <HeatTransferNode.h>
#include <fstream>
#include <elementAPI.h>


HT_TransientAnalysis::HT_TransientAnalysis(HeatTransferDomain& theDomain,
						                   TemperatureBCHandler& theHandler,
						                   HT_DOF_Numberer& theNumberer,
						                   HT_AnalysisModel& theModel,
						                   HT_SolutionAlgorithm& theSolnAlgo,		   
						                   LinearSOE& theLinSOE,
						                   HT_TransientIntegrator& theTransientIntegrator,
						                   HT_ConvergenceTest* theConvergenceTest)
:HeatTransferAnalysis(theDomain), 
 temp_bc_handler(&theHandler),
 dof_numberer(&theNumberer), 
 analysis_model(&theModel), 
 solution_algorithm(&theSolnAlgo), 
 linear_soe(&theLinSOE), 
 transient_integrator(&theTransientIntegrator), 
 convergence_test(theConvergenceTest),
 domainStamp(0)
{
    // first we set up the links needed by the elements in the 
    // aggregation
    analysis_model->setLinks(theDomain, theHandler);
    temp_bc_handler->setLinks(theDomain, theModel, theTransientIntegrator);
    dof_numberer->setLinks(theModel);
    transient_integrator->setLinks(theModel, theLinSOE, convergence_test);
    solution_algorithm->setLinks(theModel, theTransientIntegrator, theLinSOE, convergence_test);

    if (convergence_test != 0)
		solution_algorithm->setConvergenceTest(convergence_test);
    else
		convergence_test = solution_algorithm->getConvergenceTest();
}    


HT_TransientAnalysis::~HT_TransientAnalysis()
{
  // we don't invoke the destructors in case user switching
  // from a static to a direct integration analysis 
  // clearAll() must be invoked if user wishes to invoke destructor
}    


void
HT_TransientAnalysis::clearAll(void)
{
  // invoke the destructor on all the objects in the aggregation

	if (analysis_model != 0)     
		delete analysis_model;
    if (temp_bc_handler != 0) 
		delete temp_bc_handler;
    if (dof_numberer != 0)      
		delete dof_numberer;
    if (transient_integrator != 0) 
		delete transient_integrator;
    if (solution_algorithm != 0)  
		delete solution_algorithm;
    if (linear_soe != 0)
		delete linear_soe;
    if (convergence_test != 0)
		delete convergence_test;

    analysis_model = 0;
    temp_bc_handler = 0;
    dof_numberer = 0;
    transient_integrator = 0;
    solution_algorithm = 0;
    linear_soe = 0;
    convergence_test = 0;
}    


int 
HT_TransientAnalysis::initialize(void)
{
    HeatTransferDomain* the_domain = this->getDomainPtr();

    // check if domain has undergone change
    int stamp = the_domain->hasDomainChanged();

    if (stamp != domainStamp) {
		domainStamp = stamp;	
		if (this->domainChanged() < 0) {
			opserr << "HT_TransientAnalysis::initialize() - domainChanged() failed\n";
			return -1;
			}	
		}

    if (transient_integrator->initialize() < 0) {
		opserr << "HT_TransientAnalysis::initialize() - integrator initialize() failed\n";
		return -2;
		} else
			transient_integrator->commit();

    return 0;
}


int 
HT_TransientAnalysis::analyze(int numSteps, double dT,double& lastTime, double monitortime)
{
    int result = 0; int monitor = 0;  int displayTag = 1;
	int laststep = lastTime / dT;
    HeatTransferDomain* the_domain = this->getDomainPtr();
	
    for (int i = laststep; i < numSteps; i++) {
#ifdef _DEBUG
		//opserr<<"Current time: "<<the_domain->getCurrentTime()<<endln;
#endif
        if(displayTag ==1)
            opserr << "Current time: " << the_domain->getCurrentTime() << endln;

		if (analysis_model->analysisStep(dT) < 0) {
			opserr << "HT_TransientAnalysis::analyze() - the HT_AnalysisModel failed";
			opserr << " at time " << the_domain->getCurrentTime() << endln;
			the_domain->revertToLastCommit();
			return -2;
			}

		// check if domain has undergone change
		int stamp = the_domain->hasDomainChanged();
		if (stamp != domainStamp) {
			domainStamp = stamp;	
			if (this->domainChanged() < 0) {
				opserr << "HT_TransientAnalysis::analyze() - domainChanged() failed\n";
				return -1;
				}	
			}

		if (transient_integrator->newStep(dT) < 0) {
			opserr << "HT_TransientAnalysis::analyze() - the HT_TransientIntegrator failed";
			opserr << " at time " << the_domain->getCurrentTime() << endln;
			the_domain->revertToLastCommit();
			return -2;
			}

		result = solution_algorithm->solveCurrentStep();
		if (result < 0) {
			opserr << "HT_TransientAnalysis::analyze() - the HT_SolutionAlgorithm failed";
			opserr << " at time " << the_domain->getCurrentTime() << endln;
			the_domain->revertToLastCommit();	    
			transient_integrator->revertToLastStep();
			return -3;
			}    

		result = transient_integrator->commit();
		if (result < 0) {
			opserr << "DirectIntegrationAnalysis::analyze() - ";
			opserr << "the HT_TransientIntegrator failed to commit";
			opserr << " at time " << the_domain->getCurrentTime() << endln;
			the_domain->revertToLastCommit();	    
			transient_integrator->revertToLastStep();
			return -4;
			} 

        //HeatTransferDomain* the_domain = this->getDomainPtr();
		monitor=i;
		if (the_domain->getCurrentTime() == monitortime) {
			lastTime = monitortime;
			return monitor;
		}
	
	}    


#ifdef _DEBUG
	return monitor;
#else
	return result;
#endif

}

/*
int
HT_TransientAnalysis::continue (int numSteps, double dT, double monitortime)
{
	int result = 0; int monitor = monitortime/dT;

	HeatTransferDomain* the_domain = this->getDomainPtr();
	this->analyze()

	return 0;

}

*/
int
HT_TransientAnalysis::domainChanged(void)
{
    HeatTransferDomain* the_domain = this->getDomainPtr();
    int stamp = the_domain->hasDomainChanged();
    domainStamp = stamp;

    analysis_model->clearAll();    
    temp_bc_handler->clearAll();
    
    // now we invoke handle() on the constraint handler which
    // causes the creation of FE_Element and DOF_Group objects
    // and their addition to the AnalysisModel.
    temp_bc_handler->handle();

    // we now invoke number() on the numberer which causes
    // equation numbers to be assigned to all the DOFs in the
    // AnalysisModel.
    dof_numberer->numberDOF();

    // temp_bc_handler->doneNumberingDOF();  // see comments on p6 DirectIntegrationAnalysis.cpp 

    // we invoke setGraph() on the LinearSOE which
    // causes that object to determine its size
    Graph& the_graph = analysis_model->getDOFGraph();

    int result = linear_soe->setSize(the_graph);
    if (result < 0) {
		opserr << "HT_TransientAnalysis::handle() - ";
		opserr << "LinearSOE::setSize() failed";
		return -3;
		}	    

    // we invoke domainChange() on the integrator and algorithm
    transient_integrator->domainChanged();
    solution_algorithm->domainChanged();

    return 0;
}    


int 
HT_TransientAnalysis::setNumberer(HT_DOF_Numberer& new_numberer) 
{
    // invoke the destructor on the old one
    if (dof_numberer != 0)
		delete dof_numberer;

    // first set the links needed by the Algorithm
    dof_numberer = &new_numberer;
    dof_numberer->setLinks(*analysis_model);

    // invoke domainChanged() either indirectly or directly
    domainStamp = 0;
    return 0;
}


int 
HT_TransientAnalysis::setAlgorithm(HT_SolutionAlgorithm& new_algorithm) 
{
    // invoke the destructor on the old one
    if (solution_algorithm != 0)
		delete solution_algorithm;
  
    // first set the links needed by the Algorithm
    solution_algorithm = &new_algorithm;

    if (analysis_model != 0 && transient_integrator != 0 && linear_soe != 0)
		solution_algorithm->setLinks(*analysis_model, *transient_integrator, 
		*linear_soe, convergence_test);
    // invoke domainChanged() either indirectly or directly
    // domainStamp = 0;
    if (domainStamp != 0)
		solution_algorithm->domainChanged();

    return 0;
}


int 
HT_TransientAnalysis::setIntegrator(HT_TransientIntegrator& new_integrator)
{
    // set the links needed by the other objects in the aggregation
    HeatTransferDomain *the_domain = this->getDomainPtr();
    transient_integrator = &new_integrator;
    transient_integrator->setLinks(*analysis_model, *linear_soe, convergence_test);
    temp_bc_handler->setLinks(*the_domain, *analysis_model, *transient_integrator);
    solution_algorithm->setLinks(*analysis_model, *transient_integrator, *linear_soe, 
		                         convergence_test);

    // cause domainChanged to be invoked on next analyze
    //  domainStamp = 0;
    if (domainStamp != 0)
		transient_integrator->domainChanged();
   
    return 0;
}


int 
HT_TransientAnalysis::setLinearSOE(LinearSOE& new_soe)
{
    // invoke the destructor on the old one
    if (linear_soe != 0)
		delete linear_soe;

    // set the links needed by the other objects in the aggregation
    linear_soe = &new_soe;
    transient_integrator->setLinks(*analysis_model,*linear_soe, convergence_test);
    solution_algorithm->setLinks(*analysis_model, *transient_integrator, *linear_soe,
	                             convergence_test);
    // linear_soe->setLinks(*analysis_model);  // NOTE: SEEMS UNUSEFUL

    // cause domainChanged to be invoked on next analyze
    domainStamp = 0;
  
    return 0;
}


int 
HT_TransientAnalysis::setConvergenceTest(HT_ConvergenceTest& new_test)
{
    // invoke the destructor on the old one
    if (convergence_test != 0)
		delete convergence_test;
  
    // set the links needed by the other objects in the aggregation
    convergence_test = &new_test;

    if (transient_integrator != 0)
        transient_integrator->setLinks(*analysis_model, *linear_soe, convergence_test);

    if (solution_algorithm != 0)
		solution_algorithm->setConvergenceTest(convergence_test);
  
    return 0;
}


int
HT_TransientAnalysis::checkDomainChange(void)
{
    HeatTransferDomain* the_domain = this->getDomainPtr();

    // check if domain has undergone change
    int stamp = the_domain->hasDomainChanged();
    if (stamp != domainStamp) {
        domainStamp = stamp;	
        if (this->domainChanged() < 0) {
            opserr << "HT_TransientAnalysis::checkDomainChange() - domainChanged() failed\n";
            return -1;
            }	
        }

    return 0;
}


HT_SolutionAlgorithm*
HT_TransientAnalysis::getAlgorithm(void)
{
    return solution_algorithm;
}


HT_AnalysisModel*
HT_TransientAnalysis::getModel(void)
{
    return analysis_model;
}


HT_TransientIntegrator*
HT_TransientAnalysis::getIntegrator(void)
{
    return transient_integrator;
}


HT_ConvergenceTest*
HT_TransientAnalysis::getConvergenceTest(void)
{
    return convergence_test;
}




