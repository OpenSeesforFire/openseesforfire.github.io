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

#include <HT_TransientIntegrator.h>
#include <HT_FE_Element.h>
#include <LinearSOE.h>
#include <HT_AnalysisModel.h>
#include <Vector.h>
#include <HT_FE_EleIter.h>


HT_TransientIntegrator::HT_TransientIntegrator()
{

}


HT_TransientIntegrator::~HT_TransientIntegrator()
{

}

    
int
HT_TransientIntegrator::formEleResidual(HT_FE_Element* theEle)
{
    theEle->zeroResidual();
    theEle->addQtotalToResidual();
    return 0;
}    
 
