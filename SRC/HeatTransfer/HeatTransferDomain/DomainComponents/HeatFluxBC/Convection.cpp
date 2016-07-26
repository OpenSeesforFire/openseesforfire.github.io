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

#include <Convection.h>
#include <HeatTransferElement.h>

//const int Convection::type_tag;

Convection::Convection(int tag, int eleTag, int fTag, double h, double T = 20.0)
:HeatFluxBC(tag, eleTag, fTag), hc(h), Ta(T)
{

}


Convection::~Convection()
{

}


void 
Convection::applyFluxBC(double factor)
{
    if (theElement != 0)
		theElement->applyConvection(this, factor);
}


double 
Convection::getParameter(void) const
{
    return hc;
}


double 
Convection::getSurroundingTemp(void) const
{
    return Ta;
}


void
Convection::setSurroundingTemp(double T)
{
    Ta = T;
}