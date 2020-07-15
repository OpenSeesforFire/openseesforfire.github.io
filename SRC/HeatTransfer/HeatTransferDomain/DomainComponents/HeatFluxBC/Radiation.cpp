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

#include <Radiation.h>
#include <HeatTransferElement.h>

//const int Radiation::type_tag;

Radiation::Radiation(int tag, int eleTag, int fTag, double epsilon, 
					 double sigma, double alpha, double qir)
:HeatFluxBC(tag, eleTag, fTag), emissivity(epsilon), 
 absorptivity(alpha), BZM_const(sigma) ,irradiation(qir), setIndex(false)
{

}


Radiation::~Radiation()
{

}


void 
Radiation::applyFluxBC(double factor)
{
    if (!setIndex) {
        if (theElement != 0)
            theElement->applyRadiation(this, factor);
        setIndex = true;
    }
   
}


double 
Radiation::getEmissionParameter(void) const
{
    return emissivity * BZM_const;
}


double 
Radiation::getBLZMConstant(void) const
{
    return BZM_const;
}


double 
Radiation::getIrradiation(void) const
{
    return irradiation;
}


double
Radiation::getAbsorptivity(void) const
{
    return absorptivity;
}

void
Radiation::setIrradiation(double qir)
{
	irradiation = qir;
}