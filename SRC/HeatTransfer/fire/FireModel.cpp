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

#include <FireModel.h>


FireModel::FireModel(int tag, int fireTypeTag):TaggedObject(tag),
FireTypeTag(fireTypeTag),the_domain(0)
{

}


FireModel::~FireModel()
{

}


void
FireModel::setDomain(HeatTransferDomain* theDomain)
{
    the_domain = theDomain;
}


int 
FireModel::getFireTypeTag()
{
	return FireTypeTag;
}

double
FireModel::getFirePars(int parTag)
{
	return 0.0;
}

int 
FireModel::setFirePars(double time, const Vector& firePars)
{
	return 0;
}

double
FireModel::getFireOut(double time, const Vector& locs)
{
	return 0;
}

void 
FireModel::Print(OPS_Stream& s, int i)
{
	opserr << "No fire model information to print" << endln;
}