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

#include <PrescribedSurfFlux.h>
#include <HeatTransferElement.h>
#include <Vector.h>

//const int PrescribedSurfFlux::type_tag;

PrescribedSurfFlux::PrescribedSurfFlux(int tag, int eleTag, int fTag, 
									   int numNPF, const Vector& data)
:HeatFluxBC(tag, eleTag, fTag), NPF(numNPF), flux_data(0),FireType(0)
{
    int size = data.Size();
	if (numNPF != size)
		{
		opserr << "PrescribedSurfFlux::PrescribedSurfFlux() -- number of nodes on face " 
			<< fTag << " doesn't match the size of Vector.\n";
		exit(-1);	
		}

	if (flux_data != 0)
		delete flux_data;

	flux_data = new Vector(NPF);
	if (flux_data == 0) {
		opserr << "PrescribedSurfFlux::PrescribedSurfFlux() - out of memeory\n";
		exit(-1);
		}

	for (int i = 0; i < NPF; i++)
		{
		(*flux_data)(i) = data(i);
		}
}


PrescribedSurfFlux::PrescribedSurfFlux(int tag, int eleTag, int fTag, int numNPF,int fireType)
:HeatFluxBC(tag, eleTag, fTag), NPF(numNPF), flux_data(0),FireType(fireType)
{
	if (flux_data != 0)

		delete flux_data;

	flux_data = new Vector(NPF);
	if (flux_data == 0) {
		opserr << "PrescribedSurfFlux::PrescribedSurfFlux() - out of memeory\n";
		exit(-1);
		}
}


PrescribedSurfFlux::~PrescribedSurfFlux()
{

}


void 
PrescribedSurfFlux::applyFluxBC(double factor)
{
    if (theElement != 0)
		theElement->addPrecribedSurfFlux(this, factor);
}


const Vector& 
PrescribedSurfFlux::getData(void) const
{
    return *flux_data;
}


void 
PrescribedSurfFlux::setData(const Vector& theData)
{
    // first check the size
    int size1 = theData.Size();
	int size2 = flux_data->Size();
	if (size1 != size2) {
		opserr << "PrescribedSurfFlux::setData() - Vector sizes aren't matching.\n";
		exit(-1);
		} else {
			*flux_data = theData;
		}
}


int 
PrescribedSurfFlux::getFireType(void)
{
	return FireType;
}