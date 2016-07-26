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
// Modified by Liming Jiang (liming.jiang@ed.ac.uk)
//

#include <HTRecorder.h>
#include <OPS_Globals.h>

int HTRecorder::lastRecorderTag(0);


HTRecorder::HTRecorder(int tag)
:TaggedObject(tag)
{
    //lastRecorderTag++;
}





HTRecorder::~HTRecorder()
{
  
}

int
HTRecorder::domainChanged(void)
{
   return 0;
}


int 
HTRecorder::setDomain(HeatTransferDomain& theDom)
{
    return 0;
}


void
HTRecorder::Print(OPS_Stream&, int flag)
{
	return;
}
