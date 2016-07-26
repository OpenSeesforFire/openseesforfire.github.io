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

#ifndef HTRecorder_h
#define HTRecorder_h

class HeatTransferDomain;
#include <MovableObject.h>
#include <TaggedObject.h>



class HTRecorder: public TaggedObject
{
    public:
    HTRecorder(int tag);
  
    virtual ~HTRecorder();

    virtual int record(double timeStamp)=0;

    virtual int domainChanged(void);    
    virtual int setDomain(HeatTransferDomain &theDomain);
	virtual void  Print(OPS_Stream&, int flag = 0);

    protected:

    private:	
	  static int lastRecorderTag;
};


#endif