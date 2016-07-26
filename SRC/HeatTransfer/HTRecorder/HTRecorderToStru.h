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
**   Liming Jiang(liming.jiang@ed.ac.uk)                              **
** ****************************************************************** */

//
// Created by Liming Jiang (y.jiang@ed.ac.uk)
//

#ifndef HTRecorderToStru_h
#define HTRecorderToStru_h

#include <HTRecorder.h>

#include <Vector.h>
#include <ID.h>
#include <Matrix.h>


class Matrix;
class HeatTransferDomain;
class HeatTransferNode;

class HTRecorderToStru: public HTRecorder
{
    public:
    HTRecorderToStru(int tag);
	
    HTRecorderToStru(int tag, const Vector& theCrds, 
		       HeatTransferDomain& theDomain, double tolerance=0.00001); 
	
    HTRecorderToStru(int tag, const Vector& theCrds, 
		       HeatTransferDomain& theDomain, OPS_Stream &theOutputHandle,
			   double tolerance=0.00001); 
    
    ~HTRecorderToStru();

    int record(double timeStamp);
    int setDomain(HeatTransferDomain &theDomain);
	void  Print(OPS_Stream&, int = 0) {return;};

    protected:

    private:	
	int initialize(void);
	Vector response;

	ID theNodeTags;
	ID theNodeTags1;
	HeatTransferNode** theNodes;
	HeatTransferNode** theNodes1;

	HeatTransferDomain* theDomain;
	OPS_Stream *theOutputHandler;
	
	bool initializationDone;
	bool initialRecording;
	int numValidNodes;
	double Tolerance;
	double xCrd;
	Vector yCrds;

};


#endif