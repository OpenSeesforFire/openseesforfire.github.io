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
// Modified by Liming Jiang(liming.jiang@ed.ac.uk)
//

#ifndef HTNodeRecorder_h
#define HTNodeRecorder_h

#include <HTRecorder.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>
#include <PathTimeSeriesThermal.h> 

class HeatTransferDomain;
class HeatTransferNode;

class HTNodeRecorder: public HTRecorder
{
    public:
	HTNodeRecorder(int tag);
  HTNodeRecorder(int tag, const ID* theNodes,
					HeatTransferDomain& theDomain); 
					
  HTNodeRecorder(int tag, const ID* theNodes,
					HeatTransferDomain& theDomain,
					OPS_Stream &theOutputHandle);
  
  HTNodeRecorder(int tag, const ID* theNodes,
                 HeatTransferDomain& theDomain,
                 PathTimeSeriesThermal* thePathSeries);
  
  HTNodeRecorder(int tag, const ID* theNodes,
                 HeatTransferDomain& theDomain,
                 OPS_Stream &theOutputHandle,PathTimeSeriesThermal* thePathSeries);
    
    ~HTNodeRecorder();

    int record(double timeStamp);
  
    int setDomain(HeatTransferDomain &theDomain);
	void  Print(OPS_Stream&, int = 0) {return;};

    protected:

    private:	
	int initialize(void);
	
	Vector response;

	ID* theNodalTags;
	HeatTransferNode** theNodes;

	HeatTransferDomain* theDomain;

	OPS_Stream *theOutputHandler;
  
  PathTimeSeriesThermal* ThePathSeries;

	bool initializationDone;
	bool initialRecording;
	int numValidNodes;

};


#endif