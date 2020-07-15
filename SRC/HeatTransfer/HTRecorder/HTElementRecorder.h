/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
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
** ****************************************************************** */
                                                                        
// $Revision: 1.14 $
// $Date: 2009-04-14 21:14:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/HTElementRecorder.h,v $
                                                                        
                                                                        
#ifndef HTElementRecorder_h
#define HTElementRecorder_h


// Written: fmk 
// Created: 09/99
// Revision: A
//
// Description: This file contains the class definition for HTElementRecorder.
// A HTElementRecorder is used to obtain a response from an element during 
// the analysis.
//
// What: "@(#) HTElementRecorder.h, revA"

#include <HTRecorder.h>
#include <Information.h>
#include <ID.h>

class HeatTrsansferDomain;
class Vector;
class Matrix;
class HeatTransferElement;
class Response;
class FE_Datastore;

class HTElementRecorder: public HTRecorder
{
  public:
    HTElementRecorder();
    HTElementRecorder(int tag, const ID *eleID, 
		    char **argv, 
		    int argc,
		    bool echoTime, 
		    HeatTransferDomain &theDomain, 
		    OPS_Stream &theOutputHandler,
		    double deltaT = 0.0,
		    const ID *dof = 0);

    ~HTElementRecorder();

    int record( double timeStamp);
    int restart(void);    

    int setDomain(HeatTransferDomain &theDomain);
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    
  protected:

    
  private:	
    int initialize(void);
    int numDOF;
    int numEle;

    ID* eleID;
    ID* dof;

    Response **theResponses;

    HeatTransferDomain *theDomain;
    OPS_Stream *theOutputHandler;

    bool echoTimeFlag;             // flag indicating if pseudo time also printed

    double deltaT;
    double nextTimeStampToRecord;

    Vector *data;
    bool initializationDone;
    char **responseArgs;
    int numArgs;

    int addColumnInfo;

};


#endif
