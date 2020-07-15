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
                                                                        
// $Revision: 1.34 $
// $Date: 2009-04-30 23:25:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/HTElementRecorder.cpp,v $
                                                                        
// Written: fmk 
// Created: 09/99
//
// Description: This file contains the class implementatation of HTElementRecorder.
//
// What: "@(#) HTElementRecorder.C, revA"

#include <HTElementRecorder.h>
#include <HeatTransferDomain.h>
#include <HeatTransferElement.h>
#include <HT_ElementIter.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>
#include <Response.h>

#include <OPS_Globals.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <StandardStream.h>
#include <DataFileStream.h>


#include <elementAPI.h>

#include <string.h>

/*
void*
OPS_HTElementRecorder()
{
    if (OPS_GetNumRemainingInputArgs() < 5) {
        opserr << "WARING: recorder HeatTransferElement ";
        opserr << "-ele <list elements> -file <fileName> -dT <dT> reponse";
        return 0;
    }

    const char** data = 0;
    int nargrem = 0;
    OPS_Stream *theOutputStream = 0;
    const char* filename = 0;

    const int STANDARD_STREAM = 0;
    const int DATA_STREAM = 1;

    int eMode = STANDARD_STREAM;

    bool echoTimeFlag = false;
    double dT = 0.0;
    bool doScientific = false;

    int precision = 6;

    bool closeOnWrite = false;

    const char *inetAddr = 0;
    int inetPort;

    ID elements(0, 6);
    ID dofs(0, 6);

    while (OPS_GetNumRemainingInputArgs() > 0) {

        const char* option = OPS_GetString();

        if (strcmp(option, "-time") == 0) {
            echoTimeFlag = true;
        }
        else if (strcmp(option, "-file") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                filename = OPS_GetString();
            }
            eMode = DATA_STREAM;
        }
        else if (strcmp(option, "-dT") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetDoubleInput(&num, &dT) < 0) {
                    opserr << "WARNING: failed to read dT\n";
                    return 0;
                }
            }
        }
        else if (strcmp(option, "-precision") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetIntInput(&num, &precision) < 0) {
                    opserr << "WARNING: failed to read precision\n";
                    return 0;
                }
            }
        }
        else if (strcmp(option, "-ele") == 0) {
            int numEle = 0;
            while (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                int el;
                if (OPS_GetIntInput(&num, &el) < 0) {
                    OPS_ResetCurrentInputArg(-1);
                    break;
                }
                elements[numEle++] = el;
            }
        }
        else if (strcmp(option, "-eleSet") == 0) {
            int start, end;
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetIntInput(&num, &start) < 0) {
                    opserr << "WARNING: failed to read start HeatTransferElement\n";
                    return 0;
                }
            }
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetIntInput(&num, &end) < 0) {
                    opserr << "WARNING: failed to read end HeatTransferElement\n";
                    return 0;
                }
            }
            if (start > end) {
                int swap = end;
                end = start;
                start = swap;
            }
            int numEle = 0;
            for (int i = start; i <= end; i++)
                elements[numEle++] = i;
        }
        else {
            // first unknown string then is assumed to start
            // HeatTransferElement response request
            nargrem = 1 + OPS_GetNumRemainingInputArgs();
            data = new const char *[nargrem];
            data[0] = option;
            for (int i = 1; i < nargrem; i++)
                data[i] = OPS_GetString();
        }
    }

    // data handler
    if (eMode == DATA_STREAM && filename != 0)
        theOutputStream = new DataFileStream(filename, OVERWRITE, 2, 0, closeOnWrite, precision, doScientific);
    else
        theOutputStream = new StandardStream();

    theOutputStream->setPrecision(precision);

    HeatTransferDomain* domain = OPS_GetDomain();
    if (domain == 0)
        return 0;
    HTElementRecorder* recorder = new HTElementRecorder(&elements,
        data, nargrem, echoTimeFlag, *domain, *theOutputStream,
        dT, &dofs);

    return recorder;
}

*/



HTElementRecorder::HTElementRecorder()
:HTRecorder(2),
 numEle(0), numDOF(0), eleID(0), dof(0), theResponses(0), 
 theDomain(0), theOutputHandler(0),
 echoTimeFlag(true), deltaT(0), nextTimeStampToRecord(0.0), data(0), 
 initializationDone(false), responseArgs(0), numArgs(0), addColumnInfo(0)
{

}

HTElementRecorder::HTElementRecorder(int tag, const ID *ele,
				 char **argv, 
				 int argc,
				 bool echoTime, 
				 HeatTransferDomain &theDom, 
				 OPS_Stream &theOutputHandler,
				 double dT,
				 const ID *theDOFs)
:HTRecorder(2),
 numEle(0), numDOF(0), eleID(0), dof(0), theResponses(0), 
 theDomain(&theDom), theOutputHandler(&theOutputHandler),
 echoTimeFlag(echoTime), deltaT(dT), nextTimeStampToRecord(0.0), data(0),
 initializationDone(false), responseArgs(0), numArgs(0), addColumnInfo(0)
{

  if (ele != 0) {
    numEle = ele->Size();
    eleID = new ID(*ele);
    if (eleID == 0 || eleID->Size() != numEle)
      opserr << "HTElementRecorder::HTElementRecorder() - out of memory\n";
  } 

  if (theDOFs != 0) {
    dof = new ID(*theDOFs);
    numDOF = dof->Size();
  } 

  //
  // create a copy of the response request
  //

  responseArgs = new char *[argc-5];
  if (responseArgs == 0) {
    opserr << "HTElementRecorder::HTElementRecorder() - out of memory\n";
    numEle = 0;
  }
  
  for (int i=0; i<argc-5; i++) {
    responseArgs[i] = new char[strlen(argv[i+5])+1];
    if (responseArgs[i] == 0) {
      delete [] responseArgs;
      opserr << "HTElementRecorder::HTElementRecorder() - out of memory\n";
      numEle = 0;
    }
    strcpy(responseArgs[i], argv[i+5]);
  }
  
  numArgs = argc-5;
}


HTElementRecorder::~HTElementRecorder()
{
  if (theOutputHandler != 0) {
    theOutputHandler->endTag(); // Data
    delete theOutputHandler;
  }

  //
  // invoke the destructor on the response objects
  //

  if (eleID != 0)
    delete eleID;
  

  if (theResponses != 0) {
    for (int i = 0; i < numEle; i++)
      delete theResponses[i];
    delete [] theResponses;
  }

  if (data != 0)
    delete data;
  
  // 
  // invoke destructor on response args
  //

  for (int i=0; i<numArgs; i++)
    delete [] responseArgs[i];
  delete [] responseArgs;

}


int 
HTElementRecorder::record(double timeStamp)
{
  // 
  // check that initialization has been done
  //

  if (initializationDone == false) {
    if (this->initialize() != 0) {
      opserr << "HTElementRecorder::record() - failed to initialize\n";
      return -1;
    }
  }
  
  int result = 0;
  if (deltaT == 0.0 || timeStamp >= nextTimeStampToRecord) {

    if (deltaT != 0.0) 
      nextTimeStampToRecord = timeStamp + deltaT;

    int loc = 0;
    if (echoTimeFlag == true) 
      (*data)(loc++) = timeStamp;
    
    //
    // for each HeatTransferElement if responses exist, put them in response vector
    //
    for (int i=0; i< numEle; i++) {
      if (theResponses[i] != 0) {
	// ask the HeatTransferElement for the reponse
	int res;
	if (( res = theResponses[i]->getResponse()) < 0)
	  result += res;
	else {
	  Information &eleInfo = theResponses[i]->getInformation();
	  const Vector &eleData = eleInfo.getData();
	  if (numDOF == 0) {
	    for (int j=0; j<eleData.Size(); j++)
	      (*data)(loc++) = eleData(j);
	  } else {
	    int dataSize = data->Size();
	    for (int j=0; j<numDOF; j++) {
	      int index = (*dof)(j);
	      if (index >= 0 && index < dataSize)
		(*data)(loc++) = eleData(index);		
	      else
		(*data)(loc++) = 0.0;		
	    }
	  }
	}
      }
    }

    //
    // send the response vector to the output handler for o/p
    //
    theOutputHandler->write(*data);
  }
  
  // succesfull completion - return 0
  return result;
}

int
HTElementRecorder::restart(void)
{
  if (data != 0)
    data->Zero();
  return 0;
}


int 
HTElementRecorder::setDomain(HeatTransferDomain &theDom)
{
  theDomain = &theDom;
  return 0;
}

int
HTElementRecorder::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}

int 
HTElementRecorder::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
  return 0;
}

int 
HTElementRecorder::initialize(void)
{
  if (theDomain == 0)
    return 0;

  if (theResponses != 0) {
    for (int i = 0; i < numEle; i++)
      delete theResponses[i];
    delete [] theResponses;
  }

  int numDbColumns = 0;

  // Set the response objects:
  //   1. create an array of pointers for them
  //   2. iterate over the elements invoking setResponse() to get the new objects & determine size of data
  //

  int i =0;
  ID xmlOrder(0,64);
  ID responseOrder(0,64);

  if (eleID != 0) {

    //
    // if we have an eleID we know Reponse size so allocate Response holder & loop over & ask each HeatTransferElement
    //

    int eleCount = 0;
    int responseCount = 0;

    if (echoTimeFlag == true) {
      xmlOrder[0] = 0;
      responseOrder[0] = 0;
      eleCount = 1;
      responseCount =1;
    }

    // loop over ele & set Reponses
    for (i=0; i<numEle; i++) {
      HeatTransferElement *theEle = theDomain->getElement((*eleID)(i));
      if (theEle != 0) {
	xmlOrder[eleCount] = i+1;
	eleCount++;
      }
    }

    theOutputHandler->setOrder(xmlOrder);

    //
    // do time
    //

    if (echoTimeFlag == true) {
      theOutputHandler->tag("TimeOutput");
      theOutputHandler->tag("ResponseType", "time");
      theOutputHandler->endTag(); // TimeOutput
      numDbColumns += 1;
    }

    //
    // if we have an eleID we know Reponse size so allocate Response holder & loop over & ask each HeatTransferElement
    //

    // allocate memory for Reponses & set to 0
    theResponses = new Response *[numEle];
    if (theResponses == 0) {
      opserr << "HTElementRecorder::initialize() - out of memory\n";
      return -1;
    }

    for (int k=0; k<numEle; k++)
      theResponses[k] = 0;

    // loop over ele & set Reponses
    for (i=0; i<numEle; i++) {
      HeatTransferElement *theEle = theDomain->getElement((*eleID)(i));
      if (theEle == 0) {
	theResponses[i] = 0;
      } else {
	theResponses[i] = theEle->setResponse((const char **)responseArgs, numArgs, *theOutputHandler);
	if (theResponses[i] != 0) {
	  // from the response type determine no of cols for each
	  Information &eleInfo = theResponses[i]->getInformation();
	  const Vector &eleData = eleInfo.getData();
	  int dataSize = eleData.Size();
	  if (numDOF == 0)
	    numDbColumns += dataSize;
	  else
	    numDbColumns += numDOF;

	  if (addColumnInfo == 1) {
	    if (numDOF == 0)
	      for (int j=0; j<dataSize; j++)
		responseOrder[responseCount++] = i+1;
	    else
	      for (int j=0; j<numDOF; j++)
		responseOrder[responseCount++] = i+1;
	  }
	}
      }
    }

    theOutputHandler->setOrder(responseOrder);

  } else {

    if (echoTimeFlag == true) {
      theOutputHandler->tag("TimeOutput");
      theOutputHandler->tag("ResponseType", "time");
      theOutputHandler->endTag(); // TimeOutput
      numDbColumns += 1;
    }

    //
    // if no eleID we don't know response size so make initial guess & loop over & ask ele
    // if guess to small, we enlarge
    //

    // initial size & allocation
    int numResponse = 0;
    numEle = 12;
    theResponses = new Response *[numEle];

    if (theResponses == 0) {
      opserr << "HTElementRecorder::initialize() - out of memory\n";
      return -1;
    }

    for (int k=0; k<numEle; k++)
      theResponses[k] = 0;

    // loop over ele & set Reponses
    HT_ElementIter &theElements = theDomain->getElements();
    HeatTransferElement *theEle;

    while ((theEle = theElements()) != 0) {
      Response *theResponse = theEle->setResponse((const char **)responseArgs, numArgs, *theOutputHandler);
      if (theResponse != 0) {
	if (numResponse == numEle) {
	  // Why is this created locally and not used? -- MHS
	  Response **theNextResponses = new Response *[numEle*2];
	  if (theNextResponses != 0) {
	    for (int i=0; i<numEle; i++)
	      theNextResponses[i] = theResponses[i];
	    for (int j=numEle; j<2*numEle; j++)
	      theNextResponses[j] = 0;
	  }
	  numEle = 2*numEle;
	  delete [] theNextResponses;
	}
	theResponses[numResponse] = theResponse;

	// from the response type determine no of cols for each
	Information &eleInfo = theResponses[numResponse]->getInformation();
	const Vector &eleData = eleInfo.getData();
	numDbColumns += eleData.Size();

	numResponse++;

      }
    }
    numEle = numResponse;
  }

  // create the vector to hold the data
  data = new Vector(numDbColumns);

  if (data == 0) {
    opserr << "HTElementRecorder::initialize() - out of memory\n";
    return -1;
  }
  
  theOutputHandler->tag("Data");
  initializationDone = true;

  return 0;
}
