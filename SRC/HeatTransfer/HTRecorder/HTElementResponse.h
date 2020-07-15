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
                                                                        
// $Revision: 1.5 $
// $Date: 2009-12-17 23:50:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/response/HTElementResponse.h,v $
                                                                        
// Written: MHS 
// Created: Oct 2000
//
// Description: This file contains the HTElementResponse class interface

#ifndef HTElementResponse_h
#define HTElementResponse_h

#include <Response.h>
#include <Information.h>

class HeatTransferElement;

class ID;
class Vector;
class Matrix;

class HTElementResponse : public Response
{
public:
	HTElementResponse(HeatTransferElement *ele, int id);
	HTElementResponse(HeatTransferElement *ele, int id, int val);
	HTElementResponse(HeatTransferElement *ele, int id, double val);
	HTElementResponse(HeatTransferElement *ele, int id, const ID &val);
	HTElementResponse(HeatTransferElement *ele, int id, const Vector &val);
	HTElementResponse(HeatTransferElement *ele, int id, const Matrix &val);
	HTElementResponse(HeatTransferElement *ele, int id, const Vector &val1, const ID &val2);

	~HTElementResponse();

	int getResponse(void);


private:
	HeatTransferElement *theElement;
	int responseID;
};

#endif
