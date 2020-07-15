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
                                                                        
// $Revision: 1.4 $
// $Date: 2007-03-02 00:12:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/response/HTMaterialResponse.cpp,v $
                                                                        
// Written: MHS 
// Created: Oct 2000
//
// Description: This file contains the HTMaterialResponse class implementation

#include <HTMaterialResponse.h>
#include <HeatTransferMaterial.h>

HTMaterialResponse::HTMaterialResponse(HeatTransferMaterial *mat, int id):
Response(), theMaterial(mat), responseID(id)
{

}

HTMaterialResponse::HTMaterialResponse(HeatTransferMaterial *mat, int id, int val):
Response(val), theMaterial(mat), responseID(id)
{

}

HTMaterialResponse::HTMaterialResponse(HeatTransferMaterial *mat, int id, double val):
Response(val), theMaterial(mat), responseID(id)
{

}

HTMaterialResponse::HTMaterialResponse(HeatTransferMaterial *mat, int id, const ID &val):
Response(val), theMaterial(mat), responseID(id)
{

}

HTMaterialResponse::HTMaterialResponse(HeatTransferMaterial *mat, int id, const Vector &val):
Response(val), theMaterial(mat), responseID(id)
{

}

HTMaterialResponse::HTMaterialResponse(HeatTransferMaterial *mat, int id, const Matrix &val):
Response(val), theMaterial(mat), responseID(id)
{

}

HTMaterialResponse::~HTMaterialResponse()
{

}

int
HTMaterialResponse::getResponse(void)
{
	return theMaterial->getResponse(responseID, myInfo);
}

