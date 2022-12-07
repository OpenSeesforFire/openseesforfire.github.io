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
//Modified by Liming Jiang (liming.jiang@ed.ac.uk)
//
// Note: This class was adapted from Material   
#include <HeatTransferMaterial.h>
#include <Information.h>
#include <HTMaterialResponse.h>


HeatTransferMaterial::HeatTransferMaterial(int tag)
:TaggedObject(tag), k(0)
{


}


HeatTransferMaterial::~HeatTransferMaterial()
{
  // does nothing


}


const Vector&
HeatTransferMaterial::getPars() {
    opserr << "HeatTransferMaterial::getPars should not be called";
    return 0;
}


Response*
HeatTransferMaterial::setResponse(const char** argv, int argc,
    OPS_Stream& theOutput)
{
    Response* theResponse = 0;

    if((strcmp(argv[0], "phaseTag") == 0)||(strcmp(argv[0], "-phaseTag")==0)||(strcmp(argv[0], "phase") == 0)){
        theOutput.tag("ResponseType", "phaseType");
        theResponse = new HTMaterialResponse(this, 1, this->getPars()(0));

        theOutput.endTag();

    }
    else if ((strcmp(argv[0], "charTime") == 0) || (strcmp(argv[0], "-charTime") == 0) || (strcmp(argv[0], "CharTime") == 0)) {
        theOutput.tag("ResponseType", "phaseType");
        theResponse = new HTMaterialResponse(this, 2, this->getPars()(1));

        theOutput.endTag();

    }

return theResponse;

}

int
HeatTransferMaterial::getResponse(int responseID, Information& matInfo)
{

    static Vector tempData(2);  //L.jiang [SIF]
    static Information infoData(tempData);  //L.jiang [SIF]

    // each subclass must implement its own stuff  

    switch (responseID) {
    case 1:
        matInfo.setDouble(this->getPars()(0));
        return 0;
    case 2:
        matInfo.setDouble(this->getPars()(1));
        return 0;
    default:
        return -1;
    }
}


bool
HeatTransferMaterial::getIfHeatGen()
{
    return false;
}


double 
HeatTransferMaterial::getHeatGen(double par)
{
    return 0;
}