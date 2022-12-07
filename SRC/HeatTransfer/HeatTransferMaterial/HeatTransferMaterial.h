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
// Note: This class was adapted from Material   
#ifndef HeatTransferMaterial_h
#define HeatTransferMaterial_h

#include <HeatTransferDomainComponent.h>

class OPS_Stream;
class Matrix;
class Information;
class Response;

class HeatTransferMaterial: public TaggedObject
{
    public:
		HeatTransferMaterial(int tag);    
		virtual ~HeatTransferMaterial();

		// method for this material to update itself according to its new parameters
		virtual int setTrialTemperature(double T, int par=0 ) = 0;
		virtual const Matrix& getConductivity() = 0;
		virtual double getRho() = 0;
		virtual double getSpecificHeat() = 0;
		virtual double getEnthalpy() = 0;
		virtual double getEnthalpy(double temperature) = 0;
		virtual HeatTransferMaterial* getCopy() = 0;
		virtual bool getIfHeatGen();
		virtual double getHeatGen(double par =0);

		virtual int commitState() = 0;
		virtual int revertToLastCommit() = 0;
		virtual int revertToStart() = 0;
		virtual void update() = 0;
		virtual void  Print(OPS_Stream&, int = 0) {return;};

		virtual Response* setResponse(const char** argv, int argc,
			OPS_Stream& theOutputStream);
		virtual int getResponse(int responseID, Information& matInformation);

		virtual const Vector&  getPars();

    protected:
		Matrix* k;
    
    private:
};


#endif


















