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
// modified by liming jiang

#ifndef FireModel_h
#define FireModel_h

//#include <ID.h>
#include <TaggedObject.h>
#include <MovableObject.h>
#include <OPS_Globals.h>
#include <Vector.h>
class HeatTransferDomain;
class HeatFluxBC;

class FireModel: public TaggedObject
{
    public:
	  FireModel(int tag, int fireTypeTag);
	  virtual ~FireModel();
	  
	  virtual void setDomain(HeatTransferDomain* theDomain);
	  virtual void applyFluxBC(HeatFluxBC* theFlux, double time) = 0;

	  virtual int getFireTypeTag(void);
	  virtual double getFirePars(int parTag=0);
	  virtual int setFirePars(double time, const Vector& firePars = 0);
	  virtual double getFireOut(double time, const Vector& locs = 0);
	  virtual void  Print(OPS_Stream&, int = 0);

	protected:
		HeatTransferDomain* the_domain; // The domain pointer here associates a FireModel object and a HT_Domain
		int FireTypeTag;  // object thus enables a FireModel object fetch information from the 
		                                // domain.
    private:

};

#endif


//11:standard fire
//12:external fire
//13:hydroCarbon fire
//2:Parametric fire
//3:Localised fire
//5:Alpert jet fire
//6:Idealised fire
//7:Moving Localised fire
