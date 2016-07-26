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

// This class is a modified version of DOF_Numberer 
// for the heat transfer module
// Modified by: Yaqiang Jiang(y.jiang@ed.ac.uk)

// DOF_Numberer
// Written: fmk 
// Created: 9/96
// Revision: A
   
#ifndef HT_DOF_Numberer_h
#define HT_DOF_Numberer_h

#include <MovableObject.h>

class HT_AnalysisModel;
class GraphNumberer;
//class FEM_ObjectBroker;
class ID;
class Channel;

class HT_DOF_Numberer: public MovableObject
{
  public:
    //DOF_Numberer(int classTag);
    HT_DOF_Numberer(GraphNumberer& theGraphNumberer);    
    //DOF_Numberer();    
    virtual ~HT_DOF_Numberer();

    virtual void setLinks(HT_AnalysisModel& the_model);
    
    // pure virtual functions
    virtual int numberDOF(int lastDOF_Group = -1);
    //virtual int numberDOF(ID &lastDOF_Groups);    

	virtual int sendSelf(int commitTag, Channel& theChannel){return -1;};
    virtual int recvSelf(int commitTag, Channel& theChannel, 
		                 FEM_ObjectBroker& theBroker){return -1;};

  protected:
    HT_AnalysisModel* getAnalysisModelPtr(void) const;
    GraphNumberer* getGraphNumbererPtr(void) const;
    
  private:
    HT_AnalysisModel* theModel;
    GraphNumberer* theGraphNumberer;
};

#endif
