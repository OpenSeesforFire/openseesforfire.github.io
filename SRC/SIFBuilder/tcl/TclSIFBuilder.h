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
**                                                                    **                                                                 **
** ****************************************************************** */
                                                                        

/**********************************************************************
** This class is aiming to provide a Tcl interface to define models  **
** for simulating structural behaviours under fire action.           **
** Developed by:                                                     **
**   Liming Jiang (liming.jiang@ed.ac.uk)                            **
**   Praven Kamath(Praveen.Kamath@ed.ac.uk)                          **
**   Xu Dai(X.Dai@ed.ac.uk)                                          **
**   Asif Usmani(asif.usmani@ed.ac.uk)                               **
**********************************************************************/
// $Revision: 2.4.0.1 $


#ifndef TclSIFBuilder_h
#define TclSIFBuilder_h

#include <tcl.h>

#include <Domain.h>
#include <HeatTransferDomain.h>
#include <HeatTransferMaterial.h>  // for member functions adding and returning HeatTransfer Materials
#include <Simple_Entity.h>    //for member functions adding and returning Simple Entities
#include <Simple_Mesh.h>      //for member functions adding and returning Simple Meshes

#include <SIFBuilderDomain.h> //To use member functions in SIFBuilderDomain
#include <SIFMaterial.h>      //To use member functions in SIFMaterial
#include <SIFMember.h>		  //To use member functions in SIFMember
#include <SIFXWall.h>		  //To use member functions in SIFXWall
#include <SIFCompartment.h>		  //To use member functions in SIFCompartment3D
#include <SIFSlab.h>

//int inputCheck();

class TclSIFBuilder
{
public:
	TclSIFBuilder(Domain& theDomain, Tcl_Interp *interp, int builderType=1);
	~TclSIFBuilder();
    
	SIFBuilderDomain* getSIFBuilderDomain();

   // int addMaterial(SIFMaterial &theSIFMaterial);
	
	//int addHTEntity(Simple_Entity &theHTEntity);
	//int addHTMesh(Simple_Mesh &theHTMesh); 
	

    
protected:

private:
    Tcl_Interp *theInterp;
    int ndm;
	//TaggedObjectStorage *theHTEntities;
	//TaggedObjectStorage *theHTMeshes;
};

#endif







