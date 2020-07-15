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
**   Yaqiang Jiang (y.jiang@ed.ac.uk)								  **
**   Liming Jiang (Liming.Jiang@ed.ac.uk)                             **
**   Asif Usmani (asif.usmani@ed.ac.uk)                               **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 2.4.2 $


//
// Written by Liming Jiang (Liming.Jiang@ed.ac.uk)
//

#ifndef TclHeatTransferModule_h
#define TclHeatTransferModule_h

#include <tcl.h>

#include <HeatTransferDomain.h>
#include <HeatTransferMaterial.h>  // for member functions adding and returning HeatTransfer Materials
#include <Simple_Entity.h>    //for member functions adding and returning Simple Entities
#include <Simple_Mesh.h> 
#include <Simple_MeshIter.h> //for member functions adding and returning Simple Meshes
#include <HTConstants.h>  // for adding Heat transfer constants
#include <HTEleSet.h>    // for grouping elements 
#include <HTNodeSet.h>   // for grouping nodes
#include <FireModel.h>
//int inputCheck();

class TclHeatTransferModule
{
public:
	TclHeatTransferModule(int ndm, Tcl_Interp *interp);
	~TclHeatTransferModule();
    
	HeatTransferDomain* getHeatTransferDomain();

  int addHTMaterial(HeatTransferMaterial &theHTMaterial);
	int addHTEntity(Simple_Entity &theHTEntity);
	int addHTMesh(Simple_Mesh *theHTMesh);
  int addHTConstants(HTConstants *theHTConstants);
  int addHTNodeSet(HTNodeSet* theHTNodeSet);
  int addHTEleSet(HTEleSet* theHTEleSet);
  int addFireModel(FireModel* theFireModel);
  void clearAll();
	
	HeatTransferMaterial *getHTMaterial(int tag);
	Simple_Entity *getHTEntity(int tag);
	Simple_Mesh *getHTMesh(int tag);
  HTConstants *getHTConstants(int tag);
  HTEleSet *getHTEleSet(int tag);
  HTNodeSet *getHTNodeSet(int tag);
  FireModel *getFireModel(int tag);
  
  
	Simple_MeshIter &getHTMeshs();

    
protected:

private:
  Tcl_Interp *theInterp;
  
  TaggedObjectStorage *theHTMaterials;
  TaggedObjectStorage *theHTEntities;
  TaggedObjectStorage *theHTMeshes;
  TaggedObjectStorage *theHTCons;
  TaggedObjectStorage *theHTEleSets;
  TaggedObjectStorage *theHTNodeSets;
  TaggedObjectStorage *theFireModels;
  
  Simple_MeshIter *theHTMeshIter;
  
  int ndm;
  
};

#endif







