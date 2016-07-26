
#ifndef SIFSection_h
#define SIFSection_h

#include <TaggedObject.h>
#include <MovableObject.h>
#include <Vector.h>
#include <ID.h>
#include <SIFMaterial.h>
#include <SectionForceDeformation.h>

class ID;
class TaggedObject; 
class SIFSection: public TaggedObject
{
public:
	SIFSection(int tag, int SectionTypeTag, 
	 SIFMaterial* theSifMaterial);
	~SIFSection();
	
	int AssignSectionPars(const Vector& sectionPars);
	const Vector& getSectionPars(void);
	
  int getSectionTypeTag(void);
  
  
  int assignSecondMaterial(SIFMaterial* theSifMaterial);
  
  SIFMaterial* getSIFMaterialPtr(int MaterialTag =1);   //added for get pointer to SIFMaterial;
  
  SectionForceDeformation* DefineBeamSection(bool isElastic = false);
  SectionForceDeformation* DefineShellSection(bool isElastic = false);
	//const Vector& getMaterialPars(void); ,

	virtual void  Print(OPS_Stream&, int = 0) {return;};
   
private:
  int SectionTypeTag;
	Vector theSectionPars;
    SIFMaterial* SIFMaterialPtr;
	SIFMaterial* SecSIFMaterialPtr;
  
	//Vector MaterialPars(10);
 
	
};
#endif