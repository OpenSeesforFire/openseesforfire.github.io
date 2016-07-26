#ifndef SIFMember_h
#define SIFMember_h

#include <TaggedObject.h>
#include <MovableObject.h>
#include <Vector.h>
#include <ID.h>
#include <SIFSection.h>

//This class defines a base class for SIFMemeber
//Edinburgh University OpenSees Developer Group-2015

class ID;
class Vector;

class SIFMember: public TaggedObject
{
public:
	SIFMember(int tag,int memberTypeTag, int typeTag, const ID& theMemInfo);
  SIFMember();
  
	virtual ~SIFMember ();
  
  
  virtual int assignSection(SIFSection* theSection);
  virtual SIFSection* getSIFSectionPtr(void);
  virtual ID& getConnectedJoints(void);
  
  virtual int getConnectedSlab(void) = 0;
  virtual int getTypeTag(void);
  virtual int getMemberTypeTag(void);
  virtual double getGamma(void) = 0;

  virtual int setIntNodeTags(const ID& nodeTag);
  virtual const ID& getIntNodeTags(void);
  virtual int setIntEleTags(const ID& eleTag);
  virtual const ID& getIntEleTags(void);

  virtual int AddFireAction(int theFireActionID);
  virtual int WipeFireAction();
  virtual const ID& getAppliedFireActions();
  virtual const ID& getMemberInfo(void);

  virtual double getPartialDamage(Vector& DmageVec);
  virtual int setPartialDamage(const Vector& damageVec);
  
  
    //virtual int GetEntityTypeTag();
	//int cTag;
	
	
	virtual void  Print(OPS_Stream&, int = 0) {return;};

protected:
  int  TypeTag;
  int MemberTypeTag;
  Vector DamVec;
  ID theMemberInfo;

  SIFSection* SIFSectionPtr;
  ID* ConnectedJoints;
  ID IntNodeTags;
  ID IntEleTags;
  
  ID AppliedFireActions;
};
#endif