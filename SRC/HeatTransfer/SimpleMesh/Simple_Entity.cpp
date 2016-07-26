#include <Simple_Entity.h>

Simple_Entity::Simple_Entity(int tag, int EntityTypeTag):TaggedObject(tag), EntityTypeTag(EntityTypeTag),MeshTag(0)
{
	
}

Simple_Entity::~Simple_Entity()
{
	//
}

int 
Simple_Entity::getEntityTypeTag()
{
	return EntityTypeTag;
}

void
Simple_Entity::setMeshTag(int meshTag){
  MeshTag=meshTag;
}

int
Simple_Entity::getMeshTag()
{
	return MeshTag;
}

/*int 
Simple_Entity::getTag()
{

	return tag;
}*/

int
Simple_Entity::InitialMeshCtrl(Vector& MeshCtrls)
{
		return 0;
}

int
Simple_Entity::RefineSeeds(int SeedTag, const Vector& RefineSeedsInfo){
  return 0;
}