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
Simple_Entity::InitialMeshCtrl(Vector& MeshCtrls, bool numCtrl)
{
		return 0;
}

int
Simple_Entity::RefineSeeds(int SeedTag, const Vector& RefineSeedsInfo){
  return 0;
}

int Simple_Entity::GenerateNodes(HeatTransferDomain* theDomain, int ndof, const Vector& vector)
{
	return 0;
}

int Simple_Entity::GenerateEles(HeatTransferDomain* theDomain, const ID& id, HeatTransferMaterial* theHTMaterial,HeatTransferMaterial* theHTMaterial1)
{
	return 0;
}
//Entity Type:
//Line:5; Block:0;Brick:2;Isection:1;Isection3D:3; IsecProtected:11;composite 2D:6;composite3D:7