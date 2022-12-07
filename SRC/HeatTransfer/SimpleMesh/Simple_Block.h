#ifndef SimpleBlock_h
#define SimpleBlock_h

#include <Simple_Entity.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <vector>
#include <stdlib.h>

class Simple_Block : public Simple_Entity
{
public:
	Simple_Block(int tag, double a1, double b1, double a2, double b2, double a3, double b3, double a4, double b4);
	Simple_Block(int tag, double centerX, double centerY, double breadthX, double heightY);
	~Simple_Block();
   
	int InitialMeshCtrl(Vector& MeshCtrls, bool numCtrl);
	bool InitialSeeds();
	int RefineSeeds(int SeedTag, const Vector& RefinedSeedsInfo);

	const Vector& GetSeeds(int SeedTag);
  const ID& GetNumCtrlID();
	int GetNumofNodes();
	int GetNumofEles();
	int GenerateNodes(HeatTransferDomain*, int, const Vector&);
	int GenerateEles(HeatTransferDomain*, const ID&, HeatTransferMaterial*, HeatTransferMaterial* = 0);
	//virtual int GetSectionTag();

private:
    double a1,b1,a2,b2,a3,b3,a4,b4,CenterX,CenterY,BreadthX,HeightY;
	Vector Seeds1;
	Vector Seeds2;
	Vector Seeds3;
	Vector Seeds4;
	ID NumCtrlID;
};
#endif
