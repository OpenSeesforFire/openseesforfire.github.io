#ifndef Simple_Brick_h
#define Simple_Brick_h

#include <Simple_Entity.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <vector>
#include <stdlib.h>

class Simple_Brick : public Simple_Entity
{
public:
    Simple_Brick(int tag, double CenterX, double CenterY, double CenterZ, double Breadth_X, double Height_Y,double Length_Z);

	~Simple_Brick();
   
  int InitialMeshCtrl(Vector& MeshCtrls, bool numCtrl = false);
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
    double Breadth,Height,Length;
	double a1,b1,c1,a7,b7,c7;
	Vector Seeds1;
	Vector Seeds2;
	Vector Seeds3;
	ID NumCtrlID;
};
#endif
