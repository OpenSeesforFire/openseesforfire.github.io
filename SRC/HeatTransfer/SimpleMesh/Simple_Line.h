#ifndef SimpleLine_h
#define SimpleLine_h

#include <Simple_Entity.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <vector>
#include <stdlib.h>


class Simple_Line : public Simple_Entity
{
public:
	//Simple_Line(int tag, double a1, double b1, double a3, double b3,double a2, double b2, double a4, double b4);
	Simple_Line(int tag, double centerX, double lengthX);
	~Simple_Line();
   
	virtual int InitialMeshCtrl(Vector&, bool numCtrl = false);
	virtual bool InitialSeeds();
	//virtual bool RefineSeeds(int SeedTag, const Vector& RefinedSeedsInfo);

	const Vector& GetSeeds(int SeedTag);
    const ID& GetNumCtrlID();
	int GetNumofNodes();
	int GetNumofEles();
	int GenerateNodes(HeatTransferDomain*, int, const Vector&);
	int GenerateEles(HeatTransferDomain*, const ID&, HeatTransferMaterial*, HeatTransferMaterial* = 0);
	//virtual int GetSectionTag();

private:
    double a1,a2,CenterX,LengthX;
	Vector Seeds1;
	ID NumCtrlID;
};
#endif
