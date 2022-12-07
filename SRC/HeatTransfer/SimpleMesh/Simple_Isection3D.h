#ifndef SimpleIsection3D_h
#define SimpleIsection3D_h

#include <Simple_Entity.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Simple_Isection3D : public Simple_Entity
{
public:
	Simple_Isection3D(int tag, double HTI_centerX, double HTI_centerY, double HTI_centerZ,
                              double HTI_Bf, double HTI_Tf, double HTI_Hw, double HTI_Tw,
                              double HTI_UBf,double HTI_UTf, double HTI_Len);

	Simple_Isection3D(int tag, double HTI_centerX, double HTI_centerY, double HTI_centerZ,
                              double HTI_Bf, double HTI_Tf, double HTI_Hw, double HTI_Tw, double HTI_Len);
	~Simple_Isection3D();

	int InitialMeshCtrl(Vector& MeshCtrls, bool numCtrl = false);
    bool InitialSeeds(void); 
    const Vector& GetSeeds(int SeedTag);
    int GetNumofNodes(void);
	int GetNumofEles(void);
	const ID& GetNumCtrlID(void);
	int GenerateNodes(HeatTransferDomain*, int, const Vector&);
	int GenerateEles(HeatTransferDomain*, const ID&, HeatTransferMaterial*, HeatTransferMaterial* = 0);

private:
	 Vector Seeds1;
	 Vector Seeds2;
	 Vector Seeds3;

    double HTI_centerX, HTI_centerY, HTI_centerZ, HTI_Bf, HTI_Tf, HTI_Tw, HTI_Hw, HTI_UBf, HTI_UTf, HTI_Len;
	ID NumCtrlID;
};
#endif
