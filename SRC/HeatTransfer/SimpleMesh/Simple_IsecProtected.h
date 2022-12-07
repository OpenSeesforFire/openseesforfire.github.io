#ifndef SimpleIsecProtected_h
#define SimpleIsecProtected_h

#include <Simple_Entity.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Simple_IsecProtected : public Simple_Entity
{
public:
	Simple_IsecProtected(int tag, double HTI_centerX, double HTI_centerY, double HTI_Bf, double HTI_Tf, double HTI_Tw, double HTI_Hw , double HTI_coat);
	~Simple_IsecProtected();

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
	double EleX, EleY, EleX_Web, EleY_Web, Ele_coat;
 
    double HTI_centerX, HTI_centerY, HTI_Bf, HTI_Tf, HTI_Tw, HTI_Hw, HTI_UBf, HTI_UTf, HTI_Coat;
	ID NumCtrlID;
};
#endif
