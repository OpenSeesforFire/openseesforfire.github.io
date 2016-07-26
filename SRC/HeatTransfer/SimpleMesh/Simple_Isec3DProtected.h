#ifndef SimpleIsec3DProtected_h
#define SimpleIsec3DProtected_h

#include <Simple_Entity.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Simple_Isec3DProtected : public Simple_Entity
{
public:
	Simple_Isec3DProtected(int tag, double HTI_centerX, double HTI_centerY, double HTI_centerZ,
                              double HTI_Bf, double HTI_Tf, double HTI_Hw, double HTI_Tw,
                              double HTI_UBf,double HTI_UTf, double HTI_Len);

	Simple_Isec3DProtected(int tag, double HTI_centerX, double HTI_centerY, double HTI_centerZ,
                              double HTI_Bf, double HTI_Tf, double HTI_Hw, double HTI_Tw, double HTI_Len);
	~Simple_Isec3DProtected();

	virtual int InitialMeshCtrl(Vector& MeshCtrls);
    virtual bool InitialSeeds(void); 
    virtual const Vector& GetSeeds(int SeedTag);
    virtual int GetNumofNodes(void);
	virtual int GetNumofEles(void);
	virtual const ID& GetNumCtrlID(void);
private:
	 Vector Seeds1;
	 Vector Seeds2;
	 Vector Seeds3;

    double HTI_centerX, HTI_centerY, HTI_centerZ, HTI_Bf, HTI_Tf, HTI_Tw, HTI_Hw, HTI_UBf, HTI_UTf, HTI_Len;
	ID NumCtrlID;
};
#endif
