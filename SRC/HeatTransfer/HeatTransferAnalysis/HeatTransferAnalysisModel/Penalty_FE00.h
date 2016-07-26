#ifndef Penalty_FE_h
#define Penalty_FE_h

#include <HT_FE_Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>

class HeatTransferElement;
class HeatTransferIntegrator;
class HT_AnalysisModel;
class HeatTransferDomain;
class TemperatureBC;
class HeatTransferNode;

class Penalty_FE: public HT_FE_Element
{
  public:
    Penalty_FE(int tag, 
		       HeatTransferDomain& theDomain, 
		       TemperatureBC& theTempBC, 
		       double alpha = 1.0e8);    

    virtual ~Penalty_FE();    

    // public methods
    virtual int  setID(void);
    virtual const Matrix &getTangent(Integrator *theIntegrator);
    virtual const Vector &getResidual(Integrator *theIntegrator);

  protected:
    
  private:
    double alpha;
    TemperatureBC* the_temp_bc;
    HeatTransferNode* the_node;
    static Matrix tang;
    static Vector resid;
};

#endif