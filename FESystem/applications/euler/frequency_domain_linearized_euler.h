//
//  frequency_domain_linearized_euler.h
//  FESystem
//
//  Created by Manav Bhatia on 3/8/13.
//
//

#ifndef __FESystem__frequency_domain_linearized_euler__
#define __FESystem__frequency_domain_linearized_euler__

// libmesh config includes
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_USE_COMPLEX_NUMBERS

// DiffSystem framework files
#include "libmesh/fem_system.h"
#include "libmesh/auto_ptr.h"

// FEsystem includes
#include "euler/euler_elem_base.h"


using namespace libMesh;

class SurfaceMotion;


class FrequencyDomainLinearizedEuler: public FEMSystem, public EulerElemBase
{
public:
    FrequencyDomainLinearizedEuler(EquationSystems& es,
                                   const std::string& name_in,
                                   const unsigned int number_in):
    FEMSystem(es, name_in, number_in),
    EulerElemBase()
    { }
    
    void init_data();
    
    virtual void init_context(DiffContext &context);

    
    virtual bool element_time_derivative (bool request_jacobian,
                                          DiffContext &context);
    
    virtual bool side_time_derivative (bool request_jacobian,
                                       DiffContext &context);
    
    
    AutoPtr<SurfaceMotion> surface_motion;
    
protected:
    
    std::vector<unsigned int> vars;
};


#endif // LIBMESH_USE_COMPLEX_NUMBERS

#endif /* defined(__FESystem__frequency_domain_linearized_euler__) */
