//
//  assembleEuler.h
//  FESystem
//
//  Created by Manav Bhatia on 2/21/13.
//
//

#ifndef __FESystem__assembleEuler__
#define __FESystem__assembleEuler__

// DiffSystem framework files
#include "libmesh/fem_system.h"

// FESystem includes
#include "euler/euler_elem_base.h"

using namespace libMesh;


class EulerSystem : public FEMSystem, public EulerElemBase
{
public:
    // Constructor
    EulerSystem(EquationSystems& es,
                const std::string& name_in,
                const unsigned int number_in)
    : FEMSystem(es, name_in, number_in),
    EulerElemBase()
    {}
    
    void init_data();
    
    virtual void init_context(DiffContext &context);
    
    virtual void postprocess();

    virtual bool element_time_derivative (bool request_jacobian,
                                          DiffContext &context);
    
    virtual bool side_time_derivative (bool request_jacobian,
                                       DiffContext &context);
    
    virtual bool mass_residual (bool request_jacobian,
                                DiffContext& context);
    
    
protected:
    
    std::vector<unsigned int> vars;
};

#endif /* defined(__FESystem__assembleEuler__) */
