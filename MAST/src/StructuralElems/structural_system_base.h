
#ifndef __MAST__structural_system_base_h__
#define __MAST__structural_system_base_h__

// C++ includes
#include <iostream>
#include <map>
#include <memory>

// libmesh includes
#include "libmesh/libmesh_config.h"
#include "libmesh/fem_system.h"
#include "libmesh/getpot.h"
#include "libmesh/equation_systems.h"


using namespace libMesh;


enum StructuralBoundaryConditionType
{
    DIRICHLET,
    TRACTION
};


class StructuralSystemBase : public FEMSystem
{
public:
    // Constructor
    StructuralSystemBase(EquationSystems& es,
                     const std::string& name_in,
                     const unsigned int number_in)
    : FEMSystem(es, name_in, number_in),
    _infile(*es.parameters.get<GetPot*>("input_file")),
    dim(0)
    {}
    
    void init_data();
    
    virtual void init_context(DiffContext &context);

    
    virtual bool element_time_derivative (bool request_jacobian,
                                          DiffContext &context);
    
    virtual bool side_time_derivative (bool request_jacobian,
                                       DiffContext &context);
    
    virtual bool mass_residual (bool request_jacobian,
                                DiffContext& context);

    std::vector<unsigned int> vars;
    
    
protected:
    
    GetPot& _infile;
    unsigned int dim;
    std::map<unsigned int, StructuralBoundaryConditionType> _boundary_condition;
};



#endif /* defined(__MAST__structural_system_base_h__) */
