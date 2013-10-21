
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


enum StructuralAnalysisType {
    STATIC,
    DYNAMIC,
    MODAL,
    BUCKLING
};


namespace MAST
{
    
    // Forward decleration
    class ElementPropertyCardBase;
    class MaterialPropertyCardBase;
    
    
    /*!
     *    class defies the structural analysis system and implements the
     *    calculation of structural finite elements for a variety of 1, 2 and 
     *    3D elements.
     */
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
        
        /*!
         *    Calculates the A and B matrices of an eigenproblem
         *    [lambda] [A] {X} = [B] {X}
         *    The assembly procedure of the matrices depends on the
         *    nature of analysis.
         */
        //virtual bool assemble_eigenproblem_matrices();
        
        std::vector<unsigned int> vars;
        
        /*!
         *   Analysis type for this structure
         */
        StructuralAnalysisType analysis_type;
        
        /*!
         *   map of element property cards for each element
         */
        std::map<Elem*, ElementPropertyCardBase*> element_property;
        
        /*!
         *   vector of material property cards
         */
        std::vector<MaterialPropertyCardBase*> material_property;
        
        
        
    protected:
        
        GetPot& _infile;
        unsigned int dim;
        std::map<unsigned int, StructuralBoundaryConditionType> _boundary_condition;
    };
    
}

#endif /* defined(__MAST__structural_system_base_h__) */
