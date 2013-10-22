
#ifndef __MAST__structural_system_assembly_h__
#define __MAST__structural_system_assembly_h__

// C++ includes
#include <iostream>
#include <map>
#include <memory>

// libmesh includes
#include "libmesh/libmesh_config.h"
#include "libmesh/getpot.h"
#include "libmesh/equation_systems.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"


using namespace libMesh;


namespace MAST
{
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
    
    // Forward decleration
    class ElementPropertyCardBase;
    class MaterialPropertyCardBase;
    
    
    /*!
     *    class defies the structural analysis system and implements the
     *    calculation of structural finite elements for a variety of 1, 2 and 
     *    3D elements.
     */
    class StructuralSystemAssembly :
    public NonlinearImplicitSystem::ComputeResidualandJacobian
    {
    public:
        // Constructor
        StructuralSystemAssembly(System& sys,
                                 MAST::StructuralAnalysisType t,
                                 GetPot& _infile);

        /*!
         *   virtual destructor
         */
        virtual ~StructuralSystemAssembly()
        { }
        
        
        /*!
         *    function that assembles the matrices and vectors quantities for
         *    nonlinear solution
         */
        virtual void residual_and_jacobian (const NumericVector<Number>& X,
                                            NumericVector<Number>* R,
                                            SparseMatrix<Number>*  J,
                                            NonlinearImplicitSystem& S);
        
        
    protected:
        
        /*!
         *    Calculates the A and B matrices of an eigenproblem
         *    [lambda] [A] {X} = [B] {X}
         *    The assembly procedure of the matrices depends on the
         *    nature of analysis.
         */
        

        /*!
         *    System for which analysis is to be performed
         */
        System& _system;
        
        
        /*!
         *   Analysis type for this structure
         */
        MAST::StructuralAnalysisType _analysis_type;
        

        GetPot& _infile;

        /*!
         *   boundary conditions for this analysis
         */
        std::map<unsigned int, MAST::StructuralBoundaryConditionType> _boundary_condition;
        
        /*!
         *   map of element property cards for each element
         */
        std::map<const Elem*, MAST::ElementPropertyCardBase*> _element_property;
        
        
        /*!
         *   vector of material property cards
         */
        std::vector<MAST::MaterialPropertyCardBase*> _material_property;
        
    };
    
}

#endif /* defined(__MAST__structural_system_assembly_h__) */
