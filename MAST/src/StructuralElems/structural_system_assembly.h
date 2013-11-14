
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
#include "libmesh/eigen_system.h"
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
    class FunctionBase;
    class SensitivityParameters;
    class BoundaryCondition;

    /*!
     *    class defies the structural analysis system and implements the
     *    calculation of structural finite elements for a variety of 1, 2 and 
     *    3D elements.
     */
    class StructuralSystemAssembly :
    public System::Assembly,
    public NonlinearImplicitSystem::ComputeResidualandJacobian,
    public System::SensitivityAssembly,
    public EigenSystem::EigenproblemSensitivityAssembly
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
         *   adds the specified side loads for the boudnary with tag \p b_id
         */
        void add_side_load(boundary_id_type bid, MAST::BoundaryCondition& load);

        
        /*!
         *   adds the specified volume loads for the elements with
         *   subdomain tag \p s_id
         */
        void add_volume_load(subdomain_id_type bid, MAST::BoundaryCondition& load);

        
        /*!
         *    fills the set \par dof_ids with the dof ids of the Dirichlet
         *    dofs
         */
        void get_dirichlet_dofs(std::set<unsigned int>& dof_ids) const;
        
        /*!
         *    sets the same property for all cards
         */
        void set_property_for_all_elems(const MAST::ElementPropertyCardBase& prop);

        /*!
         *    sets the property for specified element
         */
        void set_property_for_elem(const Elem& e, const MAST::ElementPropertyCardBase& prop);

        /*!
         *    get the material property for the specified element
         */
        const MAST::ElementPropertyCardBase& get_property_card(const Elem& elem) const;
        
        /*!
         *   Adds the parameter and function pairing
         */
        void add_parameter(Real* par, MAST::FunctionBase* f);
        
        /*!
         *   Returns the function corresponding to a parameter
         */
        const MAST::FunctionBase* get_parameter(Real* par) const;
        
        /*!
         *    assembles the matrices for eigenproblem depending on the analysis type
         */
        virtual void assemble();

        /*!
         *    function that assembles the matrices and vectors quantities for
         *    nonlinear solution
         */
        virtual void residual_and_jacobian (const NumericVector<Number>& X,
                                            NumericVector<Number>* R,
                                            SparseMatrix<Number>*  J,
                                            NonlinearImplicitSystem& S);

        /**
         * Assembly function.  This function will be called
         * to assemble the sensitivity of system residual prior to a solve and must
         * be provided by the user in a derived class. The method provides dR/dp_i
         * for \par i ^th parameter in the vector \par parameters.
         *
         * If the routine is not able to provide sensitivity for this parameter,
         * then it should return false, and the system will attempt to use
         * finite differencing.
         */
        virtual bool sensitivity_assemble (const ParameterVector& parameters,
                                           const unsigned int i,
                                           NumericVector<Number>& sensitivity_rhs);

        /**
         * Assembly function.  This function will be called
         * to assemble the sensitivity of eigenproblem matrices.
         * The method provides dA/dp_i and dB/dpi for \par i ^th parameter
         * in the vector \par parameters.
         *
         * If the routine is not able to provide sensitivity for this parameter,
         * then it should return false, and the system will attempt to use
         * finite differencing.
         */
        virtual bool sensitivity_assemble (const ParameterVector& parameters,
                                           const unsigned int i,
                                           SparseMatrix<Number>* sensitivity_A,
                                           SparseMatrix<Number>* sensitivity_B);


    protected:
        
        /*!
         *    assembles the residual and Jacobian for static or dynamic analyses
         */
        void _assemble_residual_and_jacobian (const NumericVector<Number>& X,
                                              NumericVector<Number>* R,
                                              SparseMatrix<Number>*  J,
                                              NonlinearImplicitSystem& S,
                                              const MAST::SensitivityParameters* params);

        
        /*!
         *    Calculates the A and B matrices of a modal analysis eigenproblem
         *    [A] {X} = [B] {X} [lambda]
         */
        void _assemble_matrices_for_modal_analysis(const NumericVector<Number>& X,
                                                   SparseMatrix<Number>&  matrix_A,
                                                   SparseMatrix<Number>&  matrix_B,
                                                   const MAST::SensitivityParameters* params);
        
        /*!
         *    Calculates the A and B matrices of a buckling analysis eigenproblem
         *    [A] {X} = [B] {X} [lambda]
         */
        void _assemble_matrices_for_buckling_analysis(const NumericVector<Number>& X,
                                                      SparseMatrix<Number>&  matrix_A,
                                                      SparseMatrix<Number>&  matrix_B,
                                                      const MAST::SensitivityParameters* params);

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
         *    flag to specify if the same material card is to be used for all 
         *    elements. False by default
         */
        bool _if_same_property_for_all_elems;
        
        /*!
         *     property card if the same card is to be used
         */
        const MAST::ElementPropertyCardBase* _property;
        
        /*!
         *   boundary conditions for this analysis
         */
        std::map<unsigned int, MAST::StructuralBoundaryConditionType> _boundary_condition;
        
        /*!
         *   map of element property cards for each element
         */
        std::map<const Elem*, const MAST::ElementPropertyCardBase*> _element_property;
        
        /*!
         *   vector of material property cards
         */
        std::vector<const MAST::MaterialPropertyCardBase*> _material_property;
        
        /*!
         *   map of sensitivity parameters and the corresponding functions that
         *   are directly dependent on these parameters
         */
        std::map<Real*, const MAST::FunctionBase*> _parameter_map;
        
        /*!
         *   side boundary condition map of boundary id and load
         */
        std::multimap<boundary_id_type, MAST::BoundaryCondition*> _side_bc_map;
        
        /*!
         *   volume boundary condition map of boundary id and load
         */
        std::multimap<subdomain_id_type, MAST::BoundaryCondition*> _vol_bc_map;
    };
}

#endif /* defined(__MAST__structural_system_assembly_h__) */
