//
//  structural_system_assembly.h
//  MAST
//
//  Created by Manav Bhatia on 12/16/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

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


// MAST includes
#include "BoundaryConditions/boundary_condition.h"
#include "Numerics/constant_function.h"


using namespace libMesh;


namespace MAST
{
    enum StructuralAnalysisType {
        STATIC,
        DYNAMIC,
        MODAL,
        BUCKLING
    };
    
    // Forward decleration
    class ElementPropertyCardBase;
    class MaterialPropertyCardBase;
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
         *   returns a reference to the System object
         */
        System& get_system() {
            return _system;
        }

        
        /*!
         *    returns the analysis type
         */
        MAST::StructuralAnalysisType analysis_type() const {
            return _analysis_type;
        }
        
        /*!
         *   clear the loads and pointer to static solution system for
         *   this structural model
         */
        void clear_loads();
        
        /*!
         *    sets a pointer to the static solution system that provides
         *    deformation and its sensitivity for modal and buckling solution
         */
        void set_static_solution_system(System* sys) {
            // make sure it hasn't already not been set
            libmesh_assert(!_static_sol_system);
            _static_sol_system = sys;
        }
        
        /*!
         *    @returns a pointer to the static solution system that provides
         *    deformation and its sensitivity for modal and buckling solution
         */
        System* static_solution_system() {
            return _static_sol_system;
        }
        
        /*!
         *   adds the specified side loads for the boudnary with tag \p b_id
         */
        void add_side_load(boundary_id_type bid, MAST::BoundaryCondition& load);

        
        /*!
         *    returns a reference to the side boundary conditions
         */
        const std::multimap<boundary_id_type, MAST::BoundaryCondition*>& side_loads() const{
            return _side_bc_map;
        }
        
        
        /*!
         *   adds the specified volume loads for the elements with
         *   subdomain tag \p s_id
         */
        void add_volume_load(subdomain_id_type bid, MAST::BoundaryCondition& load);


        /*!
         *    returns a reference to the bolume boundary conditions
         */
        const std::multimap<subdomain_id_type, MAST::BoundaryCondition*>& volume_loads() const{
            return _vol_bc_map;
        }

        /*!
         *    fills the set \par dof_ids with the dof ids of the Dirichlet
         *    dofs
         */
        void get_dirichlet_dofs(std::set<unsigned int>& dof_ids) const;
        
        /*!
         *    sets the same property for all elements in the specified subdomain
         */
        void set_property_for_subdomain(const subdomain_id_type sid,
                                        const MAST::ElementPropertyCardBase& prop);
        
        /*!
         *    get the material property for the specified element
         */
        const MAST::ElementPropertyCardBase& get_property_card(const Elem& elem) const;

        /*!
         *    get the material property for the specified subdomain id \par i
         */
        const MAST::ElementPropertyCardBase& get_property_card(const unsigned int i) const;

        /*!
         *   Adds the parameter and function pairing
         */
        void add_parameter(MAST::ConstantFunction<Real>& f);
        
        /*!
         *   Returns the function corresponding to a parameter
         */
        const MAST::FieldFunctionBase* get_parameter(Real* par) const;
        
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
        
        /*!
         *    function to assemble the external forces of type specified in the 
         *    set
         */
        virtual void assemble_small_disturbance_aerodynamic_force (const NumericVector<Number>& X,
                                                                   NumericVector<Number>& F);

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

        /*!
         *    calculates the maximum element stress for the displacement provided
         *    in \par disp. This method looks at von Mises stresses for all quadrature
         *    points on the element, and returns the one with the maximum value.
         */
        virtual void calculate_max_elem_stress(const NumericVector<Number>& X,
                                               std::vector<Real>& stress,
                                               const MAST::FieldFunctionBase* params);


    protected:
        
        /*!
         *    assembles the residual and Jacobian for static or dynamic analyses
         */
        void _assemble_residual_and_jacobian (const NumericVector<Number>& X,
                                              NumericVector<Number>* R,
                                              SparseMatrix<Number>*  J,
                                              NonlinearImplicitSystem& S,
                                              const MAST::FieldFunctionBase* params);

        
        /*!
         *    Calculates the A and B matrices of a modal analysis eigenproblem
         *    [A] {X} = [B] {X} [lambda]
         *
         *    If \par params is non-NULL, then matrix sensitivities are assembled.
         *    \par static_sol is a pointer to the static deformation for 
         *    calculation of geometric stresses. If this is NULL, then no 
         *    deformation is assumed. If \par params and \par static_sol are 
         *    provided, then \par static_sol_sens should also be given.
         */
        void _assemble_matrices_for_modal_analysis(SparseMatrix<Number>&  matrix_A,
                                                   SparseMatrix<Number>&  matrix_B,
                                                   const MAST::FieldFunctionBase* params,
                                                   const NumericVector<Number>* static_sol,
                                                   const NumericVector<Number>* static_sol_sens);
        
        /*!
         *    Calculates the A and B matrices of a buckling analysis eigenproblem
         *    [A] {X} = [B] {X} [lambda]
         *
         *    If \par params is non-NULL, then matrix sensitivities are assembled.
         *    \par static_sol is a pointer to the static deformation for
         *    calculation of geometric stresses. If this is NULL, then no
         *    deformation is assumed. If \par params and \par static_sol are
         *    provided, then \par static_sol_sens should also be given.
         */
        void _assemble_matrices_for_buckling_analysis(SparseMatrix<Number>&  matrix_A,
                                                      SparseMatrix<Number>&  matrix_B,
                                                      const MAST::FieldFunctionBase* params,
                                                      const NumericVector<Number>* static_sol,
                                                      const NumericVector<Number>* static_sol_sens);

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
         *    System that provides the static solution for calculation of geometric
         *    stress. This is used in modal and buckling analysis
         */
        System *_static_sol_system;
        
        /*!
         *   map of element property cards for each element
         */
        std::map<subdomain_id_type, const MAST::ElementPropertyCardBase*> _element_property;
        
        /*!
         *   vector of material property cards
         */
        std::vector<const MAST::MaterialPropertyCardBase*> _material_property;
        
        /*!
         *   map of sensitivity parameters and the corresponding functions that
         *   are directly dependent on these parameters
         */
        std::map<Real*, const MAST::FieldFunctionBase*> _parameter_map;
        
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
