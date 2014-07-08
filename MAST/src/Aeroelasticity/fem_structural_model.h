/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2014  Manav Bhatia
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef MAST_fem_structural_model_h
#define MAST_fem_structural_model_h

// C++ includes
#include <memory>

// MAST includes
#include "Aeroelasticity/structural_model.h"
#include "BoundaryConditions/small_disturbance_motion.h"
#include "Numerics/basis_matrix.h"
#include "StructuralElems/structural_system_assembly.h"


// libMesh includes
#include "libmesh/system.h"


class FEMStructuralModel: public StructuralModel
{
public:
    FEMStructuralModel(MAST::StructuralSystemAssembly& a):
    StructuralModel(),
    assembly(a),
    structural_system(a.get_system())
    { }
    
    virtual ~FEMStructuralModel()
    { }
    
    /*!
     *    initializes the data strucutres. The eigen values must be available
     *    in both the eigen_vals vector and the modes in structural_system
     */
    void init();
    
    /*!
     *    updates the matrix to the mass matrix for the structural model.
     *    Returns false if the matrix does not exist for this model
     */
    virtual bool get_mass_matrix(RealMatrixX& m)
    {
        unsigned int n_eig = eigen_vals.size();
        m.setZero(n_eig, n_eig);
        for (unsigned i=0; i<n_eig; i++)
            m(i,i) = 1.;
        
        return true;
    }
    
    
    /*!
     *    updates the matrix to the sensitivity of mass matrix for
     *    the structural model.
     *    Returns false if the matrix does not exist for this model
     */
    virtual bool get_mass_matrix_sensitivity(const libMesh::ParameterVector& params,
                                             unsigned int p,
                                             RealMatrixX& m) {
        libMesh::SparseMatrix<Real>& mat =
        *(dynamic_cast<libMesh::ImplicitSystem&>(structural_system).matrix);
        
        const MAST::FieldFunctionBase* f = assembly.get_parameter(params[p]);
        
        // calculate the mass matrix
        assembly.assemble_mass(mat,
                               f,
                               structural_system.solution.get(),
                               &(structural_system.get_sensitivity_solution(p)));

        std::auto_ptr<libMesh::NumericVector<Real> >
        vec(structural_system.solution->zero_clone().release());
        
        // now calculate the projection
        RealVectorX projected_force;
        m.setZero(basis_matrix->n(), basis_matrix->n());
        
        for (unsigned int j_basis=0; j_basis<basis_matrix->n(); j_basis++)
        {
            projected_force.setZero(basis_matrix->n());
            
            mat.vector_mult(*vec, basis_matrix->basis(j_basis));
            // Phi^T A_FS X_F
            basis_matrix->vector_mult_transpose(projected_force, *vec);
            
            m.col(j_basis) = projected_force;
        }
        
        return true;
    }

    /*!
     *    updates the matrix to the stiffness matrix for this structural model
     *    Returns false if the matrix does not exist for this model
     */
    virtual bool get_stiffness_matrix(RealMatrixX& k)
    {
        unsigned int n_eig = eigen_vals.size();
        k.setZero(n_eig, n_eig);
        for (unsigned i=0; i<n_eig; i++)
            k(i,i) = eigen_vals(i);
        
        return true;
    }
    
    /*!
     *    updates the matrix to the sensitivity of stiffness matrix for
     *    the structural model.
     *    Returns false if the matrix does not exist for this model
     */
    virtual bool get_stiffness_matrix_sensitivity(const libMesh::ParameterVector& params,
                                                  unsigned int p,
                                                  RealMatrixX& m) {
        // create vectors for use here
        RealVectorX projected_force;
        std::auto_ptr<libMesh::NumericVector<Real> >
        vec(structural_system.solution->zero_clone().release());
        
        m.setZero(basis_matrix->n(), basis_matrix->n());

        libMesh::SparseMatrix<Real>& mat =
        *(dynamic_cast<libMesh::ImplicitSystem&>(structural_system).matrix);
        
        const MAST::FieldFunctionBase* f = assembly.get_parameter(params[p]);
        libMesh::NumericVector<Real> *sol = structural_system.solution.get(),
        *sol_sens = &(structural_system.get_sensitivity_solution(p));
        
        // calculate the mass matrix
        mat.zero();
        assembly.assemble_jacobian(mat, f, sol, sol_sens);
        
        // now calculate the projection
        for (unsigned int j_basis=0; j_basis<basis_matrix->n(); j_basis++)
        {
            projected_force.setZero(basis_matrix->n());
            
            mat.vector_mult(*vec, basis_matrix->basis(j_basis));
            // Phi^T A_FS X_F
            basis_matrix->vector_mult_transpose(projected_force, *vec);
            
            m.col(j_basis) = projected_force;
        }
        
        // do the same for the nonlinear term
        // calculate the mass matrix
        mat.zero();
        assembly.assemble_jacobian_dot_state_sensitivity(mat, sol, sol_sens);
        
        // now calculate the projection
        for (unsigned int j_basis=0; j_basis<basis_matrix->n(); j_basis++)
        {
            projected_force.setZero(basis_matrix->n());
            
            mat.vector_mult(*vec, basis_matrix->basis(j_basis));
            // Phi^T A_FS X_F
            basis_matrix->vector_mult_transpose(projected_force, *vec);
            
            m.col(j_basis) = projected_force;
        }
        
        return true;
    }


    /*!
     *    updates the matrix to the damping matrix for this structural model
     *    Returns false if the matrix does not exist for this model
     */
    virtual bool get_damping_matrix(RealMatrixX& c)
    { libmesh_assert(false); }
    
    /*!
     *    updates the matrix to the sensitivity of damping matrix for
     *    the structural model.
     *    Returns false if the matrix does not exist for this model
     */
    virtual bool get_damping_matrix_sensitivity(const libMesh::ParameterVector& params,
                                                unsigned int p,
                                                RealMatrixX& m) {
        libmesh_error(); // to be implemented
    }
    
    /*!
     *     returns the basis matrix for this structural model
     */
    virtual BasisMatrix<Real>& get_basis_matrix()
    {
        return *basis_matrix;
    }

    /*!
     *    calculates force vector for the specified unsteady pressure and
     *    motion
     */
    void assemble_force_vec(MAST::SmallDisturbanceSurfacePressure& press,
                            MAST::SurfaceMotionBase& displ,
                            libMesh::NumericVector<Real>& f_vec_real,
                            libMesh::NumericVector<Real>& f_vec_imag);

    
    /*!
     *    vector of eigenvalues
     */
    RealVectorX  eigen_vals;

    
    /*!
     *    structural assembly object
     */
    MAST::StructuralSystemAssembly& assembly;

    /*!
     *    the structural system that provides the basis of
     *    calculations for this model
     */
    libMesh::System& structural_system;
    
    /*!
     *    returns the basis matrix using the modal data in structural system
     */
    std::auto_ptr<BasisMatrix<Real> > basis_matrix;
};



inline
void
FEMStructuralModel::init()
{
    libmesh_assert(eigen_vals.size() > 0);
    
    unsigned int n_eig = eigen_vals.size();
    
    basis_matrix.reset(new BasisMatrix<Real>(structural_system.comm()));
    basis_matrix->modes.resize(n_eig);
    
    for (unsigned int i=0; i<n_eig; i++)
    {
        std::ostringstream oss;
        oss << "mode_" << i;
        basis_matrix->modes[i] = &structural_system.get_vector(oss.str());
    }
}


inline
void
FEMStructuralModel::assemble_force_vec(MAST::SmallDisturbanceSurfacePressure& press,
                                       MAST::SurfaceMotionBase& displ,
                                       libMesh::NumericVector<Real>& f_vec_real,
                                       libMesh::NumericVector<Real>& f_vec_imag) {
    MAST::SmallDisturbanceMotion load;
    load.set_deformation(displ);
    load.set_pressure(press);
    
    assembly.add_volume_load(0, load);
    assembly.assemble_small_disturbance_aerodynamic_force(*structural_system.solution,
                                                          f_vec_real,
                                                          f_vec_imag);
    assembly.clear_volume_load(0, load);
}



#endif
