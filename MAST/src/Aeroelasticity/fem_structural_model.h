//
//  fem_structural_model.h
//  MAST
//
//  Created by Manav Bhatia on 7/28/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

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
     *    updates the matrix to the damping matrix for this structural model
     *    Returns false if the matrix does not exist for this model
     */
    virtual bool get_damping_matrix(RealMatrixX& c)
    { libmesh_assert(false); }
    
    /*!
     *     returns the basis matrix for this structural model
     */
    virtual BasisMatrix<libMesh::Real>& get_basis_matrix()
    {
        return *basis_matrix;
    }

    /*!
     *    calculates force vector for the specified unsteady pressure and
     *    motion
     */
    void assemble_force_vec(MAST::SmallDisturbanceSurfacePressure& press,
                            MAST::SurfaceMotionBase& displ,
                            libMesh::NumericVector<libMesh::Real>& f_vec_real,
                            libMesh::NumericVector<libMesh::Real>& f_vec_imag);

    
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
    std::auto_ptr<BasisMatrix<libMesh::Real> > basis_matrix;
};



inline
void
FEMStructuralModel::init()
{
    libmesh_assert(eigen_vals.size() > 0);
    
    unsigned int n_eig = eigen_vals.size();
    
    basis_matrix.reset(new BasisMatrix<libMesh::Real>(structural_system.comm()));
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
                                       libMesh::NumericVector<libMesh::Real>& f_vec_real,
                                       libMesh::NumericVector<libMesh::Real>& f_vec_imag) {
    MAST::SmallDisturbanceMotion load;
    load.set_deformation(displ);
    load.set_pressure(press);
    assembly.clear_loads();
    assembly.add_volume_load(0, load);
    assembly.assemble_small_disturbance_aerodynamic_force(*structural_system.solution,
                                                          f_vec_real,
                                                          f_vec_imag);
}



#endif
