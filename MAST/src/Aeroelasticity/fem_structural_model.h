//
//  fem_structural_model.h
//  MAST
//
//  Created by Manav Bhatia on 7/28/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef MAST_fem_structural_model_h
#define MAST_fem_structural_model_h

// MAST includes
#include "Aeroelasticity/structural_model.h"
#include "StructuralElems/structural_system_base.h"
#include "Numerics/basis_matrix.h"


class FEMStructuralModel: public StructuralModel
{
public:
    FEMStructuralModel(System& system):
    StructuralModel(),
    structural_system(system)
    { }
    
    virtual ~FEMStructuralModel()
    { }
    
    
    /*!
     *    updates the matrix to the mass matrix for the structural model.
     *    Returns false if the matrix does not exist for this model
     */
    virtual bool get_mass_matrix(RealMatrixX& m)
    {
        libmesh_assert(false);
    }
    
    /*!
     *    updates the matrix to the stiffness matrix for this structural model
     *    Returns false if the matrix does not exist for this model
     */
    virtual bool get_stiffness_matrix(RealMatrixX& k)
    {
        libmesh_assert(false);
    }
    
    /*!
     *    updates the matrix to the damping matrix for this structural model
     *    Returns false if the matrix does not exist for this model
     */
    virtual bool get_damping_matrix(RealMatrixX& c)
    {
        libmesh_assert(false);
    }
    
    /*!
     *     returns the basis matrix for this structural model
     */
    virtual BasisMatrix<Number>& get_basis_matrix()
    {
        libmesh_assert(false);
    }


    /*!
     *    the structural system that provides the basis of
     *    calculations for this model
     */
    System& structural_system;
};


#endif
