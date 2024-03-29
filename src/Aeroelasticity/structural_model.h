//
//  structural_model.h
//  MAST
//
//  Created by Manav Bhatia on 7/27/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef MAST_structural_model_h
#define MAST_structural_model_h

// MAST includes
#include "Base/MAST_data_types.h"

// libMesh includes
#include "libmesh/parameter_vector.h"


class StructuralModel
{
public:
    StructuralModel()
    { }
    
    virtual ~StructuralModel()
    { }
    
    /*!
     *   @returns the number of degrees of freedom in the structural model
     */
    virtual unsigned int n_dofs() const = 0;
    
    /*!
     *    updates the matrix to the mass matrix for the structural model.
     *    Returns false if the matrix does not exist for this model
     */
    virtual bool get_mass_matrix(RealMatrixX& m) = 0;

    /*!
     *    updates the matrix to the sensitivity of mass matrix for
     *    the structural model.
     *    Returns false if the matrix does not exist for this model
     */
    virtual bool
    get_mass_matrix_sensitivity(const libMesh::ParameterVector& params,
                                unsigned int p,
                                RealMatrixX& m) = 0;

    /*!
     *    updates the matrix to the stiffness matrix for this structural model
     *    Returns false if the matrix does not exist for this model
     */
    virtual bool get_stiffness_matrix(RealMatrixX& k) = 0;
    
    /*!
     *    updates the matrix to the sensitivity of stiffness matrix for
     *    the structural model.
     *    Returns false if the matrix does not exist for this model
     */
    virtual bool
    get_stiffness_matrix_sensitivity(const libMesh::ParameterVector& params,
                                     unsigned int p,
                                     RealMatrixX& m) = 0;

    /*!
     *    updates the matrix to the damping matrix for this structural model
     *    Returns false if the matrix does not exist for this model
     */
    virtual bool get_damping_matrix(RealMatrixX& c) = 0;
    
    /*!
     *    updates the matrix to the sensitivity of damping matrix for
     *    the structural model.
     *    Returns false if the matrix does not exist for this model
     */
    virtual bool
    get_damping_matrix_sensitivity(const libMesh::ParameterVector& params,
                                   unsigned int p,
                                   RealMatrixX& m) = 0;

};


#endif
