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

#ifndef MAST_coupled_aero_structural_model_h
#define MAST_coupled_aero_structural_model_h

// C++ includes
#include <memory>

// MAST includes
#include "Base/MAST_data_types.h"
#include "Aeroelasticity/structural_model.h"
#include "Aeroelasticity/aerodynamic_model.h"


// libMesh includes
#include "libmesh/parameter_vector.h"


class CoupledAeroStructuralModel
{
public:
    CoupledAeroStructuralModel(AerodynamicModel& aero,
                               StructuralModel& structure):
    aerodynamic_model(aero),
    structural_model(structure)
    { }
    
    /*!
     *    updates the structural mass operator in \par mat.
     */
    bool get_structural_mass_matrix(RealMatrixX& mat)
    {
        return structural_model.get_mass_matrix(mat);
    }

    /*!
     *    updates the sensitivity of structural mass operator in \par mat.
     */
    bool get_structural_mass_matrix_sensitivity(const libMesh::ParameterVector& params,
                                                unsigned int p,
                                                RealMatrixX& mat)
    {
        return structural_model.get_mass_matrix_sensitivity(params,
                                                            p,
                                                            mat);
    }

    /*!
     *    updates the structural damping operator in \par mat.
     */
    bool get_structural_damping_matrix(RealMatrixX& mat)
    {
        return structural_model.get_damping_matrix(mat);
    }

    /*!
     *    updates the sensitivity of structural damping operator in \par mat.
     */
    bool get_structural_damping_matrix_sensitivity(const libMesh::ParameterVector& params,
                                                   unsigned int p,
                                                   RealMatrixX& mat)
    {
        return structural_model.get_damping_matrix_sensitivity(params,
                                                               p,
                                                               mat);
    }

    /*!
     *    updates the structural stiffness operator in \par mat.
     */
    bool get_structural_stiffness_matrix(RealMatrixX& mat)
    {
        return structural_model.get_stiffness_matrix(mat);
    }

    /*!
     *    updates the sensitivity of structural stiffness operator in \par mat.
     */
    bool get_structural_stiffness_matrix_sensitivity(const libMesh::ParameterVector& params,
                                                     unsigned int p,
                                                     RealMatrixX& mat)
    {
        return structural_model.get_stiffness_matrix_sensitivity(params,
                                                                 p,
                                                                 mat);
    }
    
    /*!
     *    updates the aerodynamic matrix operator in \par a for the
     *    reduced frequency \par k_ref. This method projects the aero
     *    matrices onto the structural degrees of freedom, so needs the
     *    coupling matrices from Structures->Fluid and Fluid->Structures.
     */
    virtual bool get_aero_operator_matrix(Real k_ref, ComplexMatrixX& a)
    {
        // needs to be implemented by the inherited class
        libmesh_assert(false);
    }

    /*!
     *    updates the aerodynamic matrix operator sensitivity in \par a for the
     *    reduced frequency \par k_ref. This method projects the aero
     *    matrices onto the structural degrees of freedom, so needs the
     *    coupling matrices from Structures->Fluid and Fluid->Structures.
     */
    virtual bool get_aero_operator_matrix_sensitivity(const libMesh::ParameterVector& params,
                                                      unsigned int p,
                                                      Real k_ref, ComplexMatrixX& a)
    {
        // needs to be implemented by the inherited class
        libmesh_assert(false);
    }

    /*!
     *    updates the aerodynamic matrix operator sensitivity wrt reduced
     *    frequency in \par a for the
     *    reduced frequency \par k_ref. This method projects the aero
     *    matrices onto the structural degrees of freedom, so needs the
     *    coupling matrices from Structures->Fluid and Fluid->Structures.
     */
    virtual bool get_aero_operator_matrix_sensitivity_for_reduced_freq(Real k_ref,
                                                                       ComplexMatrixX& a)
    {
        // needs to be implemented by the inherited class
        libmesh_assert(false);
    }

    /*!
     *    updates the aerodynamic damping matrix operator in \par a. This is
     *    applicable for only the time-domain methods. Returns \par false if
     *    the model does not have a damping term.
     */
    virtual bool get_aero_damping_matrix(RealMatrixX& a)
    {
        return aerodynamic_model.get_damping_matrix(a);
    }
    
    /*!
     *    updates the aerodynamic stiffness matrix operator in \par a. This is
     *    applicable for only the time-domain methods. Returns \par false if
     *    the model does not have a stiffness term.
     */
    virtual bool get_aero_stiffness_matrix(RealMatrixX& a)
    {
        return aerodynamic_model.get_stiffness_matrix(a);
    }

protected:

    /*!
     *    reference to the aerodynamic model
     */
    AerodynamicModel& aerodynamic_model;

    /*!
     *   reference to the structural model
     */
    StructuralModel& structural_model;
    
};

#endif
