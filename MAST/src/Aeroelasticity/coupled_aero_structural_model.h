//
//  coupled_aero_structural_model.h
//  MAST
//
//  Created by Manav Bhatia on 7/27/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef MAST_coupled_aero_structural_model_h
#define MAST_coupled_aero_structural_model_h

// C++ includes
#include <memory>

// MAST includes
#include "Base/MAST_data_types.h"
#include "Aeroelasticity/structural_model.h"
#include "Aeroelasticity/aerodynamic_model.h"


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
     *    updates the structural damping operator in \par mat.
     */
    bool get_structural_damping_matrix(RealMatrixX& mat)
    {
        return structural_model.get_damping_matrix(mat);
    }

    /*!
     *    updates the structural stiffness operator in \par mat.
     */
    bool get_structural_stiffness_matrix(RealMatrixX& mat)
    {
        return structural_model.get_stiffness_matrix(mat);
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
     *    updates the aerodynamic damping matrix operator in \par a. This is
     *    applicable for only the time-domain methods. Returns \par false if
     *    the model does not have a damping term.
     */
    virtual bool get_aero_damping_matrix(ComplexMatrixX& a)
    {
        return aerodynamic_model.get_damping_matrix(a);
    }
    
    /*!
     *    updates the aerodynamic stiffness matrix operator in \par a. This is
     *    applicable for only the time-domain methods. Returns \par false if
     *    the model does not have a stiffness term.
     */
    virtual bool get_aero_stiffness_matrix(ComplexMatrixX& a)
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
