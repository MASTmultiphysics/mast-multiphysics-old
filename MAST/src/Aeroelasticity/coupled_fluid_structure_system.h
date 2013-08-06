//
//  coupled_fluid_structure_system.h
//  MAST
//
//  Created by Manav Bhatia on 7/28/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef MAST_coupled_fluid_structure_system_h
#define MAST_coupled_fluid_structure_system_h

// MAST includes
#include "Aeroelasticity/coupled_aero_structural_model.h"
#include "Aeroelasticity/fem_structural_model.h"
#include "Aeroelasticity/cfd_aerodynamic_model.h"

// libMesh includes
#include "libmesh/equation_systems.h"

#ifdef LIBMESH_USE_COMPLEX_NUMBERS


class CoupledFluidStructureSystem: public CoupledAeroStructuralModel
{
public:
    CoupledFluidStructureSystem(CFDAerodynamicModel& aero,
                                FEMStructuralModel& structure):
    CoupledAeroStructuralModel(aero, structure)
    { }

    
    virtual ~CoupledFluidStructureSystem()
    { }
    
    
    /*!
     *    updates the aerodynamic matrix operator in \par a for the
     *    reduced frequency \par k_ref. This method projects the aero
     *    matrices onto the structural degrees of freedom, so needs the
     *    coupling matrices from Structures->Fluid and Fluid->Structures.
     */
    virtual bool get_aero_operator_matrix(Real k_ref, ComplexMatrixX& a);

protected:
    
};




bool
CoupledFluidStructureSystem::get_aero_operator_matrix(Real k_ref,
                                                      ComplexMatrixX& a)
{
    // get references to the structural and fluid models
    FEMStructuralModel& structure =
    dynamic_cast<FEMStructuralModel&> (structural_model);
    CFDAerodynamicModel& aero =
    dynamic_cast<CFDAerodynamicModel&> (aerodynamic_model);

    // get the structural basis
    BasisMatrix<Number>& structural_basis = structure.get_basis_matrix();
    
    ComplexVectorX projected_force;
    
    for (unsigned int j_basis=0; j_basis<structural_basis.n(); j_basis++)
    {
        projected_force.Zero(structural_basis.n());
        aero.solve(k_ref, structural_basis.basis(j_basis)); // J_FF^{-1} F_A A_SF X_S
        structure.project_aero_force(*aero.fluid_system.solution,
                                     projected_force); // Phi^T A_FS X_F
        a.col(j_basis) = projected_force;
    }
    
    return true;
}

#endif // LIBMESH_USE_COMPLEX_NUMBERS


#endif
