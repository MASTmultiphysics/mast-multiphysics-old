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

    ComplexVectorX projected_vec;
    
    // get the structural basis
    BasisMatrix<Number>& structural_basis = structure.get_basis_matrix();
    
    // temporary vectors for fluid and structre
    AutoPtr<NumericVector<Number> > fluid_tmp_vec =
    NumericVector<Number>::build
    (aero.fluid_system.get_equation_systems().comm());

    AutoPtr<NumericVector<Number> > structural_tmp_vec =
    NumericVector<Number>::build
    (structure.structural_system.get_equation_systems().comm());
    
    
    for (unsigned int j_basis=0; j_basis<structural_basis.n(); j_basis++)
    {
        NumericVector<Number>& basis_vec = structural_basis.basis(j_basis);
        
//        structure_to_fluid_mapping->vector_mult
//        ( *fluid_tmp_vec, basis_vec ); // F_A = A_SF X_S
        aero.solve(*fluid_tmp_vec); // X_S = J_FF^{-1} F_A
//        fluid_to_structure_mapping->vector_mult
//        ( *structural_tmp_vec, *aero.fluid_system.solution ); // F_S = A_FS X_S
        
        // now calculate the dot product
        for (unsigned int i_basis=0; i_basis<structural_basis.n(); i_basis++)
            a(i_basis, j_basis) = structural_basis.basis(i_basis).dot
            (*structural_tmp_vec); // Phi^T F_S
    }
    
    return true;
}

#endif // LIBMESH_USE_COMPLEX_NUMBERS


#endif
