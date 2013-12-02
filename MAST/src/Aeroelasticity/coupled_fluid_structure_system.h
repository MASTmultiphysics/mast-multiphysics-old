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
    CoupledAeroStructuralModel(aero, structure),
    surface_pressure(new MAST::SmallDisturbanceSurfacePressure(aero.nonlinear_fluid_system,
                                             aero.linearized_fluid_system)),
    surface_motion(new MAST::FlexibleSurfaceMotion(structure.structural_system))
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

    
    /*!
     *   this provides interpolation of the surface pressure
     */
    std::auto_ptr<MAST::SmallDisturbanceSurfacePressure> surface_pressure;
    
    
    /*!
     *   this provides the surface motion description
     */
    std::auto_ptr<MAST::FlexibleSurfaceMotion> surface_motion;
};


void assemble_force_vec(System& sys,
                        MAST::SmallDisturbanceSurfacePressure& press,
                        MAST::SurfaceMotionBase& motion,
                        NumericVector<Number>& fvec)
{ }

inline
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

    if (!structure.structural_system.have_vector("fvec"))
        structure.structural_system.add_vector("fvec");
    NumericVector<Number>& f_vec = structure.structural_system.get_vector("fvec");
    
    ComplexVectorX projected_force;
    a.setZero(structural_basis.n(), structural_basis.n());
    
    for (unsigned int j_basis=0; j_basis<structural_basis.n(); j_basis++)
    {
        projected_force.setZero(structural_basis.n());
        surface_motion->init(k_ref, 0., structural_basis.basis(j_basis));
        aero.linearized_fluid_system.perturbed_surface_motion = surface_motion.get();
        aero.linearized_fluid_system.solve(); //  X_F = J_FF^{-1} A_SF Phi

        surface_pressure->init(*aero.nonlinear_fluid_system.solution,
                               *aero.linearized_fluid_system.solution);
        
        structure.assemble_force_vec(*surface_pressure,
                                     *surface_motion,
                                     f_vec); // A_FS X_F

        structure.basis_matrix->vector_mult_transpose
        (projected_force, f_vec); // Phi^T A_FS X_F

        a.col(j_basis) = projected_force;
    }
    
    return true;
}

#endif // LIBMESH_USE_COMPLEX_NUMBERS


#endif
