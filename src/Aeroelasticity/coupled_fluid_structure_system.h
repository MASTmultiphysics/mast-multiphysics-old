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

#ifndef MAST_coupled_fluid_structure_system_h
#define MAST_coupled_fluid_structure_system_h

// MAST includes
#include "Aeroelasticity/coupled_aero_structural_model.h"
#include "Aeroelasticity/fem_structural_model.h"
#include "Aeroelasticity/cfd_aerodynamic_model.h"

// libMesh includes
#include "libmesh/equation_systems.h"


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
     *    reduced frequency \par k_red. This method projects the aero
     *    matrices onto the structural degrees of freedom, so needs the
     *    coupling matrices from Structures->Fluid and Fluid->Structures.
     */
    virtual bool get_aero_operator_matrix(const Real k_red,
                                          const Real v_ref,
                                          ComplexMatrixX& a);

    
    /*!
     *    updates the aerodynamic matrix operator sensitivity in \par a for the
     *    reduced frequency \par k_red. This method projects the aero
     *    matrices onto the structural degrees of freedom, so needs the
     *    coupling matrices from Structures->Fluid and Fluid->Structures.
     */
    virtual bool
    get_aero_operator_matrix_sensitivity(const libMesh::ParameterVector& params,
                                         unsigned int p,
                                         const Real k_red,
                                         const Real v_ref,
                                         ComplexMatrixX& a) {
        // get the structural basis
        BasisMatrix<Real>& structural_basis =
        dynamic_cast<FEMStructuralModel&>(structural_model).get_basis_matrix();
        
        a.setZero(structural_basis.n(), structural_basis.n());
        
        // this needs to be implemented
        
        return true;
    }
    
    /*!
     *    updates the aerodynamic matrix operator sensitivity wrt reduced
     *    frequency in \par a for the
     *    reduced frequency \par k_red. This method projects the aero
     *    matrices onto the structural degrees of freedom, so needs the
     *    coupling matrices from Structures->Fluid and Fluid->Structures.
     */
    virtual bool
    get_aero_operator_matrix_sensitivity_for_reduced_freq(const Real k_red,
                                                          const Real v_ref,
                                                          ComplexMatrixX& a);
    

    /*!
     *    updates the aerodynamic matrix operator sensitivity wrt reference
     *    velocity in \par a for the
     *    reduced frequency \par k_red. This method projects the aero
     *    matrices onto the structural degrees of freedom, so needs the
     *    coupling matrices from Structures->Fluid and Fluid->Structures.
     */
    virtual bool
    get_aero_operator_matrix_sensitivity_for_V_ref(const Real k_red,
                                                   const Real v_ref,
                                                   ComplexMatrixX& a);
    

    std::string nm;
    
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



#endif
