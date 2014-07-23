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
    virtual bool get_aero_operator_matrix(Real k_red, ComplexMatrixX& a);

    
    /*!
     *    updates the aerodynamic matrix operator sensitivity in \par a for the
     *    reduced frequency \par k_red. This method projects the aero
     *    matrices onto the structural degrees of freedom, so needs the
     *    coupling matrices from Structures->Fluid and Fluid->Structures.
     */
    virtual bool
    get_aero_operator_matrix_sensitivity(const libMesh::ParameterVector& params,
                                         unsigned int p,
                                         Real k_red, ComplexMatrixX& a) {
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
    get_aero_operator_matrix_sensitivity_for_reduced_freq(Real k_red,
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



//#include "libmesh/mesh_function.h"
#include "libmesh/exodusII_io.h"

inline
bool
CoupledFluidStructureSystem::get_aero_operator_matrix(Real k_red,
                                                      ComplexMatrixX& a)
{
    // get references to the structural and fluid models
    FEMStructuralModel& structure =
    dynamic_cast<FEMStructuralModel&> (structural_model);
    CFDAerodynamicModel& aero =
    dynamic_cast<CFDAerodynamicModel&> (aerodynamic_model);

    // get the structural basis
    BasisMatrix<Real>& structural_basis = structure.get_basis_matrix();

    if (!structure.structural_system.have_vector("fvec_real")) {
        structure.structural_system.add_vector("fvec_real");
        structure.structural_system.add_vector("fvec_imag");
    }
    libMesh::NumericVector<Real>
    &f_vec_real = structure.structural_system.get_vector("fvec_real"),
    &f_vec_imag = structure.structural_system.get_vector("fvec_imag");
    
    ComplexVectorX projected_force;
    a.setZero(structural_basis.n(), structural_basis.n());
    Complex iota(0., 1.);
    
    for (unsigned int j_basis=0; j_basis<structural_basis.n(); j_basis++)
    {
        projected_force.setZero(structural_basis.n());
        surface_motion->init(k_red, 0., structural_basis.basis(j_basis));
        aero.linearized_fluid_system.perturbed_surface_motion = surface_motion.get();
        aero.linearized_fluid_system.solve(); //  X_F = J_FF^{-1} A_SF Phi
        
        std::vector<unsigned int> vars(2), dval(1);
        vars[0] = aero.linearized_fluid_system.variable_number("drho_re");
        vars[1] = aero.linearized_fluid_system.variable_number("drho_im");
        libMesh::MeshFunction function( aero.linearized_fluid_system.get_equation_systems(),
                              *aero.linearized_fluid_system.solution,
                              aero.linearized_fluid_system.get_dof_map(), vars);
        function.init();
        DenseRealVector sol; sol.resize(2);
        libMesh::Point pt;
        for (unsigned int i=0; i<300; i++) {
            pt(0) = 0 + (6)*(1.*i)/299.;
            function(pt, 0., sol);
            libMesh::out
            << std::setw(15) << pt(0)
            << std::setw(15) << sol(0)
            << std::setw(15) << sol(1)
            << std::setw(15) << std::endl;
        }

        std::ostringstream file_name;
        file_name
        << nm << "_"
        << "out_"
        << std::setw(3)
        << std::setfill('0')
        << std::right
        << j_basis
        << ".exo";
        libMesh::ExodusII_IO(aero.linearized_fluid_system.get_mesh()).write_equation_systems
        (file_name.str(), aero.linearized_fluid_system.get_equation_systems());

        //libmesh_error();
        
        surface_pressure->init(*aero.nonlinear_fluid_system.solution,
                               *aero.linearized_fluid_system.solution);
        
        structure.assemble_force_vec(*surface_pressure,
                                     *surface_motion,
                                     f_vec_real,
                                     f_vec_imag); // A_FS X_F

        // Phi^T A_FS X_F
        // the real part
        structure.basis_matrix->vector_mult_transpose
        (projected_force, f_vec_real);

        a.col(j_basis) = projected_force;
        
        // the imaginary part is scaled with iota and then added
        structure.basis_matrix->vector_mult_transpose
        (projected_force, f_vec_imag);
        projected_force *= iota;
        a.col(j_basis) += projected_force;
    }
    
    return true;
}




inline
bool
CoupledFluidStructureSystem::
get_aero_operator_matrix_sensitivity_for_reduced_freq(Real k_red,
                                                      ComplexMatrixX& a) {
    
    // get references to the structural and fluid models
    FEMStructuralModel& structure =
    dynamic_cast<FEMStructuralModel&> (structural_model);
    CFDAerodynamicModel& aero =
    dynamic_cast<CFDAerodynamicModel&> (aerodynamic_model);
    
    // get the structural basis
    BasisMatrix<Real>& structural_basis = structure.get_basis_matrix();
    
    if (!structure.structural_system.have_vector("fvec_real")) {
        structure.structural_system.add_vector("fvec_real");
        structure.structural_system.add_vector("fvec_imag");
    }
    libMesh::NumericVector<Real>
    &f_vec_real = structure.structural_system.get_vector("fvec_real"),
    &f_vec_imag = structure.structural_system.get_vector("fvec_imag");
    
    ComplexVectorX projected_force;
    a.setZero(structural_basis.n(), structural_basis.n());
    Complex iota(0., 1.);
    
    for (unsigned int j_basis=0; j_basis<structural_basis.n(); j_basis++)
    {
        projected_force.setZero(structural_basis.n());
        surface_motion->init(k_red, 0., structural_basis.basis(j_basis));
        aero.linearized_fluid_system.perturbed_surface_motion = surface_motion.get();
        
        // first solve the basic solution
        aero.linearized_fluid_system.solve(); //  X_F = J_FF^{-1} A_SF Phi
        std::auto_ptr<libMesh::NumericVector<Real> >
        base_sol(aero.linearized_fluid_system.solution->clone().release());
        aero.linearized_fluid_system.base_solution = base_sol.get();
        
        // now solve the sensitivity problem
        aero.linearized_fluid_system.solution->zero();
        aero.linearized_fluid_system.if_k_red_sensitivity = true;
        aero.linearized_fluid_system.solve(); //  X_F = J_FF^{-1} A_SF Phi
        aero.linearized_fluid_system.if_k_red_sensitivity = false;
        
        std::vector<unsigned int> vars(2), dval(1);
        vars[0] = aero.linearized_fluid_system.variable_number("drho_re");
        vars[1] = aero.linearized_fluid_system.variable_number("drho_im");
        libMesh::MeshFunction function( aero.linearized_fluid_system.get_equation_systems(),
                                       *aero.linearized_fluid_system.solution,
                                       aero.linearized_fluid_system.get_dof_map(), vars);
        function.init();
        DenseRealVector sol; sol.resize(2);
        libMesh::Point pt;
        for (unsigned int i=0; i<300; i++) {
            pt(0) = 0 + (6)*(1.*i)/299.;
            function(pt, 0., sol);
            libMesh::out
            << std::setw(15) << pt(0)
            << std::setw(15) << sol(0)
            << std::setw(15) << sol(1)
            << std::setw(15) << std::endl;
        }
        
        std::ostringstream file_name;
        file_name
        << nm << "_"
        << "out_"
        << std::setw(3)
        << std::setfill('0')
        << std::right
        << j_basis
        << ".exo";
        libMesh::ExodusII_IO(aero.linearized_fluid_system.get_mesh()).write_equation_systems
        (file_name.str(), aero.linearized_fluid_system.get_equation_systems());
        
        //libmesh_error();
        
        surface_pressure->init(*aero.nonlinear_fluid_system.solution,
                               *aero.linearized_fluid_system.solution);
        
        structure.assemble_force_vec(*surface_pressure,
                                     *surface_motion,
                                     f_vec_real,
                                     f_vec_imag); // A_FS X_F
        
        // Phi^T A_FS X_F
        // the real part
        structure.basis_matrix->vector_mult_transpose
        (projected_force, f_vec_real);
        
        a.col(j_basis) = projected_force;
        
        // the imaginary part is scaled with iota and then added
        structure.basis_matrix->vector_mult_transpose
        (projected_force, f_vec_imag);
        projected_force *= iota;
        a.col(j_basis) += projected_force;
    }
    
    return true;
}




#endif
