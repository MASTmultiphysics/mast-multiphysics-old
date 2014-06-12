//
//  surface_pressure_load.cpp
//  MAST
//
//  Created by Manav Bhatia on 11/29/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//


// MAST includes
#include "BoundaryConditions/surface_pressure.h"
#include "Numerics/utility.h"


void
MAST::SmallDisturbanceSurfacePressure::init(libMesh::NumericVector<libMesh::Real>& nonlinear_sol,
                                            libMesh::NumericVector<libMesh::Real>& linearized_sol)
{
//#ifdef LIBMESH_USE_COMPLEX_NUMBERS
    libmesh_assert_equal_to(nonlinear_sol.size(),
                            nonlinear_sys.solution->size());
    libmesh_assert_equal_to(linearized_sol.size(),
                            linearized_sys.solution->size());
    
    
    // first initialize the solution to the given vector
    if (!_sol_nonlinear.get())
    {
        // vector to store the nonlinear solution
        _sol_nonlinear.reset
        (libMesh::NumericVector<libMesh::Real>::build(nonlinear_sys.comm()).release());
        _sol_nonlinear->init(nonlinear_sol.size(), true, SERIAL);
        
        // vector to store the linear solution
        _sol_linear.reset
        (libMesh::NumericVector<libMesh::Real>::build(linearized_sys.comm()).release());
        _sol_linear->init(linearized_sol.size(), true, SERIAL);
    }
    
    // now localize the give solution to this objects's vector
    nonlinear_sol.localize(*_sol_nonlinear);
    linearized_sol.localize(*_sol_linear);
    
    // if the mesh function has not been created so far, initialize it
    if (!_function_nonlinear.get())
    {
        // initialize the nonlinear solution
        unsigned int n_vars = nonlinear_sys.n_vars();
        libmesh_assert(linearized_sys.n_vars() == n_vars);
        
        std::vector<unsigned int> vars(_dim+2);
        vars[0] = nonlinear_sys.variable_number("rho");
        vars[1] = nonlinear_sys.variable_number("rhoux");
        if (_dim > 1)
            vars[2] = nonlinear_sys.variable_number("rhouy");
        if (_dim > 2)
            vars[3] = nonlinear_sys.variable_number("rhouz");
        vars[_dim+2-1] = nonlinear_sys.variable_number("rhoe");
        
        
        _function_nonlinear.reset
        (new MeshFunction( nonlinear_sys.get_equation_systems(),
                          *_sol_nonlinear, nonlinear_sys.get_dof_map(),
                          vars));
        _function_nonlinear->init();
        
        // now initialize the linearized fluid system
        vars.resize(_dim+2);
        vars[0] = linearized_sys.variable_number("drho");
        vars[1] = linearized_sys.variable_number("drhoux");
        if (_dim > 1)
            vars[2] = linearized_sys.variable_number("drhouy");
        if (_dim > 2)
            vars[3] = linearized_sys.variable_number("drhouz");
        vars[_dim +2-1] = linearized_sys.variable_number("drhoe");
        
        _function_linear.reset
        (new MeshFunction( linearized_sys.get_equation_systems(),
                          *_sol_linear,
                          linearized_sys.get_dof_map(),
                          vars));
        _function_linear->init();
    }
//#endif // LIBMESH_USE_COMPLEX_NUMBERS
}




template <typename ValType>
void
MAST::SmallDisturbanceSurfacePressure::surface_pressure(const libMesh::Real t,
                                                        const libMesh::Point& p,
                                                        libMesh::Real& cp,
                                                        ValType& dcp)
{
//#ifdef LIBMESH_USE_COMPLEX_NUMBERS
    libmesh_assert(_function_nonlinear.get()); // should be initialized before this call
    libmesh_assert(_function_linear.get()); // should be initialized before this call
    
    
    cp = 0.; dcp = 0.;
    
    // get the nonlinear and linearized solution
    libMesh::DenseVector<ValType> v_lin;
    libMesh::DenseVector<libMesh::Real> v_nonlin, v_lin_real;
    (*_function_nonlinear)(p, t, v_nonlin);
    (*_function_linear)(p, t, v_lin_real);
    
    MAST::transform_to_elem_vector(v_lin, v_lin_real);
    
    PrimitiveSolution p_sol;
    SmallPerturbationPrimitiveSolution<ValType> delta_p_sol;
    
    // now initialize the primitive variable contexts
    p_sol.init(_dim, v_nonlin,
               _flt_cond->gas_property.cp,
               _flt_cond->gas_property.cv,
               false);
    delta_p_sol.init(p_sol, v_lin);
    
    cp = p_sol.c_pressure(_flt_cond->p0(),
                          _flt_cond->q0());
    dcp = delta_p_sol.c_pressure(_flt_cond->q0());
//#endif // LIBMESH_USE_COMPLEX_NUMBERS
}



// template explicit instantiations
template
void MAST::SmallDisturbanceSurfacePressure::surface_pressure<libMesh::Real>(const libMesh::Real t,
                                                                            const libMesh::Point& p,
                                                                            libMesh::Real& cp,
                                                                            libMesh::Real& dcp);

template
void MAST::SmallDisturbanceSurfacePressure::surface_pressure<libMesh::Complex>(const libMesh::Real t,
                                                                               const libMesh::Point& p,
                                                                               libMesh::Real& cp,
                                                                               libMesh::Complex& dcp);



