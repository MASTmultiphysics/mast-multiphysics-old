//
//  surface_pressure_load.h
//  MAST
//
//  Created by Manav Bhatia on 8/1/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef MAST_surface_pressure_load_h
#define MAST_surface_pressure_load_h

// MAST includes
#include "Flight/flight_condition.h"
#include "FluidElems/fluid_elem_base.h"

// libMesh includes
#include "libmesh/mesh_function.h"
#include "libmesh/system.h"
#include "libmesh/dof_map.h"


class SurfacePressureLoad
{
public:
    SurfacePressureLoad(System& sys_nonlinar,
                        System& sys_linear):
    nonlinear_fluid_system(sys_nonlinar),
    linearized_fluid_system(sys_linear)
    { }
    
    virtual ~SurfacePressureLoad()
    { }
        
    /*!
     *   nonlinear fluid system associated with the mesh and solution vector
     */
    System& nonlinear_fluid_system;

    /*!
     *   linearized fluid system to provide the small disturbance solution
     */
    System& linearized_fluid_system;

    
    virtual void init(FlightCondition& flight,
                      NumericVector<Number>& sol_nonlinear,
                      NumericVector<Number>& sol_linear);
    
    // calculation in frequency domain
    virtual void surface_pressure(const Point& p,
                                  Number& cp, Number& dcp);
    
protected:
    
    /*!
     *   mesh function that interpolates the nonlinear solution
     */
    std::auto_ptr<MeshFunction> _function_nonlinear;

    /*!
     *   mesh function that interpolates the linearized solution
     */
    std::auto_ptr<MeshFunction> _function_linear;

    /*!
     *    numeric vector that stores the solution for nonlinear system
     */
    std::auto_ptr<NumericVector<Number> > _sol_nonlinear;

    /*!
     *    numeric vector that stores the solution for linearized system
     */
    std::auto_ptr<NumericVector<Number> > _sol_linear;
    
    /*!
     *    this provides the fluid values for calculation of cp
     */
    FlightCondition* _flt_cond;


};



void
SurfacePressureLoad::init(FlightCondition& flight,
                            NumericVector<Number>& sol_nonlinear,
                            NumericVector<Number>& sol_linear)
{
    // copy the pointer for flight condition data
    _flt_cond = &flight;
    
    // first initialize the solution to the given vector
    if (!_sol_nonlinear.get())
    {
        // vector to store the nonlinear solution
        _sol_nonlinear.reset
        (NumericVector<Number>::build(nonlinear_fluid_system.comm()).release());
        _sol_nonlinear->init(sol_nonlinear.size(), true, SERIAL);

        // vector to store the linear solution
        _sol_linear.reset
        (NumericVector<Number>::build(linearized_fluid_system.comm()).release());
        _sol_linear->init(sol_linear.size(), true, SERIAL);
    }
    
    // now localize the give solution to this objects's vector
    sol_nonlinear.localize(*_sol_nonlinear,
                           nonlinear_fluid_system.get_dof_map().get_send_list());
    sol_linear.localize(*_sol_linear,
                        linearized_fluid_system.get_dof_map().get_send_list());
    
    // if the mesh function has not been created so far, initialize it
    if (!_function_nonlinear.get())
    {
        // initialize the nonlinear solution
        std::vector<unsigned int> vars;
        nonlinear_fluid_system.get_all_variable_numbers(vars);
        _function_nonlinear.reset
        (new MeshFunction( nonlinear_fluid_system.get_equation_systems(),
                          *_sol_nonlinear, nonlinear_fluid_system.get_dof_map(),
                          vars));
        _function_nonlinear->init();
        
        // now initialize the linearized fluid system
        vars.clear();
        linearized_fluid_system.get_all_variable_numbers(vars);
        _function_linear.reset
        (new MeshFunction( linearized_fluid_system.get_equation_systems(),
                          *_sol_linear, linearized_fluid_system.get_dof_map(),
                          vars));
        _function_linear->init();
    }
}



void
SurfacePressureLoad::surface_pressure(const Point& p,
                                        Number& cp, Number& dcp)
{
    libmesh_assert(_function_nonlinear.get()); // should be initialized before this call

    // get the nonlinear and linearized solution
    DenseVector<Number> v_nonlin, v_lin;
    DenseVector<Real> v_nonlin_real;
    (*_function_nonlinear)(p, 0., v_nonlin);
    (*_function_linear)(p, 0., v_lin);
    
    PrimitiveSolution p_sol;
    SmallPerturbationPrimitiveSolution<Number> delta_p_sol;
    
    v_nonlin_real.resize(v_nonlin.size());
    for (unsigned int i=0; i<v_nonlin.size(); i++)
        v_nonlin_real(i) = real(v_nonlin(i));
    
    // now initialize the primitive variable contexts
    p_sol.init(p.size(), v_nonlin_real,
               _flt_cond->gas_property.cp,
               _flt_cond->gas_property.cv,
               false);
    delta_p_sol.init(p_sol, v_lin);
    
    cp = p_sol.c_pressure(_flt_cond->gas_property.pressure,
                          _flt_cond->q0());
    dcp = delta_p_sol.c_pressure(_flt_cond->q0());
}


#endif
