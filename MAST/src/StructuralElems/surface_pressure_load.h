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
#include "FluidElems/frequency_domain_linearized_fluid_system.h"
#include "FluidElems/fluid_elem_base.h"

// libMesh includes
#include "libmesh/mesh_function.h"
#include "libmesh/system.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh_serializer.h"


class SurfacePressureLoad
{
public:
    SurfacePressureLoad():
    _flt_cond(NULL),
    _dim(0)
    { }
    
    virtual ~SurfacePressureLoad()
    { }
    
    virtual void init(System& linearized_sys);
    
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
  
    /*!
     *    dimension of the analysis mesh
     */
    unsigned int _dim;

    /*!
     *   this serializes the mesh for use in interpolation
     */
    std::auto_ptr<MeshSerializer> _nonlinear_mesh_serializer;
    std::auto_ptr<MeshSerializer> _linear_mesh_serializer;

};



inline
void
SurfacePressureLoad::init(System& linearized_sys)
{
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
    
    MeshBase& linear_sys_mesh = linearized_sys.get_mesh();
    _linear_mesh_serializer.reset(new MeshSerializer(linear_sys_mesh, true));

    System& nonlin_sys =
    linearized_sys.get_equation_systems().get_system<System>("EulerSystem");
    MeshBase& nonlinear_sys_mesh = nonlin_sys.get_mesh();
    _nonlinear_mesh_serializer.reset(new MeshSerializer(nonlinear_sys_mesh, true));
    
    
    _dim = linearized_sys.n_vars()-2;
    
    // copy the pointer for flight condition data
    _flt_cond = dynamic_cast<FrequencyDomainLinearizedFluidSystem&>
    (linearized_sys).flight_condition;
    
    // first initialize the solution to the given vector
    if (!_sol_nonlinear.get())
    {
        // vector to store the nonlinear solution
        _sol_nonlinear.reset
        (NumericVector<Number>::build(nonlin_sys.comm()).release());
        _sol_nonlinear->init(nonlin_sys.solution->size(), true, SERIAL);

        // vector to store the linear solution
        _sol_linear.reset
        (NumericVector<Number>::build(linearized_sys.comm()).release());
        _sol_linear->init(linearized_sys.solution->size(), true, SERIAL);
    }
    
    // now localize the give solution to this objects's vector
    nonlin_sys.solution->localize
    (*_sol_nonlinear, nonlin_sys.get_dof_map().get_send_list());
    linearized_sys.solution->localize
    (*_sol_linear, linearized_sys.get_dof_map().get_send_list());
    
    // if the mesh function has not been created so far, initialize it
    if (!_function_nonlinear.get())
    {
        // initialize the nonlinear solution
        unsigned int n_vars = nonlin_sys.n_vars();
        libmesh_assert(linearized_sys.n_vars() == n_vars);
        
        std::vector<unsigned int> vars(_dim+2);
        vars[0] = nonlin_sys.variable_number("rho");
        vars[1] = nonlin_sys.variable_number("rhoux");
        if (_dim > 1)
            vars[2] = nonlin_sys.variable_number("rhouy");
        if (_dim > 2)
            vars[3] = nonlin_sys.variable_number("rhouz");
        vars[_dim+2-1] = nonlin_sys.variable_number("rhoe");
        
        
        _function_nonlinear.reset
        (new MeshFunction( nonlin_sys.get_equation_systems(),
                          *_sol_nonlinear, nonlin_sys.get_dof_map(),
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
#endif // LIBMESH_USE_COMPLEX_NUMBERS
}



inline
void
SurfacePressureLoad::surface_pressure(const Point& p,
                                        Number& cp, Number& dcp)
{
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
    libmesh_assert(_function_nonlinear.get()); // should be initialized before this call
    libmesh_assert(_function_linear.get()); // should be initialized before this call

    
    cp = 0.; dcp = 0.;
    
    // get the nonlinear and linearized solution
    DenseVector<Number> v_nonlin, v_lin;
    DenseVector<Real> v_nonlin_real;
    (*_function_nonlinear)(p, 0., v_nonlin);
    (*_function_linear)(p, 0., v_lin);
        
    PrimitiveSolution p_sol;
    SmallPerturbationPrimitiveSolution<Number> delta_p_sol;
    
    v_nonlin_real.resize(v_nonlin.size());
    for (unsigned int i=0; i<v_nonlin.size(); i++)
        v_nonlin_real(i) = std::real(v_nonlin(i));
    
    // now initialize the primitive variable contexts
    p_sol.init(_dim, v_nonlin_real,
               _flt_cond->gas_property.cp,
               _flt_cond->gas_property.cv,
               false);
    delta_p_sol.init(p_sol, v_lin);
    
    cp = p_sol.c_pressure(_flt_cond->p0(),
                          _flt_cond->q0());
    dcp = delta_p_sol.c_pressure(_flt_cond->q0());
#endif // LIBMESH_USE_COMPLEX_NUMBERS
}


#endif
