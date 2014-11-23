
#ifndef __mast_fluid_elem_initialization_h__
#define __mast_fluid_elem_initialization_h__

// C++ includes
#include <memory>

// MAST includes
#include "FluidElems/fluid_elem_base.h"
#include "Flight/flight_condition.h"

// libMesh includes
#include "libmesh/getpot.h"



struct BuildFluidElem
{
    
    BuildFluidElem() {
        
        //BOOST_TEST_MESSAGE("setup element");
        _cons_sol.reset         (new DenseRealVector(5));
        _primitive_sol.reset    (new PrimitiveSolution);
        _flight_condition.reset (new FlightCondition);
        _dummy_infile.reset     (new GetPot);
        _fluid_elem.reset       (new FluidElemBase(*_dummy_infile));
        
        // initialize the solution vector
        rho   = 1.1,
        u1    = 110.,
        u2    = 11.,
        u3    = 22.,
        T     = 330.,
        cp    = 1003.,
        cv    = 716.,
        R     = cp-cv,
        gamma = cp/cv,
        k     = .5*(u1*u1+u2*u2+u3*u3),
        tol   = 1.0e-12;
        
        // initialize the conservative solution
        (*_cons_sol)(0) =   rho;
        (*_cons_sol)(1) =   rho*u1;
        (*_cons_sol)(2) =   rho*u2;
        (*_cons_sol)(3) =   rho*u3;
        (*_cons_sol)(4) =   rho*(cv*T+k);
        
        // set the data needed for this calculation. This is a dirty
        // way of setting this data, and going forward, will need to be
        // defined in a legal way.
        
        // initialize the fluid parameter values in flight condition
        _flight_condition->gas_property.cp    = cp;
        _flight_condition->gas_property.cv    = cv;
        _flight_condition->gas_property.gamma = gamma;
        _flight_condition->gas_property.R     = R;
        
        // dimension in element
        _fluid_elem->dim = 3;
        
        
        _primitive_sol->init(3, *_cons_sol, cp, cv, false);
        _fluid_elem->flight_condition = _flight_condition.get();
    }
    
    ~BuildFluidElem() {
        
        //BOOST_TEST_MESSAGE("clearing data");
    }
    
    
    inline bool
    compare(const RealMatrixX& m0, const DenseRealMatrix& m) {
        
        unsigned int n = (unsigned int) m0.rows();
        libmesh_assert_equal_to(n,  m0.cols());
        libmesh_assert_equal_to(n,  m.m());
        libmesh_assert_equal_to(n,  m.n());
        
        
        bool pass = true;
        for (unsigned int i=0; i<n; i++) {
            for (unsigned int j=0; j<n; j++)
                if (!boost::test_tools::check_is_close(m0(i,j),
                                                       m(i,j),
                                                       boost::test_tools::percent_tolerance(tol))) {
                    std::cout
                    << "Failed comparison at (i,j) = ("
                    << i << ", " << j << ") : "
                    << "expected: " << m0(i,j) << "  , "
                    << "found: "    << m(i,j) << " : "
                    << "diff: " << m0(i,j) - m(i,j)
                    << std::endl;
                    pass = false;
                }
        }
        
        return pass;
    }
    
    
    
    // variables at which the evaluation is performed
    Real
    rho   ,
    u1    ,
    u2    ,
    u3    ,
    T     ,
    cp    ,
    cv    ,
    R     ,
    gamma ,
    k     ,
    tol   ;
    
    
    // Primitive solution for use
    std::auto_ptr<PrimitiveSolution> _primitive_sol;
    
    // Flight condition data structure for initialization of element
    std::auto_ptr<FlightCondition> _flight_condition;
    
    // Conservative solution for use
    std::auto_ptr<DenseRealVector> _cons_sol;
    
    // dummy infile
    std::auto_ptr<GetPot> _dummy_infile;
    
    // Element class for use
    std::auto_ptr<FluidElemBase> _fluid_elem;
    
};

#endif // __mast_fluid_elem_initialization_h__
