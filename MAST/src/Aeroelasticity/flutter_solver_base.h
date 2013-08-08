//
//  flutter_solver_base.h
//  MAST
//
//  Created by Manav Bhatia on 7/25/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST__flutter_solver_base__
#define __MAST__flutter_solver_base__

// C++ includes
#include <vector>
#include <memory>

// MAST includes
#include "Base/MAST_data_types.h"
#include "Flight/flight_condition.h"
#include "Aeroelasticity/coupled_aero_structural_model.h"


class FlutterRoot
{
public:
    FlutterRoot():
    V(0.),
    g(0.),
    k_ref(0.),
    root(0.)
    {}
    
    Real V, g, k_ref;
    
    Complex root;
    
    ComplexVectorX  mode;
};



class FlutterSolverBase
{
public:
    
    FlutterSolverBase()
    {}
    
    virtual ~FlutterSolverBase()
    {}
    
    
    virtual bool find_next_root() = 0;
    
    
    unsigned int n_roots_found() const
    {
        return flutter_roots.size();
    }
    
    
    FlightCondition* flight_condition;
    
    std::vector<FlutterRoot> flutter_roots;
    
    CoupledAeroStructuralModel* aero_structural_model;
};




class FrequencyDomainSolverBase: public FlutterSolverBase
{
public:
    FrequencyDomainSolverBase():
    FlutterSolverBase(),
    current_k_ref(0.)
    {}
    
    virtual ~FrequencyDomainSolverBase()
    {}
    
    /*!
     *   Current value of reduced frequency for which the aerodynamic matrices
     *   are to be evaluated.
     */
    Real current_k_ref;
};



#endif /* defined(__MAST__flutter_solver_base__) */
