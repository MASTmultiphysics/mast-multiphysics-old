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
#include <fstream>

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
    omega(0.),
    k_ref(0.),
    root(0.)
    {}
    
    Real V, g, omega, k_ref;
    
    Complex root;
    
    ComplexVectorX  mode;
    
    RealVectorX modal_participation;
};



class FlutterSolverBase
{
public:
    
    FlightCondition* flight_condition;
    
    CoupledAeroStructuralModel* aero_structural_model;

    /*!
     *    file to which the result will be written
     */
    std::ofstream _output, _mode_output;

    /*!
     *    constructor for the flutter solver base object
     */
    FlutterSolverBase()
    {}
    
    virtual ~FlutterSolverBase()
    {}
    
    
    virtual std::pair<bool, const FlutterRoot*> find_next_root() = 0;
    
    
    void set_output_file(std::string& nm)
    {
        _output.close();
        _output.open(nm.c_str(), std::ofstream::out);
        std::ostringstream oss;
        oss << "modes_" << nm;
        _mode_output.open(oss.str().c_str(), std::ofstream::out);
    }
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
