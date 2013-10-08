//
//  potential_flow_elem_base.h
//  MAST
//
//  Created by Manav Bhatia on 10/7/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST__potential_flow_elem_base__
#define __MAST__potential_flow_elem_base__


// libmesh incldues
#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"

// MAST includes
#include "FluidElems/fluid_elem_base.h"


enum PotentialFlowVars
{ RHO, PHI };


// Forward declerations
class SurfaceMotionBase;


/*!
 *    Provides the base class for potential flow equations
 */
class PotentialFlowElemBase
{
public:
    // Constructor
    PotentialFlowElemBase(GetPot& infile):
    surface_motion(NULL),
    flight_condition(NULL),
    dim(0),
    _infile(infile)
    { }
    
    virtual ~PotentialFlowElemBase()
    { }

    void init_data();

    /*!
     *   flight condition for analysis
     */
    FlightCondition* flight_condition;
    
    SurfaceMotionBase* surface_motion;
    
protected:

    /*!
     *   Input source for this class
     */
    GetPot& _infile;

    /*!
     *   flight condition for analysis
     */
    unsigned int dim;
    
    /*!
     *    map of boundary ids and boundary condition
     */
    std::multimap<unsigned int, FluidBoundaryConditionType> _boundary_condition;

};




#endif /* defined(__MAST__potential_flow_elem_base__) */
