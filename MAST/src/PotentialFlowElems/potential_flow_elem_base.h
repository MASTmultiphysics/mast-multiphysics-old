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

#ifndef __MAST__potential_flow_elem_base__
#define __MAST__potential_flow_elem_base__


// libmesh incldues
#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"

// MAST includes
#include "FluidElems/fluid_elem_base.h"


namespace MAST {
    
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
        
        MAST::SurfaceMotionBase* surface_motion;
        
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
    
}


#endif /* defined(__MAST__potential_flow_elem_base__) */
