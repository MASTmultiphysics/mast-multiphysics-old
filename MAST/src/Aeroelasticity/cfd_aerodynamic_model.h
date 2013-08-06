//
//  cfd_aerodynamic_model.h
//  MAST
//
//  Created by Manav Bhatia on 7/28/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef MAST_cfd_aerodynamic_model_h
#define MAST_cfd_aerodynamic_model_h

// C++ includes
#include <memory>

// MAST includes
#include "Aeroelasticity/aerodynamic_model.h"
#include "FluidElems/frequency_domain_linearized_fluid_system.h"
#include "FluidElems/flexible_surface_motion.h"

#ifdef LIBMESH_USE_COMPLEX_NUMBERS


class CFDAerodynamicModel: public AerodynamicModel
{
public:
    CFDAerodynamicModel(FrequencyDomainLinearizedFluidSystem& system):
    AerodynamicModel(),
    fluid_system(system)
    { }
    
    virtual ~CFDAerodynamicModel()
    { }
    
    /*!
     *    updates the aerodynamic damping matrix operator in \par a. This is
     *    applicable for only the time-domain methods. Returns \par false if
     *    the model does not have a damping term.
     */
    virtual bool get_damping_matrix(ComplexMatrixX& a)
    { libmesh_assert(false); }
    
    /*!
     *    updates the aerodynamic stiffness matrix operator in \par a. This is
     *    applicable for only the time-domain methods. Returns \par false if
     *    the model does not have a stiffness term.
     */
    virtual bool get_stiffness_matrix(ComplexMatrixX& a)
    { libmesh_assert(false); }
    

    /*!
     *   The small-disturbance fluid system object that provides the basis
     *   of calculations for this model.
     */
    FrequencyDomainLinearizedFluidSystem& fluid_system;
    
};



#endif // LIBMESH_USE_COMPLEX_NUMBERS

#endif
