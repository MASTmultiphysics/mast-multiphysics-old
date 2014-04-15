//
//  aerodynamic_model.h
//  MAST
//
//  Created by Manav Bhatia on 7/27/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef MAST_aerodynamic_model_h
#define MAST_aerodynamic_model_h

// MAST includes
#include "Base/MAST_data_types.h"


class AerodynamicModel
{
public:
    AerodynamicModel()
    { }
    
    virtual ~AerodynamicModel()
    { }
  
    /*!
     *    updates the aerodynamic damping matrix operator in \par a. This is
     *    applicable for only the time-domain methods. Returns \par false if
     *    the model does not have a damping term.
     */
    virtual bool get_damping_matrix(RealMatrixX& a)
    { libmesh_assert(false); }
    
    /*!
     *    updates the aerodynamic stiffness matrix operator in \par a. This is
     *    applicable for only the time-domain methods. Returns \par false if
     *    the model does not have a stiffness term.
     */
    virtual bool get_stiffness_matrix(RealMatrixX& a)
    { libmesh_assert(false); }
    
};




#endif
