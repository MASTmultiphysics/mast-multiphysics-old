//
//  ideal_gas.h
//  
//
//  Created by Manav Bhatia on 7/26/13.
//
//

#ifndef __ideal_gas_h__
#define __ideal_gas_h__

// MAST includes
#include "Base/MAST_data_types.h"

class GasProperty
{
public:
    GasProperty():
    rho(0.),
    T(0.),
    pressure(0.),
    gamma(0.),
    cp(0.),
    cv(0.),
    R(0.),
    a(0.),
    Pr(0.),
    k_thermal(0.),
    mu(0.),
    lambda(0.)
    {}
    
    
    void zero()
    {
        rho       = 0.;
        T         = 0.;
        pressure  = 0.;
        gamma     = 0.;
        cp        = 0.;
        cv        = 0.;
        R         = 0.;
        a         = 0.;
        Pr        = 0.;
        k_thermal = 0.;
        mu        = 0.;
        lambda    = 0.;
    }
    
    
    /*!
     *   Property values for ideal gas
     */
    Real rho, T, pressure, gamma, cp, cv, R, a;
    
    /*!
     *   Properties for viscous analysis
     */
    Real Pr, k_thermal, mu, lambda;

    /*!
     *   initializes the data
     */
    void init();
    
};



inline void
GasProperty::init()
{
    // the following data should have been set
    libmesh_assert(rho > 0.);
    libmesh_assert(T  >  0.);
    libmesh_assert(cp  > 0.);
    libmesh_assert(cv  > 0.);
    
    R        = cp-cv;
    gamma    = cp/cv;
    pressure = rho*R*T;
    a        = sqrt(gamma*R*T);
}


#endif
