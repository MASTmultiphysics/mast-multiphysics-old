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
    
};


#endif
