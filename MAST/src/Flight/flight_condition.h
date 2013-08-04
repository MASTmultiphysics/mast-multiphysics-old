//
//  flight_condition.h
//  MAST
//
//  Created by Manav Bhatia on 7/26/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef MAST_flight_condition_h
#define MAST_flight_condition_h

// MAST includes
#include "Base/MAST_data_types.h"
#include "PropertyCards/gas_property.h"

// Eigen includes
#include "Eigen/Core"


class FlightCondition
{
public:
    FlightCondition():
    ref_chord(0.),
    altitude(0.)
    {}
    
    virtual ~FlightCondition()
    {}

    /*!
     *   defines the vehicle longitudinal axis in the analysis domain. Also
     *   the body x-axis
     */
    RealVector3 body_roll_axis;
    
    /*!
     *   defines the vehicle pitch axis in the analysis domain. Also
     *   the body y-axis
     */
    RealVector3 body_pitch_axis;
    
    /*!
     *    defines the vehicle yaw axis in the analysis domain. Also
     *   the body z-axis
     */
    RealVector3 body_yaw_axis;
    
    /*!
     *   defines the three angles that define the orientation of the velocity
     *   vector with respect to the body-axis. 
     *   Three components are: roll-angle, pitch angle and yaw angle,
     *   respectively.
     */
    RealVector3 body_euler_angles;
    
    /*!
     *   defines the three angular rates angular rates of the vehicle in body
     *   axis.
     *   Three components correspond to: roll-angle, pitch angle and yaw angle,
     *   respectively.
     */
    RealVector3 body_angular_rates;

    /*!
     *  Velocity magnitude, whose direction is evaluated from the Euler angles.
     */
    Real velocity_magnitude;
    
    /*!
     *   Ambient air properties
     */
    GasProperty gas_property;
    
    /*!
     *   reference chord
     */
    Real ref_chord;
    
    /*!
     *   Flight altitude. This need not be set, but can be used to evaluate the 
     *   ambient properties if needed.
     */
    Real altitude;
    
    /*!
     *   returns the flight dynamic pressure
     */
    Real q0() const
    {
        return 0.5 * gas_property.rho * pow(velocity_magnitude, 2);
    }
    
    Real p0()
    {
        return gas_property.pressure;
    }

    Real rho() const
    {
        return gas_property.rho;
    }

    
    
protected:
    
    /*!
     *    defines the lift and drag vectors that are calculated based on the 
     *    body axis and Euler angles specified in the input
     */
    RealVector3 _flight_vector;


public:

    Real rho_u1() const
    {
        return gas_property.rho * velocity_magnitude * _flight_vector(0);
    }
    
    Real rho_u2() const
    {
        return gas_property.rho * velocity_magnitude * _flight_vector(1);
    }
    
    Real rho_u3() const
    {
        return gas_property.rho * velocity_magnitude * _flight_vector(2);
    }
    
    Real rho_e() const
    {
        return gas_property.rho * gas_property.cv * gas_property.T + q0();
    }
    
    RealVector3 lift_normal, drag_normal;

};


//    gamma = cp/cv;
//    R = cp-cv;
//    a_inf = sqrt(gamma*R*temp_inf);
//
//    u1_inf = mach_inf*a_inf*cos(aoa*pi/180.0);
//    u2_inf = mach_inf*a_inf*sin(aoa*pi/180.0);
//    u3_inf = 0.0;
//    q0_inf = 0.5*rho_inf*(u1_inf*u1_inf+u2_inf*u2_inf+u3_inf*u3_inf);
//    p_inf = R*rho_inf*temp_inf;
//
//    Real k = 0.0;
//    vars_inf(0) = rho_inf;
//    vars_inf(1) = rho_inf*u1_inf; k += u1_inf*u1_inf;
//    if (dim > 1)
//    {
//        vars_inf(2) = rho_inf*u2_inf;
//        k += u2_inf*u2_inf;
//    }
//    if (dim > 2)
//    {
//        vars_inf(3) = rho_inf*u3_inf;
//        k += u3_inf*u3_inf;
//    }
//    vars_inf(dim+2-1) = rho_inf*(cv*temp_inf+0.5*k);


#endif
