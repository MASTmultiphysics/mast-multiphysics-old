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
     *   Flight Mach number
     */
    Real mach;
    
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

    /*!
     *    initializes the data structures
     */
    void init();
    
    
    Real rho_u1() const
    {
        return gas_property.rho * velocity_magnitude * drag_normal(0);
    }
    
    Real rho_u2() const
    {
        return gas_property.rho * velocity_magnitude * drag_normal(1);
    }
    
    Real rho_u3() const
    {
        return gas_property.rho * velocity_magnitude * drag_normal(2);
    }
    
    Real rho_e() const
    {
        return gas_property.rho * gas_property.cv * gas_property.T + q0();
    }
    
    /*!
     *    defines the lift and drag vectors that are calculated based on the
     *    body axis and Euler angles specified in the input
     */
    RealVector3 lift_normal, drag_normal;
};



inline
void
FlightCondition::init()
{
    gas_property.init();
    
    // make sure that the data has been initialized
    libmesh_assert(body_roll_axis.norm()  == 1.);
    libmesh_assert(body_pitch_axis.norm() == 1.);
    libmesh_assert(body_yaw_axis.norm()   == 1.);
    libmesh_assert(mach > 0.);
    
    velocity_magnitude = mach * gas_property.a;
    
    // prepare the transformation matrix from body to inertial axis
    RealMatrix3 tmat; tmat.setZero();
    tmat.col(0) = body_roll_axis;
    tmat.col(1) = body_pitch_axis;
    tmat.col(2) = body_yaw_axis;
    
    // prepare the rotation matrices in the body coordinate system.
    // We will rotate the body vectors, and then transform the
    // resulting vector to the inertial frame.
    RealMatrix3 roll, pitch, yaw;
    roll.setZero(); pitch.setZero(); yaw.setZero();
    
    roll  << 1., 0., 0.,
             0., cos(body_euler_angles(0)), -sin(body_euler_angles(0)),
             0., sin(body_euler_angles(0)),  cos(body_euler_angles(0));
    
    pitch << cos(body_euler_angles(1)), 0., sin(body_euler_angles(1)),
             0., 1., 0.,
             -sin(body_euler_angles(1)), 0., cos(body_euler_angles(1));
    
    yaw   << cos(body_euler_angles(2)), -sin(body_euler_angles(2)), 0.,
             sin(body_euler_angles(2)),  cos(body_euler_angles(2)), 0.,
             0., 0., 1.;
    
    // drag normal is the same as flight direction, and is obtained by applying
    // the rotation matrices to the roll-axis
    RealVector3 tmp, tmpunit;
    tmpunit << 1., 0., 0.; // body roll axis
    tmp = roll * pitch * yaw * tmpunit;
    drag_normal = tmat * tmp;
    
    // lift normal is obtained by applying the rotation matrices to yaw-axis
    tmpunit << 0., 0., 1.; // body yaw axis
    tmp = roll * pitch * yaw * tmpunit;
    lift_normal = tmat * tmp;
}


#endif
