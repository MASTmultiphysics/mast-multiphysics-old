//
//  Header.h
//  MAST
//
//  Created by Manav Bhatia on 10/25/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_bernoulli_bending_operator_h__
#define __MAST_bernoulli_bending_operator_h__

// libMesh includes
#include "libmesh/point.h"
#include "libmesh/elem.h"
#include "libmesh/dense_vector.h"

// MAST includes
#include "StructuralElems/bending_operator.h"
#include "Numerics/fem_operator_matrix.h"


namespace MAST {
    class BernoulliBendingOperator : public MAST::BendingOperator1D {
        
    public:
        BernoulliBendingOperator(StructuralElementBase& elem):
        MAST::BendingOperator1D(elem),
        _length(_elem.volume())
        { }

        /*!
         *   returns true if this bending operator supports a transverse shear component
         */
        virtual bool include_transverse_shear_energy() const {
            return false;
        }

        /*!
         *   initialze the bending strain operator for Bernoulli element, withouth
         *   the y,z-location. This is useful for use with element stiffness matrix
         *   integration where the D matrix is calculated by section integration by
         *   the ElementPropertyCard1D.
         */
        virtual void initialize_bending_strain_operator (const unsigned int qp,
                                                         FEMOperatorMatrix& Bmat);
        
        /*!
         *    initializes the bending strain operator for the specified quadrature
         * point and y,z-location.
         */
        void initialize_bending_strain_operator_for_yz(const unsigned int qp,
                                                       const libMesh::Real y,
                                                       const libMesh::Real z,
                                                       FEMOperatorMatrix& Bmat_bend);

    protected:
        
        /*!
         *   element length
         */
        libMesh::Real _length;
    };
}


inline void
MAST::BernoulliBendingOperator::initialize_bending_strain_operator (const unsigned int qp,
                                                                    FEMOperatorMatrix& Bmat) {
    this->initialize_bending_strain_operator_for_yz(qp, 1., 1., Bmat);
}



inline void
MAST::BernoulliBendingOperator::initialize_bending_strain_operator_for_yz (const unsigned int qp,
                                                                           const libMesh::Real y,
                                                                           const libMesh::Real z,
                                                                           FEMOperatorMatrix& Bmat) {
    const libMesh::Real xi = _qrule.get_points()[qp](0);
    
    // shape function values
    // N1 = (length/8.0) * (4.0/length -  6.0/length*xi + 0.0 +  2.0/length*pow(xi,3));
    // N2 = (length/8.0) * (4.0/length +  6.0/length*xi + 0.0 -  2.0/length*pow(xi,3));
    // N3 = (length/8.0) * ( 1.0 - xi - pow(xi,2) + pow(xi,3));  // needs a -1.0 factor for theta_y
    // N4 = (length/8.0) * (-1.0 - xi + pow(xi,2) + pow(xi,3));  // needs a -1.0 factor for theta_y
    
    // shape function derivative
    // N1 = (1.0/4.0) * (0.0 -  6.0/length + 0.0     +  6.0/length*pow(xi,2));
    // N2 = (1.0/4.0) * (0.0 +  6.0/length + 0.0     -  6.0/length*pow(xi,2));
    // N3 = (1.0/4.0) * (0.0 -  1.0        - 2.0*xi  +         3.0*pow(xi,2));  // needs a -1.0 factor for theta_y
    // N4 = (1.0/4.0) * (0.0 -  1.0        + 2.0*xi  +         3.0*pow(xi,2));  // needs a -1.0 factor for theta_y
    
    libMesh::DenseVector<libMesh::Real> N; N.resize(2);
    
    // second order shape function derivative
    N(0) = (0.5/_length) * (  0.0     +  12.0/_length*xi);
    N(1) = (0.5/_length) * (  0.0     -  12.0/_length*xi);
    N.scale(-y);
    Bmat.set_shape_function(0, 1, N);  // v-disp

    N(0) = (0.5/_length) * (  0.0     +  12.0/_length*xi);
    N(1) = (0.5/_length) * (  0.0     -  12.0/_length*xi);
    N.scale(-z);
    Bmat.set_shape_function(1, 2, N);  // w-disp
    
    N(0) = (0.5/_length) * (- 2.0     +           6.0*xi);
    N(1) = (0.5/_length) * (  2.0     +           6.0*xi);
    N.scale(-y);
    Bmat.set_shape_function(0, 5, N); // theta-z

    N(0) = (0.5/_length) * (- 2.0     +           6.0*xi);  // needs a -1.0 factor for theta_y
    N(1) = (0.5/_length) * (  2.0     +           6.0*xi);  // needs a -1.0 factor for theta_y
    N.scale(z);
    Bmat.set_shape_function(1, 4, N); // theta-y
}

#endif // __MAST_bernoulli_bending_operator_h__
