//
//  timoshenko_bending_operator.h
//  MAST
//
//  Created by Manav Bhatia on 11/27/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_timoshenko_bending_operator_h__
#define __MAST_timoshenko_bending_operator_h__

// libMesh includes
#include "libmesh/point.h"
#include "libmesh/elem.h"
#include "libmesh/dense_vector.h"

// MAST includes
#include "StructuralElems/bending_operator.h"
#include "Numerics/fem_operator_matrix.h"


namespace MAST {
    class TimoshenkoBendingOperator : public MAST::BendingOperator {
        
    public:
        TimoshenkoBendingOperator(StructuralElementBase& elem):
        MAST::BendingOperator(elem)
        { }
        
        /*!
         *   returns true if this bending operator supports a transverse shear component
         */
        virtual bool include_transverse_shear_energy() const {
            return false;
        }

        /*!
         *   initialze the bending strain operator for DKT element
         */
        virtual void initialize_bending_strain_operator (const unsigned int qp,
                                                         FEMOperatorMatrix& Bmat) ;
        
    protected:

    };
}


inline void
MAST::TimoshenkoBendingOperator::initialize_bending_strain_operator (const unsigned int qp,
                                                                     FEMOperatorMatrix& Bmat)
{
    libmesh_error(); // to be implemented
}


#endif // __MAST_timoshenko_bending_operator_h__
