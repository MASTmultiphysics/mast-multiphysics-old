//
//  bending_operator.h
//  MAST
//
//  Created by Manav Bhatia on 11/27/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_bending_operator_h__
#define __MAST_bending_operator_h__

// C++ includes
#include <memory>

// libMesh includes
#include "libmesh/dense_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/elem.h"
#include "libmesh/quadrature.h"

// MAST includes
#include "StructuralElems/structural_elem_base.h"

// Forward declerations
class FEMOperatorMatrix;


namespace MAST {
    // Forward declerations
    class SensitivityParameters;
    
    enum BendingOperatorType {
        BERNOULLI,        // beam
        TIMOSHENKO,       // beam
        DKT,              // plate
        MINDLIN,          // plate
        DEFAULT_BENDING,
        NO_BENDING
    };
    

    class BendingOperator {
    public:
        
        BendingOperator(StructuralElementBase& elem):
        _structural_elem(elem),
        _elem(_structural_elem.local_elem()),
        _qrule(_structural_elem.quadrature_rule())
        { }
        
        virtual ~BendingOperator()
        { }
        
        /*!
         *   returns true if this bending operator supports a transverse shear component
         */
        virtual bool include_transverse_shear_energy() const = 0;

        /*!
         *   initialze the bending strain operator for the specified quadrature point
         */
        virtual void initialize_bending_strain_operator (const unsigned int qp,
                                                         FEMOperatorMatrix& Bmat) = 0;
        
        /*!
         *   calculate the transverse shear component for the element
         */
        virtual void calculate_transverse_shear_force(bool request_jacobian,
                                                      DenseVector<Real>& local_f,
                                                      DenseMatrix<Real>& local_jac,
                                                      const MAST::FieldFunctionBase* sens_params )
        { libmesh_error(); }
        
        
    protected:
        
        /*!
         *   structural element associated with this
         */
        MAST::StructuralElementBase& _structural_elem;
        
        /*!
         *    element for which bending operator is created
         */
        const Elem& _elem;
        
        /*!
         *   quadrature rule to be used. This should already be of the correct order
         */
        QBase& _qrule;

    };

    
    
    /*!
     *   Bending strain operator for 1D element
     */
    class BendingOperator1D: public MAST::BendingOperator {
    public:
        BendingOperator1D(StructuralElementBase& elem):
        MAST::BendingOperator(elem)
        { }

        /*!
         *   initialze the bending strain operator for the specified quadrature point
         */
        virtual void initialize_bending_strain_operator_for_yz (const unsigned int qp,
                                                                const Real y,
                                                                const Real z,
                                                                FEMOperatorMatrix& Bmat) = 0;

    };

    
    /*!
     *   Bending strain operator for 1D element
     */
    class BendingOperator2D: public MAST::BendingOperator {
    public:
        BendingOperator2D(StructuralElementBase& elem):
        MAST::BendingOperator(elem)
        { }
        
        /*!
         *   initialze the bending strain operator for the specified quadrature point
         */
        virtual void initialize_bending_strain_operator_for_z (const unsigned int qp,
                                                                const Real z,
                                                                FEMOperatorMatrix& Bmat) = 0;
        
    };

    
    /*!
     *   builds a bending operator and returns it in a smart-pointer
     */
    std::auto_ptr<MAST::BendingOperator>
    build_bending_operator(MAST::BendingOperatorType type,
                           StructuralElementBase& elem);
}

#endif
