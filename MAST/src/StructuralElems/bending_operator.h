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
        _elem(_structural_elem.get_elem_for_quadrature()),
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
                                                      DenseRealVector& local_f,
                                                      DenseRealMatrix& local_jac,
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
        const libMesh::Elem& _elem;
        
        /*!
         *   quadrature rule to be used. This should already be of the correct order
         */
        libMesh::QBase& _qrule;

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
