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


// MAST includes
#include "StructuralElems/bending_structural_elem.h"
#include "PropertyCards/element_property_card_1D.h"
#include "PropertyCards/element_property_card_2D.h"
#include "Numerics/fem_operator_matrix.h"
#include "BoundaryConditions/boundary_condition.h"
#include "StructuralElems/structural_element_1D.h"
#include "StructuralElems/structural_element_2D.h"


MAST::BendingStructuralElem::BendingStructuralElem(libMesh::System& sys,
                                                   const libMesh::Elem& elem,
                                                   const MAST::ElementPropertyCardBase& p):
MAST::StructuralElementBase(sys, elem, p)
{
    _local_elem.reset(MAST::build_local_elem(*this).release());
    
    _init_fe_and_qrule(_local_elem->local_elem());

    // initialize the bending operator
    MAST::BendingOperatorType bending_model =
    _property.bending_model(elem, _fe->get_fe_type());

    _bending_operator.reset(MAST::build_bending_operator(bending_model, *this).release());
}




std::auto_ptr<MAST::LocalElemBase>
MAST::build_local_elem(MAST::StructuralElementBase &elem) {
    
    std::auto_ptr<MAST::LocalElemBase> rval;
    switch (elem.elem().dim()) {
        case 1: {
            const MAST::ElementPropertyCard1D& p =
            dynamic_cast<const MAST::ElementPropertyCard1D&>(elem.elem_property());
            rval.reset(new MAST::Local1DElem(elem.elem(), p.y_vector()));
        }
            break;
            
        case 2:
            rval.reset(new MAST::Local2DElem(elem.elem()));
            break;
            
        default:
            libmesh_error();
    }
    
    return rval;
}

