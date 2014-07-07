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
#include "StructuralElems/bending_operator.h"
#include "StructuralElems/bernoulli_bending_operator.h"
#include "StructuralElems/timoshenko_bending_operator.h"
#include "StructuralElems/dkt_bending_operator.h"
#include "StructuralElems/mindlin_bending_operator.h"


std::auto_ptr<MAST::BendingOperator>
MAST::build_bending_operator(MAST::BendingOperatorType type,
                             StructuralElementBase& elem) {
    std::auto_ptr<MAST::BendingOperator> rval;
    
    switch (type) {
        case MAST::BERNOULLI:
            rval.reset(new MAST::BernoulliBendingOperator(elem));
            break;

        case MAST::TIMOSHENKO:
            rval.reset(new MAST::TimoshenkoBendingOperator(elem));
            break;

        case MAST::DKT:
            rval.reset(new MAST::DKTBendingOperator(elem));
            break;

        case MAST::MINDLIN:
            rval.reset(new MAST::MindlinBendingOperator(elem));
            break;

        case MAST::NO_BENDING:
            // nothing to be done
            break;
            
        default:
            libmesh_error(); // should not get here
            break;
    }
    
    return rval;
}


