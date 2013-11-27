//
//  bending_operator.cpp
//  MAST
//
//  Created by Manav Bhatia on 11/27/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

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


