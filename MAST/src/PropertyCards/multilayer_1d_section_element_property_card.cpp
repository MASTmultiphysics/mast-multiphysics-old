//
//  multilayer_1d_section_element_property_card.cpp
//  MAST
//
//  Created by Manav Bhatia on 2/18/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//


// MAST includes
#include "PropertyCards/multilayer_1d_section_element_property_card.h"
#include "StructuralElems/bending_structural_elem.h"


std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
MAST::Multilayer1DSectionElementPropertyCard::get_property(MAST::ElemenetPropertyMatrixType t,
                                                           const MAST::StructuralElementBase& e) const {
    
    // prepare vector of matrix functions from each layer
    std::vector<MAST::FieldFunction<DenseRealMatrix >*> layer_mats(_layers.size());
    for (unsigned int i=0; i<_layers.size(); i++)
        layer_mats[i] = _layers[i]->get_property(t, e).release();
    
    // now create the integrated object
    std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > rval
    (new MAST::Multilayer1DSectionElementPropertyCard::SectionIntegratedMatrix
     (layer_mats));
    
    return rval;
}


