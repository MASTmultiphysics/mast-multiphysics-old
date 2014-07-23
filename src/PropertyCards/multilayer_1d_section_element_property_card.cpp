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


