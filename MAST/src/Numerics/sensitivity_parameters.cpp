//
//  sensitivity_parameters.cpp
//  MAST
//
//  Created by Manav Bhatia on 11/8/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

// MAST includes
#include "Numerics/sensitivity_parameters.h"
#include "Numerics/function_base.h"


unsigned int
MAST::SensitivityParameters::parameter_order(const MAST::FunctionBase* f) const {
    SensitivityParameters::ParameterMap::const_iterator it, end;
    it = _parameters.find(f); end = _parameters.end();
    libmesh_assert(it != end);
    
    return it->second;
}



bool
MAST::SensitivityParameters::shape_sensitivity() const {
    SensitivityParameters::ParameterMap::const_iterator it, end;
    it = _parameters.begin(); end = _parameters.end();
    for ( ; it != end; it++)
        if (it->first->has_attribute(MAST::SHAPE_FUNCTION))
            return true;
    return false;
}

