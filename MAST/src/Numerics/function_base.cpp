//
//  function_base.cpp
//  MAST
//
//  Created by Manav Bhatia on 11/08/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

// MAST includes
#include "Numerics/function_base.h"
#include "Numerics/sensitivity_parameters.h"



bool
MAST::FunctionBase::depends_on(const MAST::SensitivityParameters &p) const {
    
    // only first order sensitivities are calculated at this point
    libmesh_assert_equal_to(p.total_order(), 1);
    
    const MAST::SensitivityParameters::ParameterMap& p_map = p.get_map();
    MAST::SensitivityParameters::ParameterMap::const_iterator it, end;
    it = p_map.begin(); end = p_map.end();
    
    const MAST::FunctionBase& f = *(it->first);
    return this->depends_on(f);
}



