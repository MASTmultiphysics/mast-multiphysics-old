//
//  strain_operator.h
//  MAST
//
//  Created by Manav Bhatia on 10/16/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_strain_operator_h__
#define __MAST_strain_operator_h__


// MAST includes
#include "Numerics/fem_operator_matrix.h"


namespace MAST {
    
    class StrainOperator: public FEMOperatorMatrix {
        
    public:
        
        unsigned int n_qpoints() const;
        
        const std::vector<Real>& get_JxW() const;

        void initialize_for_qp(unsigned int qp);
        
    protected:
        
    };
}


#endif // __MAST_strain_operator_h__
