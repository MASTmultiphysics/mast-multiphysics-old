//
//  StructuralDiscretization.h
//  FESystem
//
//  Created by Manav Bhatia on 4/27/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_structural_discretization_h__
#define __fesystem_structural_discretization_h__

// FESystem includes
#include "Disciplines/DiscretizationBase.h"


namespace FESystem
{
    namespace Discretization
    {
        class StructuralDiscretization: public FESystem::Discretization::DiscretizationBase
		{
		public:
            /*!
             *    Constructor
             */
			StructuralDiscretization(FESystem::Mesh::MeshBase& m, FESystem::Base::DegreeOfFreedomMap& dofs);
			
			virtual ~StructuralDiscretization();
            
            /*!
             *   Evaluates the Jacobian with respect to state given in \p vec. The Jacobian is to be returned in \p mat
             */
            virtual void evaluateJacobian(FESystemUInt iter_num, const FESystem::Numerics::VectorBase<FESystemDouble>& vec, FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            /*!
             *   Evaluates the Jacobian for a transient problem with respect to state given in \p vec. The Jacobian is to be returned in \p mat
             */
            virtual void evaluateTransientSolverData(FESystemUInt transient_iter_num, FESystemUInt sub_iter_num, FESystemDouble time,
                                                     const FESystem::Numerics::VectorBase<FESystemDouble>& x_vec,
                                                     FESystem::Numerics::VectorBase<FESystemDouble>& x_dot_vec,
                                                     FESystem::Numerics::MatrixBase<FESystemDouble>& x_dot_jac_mat);
            
		protected:
			
            
		};
    }
}



#endif  //__fesystem_structural_discretization_h__
