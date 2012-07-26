
/*
 *  DiscretizationBase.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/11/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_discretization_base_h__
#define __fesystem_discretization_base_h__


// FESystem includes
#include "Base/CountedObject.h"


namespace FESystem
{
    // Forward declerations
    namespace Mesh {class MeshBase;}
    namespace Base {class DegreeOfFreedomMap;}
    namespace Numerics {template <typename ValType> class VectorBase;}
    namespace Numerics {template <typename ValType> class MatrixBase;}
    
	namespace Discretization
	{
        /*!
         *    This provides the base class for evaluation of quantities needed for solvers. 
         */
		class DiscretizationBase//: 
		//public FESystem::Base::CountedObject<FESystem::Discretization::DiscretizationBase> 
		{
		public:
            /*!
             *    Constructor
             */
			DiscretizationBase(FESystem::Mesh::MeshBase& m, FESystem::Base::DegreeOfFreedomMap& dofs):
            mesh(m),
            dof_map(dofs)
            {}
			
			virtual ~DiscretizationBase(){}

            /*!
             *    Returns the mesh for this discretization
             */
            FESystem::Mesh::MeshBase& getMesh(){return mesh;}
			            
            /*!
             *    Returns the degree of freedom map for this discretization
             */
            FESystem::Base::DegreeOfFreedomMap& getDofMap(){return dof_map;}

            /*!
             *   Evaluates the Jacobian with respect to state given in \p vec. The Jacobian is to be returned in \p mat
             */
            virtual void evaluateJacobian(FESystemUInt iter_num, const FESystem::Numerics::VectorBase<FESystemDouble>& vec, FESystem::Numerics::MatrixBase<FESystemDouble>& mat)=0;

            /*!
             *   Evaluates the Jacobian for a transient problem with respect to state given in \p vec. The Jacobian is to be returned in \p mat
             */
            virtual void evaluateTransientSolverData(FESystemUInt transient_iter_num, FESystemUInt sub_iter_num, FESystemDouble time,
                                                     const FESystem::Numerics::VectorBase<FESystemDouble>& x_vec,
                                                     FESystem::Numerics::VectorBase<FESystemDouble>& x_dot_vec,
                                                     FESystem::Numerics::MatrixBase<FESystemDouble>& x_dot_jac_mat)=0;

		protected:
			
            
            /*!
             *   Mesh object for this discretization
             */
            FESystem::Mesh::MeshBase& mesh;

            /*!
             *   Degree of freedom map
             */
            FESystem::Base::DegreeOfFreedomMap& dof_map;

		};
		
	}
	
}


#endif // __fesystem_discretization_base_h__

