
/*
 *  Domain.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_domain_h__
#define __fesystem_domain_h__

// FESystem includes
#include "Base/FESystemTypes.h"
#include "Base/ParameterBase.h"
#include "Geom/DomainBase.h"


namespace FESystem
{
	namespace Geometry
	{		
		typedef FESystem::Base::ParameterBound<FESystemDouble> CoordinateBound;
		
		template <FESystemUInt Dim, FESystemUInt MappedDim>
		class Domain: public FESystem::Geometry::DomainBase
		{
		public: 
			Domain();
			
			virtual ~Domain();
			
			const FESystem::Geometry::CoordinateBound& getCoordinateBound(const FESystemUInt num);
			
		protected:
			
			
			
		};
	}
}


#endif // __fesystem_domain_h__
