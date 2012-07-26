
/*
 *  PolarCoordinateSystem.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_coordinate_polar_basis_system_h__
#define __fesystem_coordinate_polar_basis_system_h__

// FESystem includes
#include "Geom/CoordinateSystemBase.h"

namespace FESystem
{
	namespace Geometry
	{
		class PolarCoordinateSystem : public FESystem::Geometry::CoordinateSystemBase
		{
		public:			
			PolarCoordinateSystem();
			
			virtual ~PolarCoordinateSystem();
			
		protected:
			
			
		};
	}
}


#endif // __fesystem_coordinate_polar_basis_system_h__
