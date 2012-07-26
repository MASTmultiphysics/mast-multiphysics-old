
/*
 *  CylindricalCoordinateSystem.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_coordinate_cylindrical_basis_system_h__
#define __fesystem_coordinate_cylindrical_basis_system_h__

// FESystem includes
#include "Geom/CoordinateSystemBase.h"

namespace FESystem
{
	namespace Geometry
	{
		class CylindricalCoordinateSystem : public FESystem::Geometry::CoordinateSystemBase
		{
		public:			
			CylindricalCoordinateSystem();
			
			virtual ~CylindricalCoordinateSystem();
			
		protected:
			
			
		};
	}
}


#endif // __fesystem_coordinate_cylindrical_basis_system_h__
