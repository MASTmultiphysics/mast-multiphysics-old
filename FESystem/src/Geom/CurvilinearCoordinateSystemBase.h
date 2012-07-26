
/*
 *  CoordinateCurvilinearBasisSystemBase.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_coordinate_curvilinear_system_base_h__
#define __fesystem_coordinate_curvilinear_system_base_h__


// FESystem includes
#include "Base/FESystemTypes.h"
#include "Geom/CoordinateSystemBase.h"

namespace FESystem
{
	namespace Geometry
	{
        // Forward declerations
        class Point;
        
		template <FESystemUInt Dim, FESystemUInt MappedDim>
		class CoordinateCurvilinearBasisSystemBase: 
		public FESystem::Geometry::CoordinateSystemBase
		{
		public:
			CoordinateCurvilinearBasisSystemBase();
			
			virtual ~CoordinateCurvilinearBasisSystemBase();
			
		protected:
			
			virtual void mapCoordinates(const FESystem::Geometry::Point& p_in, 
										FESystem::Geometry::Point& p_mapped) = 0;

			virtual void inverseMapCoordinates(const FESystem::Geometry::Point& p_map, 
											   FESystem::Geometry::Point& p_inv_mapped) = 0;

		};
	}
}


#endif // __fesystem_coordinate_curvilinear_system_base_h__
