
/*
 *  RectangularCoordinateSystem.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_coordinate_rectangular_basis_system_h__
#define __fesystem_coordinate_rectangular_basis_system_h__

// FESystem includes
#include "Base/FESystemTypes.h"
#include "Geom/CoordinateSystemBase.h"


namespace FESystem
{
	namespace Geometry
	{        
		/*!
         *   This class defines a \p dim dimensional Rectangular Cartesian basis system. It derives from the base class CoordinateBasisStystemBase 
         *   and specializes the transformations for the rectangular Cartesian coordinate system. 
         */
        
		class RectangularCoordinateSystem: public FESystem::Geometry::CoordinateSystemBase
		{
		public:            
            /*! 
             *   Constructor initializes the dimension of the space and the mapped space. 
             */
			RectangularCoordinateSystem(const FESystem::Geometry::Point& origin, const FESystem::Numerics::MatrixBase<FESystemDouble>& basis);			

			virtual ~RectangularCoordinateSystem();
			
		protected:
			
            /*! 
             *   Maps the point \p p_in to the parent basis, and returns the value in \p p_mapped
             */
            virtual void mapCoordinates(const FESystem::Geometry::Point& p_in, 
										FESystem::Geometry::Point& p_mapped) const;
            
            /*! 
             *   Calculates the location of the point defined in the parent coordinate system 
             */
			virtual void inverseMapCoordinates(const FESystem::Geometry::Point& p_mapped, 
											   FESystem::Geometry::Point& p_inv_mapped) const;

		};
	}
}


#endif // __fesystem_coordinate_rectangular_basis_system_h__

