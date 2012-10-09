
/*
 *  Point.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */


#ifndef __fesystem_point_h__
#define __fesystem_point_h__

// FESystem includes
#include "Base/FESystemTypes.h"
#include "Numerics/LocalVector.h"


namespace FESystem
{
	namespace Geometry
	{
        // Forward declerations
        class CoordinateSystemBase;
        
        /*!
         *   This class defines a physical point in a given coordinate system. The coordinate system is defined by the 
         *   class \t CoordinateSystemBase, which includes a basis and origin point. 
         */
		class Point: public FESystem::Numerics::LocalVector<FESystemDouble>
		{
		public:
            /*!
             *   This mehtod can be used if the point uses the global coordinate basis. The dimension of the basis 
             *   needs to be specified and the point will assume the global basis. 
             */
			Point(FESystemUInt dim);

            /*!
             *   This mehtod can be used if the point uses the global coordinate basis. The dimension of the basis 
             *   needs to be specified and the point will assume the global basis. 
             */
			Point(const FESystem::Geometry::CoordinateSystemBase& cs);

			virtual ~Point();
            
            /*!
             *   Returns a constant reference to the coordinate system in which this point is defined
             */
            const FESystem::Geometry::CoordinateSystemBase& getCoordinateSystem() const;
                    
		protected:
            
            /*!
             *   Stores if this is the global origin, in which case it uses a global coordinate system as the
             *   coordinate system for this point
             */
            FESystemBoolean if_global_origin;
            
            /*!
             *    Coordinate system in which the coordinates of this point are defined
             */
            const FESystem::Geometry::CoordinateSystemBase* coordinate_system;
			
		};
	}
}


#endif // __fesystem_point_h__

