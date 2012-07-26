/*
 *  CoordinateSystemBase.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_coordinate_basis_system_base_h__
#define __fesystem_coordinate_basis_system_base_h__

// C++ includes
#include <memory>

// FESystem includes
#include "Base/FESystemTypes.h"
#include "Base/FESystemExceptions.h"

namespace FESystem
{

	// Forward declerations
    namespace Numerics {template <typename ValType> class VectorBase;}
    namespace Numerics {template <typename ValType> class MatrixBase;}
    namespace Functions {template <typename ValType> class FunctionMappingBase;}
    
    namespace Geometry
	{        
        // Forward declerations
        class Point;
        
        /*! 
         *   A coordinate basis is the combination of unit vectors. The coordinate basis may have a dimension \p Dim,
         *   while its unit vectors may exist in a space with dimensions \p MappedDim. Note that \f$ MappedDim >= Dim \f$.
         *   
         *   For example, a basis may have a single dimension defined by the unit vector \f$ (1,1,1)/\sqrt{3} \f$.
         *   The purpose of this class is to be able to define a basis within an upper level coordinate system, so that 
         *   a hierarchy of coordinate systems can be defined. 
         */
        class CoordinateSystemBase
		{
		public:
            /*! 
             *   Constructor initializes the dimension of the space. If the coordinate is a global basis, then there is no need to provide a parent. 
             */
			CoordinateSystemBase(const FESystem::Geometry::Point& origin);
			
			virtual ~CoordinateSystemBase();
			
            /*!
             *   Returns the parent coordinate system to which this Basis maps
             */ 
            const FESystem::Geometry::Point& getOrigin() const;
            
            /*! 
             *   Returns the dimension of this basis system
             */
            FESystemUInt getDimension() const;
                        
            /*!
             *    Maps the given vector \p vec into the current coordinate system. The method makes sure that the given vector
             *    and this coordinate system share some common parent basis. If this is true, then the transformation is performed, 
             *    otherwise, it is an error to call this method for two unrelated coordinate systems. 
             */
            void mapPointInSelf(const FESystem::Geometry::Point& p, FESystem::Numerics::VectorBase<FESystemDouble>& vec);

            /*!
             *   returns the function mapping object used for transformation of coordinate system
             */
            const FESystem::Functions::FunctionMappingBase<FESystemDouble>& getFunctionMappingObject() const;
            
        protected:
            
            /*!
             *   Pointer to the origin of the coordinate system. The parent coordinate system is inferred from the 
             *   origin, which is associated with its own basis. 
             */
            const FESystem::Geometry::Point&  origin;
                        
            /*! 
             *   defines the mapping of the point from this coordinate system to a parent basis
             */
            FESystem::Functions::FunctionMappingBase<FESystemDouble>* coordinate_mapping;
		};

        /*!
         *   this exception checks if the dimension of the coordinate system is less than or equal to its parent's dimension
         */
        DeclareException2(CoordinateBasisDimensionGreaterThanParentBasisDimension, 
						  FESystemUInt,
						  FESystemUInt,
						  << "Coordinate Basis Dimension Greater Than Parent Basis Dimension. This Basis: " << Arg1 
						  << ". Parent Basis: " << Arg2 );

        
        /*!
         *   this exception checks if the current basis is a global basis. 
         */
        DeclareException0(AttempToModifyGlobalBasis, 
						  << "Global Basis Cannot be Modified");

        /*!
         *   this exception checks if the current basis is a global basis. 
         */
        DeclareException2(CoordinateDimensionDoesNotMatchBasis,
                          FESystemUInt,
                          FESystemUInt,
						  << "Coordinate Dimension Of " << Arg1 << " Does Not Match Basis Dimension Of " << Arg2);

	}
}



#endif // __fesystem_coordinate_basis_system_base_h__

