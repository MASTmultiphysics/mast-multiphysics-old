//
//  FiniteElementBase.h
//  FESystem
//
//  Created by Manav Bhatia on 3/22/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_finite_element_base_h__
#define __fesystem_finite_element_base_h__

// C++ includes
#include <vector>

// FESystem includes
#include "Base/FESystemTypes.h"


namespace FESystem
{
    // Forward declerations
    namespace Utility {template <typename ValType> class Table;}

    namespace Geometry { class Point;}

    namespace Mesh { class ElemBase;}

    namespace Numerics { template <typename ValType> class MatrixBase;}
    namespace Numerics { template <typename ValType> class VectorBase;}

    namespace Functions {template <typename ValType> class DiscreteCompositeMappingFunctionBase;}

    
    namespace FiniteElement 
    {
        /*!
         *    Type of finite element
         */
        enum FiniteElementType
        {
            FE_LAGRANGE,
            FE_LEGENDRE
        };
        
        
        /*!
         *   This class defines the base class for finite elements. The basic purpose of this class is to 
         *   associate itself with a geometric element of 1-, 2- or 3- dimensions, and provide the shape 
         *   functions, shape function derivatives, and transformation jacobians associated at specified 
         *   points in the domain of the element. Different element types that derive from this base class 
         *   are required to provide the specializations of the shape functions. The coordinates are all 
         *   interpolated using the Lagrange functions for each node, and this is used to calculate the 
         *   transformation Jacobians. 
         */
        class FiniteElementBase
        {
        public:
            /*!
             *   Constructor. 
             */
            FiniteElementBase();
            
            virtual ~FiniteElementBase();
            
            /*!
             *   Returns the type of finite element
             */
            virtual FESystem::FiniteElement::FiniteElementType getFiniteElementType() const = 0;
            
            
            /*!
             *   Clears all the data structures of the object, before initializing it
             */
            virtual void clear();

            /*!
             *   This will initialize the data structures so that the shape functions, and the required derivatives can be calculated
             */
            virtual void reinit(const FESystem::Mesh::ElemBase& element);
            
            /*!
             *    Returns the element for which this element is initialized
             */
            const FESystem::Mesh::ElemBase& getGeometricElement() const;
            
            /*!
             *    Returns the number of shape functions for the element
             */
            virtual FESystemUInt getNShapeFunctions() const=0;
             
            /*!
             *   Returns the shape function based on the finite element initialization for the local coordinate specified in vin 
             */
            void getShapeFunction(const FESystem::Numerics::VectorBase<FESystemDouble>& vin, FESystem::Numerics::VectorBase<FESystemDouble>& vout) const;

            /*!
             *   Returns the shape function for specified boundary, based on the finite element initialization for the local coordinate specified in vin 
             */
            void getShapeFunctionForBoundary(const FESystemUInt b_id, const FESystem::Numerics::VectorBase<FESystemDouble>& vin, FESystem::Numerics::VectorBase<FESystemDouble>& vout) const;

            /*!
             *   Returns the shape function drivative wrt local coordinate of order specified in coordinate_derivative_order, based on the finite element 
             *   initialization for the local coordinate specified in vin 
             */
            void getShapeFunctionDerivativeForLocalCoordinates(const std::vector<FESystemUInt>& coordinate_derivative_order,
                                                               const FESystem::Numerics::VectorBase<FESystemDouble>& vin,
                                                               FESystem::Numerics::VectorBase<FESystemDouble>& vout) const;

            /*!
             *   Returns the shape function drivative wrt local coordinate of order specified in coordinate_derivative_order, based on the finite element 
             *   initialization for the local coordinate specified in vin. This is calculated only on the boundary \p b_id. The computational coordinate
             *   is of dimension one less than the dimension of the element, and the other coordinate is assumed from the boundary id specified. The 
             *   specified coordinates are assigned to the remaining coordinates in increasing order. 
             */
            void getShapeFunctionDerivativeForLocalCoordinatesForBoundary(const FESystemUInt b_id,
                                                                          const std::vector<FESystemUInt>& coordinate_derivative_order,
                                                                          const FESystem::Numerics::VectorBase<FESystemDouble>& vin,
                                                                          FESystem::Numerics::VectorBase<FESystemDouble>& vout) const;

            
            /*!
             *   Returns the shape function drivative wrt physical coordinate of order specified in coordinate_derivative_order, based on the finite element 
             *   initialization for the local coordinate specified in vin. It is important to note that the physical coordinate for which the derivative is 
             *   calculated is the same as the coordinate system of the element. 
             */
            void getShapeFunctionDerivativeForPhysicalCoordinates(const std::vector<FESystemUInt>& coordinate_derivative_order,
                                                                  const FESystem::Numerics::VectorBase<FESystemDouble>& vin,
                                                                  FESystem::Numerics::VectorBase<FESystemDouble>& vout) const;

            /*!
             *   Returns the shape function drivative wrt physical coordinate of order specified in coordinate_derivative_order, based on the finite element 
             *   initialization for the local coordinate specified in vin. It is important to note that the physical coordinate for which the derivative is 
             *   calculated is the same as the coordinate system of the element. 
             */
            void getShapeFunctionDerivativeForPhysicalCoordinatesForBoundary(const FESystemUInt b_id,
                                                                             const std::vector<FESystemUInt>& coordinate_derivative_order,
                                                                             const FESystem::Numerics::VectorBase<FESystemDouble>& vin,
                                                                             FESystem::Numerics::VectorBase<FESystemDouble>& vout) const;

            /*!
             *   Returns the Jacobiain in the matrix \p jac. The Jacobian is defined for the mapping from the computational to the physical coordinate
             *   \f$ J_{ij} [ \partial X_i / \partial \xi_j ] \f$, where X_i is the physical coordinate and \xi_j is the computational coordinate
             */
            void getJacobianMatrix(const FESystem::Numerics::VectorBase<FESystemDouble>& vin, FESystem::Numerics::MatrixBase<FESystemDouble>& jac) const;

            /*!
             *   Returns determinant of the Jacobian matrix at the specified point. For a 1-D element this is the differential length, for 2-D this is the 
             *   differential area and the differential volume in the case of a 3-D element. 
             */
            FESystemDouble getJacobianValue(const FESystem::Numerics::VectorBase<FESystemDouble>& vin) const;

            /*!
             *   Returns the boundary determinant of the Jacobian matrix at the specified point. For a 2-D problem this is the differential length of the boundary, 
             *   for 3-D element this is the differential area of the specified boundary. For a 1-D element, a unit value is returned. 
             */
            FESystemDouble getJacobianValueForBoundary(const FESystemUInt b_id, const FESystem::Numerics::VectorBase<FESystemDouble>& vin) const;

        protected:
            
            
            /*!
             *   Virtual method to initialize the shape functions in the matrices that have been set up by reinit. This method is to be called 
             *   by the method reinit. 
             */
            virtual void initializeMaps() = 0;
                        
            /*!
             *   Sotres the initialization state of the element
             */
            FESystemBoolean if_initialized;
            
            /*!
             *   Pointer to the geometric element
             */
            const FESystem::Mesh::ElemBase*  geom_elem;
            
            /*!
             *   Pointer to the composite map that provides the mapping from computational to physical coordinate and solution variable
             *   interpolations
             */
            FESystem::Functions::DiscreteCompositeMappingFunctionBase<FESystemDouble>* composite_map;
        };
    }
}



#endif  // __fesystem_finite_element_base_h__

