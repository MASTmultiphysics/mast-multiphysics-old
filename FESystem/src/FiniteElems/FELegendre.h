//
//  FELegendre.h
//  FESystem
//
//  Created by Manav Bhatia on 1/25/13.
//
//

#ifndef __FESystem__FELegendre__
#define __FESystem__FELegendre__

// FESystem includes
#include "FiniteElems/FiniteElementBase.h"


namespace FESystem
{
    // Forward declerations
    namespace Functions {template <typename ValType> class LagrangeFunction;}
    namespace Functions {template <typename ValType> class LegendreFunction;}
    namespace Functions {template <typename ValType> class LagrangeFunctionMapping;}
    namespace Functions {template <typename ValType> class LegendreFunctionMapping;}
    
    
    namespace FiniteElement
    {
        /*!
         *   This derives from the base class FiniteElementBase and implements the Lagrange shape functions.
         */
        class FELegendre: public FESystem::FiniteElement::FiniteElementBase
        {
        public:
            /*!
             *   Constructor
             */
            FELegendre();
            
            virtual ~FELegendre();
            
            /*!
             *   Returns the type of finite element
             */
            virtual FESystem::FiniteElement::FiniteElementType getFiniteElementType() const;
            
            /*!
             *   Method to clear the shape functions in the matrices that have been set up by reinit. This method is to be called
             *   by the method reinit.
             */
            virtual void clear();
            
            /*!
             *   This will initialize the data structures so that the shape functions, and the required derivatives can be calculated
             */
            virtual void reinit(const FESystem::Mesh::ElemBase& element, const FESystemUInt order);
            
            /*!
             *   Returns the number of shape functions, which for a Legendre element is based on the order of interpolation, and a combination of 
             *   polynomial terms from each interpolation order. Thus, the final number = (order+1)^dim
             */
            virtual FESystemUInt getNShapeFunctions() const;
            
        protected:
            
            /*!
             *   Method to initialize the shape functions in the matrices that have been set up by reinit. This method is to be called
             *   by the method reinit.
             */
            virtual void initializeMaps();
            
            /*!
             *   highest order of the Legendre function to be used
             */
            FESystemUInt legendre_order;
            
            /*!
             *   Lagrange functions for each dimension, used to interpolate both geometry
             */
            std::vector<FESystem::Functions::LagrangeFunction<FESystemDouble>*> lagrange_functions;

            /*!
             *   Legendre function object pointer. Only one function is used since it does not need to be initialized like the LagrangeFunction and 
             *   the same object can be used for each dimension
             */
            FESystem::Functions::LegendreFunction<FESystemDouble>* legendre_function;

            /*!
             *   Lagrange mapping based on the Lagrange Functions to perform interpolations and calculate shape funcitons values for geometry interpolation.
             */
            FESystem::Functions::LagrangeFunctionMapping<FESystemDouble>* lagrange_map;

            /*!
             *   Mapping based on the Legendre Functions to perform interpolations and calculate shape funcitons values for variable interpolation.
             */
            FESystem::Functions::LegendreFunctionMapping<FESystemDouble>* legendre_map;

        };
    }
}

#endif /* defined(__FESystem__FELegendre__) */
