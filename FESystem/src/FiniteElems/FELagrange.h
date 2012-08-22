//
//  FELagrange.h
//  FESystem
//
//  Created by Manav Bhatia on 3/23/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_fe_lagrange_h__
#define __fesystem_fe_lagrange_h__


// FESystem includes
#include "FiniteElems/FiniteElementBase.h"


namespace FESystem
{
    // Forward declerations
    namespace Functions {template <typename ValType> class LagrangeFunction;}
    namespace Functions {template <typename ValType> class LagrangeFunctionMapping;}
    
    
    namespace FiniteElement
    {
        /*!
         *   This derives from the base class FiniteElementBase and implements the Lagrange shape functions.
         */
        class FELagrange: public FESystem::FiniteElement::FiniteElementBase
        {
        public:            
            /*!
             *   Constructor
             */
            FELagrange();
            
            ~FELagrange();
            
            /*!
             *   Method to clear the shape functions in the matrices that have been set up by reinit. This method is to be called 
             *   by the method reinit. 
             */
            virtual void clear();

            /*!
             *   Returns the number of shape functions, which for a Lagrange element is the number of nodes
             */
            virtual FESystemUInt getNShapeFunctions() const;

        protected:
            
            /*!
             *   Method to initialize the shape functions in the matrices that have been set up by reinit. This method is to be called 
             *   by the method reinit. 
             */
            virtual void initializeMaps();
        
            /*!
             *   Lagrange functions for each dimension, used to interpolate both geometry and variables
             */
            std::vector<FESystem::Functions::LagrangeFunction<FESystemDouble>*> lagrange_functions;
            
            /*!
             *   Lagrange mapping based on the Lagrange Functions to perform interpolations and calculate shape funcitons values.
             */
            FESystem::Functions::LagrangeFunctionMapping<FESystemDouble>* lagrange_map;
            
        };
    }
}


#endif  // __fesystem_fe_lagrange_h__
