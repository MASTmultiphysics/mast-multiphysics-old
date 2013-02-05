//
//  QuadratureBase.h
//  FESystem
//
//  Created by Manav Bhatia on 4/9/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_quadrature_base_h__
#define __fesystem_quadrature_base_h__

// C++ includes
#include <vector>

// FESystem includes
#include "Base/FESystemTypes.h"


namespace FESystem
{
    // Forward declerations
    namespace Geometry {class Point;}
    namespace Mesh {class ElemBase;}
    namespace FiniteElement {class FiniteElementBase;}
    
    namespace Quadrature
    {        
        /*!
         *    This defines the base class of quadrature rules. The basic interface is specified here, and the 
         *    rule specific routines will be implemented by the derived classes. The integration domains are assumed 
         *    to be between -1 and 1. 
         */
        class QuadratureBase
        {
        public:
            /*!
             *   Default constructor
             */ 
            QuadratureBase();
            
            virtual ~QuadratureBase();
            
            /*!
             *    clears the data structures for reinitialization
             */
            virtual void clear();
            
            /*!
             *   Returns the dimension of this quadrature rule
             */ 
            FESystemUInt getDimention() const;

            /*!
             *   Returns the order of this quadrature rule
             */ 
            FESystemUInt getOrder() const;

            /*!
             *   Initializes the quadrature rule for the specified dimension and order of the function that is exactly integrated
             */ 
            virtual void init(FESystemUInt dim, FESystemUInt order)=0;
            
            /*!
             *   Returns a constant reference to the vector of quadrature points
             */ 
            const std::vector<FESystem::Geometry::Point*>& getQuadraturePoints() const;

            /*!
             *   Returns a constant reference to the vector of quadrature point weights
             */ 
            const std::vector<FESystemDouble>& getQuadraturePointWeights() const;
            
        protected:
            
            /*!
             *   initializes the origin for the specified dimension
             */
            void initializeOrigin(FESystemUInt dim);
            
            /*!
             *   Flag about whether the object is initialized
             */
            FESystemBoolean if_initialized;
            
            /*!
             *   Dimension of the space for which this quadrature rule is initialized
             */ 
            FESystemUInt dimension;
            
            /*!
             *   Order of the function for which this quadrature rule is required to be accurate
             */ 
            FESystemUInt order;
            
            /*!
             *   origin for the coordinates of the quadature point
             */
            FESystem::Geometry::Point* origin;

            /*!
             *   vector of quadrature points 
             */ 
            std::vector<FESystem::Geometry::Point*> quadrature_points;

            /*!
             *   vector of quadrature point weights
             */ 
            std::vector<FESystemDouble> quadrature_point_weights;
        };
        
        /*!
         *    Initializes the qrule for the given geometric element and the finite element
         */
        void initializeQRuleForElem(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, FESystem::Quadrature::QuadratureBase& qrule);

        /*!
         *    Initializes the qrule for use on boundary of the given geometric element and the finite element
         */
        void initializeQRuleForElemBoundary(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, FESystem::Quadrature::QuadratureBase& qrule);

    }
}



#endif //__fesystem_quadrature_base_h__

