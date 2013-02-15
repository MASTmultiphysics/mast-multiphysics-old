//
//  DegreeOfFreedomObject.h
//  FESystem
//
//  Created by Manav Bhatia on 4/6/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_degree_of_freedom_object_h__
#define __fesystem_degree_of_freedom_object_h__

// C++ includes
#include <vector>
#include <set>

// FESystem includes
#include "Base/FESystemTypes.h"


namespace FESystem
{
    // Forward decleration
    namespace Utility {template <typename ValType> class AutoPtrTable;}
    namespace Mesh {class ElemBase;}
    namespace Mesh {class ElemBoundaryBase;}
    namespace Geometry {class Point;}
    
    namespace DegreeOfFreedom
    {
        // Forward declerations
        class DegreeOfFreedomUnit;
        
        /*!
         *   This class defines contains information about the degree of freedoms. The mesh entities like nodes and elements all 
         *   derive from this class. 
         */
        class DegreeOfFreedomObject
        {
        public:
            /*!
             *   Constructor
             */
            DegreeOfFreedomObject();
            
            /*!
             *   Destructor
             */
            virtual ~DegreeOfFreedomObject();
            
            /*!
             *    adds the association to an element and sets the point association
             */
            void setElement(FESystem::Mesh::ElemBase& elem, FESystemBoolean if_point_associated, const FESystem::Geometry::Point* pt = NULL);

            /*!
             *    adds the association to a boundary entity, be it a vertex, edge or face. The element checks internally that this is indeed the case
             *    based on the point association. Hence, this method can only be called after the setElement method has been called.
             */
            void setBoundaryAssociation(FESystem::Mesh::ElemBoundaryBase& boundary);
            
            /*!
             *    tells the object that it is a dependent degree of freedom object and sets the coefficients and parent dof objects
             */
            void setDegreeOfFreedomConstraints(std::vector<FESystemDouble, FESystem::DegreeOfFreedom::DegreeOfFreedomObject*>& parent_dofs);
            
            /*!
             *    tells the object that the specified dof object is dependent on this object
             */
            void setDependentDegreeOfFreedomObject(FESystem::DegreeOfFreedom::DegreeOfFreedomObject& dof);

            /*!
             *    tells the object that the specified dof object is dependent on this object
             */
            void setDependentDegreeOfFreedomObject(std::set<FESystem::DegreeOfFreedom::DegreeOfFreedomObject>& dof);

            /*!
             *    Initializes the table to the required size of number of variables \p n_vars
             */
            void init(FESystemUInt n_vars);
            
            /*!
             *   returns the number of degree of freedom units in this degree of freedom object
             */
            FESystemUInt getNDegreeOfFreedomUnits() const;
            
            /*!
             *   Adds a new DegreeOfFreedomUnit for variable \p var_id to this DegreeOfFreedomObject. 
             *   It is an error to add new variables to the same system and variable ids.
             */
            FESystem::DegreeOfFreedom::DegreeOfFreedomUnit& addDegreeOfFreedomUnit(FESystemUInt var_id);
            
            /*!
             *   Returns the DegreeOfFreedomUnit for variable \p var_id in this DegreeOfFreedomObject. 
             *   An exception is thrown if the variable does not exist
             */
            FESystem::DegreeOfFreedom::DegreeOfFreedomUnit& getDegreeOfFreedomUnit(FESystemUInt var_id);

            /*!
             *   Returns the DegreeOfFreedomUnit for variable \p var_id to this DegreeOfFreedomObject. 
             *   An exception is thrown if the variable does not exist
             */
            const FESystem::DegreeOfFreedom::DegreeOfFreedomUnit& getDegreeOfFreedomUnit(FESystemUInt var_id) const;
            
        protected:
            
            /*!
             *   If the element has been initialized
             */
            FESystemBoolean if_initialized;
            
            /*!
             *   Element to which this degree of freedom belongs
             */
            FESystem::Mesh::ElemBase* elem;
            
            /*!
             *   If the object is associated with a computational coordinate location, or if it is independent of a coordinate
             */
            FESystemBoolean if_associated_with_computational_location;
            
            /*!
             *   If the object is associated with a coordinate on the element, then this defines the value
             */
            FESystem::Geometry::Point* location;
            
            /*!
             *   Set of boundaries that this degree of freedom lies on: this will include the vertices, edges and faces. The
             *   first element of the vector has boundaries of dimension 0, followed by those of dimension 1 and so on. Note that
             *   this needs to be less than the dimensionality of the element.
             */
            std::vector<std::set<FESystem::Mesh::ElemBoundaryBase*> > associated_boundaries;
            
            /*!
             *   If the object is an independent degree of freedom, or if it is constrained to other degree of freedom objects
             */
            FESystemBoolean if_independent_dof;
            
            /*!
             *   The degrees of freedom that this is dependent on. The first element of the pair defines the coefficient. The sum of all the coefficients should
             *   equal 1.0
             */
            std::vector<std::pair<FESystemDouble, FESystem::DegreeOfFreedom::DegreeOfFreedomObject*> > parent_dof_objets;
            
            /*!
             *   Dof objects that may depend on this object. It is imporant that if this object is going out of scope, it tells these
             *   objects that they are now independent
             */
            std::set<FESystem::DegreeOfFreedom::DegreeOfFreedomObject*> dependent_dof_objects;
            
            /*!
             *   Table of variables for this degree of freedom object. The table is two dimensional, where the
             *   first dimension is the system id, and the second dimension is the variable number for the system in
             *   question
             */
            FESystem::Utility::AutoPtrTable<FESystem::DegreeOfFreedom::DegreeOfFreedomUnit>* variable_table;
        };
    }
}



#endif  // __fesystem_degree_of_freedom_object_h__
