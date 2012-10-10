
/*
 *  Node.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */


#ifndef __fesystem_node_h__
#define __fesystem_mesh_h__

#include <set>

// FESystem includes
#include "Geom/Point.h"
#include "Base/DegreeOfFreedomObject.h"
#include "Base/IdObject.h"


namespace FESystem
{

    // Forward declerations
    namespace Geometry {class CoordinateSystemBase;}
    
	namespace Mesh
	{
        // Forward declerations
        class MeshBase;
        class ElemBase;
        
        /*!
         *   A node is a basic entity of a computational mesh, which has a physical location in a domain, and may or may not have 
         *   computable degrees of freedom associated with it. The physical location is obtained though the class \t Point, whilt the 
         *   degree of freedom is obtained by the class \t DegreeOfFreedom. 
         */ 
		class Node: public FESystem::Geometry::Point, public FESystem::Base::DegreeOfFreedomObject, public FESystem::Base::IdObject
		{
		public:
            /*!
             *  The constructor requires the \t MeshBase object to which this node belongs, and the underlying coordinate system in 
             *  which the location of this node is defined
             */
			Node(const FESystem::Geometry::CoordinateSystemBase& cs);
			
			virtual ~Node();
            
            /*!
             *   Tells the node that it is connected to the element \el. This is needed in preparation of degree of freedom connectivity
             */
            void addElementConnection(FESystem::Mesh::ElemBase& el);

            /*!
             *   Returns true if the node is connected to the element \p el
             */
            FESystemBoolean ifConnectedToElement(FESystem::Mesh::ElemBase& el);
            
            /*!
             *   Returns a constant reference to the set of pointers to which this element is connected
             */
            const std::set<FESystem::Mesh::ElemBase*>& getElementConnectivitySet() const;

            /*!
             *   Returns the number of elements that share this node
             */
            FESystemUInt getNConnectedElements() const;
            
            /*!
             *   Returns the number of elements of specified dimension connected to this node
             */
            FESystemUInt getNConnectedElementsOfDim(FESystemUInt dim) const;
            
            
		protected:
			
            /*!
             *   Set of elements that this node is connected to
             */
            std::set<FESystem::Mesh::ElemBase*> element_connections;
            
		};
	}
}

#endif // __fesystem_mesh_h__
