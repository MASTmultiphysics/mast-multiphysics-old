
/*
 *  MeshBase.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/15/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_mesh_base_h__
#define __fesystem_mesh_base_h__

// C++ includes
#include <vector>
#include <map>
#include <memory>

// FESystem includes
#include "Base/FESystemTypes.h"
#include "Mesh/ElementType.h"

namespace FESystem
{
    // Forward declerations
    namespace Geometry {class CoordinateSystemBase;}
    
    namespace Mesh
    {
        // Forward decleration
        class Node;
        class ElemBase;
        
        /*!
         *   This defines the base mesh class that stores the nodes and elements 
         */
        
        class MeshBase 
        {
        public:
            /*!
             *   Default constructor. It initializes the data structures. 
             */
            MeshBase();
            
            virtual ~MeshBase();
            
            
            /*!
             *   Reinitializes the data structures
             */
            void reinit();
            
            /*!
             *   Returns the number of nodes in this mesh
             */
            virtual FESystemUInt getNNodes() const;
            
            
            /*!
             *   Returns the number of elements in this mesh
             */
            virtual FESystemUInt getNElements() const;


            /*!
             *   Returns the number of elements of type \p t in the mesh
             */
            FESystemUInt getNElemsOfType(FESystem::Mesh::ElementType t) const;
            
            
            /*!
             *   Returns a reference to the vector of Nodes in this mesh
             */
            std::vector<FESystem::Mesh::Node*>& getNodes();

            
            /*!
             *   Returns a constant reference to the vector of Nodes in this mesh
             */
            const std::vector<FESystem::Mesh::Node*>& getNodes() const;

            
            /*!
             *   Returns a reference to the vector of elements in this mesh
             */
            std::vector<FESystem::Mesh::ElemBase*>& getElements();
            
            /*!
             *   Returns a constant reference to the vector of elements in this mesh
             */
            const std::vector<FESystem::Mesh::ElemBase*>& getElements() const;
            
            
            /*!
             *    Returns the node corresponding to the internal ID
             */
            FESystem::Mesh::Node& getNodeFromInternalID(FESystemUInt i); 

            
            /*!
             *    Returns the node corresponding to the internal ID
             */
            FESystem::Mesh::Node& getNodeFromExternalID(FESystemUInt i); 

            /*!
             *    Returns the node corresponding to the internal ID
             */
            const FESystem::Mesh::Node& getNodeFromInternalID(FESystemUInt i) const; 

            
            /*!
             *    Returns the node corresponding to the internal ID
             */
            const FESystem::Mesh::Node& getNodeFromExternalID(FESystemUInt i) const; 

            /*!
             *    Returns the element corresponding to the internal ID
             */
            FESystem::Mesh::ElemBase& getElemFromInternalID(FESystemUInt i); 
            
            
            /*!
             *    Returns the element corresponding to the internal ID
             */
            FESystem::Mesh::ElemBase& getElemFromExternalID(FESystemUInt i); 

            /*!
             *    Returns the element corresponding to the internal ID
             */
            const FESystem::Mesh::ElemBase& getElemFromInternalID(FESystemUInt i) const; 
            
            
            /*!
             *    Returns the element corresponding to the internal ID
             */
            const FESystem::Mesh::ElemBase& getElemFromExternalID(FESystemUInt i) const; 

            /*!
             *   Creates a new node in this mesh and return the pointer. The node will be defined in the given coordinate system \p cs. 
             *   The node will be associated to this mesh. The node location can be reset. 
             */ 
            FESystem::Mesh::Node& createNode(const FESystem::Geometry::CoordinateSystemBase& cs);
            

            /*!
             *   Creates \p n_nodes new node in this mesh and returns the pointers to them through the smart pointer array. The nodes are all defined in the 
             *   given coordinate system \p cs, and are associated with this mesh. The node locations can be reset. 
             */ 
            std::auto_ptr<std::vector<FESystem::Mesh::Node*> > createNodes(FESystemUInt n_nodes, const FESystem::Geometry::CoordinateSystemBase& cs);

            /*!
             *   Creates \p n_elems new elements of type \p elem_type in this mesh and returns the pointers to them through the smart pointer array. 
             *   The element node definitions can be defined using these pointers
             */ 
            std::auto_ptr<std::vector<FESystem::Mesh::ElemBase*> > createElements(FESystemUInt n_elems, FESystem::Mesh::ElementType elem_type);

            
        protected:
            
            /*!
             *   Stores if the mesh is initialized
             */
            FESystemBoolean if_initialized;
            
            /*!
             *   Vector of pointers to all nodes in this mesh
             */
            std::vector<FESystem::Mesh::Node*> nodes;
            
            /*!
             *   Vector of pointers to all elements in this mesh
             */
            std::vector<FESystem::Mesh::ElemBase*> elements;
            
            /*!
             *    Map of internal ID to node pointer
             */
            std::map<FESystemUInt, FESystem::Mesh::Node*> internal_id_to_node_map;

            /*!
             *    Map of external ID to node pointer
             */
            std::map<FESystemUInt, FESystem::Mesh::Node*> external_id_to_node_map;

            /*!
             *    Map of internal ID to elem pointer
             */
            std::map<FESystemUInt, FESystem::Mesh::ElemBase*> internal_id_to_elem_map;

            /*!
             *    Map of external ID to elem pointer
             */
            std::map<FESystemUInt, FESystem::Mesh::ElemBase*> external_id_to_elem_map;

                        
        };
        
        
    }
}

#endif // __fesystem_mesh_base_h__

