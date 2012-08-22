/*
 *  ElemBase.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */


#ifndef __fesystem_elem_base_h__
#define __fesystem_elem_base_h__

// C++ includes
#include <vector>
#include <set>
#include <map>

// FESystem includes
#include "Base/FESystemTypes.h"
#include "Base/FESystemExceptions.h"
#include "Mesh/ElementType.h"
#include "Base/DegreeOfFreedomObject.h"
#include "Base/IdObject.h"


namespace FESystem
{
    // Forward declerations
    namespace Geometry {class Point;}
    namespace Geometry {class CoordinateSystemBase;}
    namespace Utility {template <typename ValType> class Table;}
    namespace Numerics {template <typename ValType> class VectorBase;}
    namespace Numerics {template <typename ValType> class MatrixBase;}
    namespace FiniteElement {class FiniteElementBase;}
    namespace Quadrature {class QuadratureBase;}

	namespace Mesh
	{
        // Forward declerations
        class MeshBase;
        class Node;

        /*!
         *   This provides the base class for all geometric element types. An element is defined by a shape, 1-D, 2-D or 3-D, and 
         *   the nodes that it is connected by. For each element, the sequence of nodes is important as it defined the face normal directions. 
         *   The requirement is specifically defined in each element. 
         */
		class ElemBase: public FESystem::Base::DegreeOfFreedomObject, public FESystem::Base::IdObject
		{
		public:
            /**
             *   The constructor takes the \t MeshBase class to which this element is associated. 
             */ 
			ElemBase(FESystemUInt n_nodes, FESystem::Mesh::ElementType type);
			
			virtual ~ElemBase();
            
            /*!
             *  returns the element type 
             */
            FESystem::Mesh::ElementType getElementType() const;
            
            /*!
             *   Returns the number of nodes that belong to this element
             */
            FESystemUInt getNNodes() const; 
                        
            /*!
             *   Sets the \p i th node for this element to \p n
             */
            void setNode(FESystemUInt i, FESystem::Mesh::Node& n); 

            /*!
             *   Gets the \p i th node for this element
             */
            FESystem::Mesh::Node& getNode(FESystemUInt i); 

            /*!
             *   Gets the \p i th node for this element
             */
            const FESystem::Mesh::Node& getNode(FESystemUInt i) const; 

            /*!
             *   Returns the number of dimensions of this element
             */
            virtual FESystemUInt getDimension() const = 0; 

            /*!
             *   Returns the geometry order for this element. 
             */
            virtual FESystemUInt getGeometryOrder() const = 0;
            
            /*!
             *   Returns the length of the element in the dimensionality of the element: length for 1-D element, area for 2-D and volume for 3-D. 
             *   This requires a finite element and quadrature base initialized for this element.
             */
            virtual FESystemDouble getElementSize(const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule) const;

            /*!
             *   returns the location of the \p i_node^th node in the local physical coordinate
             */
            void getNodeLocationInLocalPhysicalCoordinate(FESystemUInt i_node, FESystem::Numerics::VectorBase<FESystemDouble>& vec) const;
        

            /*!
             *   Returns the local coordinate for \p n^th node along the \p d^th dimension. Note that this is to support the tensorial 
             *   products for creation of shape functions.
             */
            virtual FESystemDouble getLocalComputationalCoordinateForNodeAlongDim(FESystemUInt d, FESystemUInt n) const;
            
            /*!
             *   Returns the local coordinate for \p id^th node in the element.
             */
            virtual void getLocalComputationalCoordinateForLocalNodeID(FESystemUInt id, FESystem::Geometry::Point& p) const = 0;

            /*!
             *   Returns a constant reference to the table that identifies the ID of the node along the local dimension
             *   This is useful for Lagrange function mapping for this element. The table is 2-D, with number of rows 
             *   as the number of nodes, and the number of columns as the number of local dimensions. 
             */   
            virtual const FESystem::Utility::Table<FESystemUInt>& getPointDimensionIDTable() const = 0;
            
            /*! 
             *   This defines the number of nodes along the \p n^th dimension of the element. This
             *   is useful for tasks like tensorial products of unidimensional shape functions to create 
             *   shape functions in multiple dimensions. 
             */
            virtual FESystemUInt getNNodesAlongDimension(FESystemUInt n) const = 0;

            /*!
             *   This method supports the tensor prouct for creation of shape functions. It returns the number of the node with 
             *   internal \p id along the \p d dimension. 
             *   
             */
            virtual FESystemUInt getLocalNodeIDALongDim(FESystemUInt n_id, FESystemUInt d) const = 0;
            

            /*!
             *   Returns a constant referenc to the local physical coordinate system for the element. This is a 
             *   rectangular coordinate system with the origin at node zero. The dimension of this coordinate system is 
             *   same as the dimensions of this element. 
             */
            const FESystem::Geometry::CoordinateSystemBase& getLocalPhysicalCoordinateSystem() const;
            
            /*!
             *   Returns true if the element is a degerate element from Quad or Hex element. This is used for 
             *   elements like Triangle, Tetrahedron, Prisms, Pyramids, etc. 
             */
            virtual bool ifDegerateElement() const = 0;
            
            /*!
             *   Returns the parent nondegerate element from which this element is derived. 
             */
            virtual const FESystem::Mesh::ElemBase& getParentNondegenerateElem() const=0;
            
            /*!
             *   Returns the matrix that maps the shape functions from the parent nondegenrate to this 
             *   element
             */
            virtual const FESystem::Numerics::MatrixBase<FESystemDouble>& getParentToDegenerateElemMappingMatrix() const=0;

            /*!
             *    returns the internal numbering of the node in this element
             */
            FESystemUInt getInternalIDForNode(const FESystem::Mesh::Node& node) const;
            
            /*!
             *   Returns the id of the boundary that contains the nodes specified in the set
             */  
            virtual FESystemUInt getIDOfBoundaryWithNodes(const std::set<FESystem::Mesh::Node*>& b_nodes) const;

            /*!
             *   Returns the id of the boundary that contains the node ids specified in the set.  
             */  
            virtual FESystemUInt getIDOfBoundaryWithNodes(const std::set<FESystemUInt>& b_nodes) const;

            /*!
             *   Returns the map of boundary id and nodes defining the boundary
             */
            virtual const std::map<FESystemUInt, std::vector<FESystemUInt> >& getBoundaryIDAndBoundaryNodeMap() const = 0;
            
            /*!
             *   Each bounadry in FE is defined by one of the computational coordinates being constant. This method returns the 
             *   number of that computational coordinate and its value for the boundary, in the provided \p coord_id and \p coord_val variables. 
             *   This is useful for doing computations on the boundaries, for example, application of boundary conditions
             */
            virtual void getConstantCoordinateIDAndValueForBoundary(const FESystemUInt b_id, FESystemUInt& coord_id, FESystemDouble& coord_val) const = 0; 

            /*!
             *   Calculates the normal of the specified boundary. This method is to be used for only the linear elements that have straight edges. 
             *   For curved elements, the other method should be used, since it requires an interpolation method. 
             */ 
            void calculateBoundaryNormal(const FESystemUInt b_id, FESystem::Numerics::VectorBase<FESystemDouble>& n_vec) const;
            
            /*!
             *   Calaultes the surface normal for higher order elements that require a function mapping for calculation the element geometry from specified points
             */
            // void calculateBoundaryNormal(const FESystemUInt b_id, const FESystem::Functions::FunctionMappingBase<FESystemDouble>& geom_mapping_func, 
            //                             const FESystem::Numerics::VectorBase<FESystemDouble>& b_pt, FESystem::Numerics::VectorBase<FESystemDouble>& n_vec ) const;
            
            

        protected:

            /*!
             *   Initializes the local coordinate system
             */
            virtual void initializeLocalPhysicalCoordinateSystem()=0;
                        
            /*!
             *   Initialize parent nondegenerate element. Is defined for each inherited element
             */
            virtual void initializeParentNondegenerateElement()=0;
                        
            /*!
             *   Type of this element
             */
            FESystem::Mesh::ElementType element_type;
            
 			/*!
             *  This vector container stores the nodes that comprise this element. This stores the physical coordinates of the nodes
             */
            std::vector<FESystem::Mesh::Node*> physical_nodes;
            
            /*!
             *   Coordinate system in the element plane: this is to calculate the physical coordinates and not computational coordinates
             */
            FESystem::Geometry::CoordinateSystemBase* local_coordinate_system;                        

            /*!
             *  parent nondegenrate element
             */
            FESystem::Mesh::ElemBase* parent_nondegenerate_elem;
        };
        
        /*!
         *  Exception for resetting element node, which has already been defined
         */
        DeclareException1(CannotResetElementNode, 
                          FESystemUInt,
						  << "Cannot Reset Element Node Number: " << Arg1);
        
        /*!
         *  Exception for accessing a node that has not been set
         */
        DeclareException1(ElementNodeNotDefined, 
                          FESystemUInt,
						  << "Element Node Number: " << Arg1 << " : Not Defined.");

        /*!
         *   Function to create elements of a given type \p elem_type, associated with the given \p mesh.
         */
        std::auto_ptr<FESystem::Mesh::ElemBase> 
        ElementCreate(FESystem::Mesh::MeshBase& mesh, FESystem::Mesh::ElementType elem_type);

	}
}


#endif // __fesystem_elem_base_h__
