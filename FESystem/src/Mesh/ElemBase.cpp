
/*
 *  ElemBase.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

// FESystem includes
#include "Mesh/ElemBase.h"
#include "Mesh/Edge2.h"
#include "Mesh/Edge3.h"
#include "Mesh/Quad4.h"
#include "Mesh/Quad8.h"
#include "Mesh/Quad9.h"
#include "Mesh/Tri3.h"
#include "Mesh/Tri6.h"
#include "Mesh/Tri7.h"
#include "Mesh/Node.h"
#include "Geom/CoordinateSystemBase.h"
#include "Functions/FunctionMappingBase.h"
#include "Numerics/DenseMatrix.h"
#include "FiniteElems/FiniteElementBase.h"
#include "Quadrature/QuadratureBase.h"

FESystem::Mesh::ElemBase::ElemBase(FESystemUInt n_nodes, FESystem::Mesh::ElementType type, FESystemBoolean local_cs_same_as_global):
FESystem::Base::DegreeOfFreedomObject(),
element_type(type),
if_local_physical_cs_same_as_global(local_cs_same_as_global),
physical_nodes(n_nodes),
local_coordinate_system(NULL),
parent_nondegenerate_elem(NULL)
{
    // initialize all nodes to NULL
    for (FESystemUInt i=0; i<n_nodes; i++)
        this->physical_nodes[i] = NULL;
}



FESystem::Mesh::ElemBase::~ElemBase()
{
    if ((!this->if_local_physical_cs_same_as_global) && (this->local_coordinate_system != NULL))
        delete this->local_coordinate_system;
}
    

FESystem::Mesh::ElementType 
FESystem::Mesh::ElemBase::getElementType() const
{
    return this->element_type;
}



FESystemUInt 
FESystem::Mesh::ElemBase::getNNodes() const
{
    return this->physical_nodes.size();
}
    


void
FESystem::Mesh::ElemBase::updateAfterMeshDeformation()
{
    // if all the nodes for this element have been set, initialize the local coordinate system
    for (FESystemUInt i_node=0; i_node<this->getNNodes(); i_node++)
        FESystemAssert0(this->physical_nodes[i_node] != NULL, FESystem::Exception::NULLQuantity);
    
    this->clearLocalPhysicalCoordinateSystem();
    this->clearParentNondegenerateElement();
    
    this->initializeLocalPhysicalCoordinateSystem();
    if (this->ifDegerateElement()) this->initializeParentNondegenerateElement();
}


void
FESystem::Mesh::ElemBase::setNode(FESystemUInt i, FESystem::Mesh::Node& n)
{
    FESystemAssert2(i < this->physical_nodes.size(), FESystem::Exception::IndexOutOfBound, i, this->physical_nodes.size()-1);
    FESystemAssert1(this->physical_nodes[i] == NULL, FESystem::Mesh::CannotResetElementNode, i);
    
    // set the pointer to this node
    this->physical_nodes[i] = &n;
    
    // tell the node that it is connected to this element
    n.addElementConnection(*this);
    
    // if all the nodes for this element have been set, initialize the local coordinate system
    FESystemUInt counter = 0;
    for (FESystemUInt i_node=0; i_node<this->getNNodes(); i_node++)
        if (this->physical_nodes[i_node] != NULL)
            counter++;
    if (counter == this->getNNodes())
    {
        if (!this->if_local_physical_cs_same_as_global)
            this->initializeLocalPhysicalCoordinateSystem();
        else
            this->local_coordinate_system = &this->getNode(0).getCoordinateSystem();
        if (this->ifDegerateElement()) this->initializeParentNondegenerateElement();
    }
}
    

FESystem::Mesh::Node& 
FESystem::Mesh::ElemBase::getNode(FESystemUInt i)
{
    FESystemAssert2(i < this->physical_nodes.size(), FESystem::Exception::IndexOutOfBound, i, this->physical_nodes.size()-1);
    FESystemAssert1(this->physical_nodes[i] != NULL, FESystem::Mesh::ElementNodeNotDefined, i);

    return *(this->physical_nodes[i]);
}



const FESystem::Mesh::Node& 
FESystem::Mesh::ElemBase::getNode(FESystemUInt i) const
{
    FESystemAssert2(i < this->physical_nodes.size(), FESystem::Exception::IndexOutOfBound, i, this->physical_nodes.size()-1);
    FESystemAssert1(this->physical_nodes[i] != NULL, FESystem::Mesh::ElementNodeNotDefined, i);
    
    return *(this->physical_nodes[i]);
}


void 
FESystem::Mesh::ElemBase::getNodeLocationInLocalPhysicalCoordinate(FESystemUInt i_node, FESystem::Numerics::VectorBase<FESystemDouble>& vec) const
{
    // make sure that the dimensions of the vector are accurate
    FESystemAssert2(vec.getSize() == this->getNode(i_node).getSize(), FESystem::Exception::DimensionsDoNotMatch, vec.getSize(), this->getNode(i_node).getSize());
    
    if (!this->if_local_physical_cs_same_as_global)
        this->local_coordinate_system->mapPointInSelf(this->getNode(i_node), vec);
    else
        vec.copyVector(this->getNode(i_node));
}




FESystemDouble
FESystem::Mesh::ElemBase::getElementSize(const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule) const
{
    // make sure that the dimensionality of the element and quadrature rule is appropriate
    FESystemAssert0(&(fe.getGeometricElement()) == this, FESystem::Exception::InvalidValue);
    FESystemAssert0(q_rule.getDimention() == this->getDimension(), FESystem::Exception::InvalidValue);
    
    FESystemDouble val=0.0;
        
    const std::vector<FESystem::Geometry::Point*>& q_pts = q_rule.getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = q_rule.getQuadraturePointWeights();

    for (FESystemUInt i=0; i<q_pts.size(); i++)
        val += fe.getJacobianValue(*(q_pts[i])) * q_weight[i];
    
    return val;
}



const FESystem::Geometry::CoordinateSystemBase&
FESystem::Mesh::ElemBase::getLocalPhysicalCoordinateSystem() const
{
    FESystemAssert0(this->local_coordinate_system != NULL, FESystem::Exception::InvalidState);
    return *(this->local_coordinate_system);
}



FESystemDouble
FESystem::Mesh::ElemBase::getLocalComputationalCoordinateForNodeAlongDim(FESystemUInt d, FESystemUInt n) const
{
    FESystemDouble val = 0;
    
    switch (d) {
        case 0:  // xi  for 1D, 2D and 3D elements
        case 1:  // eta for 2D and 3D elements
        case 2:  // phi for 3D elements
            switch (n) {
                case 0:  // left edge node 
                    val = -1.0;
                    break;
                    
                case 1:  // right edge node 
                    val = +1.0;
                    break;
                    
                case 2:  // middle node 
                    val =  0.0;
                    break;

                case 3:  // quarter point node
                    val =  -0.5;
                    break;

                case 4:  // three quarter point node
                    val =  0.5;
                    break;

                default:
                    FESystemAssert0(false, FESystem::Exception::InvalidValue);
                    break;
            }
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
    
    return val;
}



FESystemUInt 
FESystem::Mesh::ElemBase::getIDOfBoundaryWithNodes(const std::set<FESystemUInt>& b_nodes) const
{
    // get the map of boundary id and node ids for each boundary
    const std::map<FESystemUInt, std::vector<FESystemUInt> >& b_n_map = this->getBoundaryIDAndBoundaryNodeMap();
    
    FESystemUInt n_nodes_found = 0;

    // iterate over each boundary and try to find a match
    std::map<FESystemUInt, std::vector<FESystemUInt> >::const_iterator it=b_n_map.begin(), end=b_n_map.end();
    std::vector<FESystemUInt>::const_iterator v_it, v_end;
    
    for ( ; it != end; it++)
    {
        v_it = it->second.begin();
        v_end = it->second.end();
        n_nodes_found = 0;
        for ( ; v_it != v_end; v_it++)
            if (b_nodes.count(*v_it))
                n_nodes_found++;
        
        if (n_nodes_found == it->second.size())
            return it->first;
    }
    
    // if the code gets here, then the boundary could not be identified. This is a logic error
    FESystemAssert0(false, FESystem::Exception::InvalidValue);
}



FESystemUInt 
FESystem::Mesh::ElemBase::getInternalIDForNode(const FESystem::Mesh::Node& node) const
{
    for (FESystemUInt i=0; i<this->physical_nodes.size(); i++)
        if (this->physical_nodes[i] == &node)
            return i;
    
    // if the code gets here, then this element does not contain this node, and it is a logic error
    FESystemAssert0(false, FESystem::Exception::InvalidValue);
}


FESystemUInt 
FESystem::Mesh::ElemBase::getIDOfBoundaryWithNodes(const std::set<FESystem::Mesh::Node*>& b_nodes) const
{
    std::set<FESystem::Mesh::Node*>::const_iterator n_it=b_nodes.begin(), n_end=b_nodes.end();
    std::set<FESystemUInt> node_ids;
    
    for ( ; n_it != n_end; n_it++)
        node_ids.insert(this->getInternalIDForNode(**n_it));
    
    return this->ElemBase::getIDOfBoundaryWithNodes(node_ids);
}



void 
FESystem::Mesh::ElemBase::calculateBoundaryNormal(const FESystemUInt b_id, FESystem::Numerics::VectorBase<FESystemDouble>& n_vec) const
{
    // only to be used for first order elements
    //FESystemAssert0(this->getGeometryOrder() == 1, FESystem::Exception::InvalidFunctionCall);
    
    FESystem::Numerics::LocalVector<FESystemDouble> v1, v2, v3;
    FESystem::Numerics::DenseMatrix<FESystemDouble> jac;
    v1.resize(3);
    v2.resize(3);
    v3.resize(3);
    jac.resize(3, 3);
    
    switch (this->getDimension())
    {
        case 1:
        {
            n_vec.copyVector(this->getNode(1));
            n_vec.add(-1.0, this->getNode(0));
            n_vec.scaleToUnitLength();
            if (b_id == 0)
                n_vec.scale(-1.0);
        }
            break;

        case 2:
        {
            // first calculate the z-vector
            v1.copyVector(this->getNode(1));
            v1.add(-1.0, this->getNode(0));
            
            v2.copyVector(this->getNode(2));        
            v2.add(-1.0, this->getNode(0));
            
            v1.crossProduct(v2, v3);
            v3.scaleToUnitLength(); // z-axis
            
            // get the nodes for specified boundary ID
            const std::map<FESystemUInt, std::vector<FESystemUInt> >& b_n_map = this->getBoundaryIDAndBoundaryNodeMap();
            std::map<FESystemUInt, std::vector<FESystemUInt> >::const_iterator it = b_n_map.find(b_id);
            FESystemAssert1( it != b_n_map.end(), FESystem::Exception::InvalidID, b_id);
            FESystemAssert2(it->second.size() == 2, FESystem::Exception::DimensionsDoNotMatch, it->second.size(), 2);
            std::vector<FESystemUInt>::const_reverse_iterator s_it = it->second.rbegin();
            
            v2.copyVector(this->getNode(*s_it++));
            v2.add(-1.0, this->getNode(*s_it));
            v2.crossProduct(v3, n_vec); // cross-product with the z-vector gives the outward normal
            n_vec.scaleToUnitLength(); // scale to unit length
        }
            break;

        case 3:
        {
            // get the nodes for specified boundary ID
            const std::map<FESystemUInt, std::vector<FESystemUInt> >& b_n_map = this->getBoundaryIDAndBoundaryNodeMap();
            std::map<FESystemUInt, std::vector<FESystemUInt> >::const_iterator it = b_n_map.find(b_id);
            FESystemAssert1( it != b_n_map.end(), FESystem::Exception::InvalidID, b_id);
            FESystemAssert2(it->second.size() >= 3, FESystem::Exception::DimensionsDoNotMatch, it->second.size(), 3);
            std::vector<FESystemUInt>::const_iterator s_it = it->second.begin();
            std::vector<FESystemUInt>::const_reverse_iterator s_it_reverse = it->second.rbegin();
            
            v2.copyVector(this->getNode(*s_it_reverse++));
            v2.add(-1.0, this->getNode(*s_it));
            v1.copyVector(this->getNode(*s_it_reverse++));
            v1.add(-1.0, this->getNode(*s_it));

            v1.crossProduct(v2, n_vec); // cross-product with the z-vector gives the outward normal
            n_vec.scaleToUnitLength(); // scale to unit length
        }
            break;

        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
}



std::auto_ptr<FESystem::Mesh::ElemBase> 
FESystem::Mesh::ElementCreate(FESystem::Mesh::MeshBase& mesh, FESystem::Mesh::ElementType elem_type, FESystemBoolean if_local_physical_cs_same_as_global)
{
    std::auto_ptr<FESystem::Mesh::ElemBase> elem;
    
    switch (elem_type)
    {
        case FESystem::Mesh::EDGE2:
            elem.reset(new FESystem::Mesh::Edge2(if_local_physical_cs_same_as_global));
            break;

        case FESystem::Mesh::EDGE3:
            elem.reset(new FESystem::Mesh::Edge3(if_local_physical_cs_same_as_global));
            break;

        case FESystem::Mesh::TRI3:
            elem.reset(new FESystem::Mesh::Tri3(if_local_physical_cs_same_as_global));
            break;

        case FESystem::Mesh::TRI6:
            elem.reset(new FESystem::Mesh::Tri6(if_local_physical_cs_same_as_global));
            break;

        case FESystem::Mesh::TRI7:
            elem.reset(new FESystem::Mesh::Tri7(if_local_physical_cs_same_as_global));
            break;

        case FESystem::Mesh::QUAD4:
            elem.reset(new FESystem::Mesh::Quad4(if_local_physical_cs_same_as_global));
            break;

        case FESystem::Mesh::QUAD8:
            elem.reset(new FESystem::Mesh::Quad8(if_local_physical_cs_same_as_global));
            break;

        case FESystem::Mesh::QUAD9:
            elem.reset(new FESystem::Mesh::Quad9(if_local_physical_cs_same_as_global));
            break;

        default:
            FESystemAssert1(false, FESystem::Exception::EnumerationNotHandled, elem_type);
            break;
    }

    return elem;
}


