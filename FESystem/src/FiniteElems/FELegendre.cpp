//
//  FELegendre.cpp
//  FESystem
//
//  Created by Manav Bhatia on 1/25/13.
//
//

// FESystem includes
#include "FiniteElems/FELegendre.h"
#include "Functions/LegendreFunction.h"
#include "Functions/LagrangeFunction.h"
#include "Functions/LegendreFunctionMapping.h"
#include "Functions/LagrangeFunctionMapping.h"
#include "Functions/CompositeFunctionMappingBase.h"
#include "Mesh/ElemBase.h"
#include "Mesh/Node.h"
#include "Geom/Point.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "Utils/Table.h"


FESystem::FiniteElement::FELegendre::FELegendre():
FESystem::FiniteElement::FiniteElementBase(),
legendre_order(0)
{
    this->composite_map = new FESystem::Functions::DiscreteCompositeMappingFunctionBase<FESystemDouble>;
    this->lagrange_map = new FESystem::Functions::LagrangeFunctionMapping<FESystemDouble>;
    this->legendre_map = new FESystem::Functions::LegendreFunctionMapping<FESystemDouble>;
    this->lagrange_functions.resize(3); // keep a 3-D by default, and if the element needs less based on dimensionality, it will use less
    for (FESystemUInt i=0; i<this->lagrange_functions.size(); i++)
        this->lagrange_functions[i] = new FESystem::Functions::LagrangeFunction<FESystemDouble>;
    this->legendre_function = new FESystem::Functions::LegendreFunction<FESystemDouble>;
}


FESystem::FiniteElement::FELegendre::~FELegendre()
{
    this->clear();
    delete this->composite_map;
    delete this->lagrange_map;
    delete this->legendre_map;
    for (FESystemUInt i=0; i<this->lagrange_functions.size(); i++)
        delete this->lagrange_functions[i];
    delete this->legendre_function;
}



FESystem::FiniteElement::FiniteElementType
FESystem::FiniteElement::FELegendre::getFiniteElementType() const
{
    return FESystem::FiniteElement::FE_LEGENDRE;
}


void
FESystem::FiniteElement::FELegendre::clear()
{
    // iterate over the functions and clear them
    for (FESystemUInt i=0; i<this->lagrange_functions.size(); i++)
        this->lagrange_functions[i]->clear();
    this->composite_map->clear();
    this->lagrange_map->clear();
    this->legendre_map->clear();
    
    FESystem::FiniteElement::FiniteElementBase::clear();
}




void
FESystem::FiniteElement::FELegendre::reinit(const FESystem::Mesh::ElemBase &element, const FESystemUInt order)
{
    FESystem::FiniteElement::FiniteElementBase::reinit(element);
    this->legendre_order = order;
}


FESystemUInt
FESystem::FiniteElement::FELegendre::getNShapeFunctions() const
{
    return pow(this->legendre_order+1, this->geom_elem->getParentNondegenerateElem().getDimension());
}


void
FESystem::FiniteElement::FELegendre::initializeMaps()
{
    // make sure that the element has been set
    FESystemAssert0(this->geom_elem != NULL, FESystem::Exception::NULLQuantity);
    
    FESystemUInt dim = this->geom_elem->getParentNondegenerateElem().getDimension(), n_nodes_dim=0;
    
    // initialize the Lagrange polynomials with the point locations
    std::vector<FESystemDouble> vals(10), physical_loc(150); // array to store node locaitons
    FESystem::Numerics::LocalVector<FESystemDouble> node_locations, vec;
    FESystem::Numerics::DenseMatrix<FESystemDouble> node_physical_locations;
    if (vec.getSize() != 3)
        vec.resize(3);
    if (physical_loc.size() < dim*this->geom_elem->getParentNondegenerateElem().getNNodes())
        physical_loc.resize(dim*this->geom_elem->getParentNondegenerateElem().getNNodes());
    node_physical_locations.resize(&physical_loc[0], dim, this->geom_elem->getParentNondegenerateElem().getNNodes());
    
    for (FESystemUInt i=0; i<dim; i++)
    {
        n_nodes_dim = this->geom_elem->getParentNondegenerateElem().getNNodesAlongDimension(i);
        
        if (n_nodes_dim > vals.size())
            vals.resize(n_nodes_dim); // this is to minimize reallocation of memory, and to ensure that the code still runs
        node_locations.resize(n_nodes_dim, &vals[0]);
        
        for (FESystemUInt j=0; j<n_nodes_dim; j++)
            node_locations.setVal(j, this->geom_elem->getParentNondegenerateElem().getLocalComputationalCoordinateForNodeAlongDim(i, j));
        
        this->lagrange_functions[i]->reinit(node_locations);
    }
    
    // create the mapping routine
    this->lagrange_map->reinit(dim, this->lagrange_functions, this->geom_elem->getParentNondegenerateElem().getPointDimensionIDTable());
    this->legendre_map->reinit(dim, this->legendre_order, *this->legendre_function);
    
    // set the values of x-y-z locations in the map
    for (FESystemUInt i_node=0; i_node<this->geom_elem->getParentNondegenerateElem().getNNodes(); i_node++)
    {
        this->geom_elem->getParentNondegenerateElem().getNodeLocationInLocalPhysicalCoordinate(i_node, vec);
        
        for (FESystemUInt i_dim=0; i_dim<dim; i_dim++)
            node_physical_locations.setVal(i_dim,i_node,vec.getVal(i_dim));
    }
    
    this->lagrange_map->setDiscreteFunctionValues(node_physical_locations);
    
    // the Lagrange finite element uses the same mapping for both the variable and geometry mapping
    this->composite_map->reinit(*(this->legendre_map), *(this->lagrange_map));
    
}


