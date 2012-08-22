//
//  FiniteElementBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 3/22/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem includes
#include "FiniteElems/FiniteElementBase.h"
#include "Numerics/MatrixBase.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "Functions/CompositeFunctionMappingBase.h"
#include "Base/FESystemExceptions.h"
#include "Mesh/ElemBase.h"


FESystem::FiniteElement::FiniteElementBase::FiniteElementBase():
if_initialized(false),
geom_elem(NULL),
composite_map(NULL)
{
    
}

FESystem::FiniteElement::FiniteElementBase::~FiniteElementBase()
{
    
}
        

void 
FESystem::FiniteElement::FiniteElementBase::clear()
{
    this->if_initialized = false;
    this->geom_elem = NULL;
}


void
FESystem::FiniteElement::FiniteElementBase::reinit(const FESystem::Mesh::ElemBase& element)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    
    this->geom_elem = &element;
    
    this->initializeMaps();
            
    // reset the initialization flag once this is done. 
    this->if_initialized = true;

}


const FESystem::Mesh::ElemBase&
FESystem::FiniteElement::FiniteElementBase::getGeometricElement() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return *(this->geom_elem);
}



void
FESystem::FiniteElement::FiniteElementBase::getShapeFunction(const FESystem::Numerics::VectorBase<FESystemDouble>& vin,
                                                             FESystem::Numerics::VectorBase<FESystemDouble>& vout) const
{
    // make sure that the element is initialized before access to any of these data structures is requested
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);

    static FESystem::Numerics::LocalVector<FESystemDouble> vec1;
    vec1.resize(this->geom_elem->getParentNondegenerateElem().getNNodes());
    
    if (this->geom_elem->ifDegerateElement())
    {
        vec1.zero();
        this->composite_map->getDiscreteFunction(vin, vec1);
        this->geom_elem->getParentToDegenerateElemMappingMatrix().rightVectorMultiply(vec1, vout);
    }
    else
        this->composite_map->getDiscreteFunction(vin, vout);
}
        

void 
FESystem::FiniteElement::FiniteElementBase::getShapeFunctionForBoundary(const FESystemUInt b_id, const FESystem::Numerics::VectorBase<FESystemDouble>& vin, FESystem::Numerics::VectorBase<FESystemDouble>& vout) const
{
    const FESystemUInt elem_dim=this->geom_elem->getDimension();
    
    static FESystem::Numerics::LocalVector<FESystemDouble> vin_b, v2, v3;
    vin_b.resize(elem_dim);
    vin_b.zero();
    v2.resize(elem_dim);
    v3.resize(elem_dim);
    
    FESystemUInt dof=0, n_assigned=0;
    FESystemDouble dof_val=0.0;
    
    // get the coordinate that remains constant for the boundary
    this->geom_elem->getConstantCoordinateIDAndValueForBoundary(b_id, dof, dof_val);
    
    // set the coordinate value: the constant coordinate is borrowed from the boundary, and the other are set from the given vector
    for (FESystemUInt i=0; i<elem_dim; i++)
        if (i == dof)
            vin_b.setVal(i, dof_val);
        else 
            vin_b.setVal(i, vin.getVal(n_assigned++));
    
    // get the shape function
    this->getShapeFunction(vin_b, vout); 
}


void
FESystem::FiniteElement::FiniteElementBase::getShapeFunctionDerivativeForLocalCoordinates(const std::vector<FESystemUInt>& coordinate_derivative_order,
                                                                                          const FESystem::Numerics::VectorBase<FESystemDouble>& vin, 
                                                                                          FESystem::Numerics::VectorBase<FESystemDouble>& vout) const
{
    // make sure that the element is initialized before access to any of these data structures is requested
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    static FESystem::Numerics::LocalVector<FESystemDouble> vec1;
    vec1.resize(this->geom_elem->getParentNondegenerateElem().getNNodes());
    
    if (this->geom_elem->ifDegerateElement())
    {
        vec1.zero();
        this->composite_map->getOriginalFunctionMapping().getDiscreteFunctionDerivative(coordinate_derivative_order, vin, vec1);
        this->geom_elem->getParentToDegenerateElemMappingMatrix().rightVectorMultiply(vec1, vout);
    }
    else
        this->composite_map->getOriginalFunctionMapping().getDiscreteFunctionDerivative(coordinate_derivative_order, vin, vout);
}



void
FESystem::FiniteElement::FiniteElementBase::getShapeFunctionDerivativeForLocalCoordinatesForBoundary(const FESystemUInt b_id,
                                                                                                     const std::vector<FESystemUInt>& coordinate_derivative_order,
                                                                                                     const FESystem::Numerics::VectorBase<FESystemDouble>& vin,
                                                                                                     FESystem::Numerics::VectorBase<FESystemDouble>& vout) const
{
    const FESystemUInt elem_dim=this->geom_elem->getDimension();
    
    static FESystem::Numerics::LocalVector<FESystemDouble> vin_b, v2, v3;
    vin_b.resize(elem_dim);
    vin_b.zero();
    v2.resize(elem_dim);
    v3.resize(elem_dim);
    
    FESystemUInt dof=0, n_assigned=0;
    FESystemDouble dof_val=0.0;
    
    // get the coordinate that remains constant for the boundary
    this->geom_elem->getConstantCoordinateIDAndValueForBoundary(b_id, dof, dof_val);
    
    // set the coordinate value: the constant coordinate is borrowed from the boundary, and the other are set from the given vector
    for (FESystemUInt i=0; i<elem_dim; i++)
        if (i == dof)
            vin_b.setVal(i, dof_val);
        else 
            vin_b.setVal(i, vin.getVal(n_assigned++));
    
    // get the shape function derivative
    this->getShapeFunctionDerivativeForLocalCoordinates(coordinate_derivative_order, vin_b, vout); 
}


void
FESystem::FiniteElement::FiniteElementBase::getShapeFunctionDerivativeForPhysicalCoordinates(const std::vector<FESystemUInt>& coordinate_derivative_order,
                                                                                             const FESystem::Numerics::VectorBase<FESystemDouble>& vin, 
                                                                                             FESystem::Numerics::VectorBase<FESystemDouble>& vout) const
{
    // make sure that the element is initialized before access to any of these data structures is requested
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    static FESystem::Numerics::LocalVector<FESystemDouble> vec1;
    vec1.resize(this->geom_elem->getParentNondegenerateElem().getNNodes());
    
    if (this->geom_elem->ifDegerateElement())
    {
        vec1.zero();
        this->composite_map->getDiscreteFunctionDerivative(coordinate_derivative_order, vin, vec1);
        this->geom_elem->getParentToDegenerateElemMappingMatrix().rightVectorMultiply(vec1, vout);
    }
    else
        this->composite_map->getDiscreteFunctionDerivative(coordinate_derivative_order, vin, vout);
}

        
void 
FESystem::FiniteElement::FiniteElementBase::getShapeFunctionDerivativeForPhysicalCoordinatesForBoundary(const FESystemUInt b_id,
                                                                                                        const std::vector<FESystemUInt>& coordinate_derivative_order,
                                                                                                        const FESystem::Numerics::VectorBase<FESystemDouble>& vin,
                                                                                                        FESystem::Numerics::VectorBase<FESystemDouble>& vout) const
{
    const FESystemUInt elem_dim=this->geom_elem->getDimension();
    
    static FESystem::Numerics::LocalVector<FESystemDouble> vin_b, v2, v3;
    vin_b.resize(elem_dim);
    vin_b.zero();
    v2.resize(elem_dim);
    v3.resize(elem_dim);
    
    FESystemUInt dof=0, n_assigned=0;
    FESystemDouble dof_val=0.0;
    
    // get the coordinate that remains constant for the boundary
    this->geom_elem->getConstantCoordinateIDAndValueForBoundary(b_id, dof, dof_val);
    
    // set the coordinate value: the constant coordinate is borrowed from the boundary, and the other are set from the given vector
    for (FESystemUInt i=0; i<elem_dim; i++)
        if (i == dof)
            vin_b.setVal(i, dof_val);
        else 
            vin_b.setVal(i, vin.getVal(n_assigned++));
    
    // get the shape function derivative
    this->getShapeFunctionDerivativeForPhysicalCoordinates(coordinate_derivative_order, vin_b, vout); 
}


void
FESystem::FiniteElement::FiniteElementBase::getJacobianMatrix(const FESystem::Numerics::VectorBase<FESystemDouble>& vin, FESystem::Numerics::MatrixBase<FESystemDouble>& jac) const
{
    // make sure that the element is initialized before access to any of these data structures is requested
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    this->composite_map->getSubFunctionMapping().getMappingJacobian(vin, jac);
}
        


FESystemDouble 
FESystem::FiniteElement::FiniteElementBase::getJacobianValue(const FESystem::Numerics::VectorBase<FESystemDouble>& vin) const
{
    static FESystem::Numerics::DenseMatrix<FESystemDouble> jac;
    jac.resize(this->geom_elem->getDimension(), this->geom_elem->getDimension());
    jac.zero();
    
    this->getJacobianMatrix(vin, jac);
    return jac.getDeterminant();
}


FESystemDouble 
FESystem::FiniteElement::FiniteElementBase::getJacobianValueForBoundary(const FESystemUInt b_id, const FESystem::Numerics::VectorBase<FESystemDouble>& vin) const
{
    // return 1 for the 1-D element 
    const FESystemUInt elem_dim=this->geom_elem->getDimension();
    if (elem_dim == 1)
        return 1.0;
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> jac;
    static FESystem::Numerics::LocalVector<FESystemDouble> vin_b, v2, v3;
    jac.resize(elem_dim, elem_dim);
    jac.zero();
    vin_b.resize(elem_dim);
    vin_b.zero();
    v2.resize(elem_dim);
    v3.resize(elem_dim);
    
    FESystemUInt dof=0, n_assigned=0, id1=0,id2=0;
    FESystemDouble dof_val=0.0;

    // get the coordinate that remains constant for the boundary
    this->geom_elem->getConstantCoordinateIDAndValueForBoundary(b_id, dof, dof_val);

    // set the coordinate value: the constant coordinate is borrowed from the boundary, and the other are set from the given vector
    for (FESystemUInt i=0; i<elem_dim; i++)
        if (i == dof)
            vin_b.setVal(i, dof_val);
        else 
            vin_b.setVal(i, vin.getVal(n_assigned++));
    
    // get the jacobian matrix 
    this->getJacobianMatrix(vin_b, jac);
    
    // the Jacobian is d X_i / d xi_j . The element length for the \xi_j ^th coordinate is the j^th column
    if (this->geom_elem->getDimension() == 2) // calculate the norm of the colum, which is the length of the element length
    {
        id1 = 1-dof; // if dof=0, then the second column is of interest, else, the first column is of interest
        jac.getColumnVals(id1, 0, elem_dim-1, vin_b);
        return vin_b.getL2Norm();
    }
    else // for a 3-D element, calculate the cross product of the element
    {
        if (dof==0)
        {
            id1=1; id2=2;
        }
        else if (dof==1)
        {
            id1=0; id2=2;
        }
        else  // dof == 2
        {
            id1=0; id2=1;
        }

        jac.getColumnVals(id1, 0, elem_dim, vin_b);
        jac.getColumnVals(id2, 0, elem_dim, v2);
        vin_b.crossProduct(v2, v3); // this calculates the differential area
        
        return v3.getL2Norm();
    }
    
}

