//
//  DegreeOfFreedomMap.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/7/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem includes
#include "DegreesOfFreedom/DegreeOfFreedomMap.h"
#include "DegreesOfFreedom/DegreeOfFreedomUnit.h"
#include "Mesh/MeshBase.h"
#include "Mesh/Node.h"
#include "Mesh/ElemBase.h"
#include "Base/FESystemExceptions.h"
#include "Numerics/MatrixBase.h"
#include "Numerics/VectorBase.h"
#include "Numerics/SparsityPattern.h"


FESystem::DegreeOfFreedom::DegreeOfFreedomMap::DegreeOfFreedomMap(FESystem::Mesh::MeshBase& m):
mesh(m),
sparsity_pattern(NULL),
if_h_adaptable(false),
if_p_adaptable(false),
n_dofs(0)
{
    
}



FESystem::DegreeOfFreedom::DegreeOfFreedomMap::~DegreeOfFreedomMap()
{
    delete this->sparsity_pattern;
}


FESystemUInt
FESystem::DegreeOfFreedom::DegreeOfFreedomMap::getNVariables() const
{
    return this->variable_vec.size();
}



FESystemUInt 
FESystem::DegreeOfFreedom::DegreeOfFreedomMap::addVariable(std::string& var_name, FESystemUInt order)
{
    // make sure the name does not already exist in the system. 
    for (FESystemUInt i=0; i<this->variable_vec.size(); i++)
        FESystemAssert1(this->variable_vec[i].name != var_name, FESystem::DegreeOfFreedom::VariableNameAlreadyUsed, var_name);

    // create the variable and set the appropriate values
    FESystem::DegreeOfFreedom::Variable var;
    var.name = var_name;
    var.order = order;
    var.variable_id = this->variable_vec.size();
    var.system_id = 0; // so far this is unused, but may be used in future
    
    this->variable_vec.push_back(var);
    
    // return the id of the variable
    return var.variable_id;
}


FESystemUInt
FESystem::DegreeOfFreedom::DegreeOfFreedomMap::getNDofs() const
{
    return this->n_dofs;
}


const FESystem::Numerics::SparsityPattern&
FESystem::DegreeOfFreedom::DegreeOfFreedomMap::getSparsityPattern() const
{
    return *(this->sparsity_pattern);
}



void 
FESystem::DegreeOfFreedom::DegreeOfFreedomMap::reinit()
{
    // get nodes and elements
    std::vector<FESystem::Mesh::Node*>& nodes = this->mesh.getNodes();
    std::vector<FESystem::Mesh::ElemBase*>& elems = this->mesh.getElements();
    
    // this is the counter of the global dof id
    FESystemUInt global_dof_id = 0;
    
    // iterate over each node and element and initialize their dofs and the connectivity model
    for (FESystemUInt i_node=0; i_node<nodes.size(); i_node++)
    {
        nodes[i_node]->init(this->getNVariables());
        for (FESystemUInt i_var=0; i_var<this->getNVariables(); i_var++)
        {
            FESystem::DegreeOfFreedom::DegreeOfFreedomUnit& dof = nodes[i_node]->addDegreeOfFreedomUnit(i_var);
            dof.system_id = 0;
            dof.variable_id = i_var;
            dof.p_level = this->variable_vec[i_var].order;
            dof.if_p_adaptable = false;   // none of the nodes dofs are p-adaptable. Only the element interior dofs can be incremented
            dof.global_dof_id.push_back(global_dof_id++);  
            this->dof_object_vector.push_back(nodes[i_node]);
        }
    }
    
    // this needs to be revisited, as independent dofs will need to be added to the sub-dimensional elements, in addition to the
    // elem itself. That is, for a 3-D element, element will get its bubble function, in addition to bubble functions of the 
    // faces and bubble function for the edges. 
//    // now iterate over all the elements and add the degrees of freedom for the variables
//    for (FESystemUInt i_elem=0; i_elem<elems.size(); i_elem++)
//    {
//        elems[i_elem]->init(this->getNVariables());
//        for (FESystemUInt i_var=0; i_var<this->getNVariables(); i_var++)
//        {
//            FESystem::DegreeOfFreedom::DegreeOfFreedomUnit& dof = elems[i_elem]->addDegreeOfFreedomUnit(i_var);
//            dof.system_id = 0;
//            dof.variable_id = i_var;
//            dof.p_level = this->variable_vec[i_var].order;
//            dof.if_p_adaptable = this->if_p_adaptable;   // all node dofs are not p-adaptable. Only the element interior dofs can be incremented
//            // currently the p adaptability is not fully implemented
//            if (if_p_adaptable)
//            {
//                FESystemAssert0(false, FESystem::Exception::InvalidState);
//            }
//            else
//                // add the number of variables
//                dof.global_dof_id.push_back(global_dof_id++);  
//        }   
//    }
    
    this->n_dofs = global_dof_id;
    
    //    // calculate the fill reducing ordering
    //    if (this->sparsity_pattern == NULL)
    //        this->sparsity_pattern = new FESystem::Numerics::SparsityPattern;
    //    this->sparsity_pattern->clear();
    //    
    //    this->sparsity_pattern->nonzero_column_ids_per_row.resize(this->getNDofs());
    //    
    //    // iterate over the elements
    //    std::vector<FESystem::Mesh::ElemBase*>::const_iterator el_it=elems.begin();
    //    
    //    // this is the counter of the global dof id
    //    std::set<FESystemUInt> ids;
    //    std::set<FESystemUInt>::const_iterator s_it1;
    //    
    //    // iterate over each node and element and initialize their dofs and the connectivity model
    //    for ( ; el_it != elems.end(); el_it++)
    //    {
    //        ids.clear();
    //        for (FESystemUInt i_nodes=0; i_nodes<(*el_it)->getNNodes(); i_nodes++)
    //            for (FESystemUInt i_dofs=0; i_dofs<(*el_it)->getNode(i_nodes).getNDegreeOfFreedomUnits(); i_dofs++)
    //                // now add the connectivity for the respective degree of freedom units
    //                ids.insert((*el_it)->getNode(i_nodes).getDegreeOfFreedomUnit(i_dofs).global_dof_id[0]); //  currently works only for a single p-level
    //        
    //        // now that the dofs ids for this element are known, tell each id that it is connected to the others; this requires a double loop
    //        for (s_it1 = ids.begin(); s_it1 != ids.end(); s_it1++)
    //            this->sparsity_pattern->addNonZeroColumnsForRow(*s_it1, ids);
    //    }
    //    
    //    // calculate the reordering
    //    std::vector<FESystemUInt> reordered_dofs(this->n_dofs);
    //    this->sparsity_pattern->calculateFillReducingOrdering(reordered_dofs);
    //    
    //    // now reassign the dofs and clear the sparsity pattern
    //    for (FESystemUInt i=0; i<this->dof_object_vector.size(); i++)
    //        for (FESystemUInt j=0; j<this->dof_object_vector[i]->getNDegreeOfFreedomUnits(); j++)
    //        {
    //            std::vector<FESystemUInt>& global_dof_id = this->dof_object_vector[i]->getDegreeOfFreedomUnit(j).global_dof_id;
    //            global_dof_id[0] = reordered_dofs[global_dof_id[0]];
    //        }
    //        

    // now that the dofs have been distributed, initialize the sparsity pattern
    this->initializeSparsityPattern();

}

    
    
template <typename ValType>
void
FESystem::DegreeOfFreedom::DegreeOfFreedomMap::addToGlobalMatrix(const FESystem::Mesh::ElemBase& elem, const FESystem::Numerics::MatrixBase<ValType>& elem_mat, 
                                                      FESystem::Numerics::MatrixBase<ValType>& global_mat) const
{
    FESystemUInt n_dofs = this->getNDofsForElem(elem);

    // create the sequence of the dofs    
    std::vector<FESystemUInt> dof_ids;
    if (dof_ids.size() != n_dofs)
        dof_ids.resize(n_dofs);
    
    n_dofs=0;
    for (FESystemUInt i_var=0; i_var<this->getNVariables(); i_var++)
        for (FESystemUInt i_node=0; i_node<elem.getNNodes(); i_node++)
            dof_ids[n_dofs++]=elem.getNode(i_node).getDegreeOfFreedomUnit(i_var).global_dof_id[0]; // this will change when p-adaptivity is added. 

    // make sure that the matrix sizes are accurate
    std::pair<FESystemUInt, FESystemUInt> s_elem = elem_mat.getSize(), 
    s_global = global_mat.getSize();

    FESystemAssert0(s_elem.first == s_elem.second, FESystem::Numerics::MatrixMustBeSquareForOperation);
    FESystemAssert0(s_global.first == s_global.second, FESystem::Numerics::MatrixMustBeSquareForOperation);
    FESystemAssert4(s_elem.first == dof_ids.size(), FESystem::Numerics::MatrixSizeMismatch, dof_ids.size(), dof_ids.size(), s_elem.first, s_elem.second);
    FESystemAssert4(s_global.first == this->getNDofs(), FESystem::Numerics::MatrixSizeMismatch, this->getNDofs(), this->getNDofs(), s_global.first, s_global.second);
    
    
    // now that the ids are available, iterate over all and assemble the matrix
    global_mat.addVal(dof_ids, dof_ids, elem_mat);
}



FESystemUInt
FESystem::DegreeOfFreedom::DegreeOfFreedomMap::getNDofsForElem(const FESystem::Mesh::ElemBase& elem) const
{
    FESystemUInt n_dofs=0;
    for (FESystemUInt i_var=0; i_var<this->getNVariables(); i_var++)
        for (FESystemUInt i_node=0; i_node<elem.getNNodes(); i_node++)
            n_dofs++; // this will change when p-adaptivity is added. 
    return n_dofs;
}
    
    
template <typename ValType>
void
FESystem::DegreeOfFreedom::DegreeOfFreedomMap::addToGlobalVector(const FESystem::Mesh::ElemBase& elem, const FESystem::Numerics::VectorBase<ValType>& elem_vec, 
                                                      FESystem::Numerics::VectorBase<ValType>& global_vec) const
{
    FESystemUInt n_dofs = this->getNDofsForElem(elem);
    
    // create the sequence of the dofs    
    std::vector<FESystemUInt> dof_ids;
    if (dof_ids.size() != n_dofs)
        dof_ids.resize(n_dofs);
    
    n_dofs=0;
    for (FESystemUInt i_var=0; i_var<this->getNVariables(); i_var++)
        for (FESystemUInt i_node=0; i_node<elem.getNNodes(); i_node++)
            dof_ids[n_dofs++]=elem.getNode(i_node).getDegreeOfFreedomUnit(i_var).global_dof_id[0]; // this will change when p-adaptivity is added. 
    
    // make sure that the matrix sizes are accurate
    FESystemUInt s_elem = elem_vec.getSize(),  s_global = global_vec.getSize();
    
    FESystemAssert2(s_elem == dof_ids.size(), FESystem::Exception::DimensionsDoNotMatch, s_elem, dof_ids.size());
    FESystemAssert2(s_global == this->getNDofs(), FESystem::Exception::DimensionsDoNotMatch, s_global, this->getNDofs());
    
    
    // now that the ids are available, iterate over all and assemble the vector
    for (FESystemUInt i=0; i<dof_ids.size(); i++)
        global_vec.addVal(dof_ids[i], elem_vec.getVal(i));
}




template <typename ValType>
void
FESystem::DegreeOfFreedom::DegreeOfFreedomMap::getFromGlobalVector(const FESystem::Mesh::ElemBase& elem, const FESystem::Numerics::VectorBase<ValType>& global_vec,
                                                        FESystem::Numerics::VectorBase<ValType>& elem_vec) const
{
    // create the sequence of the dofs    
    std::vector<FESystemUInt> dof_ids;
    
    for (FESystemUInt i_var=0; i_var<this->getNVariables(); i_var++)
        for (FESystemUInt i_node=0; i_node<elem.getNNodes(); i_node++)
            dof_ids.push_back(elem.getNode(i_node).getDegreeOfFreedomUnit(i_var).global_dof_id[0]); // this will change when p-adaptivity is added. 
    
    // make sure that the matrix sizes are accurate
    FESystemUInt s_elem = elem_vec.getSize(),  s_global = global_vec.getSize();
    
    FESystemAssert2(s_elem == dof_ids.size(), FESystem::Exception::DimensionsDoNotMatch, s_elem, dof_ids.size());
    FESystemAssert2(s_global == this->getNDofs(), FESystem::Exception::DimensionsDoNotMatch, s_global, this->getNDofs());
    
    
    // now that the ids are available, and extract the vector
    global_vec.getSubVectorValsFromIndices(dof_ids, elem_vec);
}





void
FESystem::DegreeOfFreedom::DegreeOfFreedomMap::initializeSparsityPattern()
{
    if (this->sparsity_pattern == NULL)
        this->sparsity_pattern = new FESystem::Numerics::SparsityPattern;
    this->sparsity_pattern->clear();
    
    this->sparsity_pattern->nonzero_column_ids_per_row.resize(this->getNDofs());
    
    // iterate over the elements
    const std::vector<FESystem::Mesh::ElemBase*>& elems = this->mesh.getElements();
    std::vector<FESystem::Mesh::ElemBase*>::const_iterator el_it=elems.begin();
    
    // this is the counter of the global dof id
    std::set<FESystemUInt> ids;
    std::set<FESystemUInt>::const_iterator s_it1;
        
    // iterate over each node and element and initialize their dofs and the connectivity model
    for ( ; el_it != elems.end(); el_it++)
    {
        ids.clear();
        for (FESystemUInt i_nodes=0; i_nodes<(*el_it)->getNNodes(); i_nodes++)
            for (FESystemUInt i_dofs=0; i_dofs<(*el_it)->getNode(i_nodes).getNDegreeOfFreedomUnits(); i_dofs++)
                // now add the connectivity for the respective degree of freedom units
                ids.insert((*el_it)->getNode(i_nodes).getDegreeOfFreedomUnit(i_dofs).global_dof_id[0]); //  currently works only for a single p-level
        
        // now that the dofs ids for this element are known, tell each id that it is connected to the others; this requires a double loop
        for (s_it1 = ids.begin(); s_it1 != ids.end(); s_it1++)
            this->sparsity_pattern->addNonZeroColumnsForRow(*s_it1, ids);
    }
    
    std::vector<FESystemUInt> reordered_ids(this->n_dofs);
    this->sparsity_pattern->calculateFillReducingOrdering(reordered_ids);
    

    // now redefine the ordeing
    std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
    std::vector<FESystem::Mesh::Node*>::iterator it = nodes.begin(), end = nodes.end();
    
    for ( ; it != end; it++)
        for (FESystemUInt i_var=0; i_var<this->getNVariables(); i_var++)
            for (FESystemUInt i_order=0; i_order<(*it)->getDegreeOfFreedomUnit(i_var).global_dof_id.size(); i_order++)
                (*it)->getDegreeOfFreedomUnit(i_var).global_dof_id[i_order] = reordered_ids[(*it)->getDegreeOfFreedomUnit(i_var).global_dof_id[i_order]];
    

    
    // now reset the data structure of sparsity pattern with the updated orderin
    this->sparsity_pattern->clear();
    
    this->sparsity_pattern->nonzero_column_ids_per_row.resize(this->getNDofs());
    
    // iterate over the elements
    el_it=elems.begin();
    
    // iterate over each node and element and initialize their dofs and the connectivity model
    for ( ; el_it != elems.end(); el_it++)
    {
        ids.clear();
        for (FESystemUInt i_nodes=0; i_nodes<(*el_it)->getNNodes(); i_nodes++)
            for (FESystemUInt i_dofs=0; i_dofs<(*el_it)->getNode(i_nodes).getNDegreeOfFreedomUnits(); i_dofs++)
                // now add the connectivity for the respective degree of freedom units
                ids.insert((*el_it)->getNode(i_nodes).getDegreeOfFreedomUnit(i_dofs).global_dof_id[0]); //  currently works only for a single p-level
        
        // now that the dofs ids for this element are known, tell each id that it is connected to the others; this requires a double loop
        for (s_it1 = ids.begin(); s_it1 != ids.end(); s_it1++)
            this->sparsity_pattern->addNonZeroColumnsForRow(*s_it1, ids);
    }
    
    
    // now tell the patten to reinit its data structures
    this->sparsity_pattern->reinit();
}




//********************************************************************************************************************************
// template instantiation

template void FESystem::DegreeOfFreedom::DegreeOfFreedomMap::addToGlobalMatrix(const FESystem::Mesh::ElemBase& elem, const FESystem::Numerics::MatrixBase<FESystemDouble>& elem_mat,
                                                                    FESystem::Numerics::MatrixBase<FESystemDouble>& global_mat) const;
template void FESystem::DegreeOfFreedom::DegreeOfFreedomMap::addToGlobalVector(const FESystem::Mesh::ElemBase& elem, const FESystem::Numerics::VectorBase<FESystemDouble>& elem_vec,
                                                                    FESystem::Numerics::VectorBase<FESystemDouble>& global_vec) const;
template void FESystem::DegreeOfFreedom::DegreeOfFreedomMap::getFromGlobalVector(const FESystem::Mesh::ElemBase& elem, const FESystem::Numerics::VectorBase<FESystemDouble>& global_vec,
                                                                      FESystem::Numerics::VectorBase<FESystemDouble>& elem_vec) const;

//********************************************************************************************************************************

