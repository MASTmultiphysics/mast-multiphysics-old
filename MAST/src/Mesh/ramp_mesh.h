//
//  ramp_mesh.h
//  MAST
//
//  Created by Manav Bhatia on 10/10/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_ramp_mesh_h__
#define __MAST_ramp_mesh_h__


// MAST includes
#include "Mesh/mesh_initializer.h"


class RampMesh2D: public MeshInitializer
{
public:
    RampMesh2D():
    _tc_ratio(0.),
    _x0(0.),
    _y0(0.), _y1(0.),
    _panel_bc_id(0),
    _symmetry_bc_id(0)
    { }
    
    /*!
     *   initializes the object with the division for each dimension. Sets-up the
     *   mesh for the provided information, and then uses the provided funciton
     *   move the mesh points.
     */
    virtual void init (const libMesh::Real tc,
                       const unsigned int panel_bc_id,
                       const unsigned int symmetry_bc_id,
                       const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
                       UnstructuredMesh& mesh, ElemType t);
    
protected:
    
    /*!
     *   move the mesh points to account for a bump, and apply boudnary condition
     */
    virtual void process_mesh( );
    
    /*!
     *   t/c ratio of the panel
     */
    libMesh::Real _tc_ratio;
    
    libMesh::Real _x0;
    
    libMesh::Real _y0, _y1;
    
    unsigned int _panel_bc_id, _symmetry_bc_id;
    
};




inline void
RampMesh2D::init (const libMesh::Real tc_ratio,
                   const unsigned int panel_bc_id,
                   const unsigned int symmetry_bc_id,
                   const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
                   UnstructuredMesh& mesh, ElemType t)
{
    libmesh_assert(divs.size() == 2);
    libmesh_assert(divs[0]->n_divs() == 2);
    libmesh_assert(divs[1]->n_divs() >= 1);
    libmesh_assert(panel_bc_id > 3);
    libmesh_assert(symmetry_bc_id > 3);
    
    _tc_ratio = tc_ratio;
    _x0 = divs[0]->div_location(1);
    _y0 = divs[1]->div_location(0);
    _y1 = divs[1]->div_location(1);
    _panel_bc_id = panel_bc_id;
    _symmetry_bc_id = symmetry_bc_id;
    MeshInitializer::init(divs, mesh, t);
}



inline void
RampMesh2D::process_mesh( )
{
    // note that we apply the boundary conditions before moving the mesh since
    // the application of boudnary conditions is contingent upon the panel surface
    // points lying on the y = _y0, and the mesh movement ends up altering that
    
    //march over all the elmeents and tag the sides that all lie on the panel suface
    MeshBase::element_iterator e_it = _mesh->elements_begin();
    
    for ( ; e_it != _mesh->elements_end(); e_it++)
    {
        // iterate over the sides of each element and check if all
        // nodes satisfy the requirement
        
        for (unsigned int i_side=0; i_side<(*e_it)->n_sides(); i_side++)
        {
            AutoPtr<Elem> side_elem ((*e_it)->side(i_side).release());
            std::vector<bool> side_on_panel(side_elem->n_nodes()),
            side_on_slip_wall(side_elem->n_nodes());
            std::fill(side_on_panel.begin(), side_on_panel.end(), false);
            std::fill(side_on_slip_wall.begin(), side_on_slip_wall.end(), false);
            
            for (unsigned int i_node=0; i_node<side_elem->n_nodes(); i_node++)
            {
                const Node& n = *(side_elem->get_node(i_node));
                if ((n(1)==_y0) && (n(0) >= _x0-1.0e-6))
                    side_on_panel[i_node] = true;
                
                if ((n(1)==_y0) && (n(0) <= _x0+1.0e-6))
                    side_on_slip_wall[i_node] = true;
            }
            
            // check for side on panel
            bool if_apply_bc = true;
            for (unsigned int i_node=0; i_node<side_elem->n_nodes(); i_node++)
                if_apply_bc = side_on_panel[i_node] && if_apply_bc;
            if (if_apply_bc)
                _mesh->boundary_info->add_side(*e_it, i_side, _panel_bc_id);
            
            // now check for the slip wall
            if_apply_bc = true;
            for (unsigned int i_node=0; i_node<side_elem->n_nodes(); i_node++)
                if_apply_bc = side_on_slip_wall[i_node] && if_apply_bc;
            if (if_apply_bc)
                _mesh->boundary_info->add_side(*e_it, i_side, _symmetry_bc_id);
        }
    }
    
    // set the boudnary id names
    _mesh->boundary_info->sideset_name(_panel_bc_id) = "Panel";
    _mesh->boundary_info->sideset_name(_symmetry_bc_id) = "Symmetry";
    
    // now move the mesh points
    MeshBase::node_iterator   n_it  = _mesh->nodes_begin();
    const Mesh::node_iterator n_end = _mesh->nodes_end();
    
    const libMesh::Real pi = acos(-1.);
    libMesh::Real x, y;
    
    for (; n_it != n_end; n_it++)
    {
        Node& n =  **n_it;
        
        // this is for sine bump
        if (n(0) >= _x0)
        {
            x = n(0);
            y = n(1);
            
            n(1) += _tc_ratio * ( 1. -(y-_y0)/(_y1-_y0) ) * (x-_x0);
        }
    }
}




#endif //__MAST_ramp_mesh_h__
