/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2014  Manav Bhatia
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef __MAST_stiffened_panel_h__
#define __MAST_stiffened_panel_h__

// MAST includes
#include "Mesh/mesh_initializer.h"

// libMesh includes
#include "libmesh/serial_mesh.h"


namespace MAST {
    class StiffenedPanel {
    public:
        StiffenedPanel() { }
        
        ~StiffenedPanel() { }
        
        
        void init (const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
                   libMesh::UnstructuredMesh& mesh, libMesh::ElemType t, bool beam_stiffeners);
        
    protected:
        
        enum Component {
            PANEL,
            STIFFENER_X,
            STIFFENER_Y
        };
        
        void _combine_mesh(libMesh::UnstructuredMesh& panel,
                           libMesh::UnstructuredMesh& stiffener,
                           MAST::StiffenedPanel::Component c,
                           Real stiff_offset,
                           libMesh::subdomain_id_type sid);

    };
}


inline
void
MAST::StiffenedPanel::init(const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
                           libMesh::UnstructuredMesh& mesh, libMesh::ElemType t, bool beam_stiffeners) {
    
    libmesh_assert_equal_to(divs.size(), 3);
    
    mesh.set_mesh_dimension(2);
    
    std::vector<MeshInitializer::CoordinateDivisions*> panel_divs(2), stiff_divs(2);
    MeshInitializer init;
    libMesh::ElemType stiff_t = t;
    panel_divs[0] = divs[0];
    panel_divs[1] = divs[1];
    if (beam_stiffeners) {
        stiff_divs.resize(1);
        switch (t) {
            case libMesh::TRI3:
            case libMesh::QUAD4:
                stiff_t = libMesh::EDGE2;
                break;
                
            case libMesh::TRI6:
            case libMesh::QUAD8:
            case libMesh::QUAD9:
                stiff_t = libMesh::EDGE3;
                break;
                
            default:
                libmesh_error();
        }
    }
    
    // initialize the main mesh
    {
        libMesh::SerialMesh panel(mesh.comm());
        init.init(panel_divs, panel, t);
        // use subdomain id for panel as 0
        _combine_mesh(mesh, panel, MAST::StiffenedPanel::PANEL, 0., 0);
    }
    
    const unsigned int n_x_stiff = divs[1]->n_divs()-1,
    n_y_stiff = divs[0]->n_divs()-1;
    
    
    // now iterate over the stiffeners and create them
    for (unsigned int i=0; i<n_x_stiff; i++) {
        libMesh::SerialMesh stiff(mesh.comm());
        stiff_divs[0] = divs[0];
        if (!beam_stiffeners)
            stiff_divs[1] = divs[2];
        
        init.init(stiff_divs, stiff, stiff_t);
        
        // add the elements and nodes to the panel mesh
        // the y-coordinate for this mesh will be the x-coordinate
        // subdomain id for each x-stiffener is 1+i
        _combine_mesh(mesh, stiff, MAST::StiffenedPanel::STIFFENER_X,
                     divs[1]->div_location(i+1), i+1);
    }
    
    for (unsigned int i=0; i<n_y_stiff; i++) {
        libMesh::SerialMesh stiff(mesh.comm());
        stiff_divs[0] = divs[1];
        if (!beam_stiffeners)
            stiff_divs[1] = divs[2];
        
        init.init(stiff_divs, stiff, stiff_t);
        
        // add the elements and nodes to the panel mesh
        // the y-coordinate for this mesh will be the x-coordinate
        // subdomain id for each y-stiffener is n_x_stiff + i
        _combine_mesh(mesh, stiff, MAST::StiffenedPanel::STIFFENER_Y,
                     divs[0]->div_location(i+1), n_x_stiff+i+1);
    }
    
    mesh.prepare_for_use();
}



inline
void
MAST::StiffenedPanel::_combine_mesh(libMesh::UnstructuredMesh& panel,
                                    libMesh::UnstructuredMesh& component,
                                    MAST::StiffenedPanel::Component c,
                                    Real stiff_offset,
                                    libMesh::subdomain_id_type sid) {
    libMesh::BoundaryInfo &panel_binfo = *panel.boundary_info,
    &component_binfo = *component.boundary_info;

    panel.reserve_nodes(component.n_nodes());
    panel.reserve_elem(component.n_elem());
    
    libMesh::MeshBase::const_element_iterator
    el_it = component.elements_begin(),
    el_end = component.elements_end();
    
    std::map<libMesh::Node*, libMesh::Node*> old_to_new;
    libMesh::Node *old_node, *new_node;
    
    for ( ; el_it != el_end; el_it++ ) {
        libMesh::Elem* old_elem = *el_it;
        
        libMesh::Elem *new_elem = panel.add_elem(libMesh::Elem::build(old_elem->type()).release());
        new_elem->subdomain_id() = sid;
        
        // add boundary condition tags for the panel boundary
        if (c == MAST::StiffenedPanel::PANEL)
            for (unsigned short int n=0; n<old_elem->n_sides(); n++)
                if (component_binfo.n_boundary_ids(old_elem, n)) {
                    // add the boundary tags to the panel mesh
                    std::vector<libMesh::boundary_id_type> bc_ids = component_binfo.boundary_ids(old_elem, n);
                    for ( unsigned int bid=0; bid < bc_ids.size(); bid++)
                        panel_binfo.add_side(new_elem, n, bc_ids[bid]);
                }
        
        
        for (unsigned int n=0; n<old_elem->n_nodes(); n++) {
            old_node = old_elem->get_node(n);
            
            if (!old_to_new.count(old_node)) {
                libMesh::Point p;
                switch (c) {
                    case MAST::StiffenedPanel::PANEL:
                        p = (*old_node);
                        break;
                        
                    case MAST::StiffenedPanel::STIFFENER_X:
                        p(0) = (*old_node)(0);
                        p(1) = stiff_offset;
                        p(2) = (*old_node)(1);
                        break;
                        
                    case MAST::StiffenedPanel::STIFFENER_Y:
                        p(0) = stiff_offset;
                        p(1) = (*old_node)(0);
                        p(2) = (*old_node)(1);
                        break;
                        
                    default:
                        libmesh_error();
                }

                // first see if the node can be found on the panel
                // this will be done only when z-coordinate is zero
                if ((c != MAST::StiffenedPanel::PANEL) &&
                    p(2) == 0.) {
                    libMesh::MeshBase::node_iterator
                    n_it = panel.nodes_begin(),
                    n_end = panel.nodes_end();
                    for ( ; n_it != n_end; n_it++)
                        if ((**n_it)(0) == p(0) &&
                            (**n_it)(1) == p(1)) {
                            old_to_new[old_node] = *n_it;
                            break;
                        }
                }
                else
                    old_to_new[old_node] = panel.add_point(p);
            }

            new_node = old_to_new[old_node];
            new_elem->set_node(n) = new_node;
        }
    }
}





#endif // __MAST_stiffened_panel_h__
