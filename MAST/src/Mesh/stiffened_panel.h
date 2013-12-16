//
//  stiffened_panel.h
//  MAST
//
//  Created by Manav Bhatia on 12/9/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_stiffened_panel_h__
#define __MAST_stiffened_panel_h__

// MAST includes
#include "Mesh/mesh_initializer.h"


namespace MAST {
    class StiffenedPanel {
    public:
        StiffenedPanel() { }
        
        ~StiffenedPanel() { }
        
        
        void init (const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
                   UnstructuredMesh& mesh, ElemType t);
        
    protected:
        
        enum Component {
            PANEL,
            STIFFENER_X,
            STIFFENER_Y
        };
        
        void _combine_mesh(UnstructuredMesh& panel,
                           UnstructuredMesh& stiffener,
                           MAST::StiffenedPanel::Component c,
                           Real stiff_offset,
                           subdomain_id_type sid);

    };
}


void
MAST::StiffenedPanel::init(const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
                           UnstructuredMesh& mesh, ElemType t) {
    
    libmesh_assert_equal_to(divs.size(), 3);
    
    std::vector<MeshInitializer::CoordinateDivisions*> stiff_divs(2);
    MeshInitializer init;
    stiff_divs[0] = divs[0];
    stiff_divs[1] = divs[1];

    // initialize the main mesh
    {
        SerialMesh panel(mesh.comm());
        init.init(stiff_divs, panel, t);
        // use subdomain id for panel as 0
        _combine_mesh(mesh, panel, MAST::StiffenedPanel::PANEL, 0., 0);
    }
    
    const unsigned int n_x_stiff = divs[1]->n_divs()-1,
    n_y_stiff = divs[0]->n_divs()-1;
    
    
    // now iterate over the stiffeners and create them
    for (unsigned int i=0; i<n_x_stiff; i++) {
        SerialMesh stiff(mesh.comm());
        stiff_divs[0] = divs[0];
        stiff_divs[1] = divs[2];
        
        init.init(stiff_divs, stiff, t);
        
        // add the elements and nodes to the panel mesh
        // the y-coordinate for this mesh will be the x-coordinate
        // subdomain id for each x-stiffener is 1+i
        _combine_mesh(mesh, stiff, MAST::StiffenedPanel::STIFFENER_X,
                     divs[1]->div_location(i+1), i+1);
    }
    
    for (unsigned int i=0; i<n_y_stiff; i++) {
        SerialMesh stiff(mesh.comm());
        stiff_divs[0] = divs[1];
        stiff_divs[1] = divs[2];
        
        init.init(stiff_divs, stiff, t);
        
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
MAST::StiffenedPanel::_combine_mesh(UnstructuredMesh& panel,
                                    UnstructuredMesh& component,
                                    MAST::StiffenedPanel::Component c,
                                    Real stiff_offset,
                                    subdomain_id_type sid) {
    BoundaryInfo &panel_binfo = *panel.boundary_info,
    &component_binfo = *component.boundary_info;

    panel.reserve_nodes(component.n_nodes());
    panel.reserve_elem(component.n_elem());
    
    MeshBase::const_element_iterator
    el_it = component.elements_begin(),
    el_end = component.elements_end();
    
    std::map<Node*, Node*> old_to_new;
    Node *old_node, *new_node;
    
    for ( ; el_it != el_end; el_it++ ) {
        Elem* old_elem = *el_it,
        *new_elem = panel.add_elem(Elem::build(old_elem->type()).release());
        new_elem->subdomain_id() = sid;
        
        // add boundary condition tags for the panel boundary
        if (c == MAST::StiffenedPanel::PANEL)
            for (unsigned short int n=0; n<old_elem->n_sides(); n++)
                if (component_binfo.n_boundary_ids(old_elem, n)) {
                    // add the boundary tags to the panel mesh
                    std::vector<boundary_id_type> bc_ids = component_binfo.boundary_ids(old_elem, n);
                    for ( unsigned int bid=0; bid < bc_ids.size(); bid++)
                        panel_binfo.add_side(new_elem, n, bc_ids[bid]);
                }
        
        
        for (unsigned int n=0; n<old_elem->n_nodes(); n++) {
            old_node = old_elem->get_node(n);
            
            if (!old_to_new.count(old_node)) {
                Point p;
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
                    MeshBase::node_iterator
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
