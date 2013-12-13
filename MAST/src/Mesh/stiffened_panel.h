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
        
        void combine_mesh(UnstructuredMesh& panel,
                          UnstructuredMesh& stiffener,
                          MAST::StiffenedPanel::Component c,
                          Real stiff_offset);
        
    };
}


void
MAST::StiffenedPanel::init(const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
                           UnstructuredMesh& mesh, ElemType t) {
    
    MeshInitializer init;
    // initialize the main mesh
    {
        SerialMesh panel(mesh.comm());
        init.init(divs, panel, t);
        combine_mesh(mesh, panel, MAST::StiffenedPanel::PANEL, 0.);
    }
    
    const unsigned int n_x_stiff = divs[1]->n_divs()-1,
    n_y_stiff = divs[0]->n_divs()-1;
    
    std::vector<MeshInitializer::CoordinateDivisions*> stiff_divs(2);
    MeshInitializer::CoordinateDivisions height;
    std::vector<Real> div_loc(2); div_loc[0] = 0.; div_loc[1] = 0.1;
    std::vector<Real> dx_rel(2); dx_rel[0] = 1.; dx_rel[1] = 1.;
    std::vector<unsigned int> n_dx(1); n_dx[0] = 4;
    height.init(1, div_loc, dx_rel, n_dx);
    
    // now iterate over the stiffeners and create them
    for (unsigned int i=0; i<n_x_stiff; i++) {
        SerialMesh stiff(mesh.comm());
        stiff_divs[0] = divs[0];
        stiff_divs[1] = &height;
        
        init.init(stiff_divs, stiff, t);
        
        // add the elements and nodes to the panel mesh
        // the y-coordinate for this mesh will be the x-coordinate
        combine_mesh(mesh, stiff, MAST::StiffenedPanel::STIFFENER_X,
                     divs[1]->div_location(i+1));
    }
    
    for (unsigned int i=0; i<n_y_stiff; i++) {
        SerialMesh stiff(mesh.comm());
        stiff_divs[0] = divs[1];
        stiff_divs[1] = &height;
        
        init.init(stiff_divs, stiff, t);
        
        // add the elements and nodes to the panel mesh
        // the y-coordinate for this mesh will be the x-coordinate
        combine_mesh(mesh, stiff, MAST::StiffenedPanel::STIFFENER_Y,
                     divs[0]->div_location(i+1));
    }
}


void
MAST::StiffenedPanel::combine_mesh(UnstructuredMesh& panel,
                                   UnstructuredMesh& stiffener,
                                   MAST::StiffenedPanel::Component c,
                                   Real stiff_offset) {
    panel.reserve_nodes(stiffener.n_nodes());
    panel.reserve_elem(stiffener.n_elem());
    
    MeshBase::const_element_iterator
    el_it = stiffener.elements_begin(),
    el_end = stiffener.elements_end();
    
    std::map<Node*, Node*> old_to_new;
    Node *old_node, *new_node;
    
    for ( ; el_it != el_end; el_it++ ) {
        Elem* old_elem = *el_it,
        *new_elem = panel.add_elem(Elem::build(old_elem->type()).release());
        
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