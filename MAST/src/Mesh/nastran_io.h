//
//  nastran_io.h
//  MAST
//
//  Created by Manav Bhatia on 12/16/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_nastran_io_h__
#define __MAST_nastran_io_h__


// libMesh includes
#include "libmesh/mesh_input.h"
#include "libmesh/mesh_output.h"


// MAST includes
#include "StructuralElems/structural_system_assembly.h"
#include "PropertyCards/element_property_card_1D.h"
#include "PropertyCards/element_property_card_2D.h"



namespace MAST {
    class NastranIO:  public MeshInput<MeshBase>,
    public MeshOutput<MeshBase> {
      
    public:
        
        NastranIO(MAST::StructuralSystemAssembly& assembly):
        MeshInput<MeshBase>(assembly.get_system().get_mesh()),
        MeshOutput<MeshBase>(assembly.get_system().get_mesh()),
        _assembly(assembly)
        { }
        
        virtual ~NastranIO() { }
        
        /**
         * This method implements reading a mesh from a NASTRAN input file
         */
        virtual void read (const std::string&)
        { libmesh_error(); } // to be implemented
        
        /**
         * This method implements writing a mesh to a specified file.
         */
        virtual void write (const std::string&);


    protected:
        
        void write_mesh (std::ostream& out_stream);
      
        
        void _nastran_element(const ElemType t,
                              std::string& nm,
                              unsigned int& n_nodes);
        
        void _write_property_cards(std::ostream& out_stream,
                                   std::set<const MAST::ElementPropertyCardBase*>& pcards);
        

        /*!
         *   reference to the structural system assembly for which the
         *   data is to be read/written
         */
        MAST::StructuralSystemAssembly& _assembly;
    };
}





inline
void
MAST::NastranIO::write(const std::string& name) {

    if (MeshOutput<MeshBase>::mesh().processor_id() == 0)
    {
        // Open the output file stream
        std::ofstream out_stream (name.c_str());
        
        // Make sure it opened correctly
        if (!out_stream.good())
            libmesh_file_error(name.c_str());
        
        this->write_mesh (out_stream);
    }
}





inline
void
MAST::NastranIO::write_mesh (std::ostream& out_stream)
{
    // Be sure that the stream is valid.
    libmesh_assert (out_stream.good());
    
    // Get a const reference to the mesh
    const MeshBase& mesh = MeshOutput<MeshBase>::mesh();
    
    // Note: we are using version 2.0 of the gmsh output format.
    
    {
        unsigned int sol_num = 0;
        switch (_assembly.analysis_type()) {
            case MAST::STATIC: {
                sol_num = 101;
            }
                break;
                
            case MAST::MODAL: {
                sol_num = 103;
            }
                break;
                
            case MAST::BUCKLING: {
                sol_num = 105;
            }
                break;
                
            default:
                libmesh_error();
        }
        
        // Write the file header.
        out_stream
        << "SOL " << sol_num << std::endl
        << "CEND" << std::endl
        << "POST TOFILE 12 DISPLACE" << std::endl
        << "SPC = 1" << std::endl
        << "DISPLACEMENT = all" << std::endl
        << "METHOD = 1" << std::endl
        << "BEGIN BULK" << std::endl
        << "PARAM, POST, 0" << std::endl
        //<< "PARAM, LMODES, 20" << std::endl
        << "EIGRL,1,0.0,,20" << std::endl
        << "$" << std::endl;
    }
    
    {
        // write the nodes in (n x y z) format
        for (unsigned int v=0; v<mesh.n_nodes(); v++)
            out_stream
            << std::setw(8)  << "GRID*   "
            << std::setw(16) << mesh.node(v).id()+1
            << std::setw(16) << " "
            << std::setw(16) << std::showpoint << mesh.node(v)(0)
            << std::setw(16) << std::showpoint << mesh.node(v)(1)
            << std::setw(8)  << "      ++" << std::endl
            << std::setw(8)  << "*++     "
            << std::setw(16) << std::showpoint << mesh.node(v)(2)
            << std::endl;
    }
    
    std::set<const MAST::ElementPropertyCardBase*> pcards;
    
    {
        // write the connectivity
        MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
        const MeshBase::const_element_iterator end = mesh.active_elements_end();
        
        // loop over the elements
        for ( ; it != end; ++it)
        {
            const Elem* elem = *it;
            
            // consult the export element table
            std::string nm;
            unsigned int n_nastran_nodes;
            _nastran_element(elem->type(), nm, n_nastran_nodes);
            nm += "*";
            
            
            const MAST::ElementPropertyCardBase& prop = _assembly.get_property_card(*elem);
            
            unsigned int count=0, pid = prop.id();
            
            // this property card needs to be written
            pcards.insert(&prop);
            
            out_stream
            << std::setw(8) << nm
            << std::setw(16) << elem->id()+1
            << std::setw(16) << pid+1;
            // now write the nodes
            
            count = 2; // two entries have already been written
            for (unsigned int i=0; i<n_nastran_nodes; i++) {
                out_stream
                << std::setw(16) << elem->node(i)+1;
                count++;
                if (count == 4) {
                    out_stream
                    << std::setw(8)  << "      ++" << std::endl
                    << std::setw(8)  << "*++     ";
                    count = 0;
                }
            }
            out_stream << std::endl;
        } // element loop
    }
    
    
    // write the element property and material cards
    _write_property_cards(out_stream, pcards);
    
    
    
    {
        // iterate over the loads in the assembly object and write the
        // displacement boundary conditions

        // first prepare a map of boundary ids and the constrained vars on that
        // boundary
        std::map<boundary_id_type, std::vector<unsigned int> >  constrained_vars_map;
        
        // now populate the map for all the give boundaries
        std::multimap<boundary_id_type, MAST::BoundaryCondition*>::const_iterator
        it = _assembly.side_loads().begin(), end = _assembly.side_loads().end();
        
        for ( ; it != end; it++)
            if (it->second->type() == MAST::DISPLACEMENT_DIRICHLET) {
                // get the displacement dirichlet condition
                DirichletBoundary& dirichlet_b =
                (dynamic_cast<MAST::DisplacementDirichletBoundaryCondition*>(it->second))->dirichlet_boundary();
                
                constrained_vars_map[it->first] = dirichlet_b.variables;
            }
        
        // iterate over elements to write the displacement boudnary conditions
        MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
        
        for ( ; el != end_el; ++el) {
            const Elem* elem = *el;
            
            // boundary condition is applied only on sides with no neighbors
            // and if the side's boundary id has a boundary condition tag on it
            for (unsigned int s=0; s<elem->n_sides(); s++)
                if ((*el)->neighbor(s) == NULL &&
                    mesh.boundary_info->n_boundary_ids(elem, s)) {
                    
                    std::vector<boundary_id_type> bc_ids = mesh.boundary_info->boundary_ids(elem, s);
                    
                    // the boundary element
                    AutoPtr<Elem> side(elem->build_side(s).release());
                    // variables that are constrainted
                    std::stringstream oss;
                    
                    for (unsigned int i_bid=0; i_bid<bc_ids.size(); i_bid++)
                        if (constrained_vars_map.count(bc_ids[i_bid])) {
                            
                            const std::vector<unsigned int>& vars = constrained_vars_map[bc_ids[i_bid]];
                            for (unsigned int i_var=0; i_var<vars.size(); i_var++)
                                oss << std::setw(1) << vars[i_var]+1;
                        } // end of boundary id loop
                    
                    // write the SPC entry for each
                    for (unsigned int i_nd=0; i_nd<side->n_nodes(); i_nd++)
                        out_stream
                        << std::setw(8) << "SPC     "
                        << std::setw(8) << 1               // SPC set id
                        << std::setw(8) << side->get_node(i_nd)->id()+1
                        << std::setw(8) << oss.str() <<  std::endl;
                } // end of side loop
        } // end of element loop
    }
    
    out_stream
    << "ENDDATA" << std::endl;
}



inline
void
MAST::NastranIO::_nastran_element(const ElemType t,
                                  std::string& nm,
                                  unsigned int& n_nodes) {
    switch (t) {
            
        case EDGE2:
            nm = "CBEAM";
            n_nodes = 2;
            break;
            
        case EDGE3:
            nm = "CBEAM3";
            n_nodes = 2;
            break;

        case TRI3:
            nm = "CTRIA3";
            n_nodes = 3;
            break;

        case TRI6:
            nm = "CTRIA6";
            n_nodes = 6;
            break;

        case QUAD4:
            nm = "CQUAD4";
            n_nodes = 4;
            break;
            
        case QUAD8:
        case QUAD9:
            nm = "CQUAD8";
            n_nodes = 8;
            break;

        case TET4:
            nm = "CTETRA";
            n_nodes = 4;
            break;

        case TET10:
            nm = "CTETRA";
            n_nodes = 10;
            break;

        case HEX8:
            nm = "CHEXA";
            n_nodes = 8;
            break;

        case HEX20:
        case HEX27:
            nm = "CHEXA";
            n_nodes = 20;
            break;
            
        default:
            libmesh_error(); // should not get here
    }
}



inline
void
MAST::NastranIO::_write_property_cards(std::ostream& out_stream,
                                       std::set<const MAST::ElementPropertyCardBase *> &pcards) {
    
    // set of material cards to be written
    std::set<const MAST::MaterialPropertyCardBase*> mcards;
    
    {
        std::set<const MAST::ElementPropertyCardBase*>::const_iterator
        it = pcards.begin(), end = pcards.end();
        
        for ( ; it != end; it++) {
            const MAST::ElementPropertyCardBase& prop = **it;
            libmesh_assert(prop.if_isotropic());
            mcards.insert(&prop.get_material());
            
            // currently only isotropic material sections are written
            switch (prop.dim()) {
                case 1: {
                    const MAST::ElementPropertyCard1D& prop1d =
                    dynamic_cast<const MAST::ElementPropertyCard1D&>(prop);
                    out_stream
                    << std::setw(8)  << "PBEAM*  "
                    << std::setw(16) << prop1d.id() + 1
                    << std::setw(16) << prop1d.get_material().id()+1
                    << std::setw(16) << prop1d.value("A")
                    << std::setw(16) << prop1d.value("IZZ")
                    << std::setw(8)  << "      ++" << std::endl
                    << std::setw(8)  << "*++     "
                    << std::setw(16) << prop1d.value("IYY")
                    << std::setw(16) << prop1d.value("IYZ")
                    << std::setw(16) << prop1d.value("J")
                    << std::endl;
                }
                    break;
                    
                case 2: {
                    out_stream
                    << std::setw(8)  << "PSHELL* "
                    << std::setw(16) << prop.id() + 1
                    << std::setw(16) << prop.get_material().id()+1
                    << std::setw(16) << prop.get<Real>("h")()
                    << std::setw(16) << prop.get_material().id()+1
                    << std::setw(8)  << "      ++" << std::endl
                    << std::setw(8)  << "*++     "
                    << std::setw(16) << " "
                    << std::setw(16) << prop.get_material().id()+1;
                    
                    if (prop.get_material().contains("kappa"))
                        out_stream
                        << std::setw(16) << prop.get_material().get<Real>("kappa")()
                        << std::endl;
                    else
                        out_stream
                        << std::setw(16) << 0. << std::endl;
                    
                }
                    break;
                    
                case 3: {
                    out_stream
                    << std::setw(8)  << "PSOLID* "
                    << std::setw(16) << prop.id() + 1
                    << std::setw(16) << prop.get_material().id()+1
                    << std::endl;
                }
                    break;
                    
                default:
                    // should not get here
                    libmesh_error();
                    
            }
        }
    }
    
    {
        std::set<const MAST::MaterialPropertyCardBase*>::const_iterator
        it = mcards.begin(), end = mcards.end();
        
        for ( ; it != end; it++) {
            const MAST::MaterialPropertyCardBase& prop = **it;
            
            Real E = prop.get<Real>("E")(),
            nu = prop.get<Real>("nu")(),
            G = E/2./(1.+nu);

            out_stream
            << std::setw(8)  << "MAT1*   "
            << std::setw(16) << prop.id() + 1
            << std::setw(16) << E
            << std::setw(16) << G
            << std::setw(16) << nu
            << std::setw(8)  << "      ++" << std::endl
            << std::setw(8)  << "*++     "
            << std::setw(16) << prop.get<Real>("rho")()
            << std::endl;
        }
    }
}




#endif
