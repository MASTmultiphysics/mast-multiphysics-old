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
        // Write the file header.
        out_stream
        << "SOL 145" << std::endl
        << "CEND" << std::endl
        << "SPC = 1" << std::endl
        << "DISPLACEMENT = all" << std::endl
        << "SDISPLACEMENT = all" << std::endl
        << "SVECTOR = all" << std::endl
        << "METHOD = 1" << std::endl
        << "FMETHOD = 1" << std::endl
        << "CMETHOD = 1" << std::endl
        << "BEGIN BULK" << std::endl
        << "PARAM, POST, 0" << std::endl
        << "PARAM, LMODES, 20" << std::endl
        << "EIGRL,1,0.0,,20" << std::endl
        << "EIGC,1,CLAN,,,,,20" << std::endl
        << "$" << std::endl
        << "$TABDMP1 , 1, G," << std::endl
        << "$,0.0,0.02,10.0,0.02,ENDT" << std::endl
        << "$PARAM, KDAMP, -1" << std::endl;
    }
    
    {
        // write the nodes in (n x y z) format
        for (unsigned int v=0; v<mesh.n_nodes(); v++)
            out_stream
            << std::setw(8)  << "GRID*   "
            << std::setw(16) << mesh.node(v).id()+1
            << std::setw(16) << " "
            << std::setw(16) << mesh.node(v)(0)
            << std::setw(16) << mesh.node(v)(1)
            << std::setw(8)  << "   ++   " << std::endl
            << std::setw(8)  << "   ++   "
            << std::setw(16) << mesh.node(v)(2) << std::endl;
    }
    
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
            
            unsigned int pid, count=0;
            
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
                    << std::setw(8)  << "   ++   " << std::endl
                    << std::setw(8)  << "   ++   ";
                    count = 0;
                }
            }
            out_stream << std::endl;
        } // element loop
    }
    
    
    {
        // iterate over the loads in the assembly object and write the
        // displacement boundary conditions
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
            nm = "CTRI3";
            n_nodes = 3;
            break;

        case TRI6:
            nm = "CTRI6";
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


//{
//    do i=1,n_elems
//        
//    do i=1,n_constrained_nodes
//    write(1,*) "SPC,1,",constrained_nodes(i),",123456,0.0"
//    end do
//    do i=1,n_materials
//    write (1,"(A5,I4,A1,D10.5,A3,F10.5,A3,F10.5)") "MAT1,",i,","
//    $      ,material_prop(i,2),",,",material_prop(i,3),","
//    $      ,material_prop(i,1)
//    end do
//    do i=1,n_elems
//    n1 = elem_nodes(i,1)
//    n2 = elem_nodes(i,2)
//    val = 0.0D0
//    do j =1,3
//    v(j) = node_coords(n2,j) - node_coords(n1,j)
//    val = val + v(j)**2
//    end do
//    val = val ** 0.5
//    do j=1,3
//    v(j) = v(j)/val
//    end do
//    write (1,"(A8,I4,A14,F10.5,A2,F10.5)") "cord2r,",i
//    $         ,",,0.,0.,0.,0.,",-v(3),",",v(2)
//    write (1,*) ",1.,",-v(3),",",v(2)
//    area = 2.0*elem_chord(i)*(elem_tc(i)*elem_skin_thickness(i,2)
//                              $         + elem_structural_to_aero_chord_ratio(i)*
//                              $         elem_skin_thickness(i,1))
//    Itr_l = Itr(elem_skin_thickness(i,1), elem_skin_thickness(i
//                                                              $         ,2), elem_tc(i), elem_chord(i),
//                $         elem_structural_to_aero_chord_ratio(i))
//    Ich_l = Ich(elem_skin_thickness(i,1), elem_skin_thickness(i
//                                                              $         ,2), elem_tc(i), elem_chord(i),
//                $         elem_structural_to_aero_chord_ratio(i))
//    Jfactor = GJfactor(elem_skin_thickness(i,1),
//                       $         elem_skin_thickness(i,2), elem_tc(i), elem_chord(i)
//                       $         ,elem_structural_to_aero_chord_ratio(i), 1.0D0, 1.0D0)
//    write (1,"(A6,I4,A2,I4,A2,F10.5,A3,F10.5,A3,F10.5,A7,F10.5,A4)")
//    $         "pbeam,",i,",",elem_property(i,1),",",area,",",Itr_l,",",
//    $         Ich_l,",0.0,",Jfactor,","
//    write (1,*) ",,"
//    write (1,*) ",0.0,0.0,"
//    write (1,*) ",,"
//    
//    c          write (1,"(A7,I5,A2,I5,A7)") "PBEAML,",i,",",1,",,BOX"
//    c          write (1,"(4(A1,F10.5))") ",",elem_chord(i)
//    c     $         *elem_structural_to_aero_chord_ratio(i),",",elem_tc(i)
//    c     $         *elem_chord(i),",",elem_skin_thickness(i,1),","
//    c     $         ,elem_skin_thickness(i,2)
//    if ((abs(v(1)) .lt. 1.0e-8) .and. (abs(v(1)) .lt. 1.0e-8))
//    $         then
//    v(2) = 1.0D0
//    end if
//    
//    write(1,"(A5,I4,A1,I4,A1,I4,A2,I4,A2,F10.5,A2,F10.5,A2F10.5)")
//    $         "CBEAM,",i,",",i,",",n1,",",n2,",",0.0D0,",",-v(3),",",
//    $         -v(2)
//    
//}




#endif
