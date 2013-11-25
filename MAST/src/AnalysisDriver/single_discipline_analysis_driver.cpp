//
//  analysis_driver.cpp
//  MAST
//
//  Created by Manav Bhatia on 11/25/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

// MAST includes
#include "AnalysisDriver/analysis_driver.h"




MAST::AnalysisDriver::AnalysisDriver() {
    
}



MAST::AnalysisDriver::~AnalysisDriver() {
    
}



void
MAST::AnalysisDriver::init(LibMeshInit& init, GetPot& infile,
                           int argc, char* const argv[]) {
    
    
}


void
MAST::AnalysisDriver::init_geometry() {
    if (if_panel_mesh)
    {
        const unsigned int nx_divs = infile("nx_divs",0),
        ny_divs = infile("ny_divs",0),
        nz_divs = infile("nz_divs",0);
        const Real t_by_c =  infile("t_by_c", 0.0);
        ElemType elem_type =
        Utility::string_to_enum<ElemType>(infile("elem_type", "QUAD4"));
        
        std::vector<Real> x_div_loc(nx_divs+1), x_relative_dx(nx_divs+1),
        y_div_loc(ny_divs+1), y_relative_dx(ny_divs+1),
        z_div_loc(nz_divs+1), z_relative_dx(nz_divs+1);
        std::vector<unsigned int> x_divs(nx_divs), y_divs(ny_divs), z_divs(nz_divs);
        std::auto_ptr<MeshInitializer::CoordinateDivisions>
        x_coord_divs (new MeshInitializer::CoordinateDivisions),
        y_coord_divs (new MeshInitializer::CoordinateDivisions),
        z_coord_divs (new MeshInitializer::CoordinateDivisions);
        std::vector<MeshInitializer::CoordinateDivisions*> divs(dim);
        
        // now read in the values: x-coord
        if (nx_divs > 0)
        {
            for (unsigned int i_div=0; i_div<nx_divs+1; i_div++)
            {
                x_div_loc[i_div]     = infile("x_div_loc", 0., i_div);
                x_relative_dx[i_div] = infile( "x_rel_dx", 0., i_div);
                if (i_div < nx_divs) //  this is only till nx_divs
                    x_divs[i_div] = infile( "x_div_nelem", 0, i_div);
            }
            divs[0] = x_coord_divs.get();
            x_coord_divs->init(nx_divs, x_div_loc, x_relative_dx, x_divs);
        }
        
        // now read in the values: y-coord
        if ((dim > 1) && (ny_divs > 0))
        {
            for (unsigned int i_div=0; i_div<ny_divs+1; i_div++)
            {
                y_div_loc[i_div]     = infile("y_div_loc", 0., i_div);
                y_relative_dx[i_div] = infile( "y_rel_dx", 0., i_div);
                if (i_div < ny_divs) //  this is only till ny_divs
                    y_divs[i_div] = infile( "y_div_nelem", 0, i_div);
            }
            divs[1] = y_coord_divs.get();
            y_coord_divs->init(ny_divs, y_div_loc, y_relative_dx, y_divs);
        }
        
        // now read in the values: z-coord
        if ((dim == 3) && (nz_divs > 0))
        {
            for (unsigned int i_div=0; i_div<nz_divs+1; i_div++)
            {
                z_div_loc[i_div]     = infile("z_div_loc", 0., i_div);
                z_relative_dx[i_div] = infile( "z_rel_dx", 0., i_div);
                if (i_div < nz_divs) //  this is only till nz_divs
                    z_divs[i_div] = infile( "z_div_nelem", 0, i_div);
            }
            divs[2] = z_coord_divs.get();
            z_coord_divs->init(nz_divs, z_div_loc, z_relative_dx, z_divs);
        }
        
        const std::string mesh_type = infile("mesh_type", std::string(""));
        if (mesh_type == "panel")
        {
            const bool if_cos_bump = infile("if_cos_bump", false);
            const unsigned int n_max_bumps_x = infile("n_max_bumps_x", 1),
            n_max_bumps_y = infile("n_max_bumps_x", 1),
            panel_bc_id = infile("panel_bc_id", 10),
            symmetry_bc_id = infile("symmetry_bc_id", 11);
            
            
            if (dim == 2)
                PanelMesh2D().init(t_by_c, if_cos_bump, n_max_bumps_x,
                                   panel_bc_id, symmetry_bc_id,
                                   divs, mesh, elem_type);
            else if (dim == 3)
                PanelMesh3D().init(t_by_c, if_cos_bump,
                                   n_max_bumps_x, n_max_bumps_y,
                                   panel_bc_id, symmetry_bc_id,
                                   divs, mesh, elem_type);
            else
                libmesh_error();
        }
        else if (mesh_type == "gaussian")
        {
            if (dim == 2)
                GaussianBumpMesh2D().init(t_by_c, divs, mesh, elem_type);
            else if (dim == 3)
                GaussianBumpMesh3D().init(t_by_c, divs, mesh, elem_type);
            else
                libmesh_error();
        }
        else if (mesh_type == "ringleb")
            RinglebMesh().init(nx_divs, ny_divs, mesh, elem_type);
        else
            libmesh_error();
    }
    else
    {
        mesh.set_mesh_dimension(dim);
        const std::string
        mesh_input_file = infile("mesh_input", std::string("")),
        mesh_type = infile("mesh_type", std::string(""));
        
        std::auto_ptr<MeshInput<UnstructuredMesh> > mesh_io;
        if (mesh_type == "gmsh")
            GmshIO(mesh).read(mesh_input_file);
        else if (mesh_type == "exodus")
            ExodusII_IO(mesh).read(mesh_input_file);
        else if (mesh_type == "exodus_parallel")
            ExodusII_IO(mesh).read_parallel(mesh_input_file);
        else if (mesh_type == "nemesis")
            Nemesis_IO(mesh).read(mesh_input_file);
        else
            libmesh_error();
        
        mesh.prepare_for_use();
    }
#else
    
    mesh.read("saved_mesh.xdr");
    
#endif // LIBMESH_USE_COMPLEX_NUMBERS
    
    // uniformly refine the mesh
    for (unsigned int i=0; i<n_uniform_refine; i++)
        mesh_refinement.uniformly_refine();
    
    // Print information about the mesh to the screen.
    mesh.print_info();

}


void
MAST::AnalysisDriver::init_solver() {
    
}


void
MAST::AnalysisDriver::perform_amr() {
    
}


void
MAST::AnalysisDriver::output() {
    
}


