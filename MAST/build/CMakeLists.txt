project (MAST)

cmake_minimum_required(VERSION 2.8)

enable_language (CXX Fortran)

set (MAST_VERSION_MAJOR 0)
set (MAST_VERSION_MINOR 9)

set (libmesh_dir       "libmesh-dir" CACHE STRING "Directory containing libMesh header and library")
set (vtk_include_dir   "vtk-include-dir" CACHE STRING "Directory containing VTK header files")
set (vtk_lib_dir       "vtk-lib-dir" CACHE STRING "Directory containing  VTK library files" )
set (vtk_lib_file      "vtkIO vtkCommon vtkFiltering" CACHE STRING "VTK library files" )
separate_arguments (vtk_lib_file)
set (tbb_include_dir   "tbb-include-dir" CACHE STRING "Directory containing TBB header files")
set (tbb_lib_dir       "tbb-lib-dir" CACHE STRING "Directory containing TBB library files")
set (tbb_lib_file      "tbb tbbmalloc" CACHE STRING "TBB library files")
separate_arguments (tbb_lib_file)
set (boost_include_dir "boost-include-dir" CACHE STRING "Directory containing boost header files")
set (mpi_include_dir   "mpi-include-dir" CACHE STRING "Directory containing MPI header files")
set (mpi_lib_dir       "mpi-lib-dir" CACHE STRING "Directory containing MPI library files")
set (lapack_lib_file   "lapack" CACHE STRING "LAPACK library")
separate_arguments (lapack_lib_file)
set (blas_lib_file     "blas" CACHE STRING "BLAS library")
separate_arguments (blas_lib_file)
set (lapack_lib_dir    "lapack-lib-dir" CACHE STRING "Directory containing LAPACK library")
set (blas_lib_dir      "blas-lib-dir" CACHE STRING "Directory containing BLAS library")
set (CMAKE_Fortran_COMPILER "fortran-compiler" CACHE STRING "Fortran compiler")


file (GLOB gcmma_source_files ../src/Optimization/gcmma/*.f)
add_library (gcmma ${gcmma_source_files})


file (GLOB_RECURSE mast_source_files ../src/*.cpp)
file (GLOB_RECURSE mast_header_files ../src/*.h)



#set (EXTRA_LIBS ${EXTRA_LIBS} MathFunctions)


include_directories (../src 
		    ${vtk_include_dir}
		    ${tbb_include_dir}
		    ${boost_include_dir})

include_directories (SYSTEM ${mpi_include_dir})

add_executable (mast-real ${mast_source_files} ${mast_header_files})
add_executable (mast-complex ${mast_source_files} ${mast_header_files})

set_property (TARGET mast-real APPEND 
	     PROPERTY INCLUDE_DIRECTORIES 
	     ${libmesh_dir}/real/include) 

set_property (TARGET mast-complex APPEND 
	     PROPERTY INCLUDE_DIRECTORIES 
	     ${libmesh_dir}/complex/include) 


find_library (mesh_dbg_real mesh_dbg ${libmesh_dir}/real/lib)
find_library (mesh_opt_real mesh_opt ${libmesh_dir}/real/lib)
find_library (mesh_dbg_complex mesh_dbg ${libmesh_dir}/complex/lib)
find_library (mesh_opt_complex mesh_opt ${libmesh_dir}/complex/lib)

find_library (mpicxx_lib    mpi_cxx              ${mpi_lib_dir})
find_library (mpi_lib       mpi                  ${mpi_lib_dir})

foreach (f ${tbb_lib_file})
   find_library (tbb_lib    ${f}   ${tbb_lib_dir})
endforeach (f)

foreach (f ${vtk_lib_file})
   find_library (vtk_lib    ${f}   ${vtk_lib_dir})
endforeach (f)

foreach (f ${lapack_lib_file})
   find_library (lapack_lib    ${f}   ${lapack_lib_dir})
endforeach (f)

foreach (f ${blas_lib_file})
find_library (blas_lib      ${f}     ${blas_lib_dir})
endforeach (f)


target_link_libraries (mast-real  
		      debug ${mesh_dbg_real}
		      optimized ${mesh_opt_real}
		      gcmma
		      ${lapack_lib}
		      ${blas_lib}
		      ${vtk_lib}
		      ${tbb_lib}
		      ${mpicxx_lib}
		      ${mpi_lib})

target_link_libraries (mast-complex  
		      debug ${mesh_dbg_complex}
		      optimized ${mesh_opt_complex}
		      gcmma
		      ${lapack_lib}
		      ${blas_lib}
	              ${vtk_lib}
	              ${tbb_lib}
		      ${mpicxx_lib}
		      ${mpi_lib})

