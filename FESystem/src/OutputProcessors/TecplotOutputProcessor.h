//
//  TecplotOutputProcessor.h
//  FESystem
//
//  Created by Manav Bhatia on 4/18/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_tecplot_output_processor_h__
#define __fesystem_tecplot_output_processor_h__

// FESystem includes
#include "OutputProcessors/OutputProcessorBase.h"
#include "Mesh/ElementType.h"

namespace FESystem
{
    namespace OutputProcessor
    {
        /*!
         *   Inherits from the OutputProcessorBase class to implement methods used for output of data to 
         *   Gmsh format
         */
        class TecplotOutputProcessor: public FESystem::OutputProcessor::OutputProcessorBase
        {
        public:
            /*!
             *    constructor
             */
            TecplotOutputProcessor();
            
            virtual ~TecplotOutputProcessor();
            
            /*!
             *   Write \p mesh to the output stream \p output in a Gmsh readable format
             */
            virtual void writeMesh(std::ostream& output, const FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map);
            
            /*!
             *   Writes the variables specified in \p variables_to_write to the output stream. The variables should be 
             *   corresponding to the \p mesh and associated \p dof_map. The values of which are stored in the vector \p vec. 
             */
            virtual void writeSolution(std::ostream& output, 
                                       const std::string& data_name,
                                       const FESystem::Mesh::MeshBase& mesh, 
                                       const FESystem::Base::DegreeOfFreedomMap& dof_map,
                                       const std::vector<FESystemUInt>& variables_to_write,
                                       const FESystem::Numerics::VectorBase<FESystemDouble>& vec);
            
            
        protected:
            
            /*!
             *   value of previous variable zone, used for reference
             */
            FESystemUInt previous_var_zone_number;
            
            /*!
             *   Returns the preferred zone name and element type name for the zone based on the element type for which data is to be written
             */
            void getFiniteElementZoneName(FESystem::Mesh::ElementType elem_type,  std::string& zone_name, std::string& elem_name);
            
            /*!
             *   Populates the vector with the types of elements used by GMSH
             */
            void populateElemTypeVector(const FESystem::Mesh::MeshBase& mesh, std::vector<FESystem::Mesh::ElementType>& elem_type_vec);
            
            /*!
             *    Writes the header for a view with the specified elem_type_vec used by GMSH, and the number of 
             *    scalar, vector and tensor data types. The vector is 3-D and tensor has nine components. 
             */
            void writeViewHeader(std::ostream& output,
                                 const std::vector<FESystem::Mesh::ElementType>& elem_type_vec,
                                 const FESystem::Mesh::MeshBase& mesh,
                                 FESystemUInt n_scalars,
                                 FESystemUInt n_vectors,
                                 FESystemUInt n_tensors);
            
            /*!
             *    Writes the solution to the output stream. If the solution is to be written (specified by the flag), then 
             *    the three pointer data should be provided. Otherwise, no data need be provided, and a scalar of 0 value is automatically
             *    written. 
             */
            void writeSolutionVector(std::ostream& output,
                                     const std::vector<FESystem::Mesh::ElementType>& elem_type_vec,
                                     const FESystem::Mesh::MeshBase& mesh,
                                     FESystemBoolean if_write_solution,
                                     const FESystem::Base::DegreeOfFreedomMap* dof_map = NULL,
                                     const std::vector<FESystemUInt>* variables_to_write = NULL,
                                     const FESystem::Numerics::VectorBase<FESystemDouble>* vec = NULL);
        };
    }
    
}



#endif  //__fesystem_tecplot_output_processor_h__

