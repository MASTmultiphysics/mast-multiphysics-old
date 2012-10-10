//
//  GmshOutputProcessor.h
//  FESystem
//
//  Created by Manav Bhatia on 4/18/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_gmsh_output_processor_h__
#define __fesystem_gmsh_output_processor_h__

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
        class GmshOutputProcessor: public FESystem::OutputProcessor::OutputProcessorBase
        {
        public:
            /*!
             *    constructor
             */
            GmshOutputProcessor();
            
            virtual ~GmshOutputProcessor();
            
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
             *   Returns the Gmsh cell type number for the given FESystem element \p type. 
             */
            void getGmshElemTypeNum(FESystem::Mesh::ElementType type, FESystemUInt& elem_type_num, FESystemUInt& n_nodes_to_write);
            
        };
    }
    
}



#endif  //__fesystem_gmsh_output_processor_h__
