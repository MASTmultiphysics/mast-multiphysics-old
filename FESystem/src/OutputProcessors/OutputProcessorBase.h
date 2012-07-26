//
//  OutputProcessorBase.h
//  FESystem
//
//  Created by Manav Bhatia on 4/18/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_output_processor_base_h__
#define __fesystem_output_processor_base_h__

// C++ includes
#include <ostream>
#include <vector>
#include <string>

// FESystem includes
#include "Base/FESystemTypes.h"


namespace FESystem
{
    // Forward declerations
    namespace Mesh {class MeshBase;}
    namespace Base {class DegreeOfFreedomMap;}
    namespace Numerics {template <typename ValType> class VectorBase;}
    
    namespace OutputProcessor
    {
        /*!
         *   Provides an abstract class to write data to output stream. 
         */
        class OutputProcessorBase
        {
        public:
                        
            /*!
             *   Provides the interface to write the mesh to the output stream
             */
            /*!
             *   Write \p mesh to the output stream \p output stream
             */
            virtual void writeMesh(std::ostream& output, const FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map)=0;
            
            /*!
             *   Writes the variables specified in \p variables_to_write to the output stream. The variables should be 
             *   corresponding to the \p mesh and associated \p dof_map. The values of which are stored in the vector \p vec. 
             */
            virtual void writeSolution(std::ostream& output, 
                                       const std::string& data_name,
                                       const FESystem::Mesh::MeshBase& mesh, 
                                       const FESystem::Base::DegreeOfFreedomMap& dof_map,
                                       const std::vector<FESystemUInt>& variables_to_write,
                                       const FESystem::Numerics::VectorBase<FESystemDouble>& vec)=0;
            
        protected:
            
        };
    }
}


#endif  // __fesystem_output_processor_base_h__
