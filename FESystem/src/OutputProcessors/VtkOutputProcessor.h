//
//  VtkOutputProcessor.h
//  FESystem
//
//  Created by Manav Bhatia on 4/18/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_vtk_output_processor_h__
#define __fesystem_vtk_output_processor_h__


// FESystem includes
#include "OutputProcessors/OutputProcessorBase.h"
#include "Mesh/ElementType.h"

namespace FESystem
{
    namespace OutputProcessor
    {
        /*!
         *   Inherits from the OutputProcessorBase class to implement methods used for output of data to 
         *   Vtk format
         */
        class VtkOutputProcessor: public FESystem::OutputProcessor::OutputProcessorBase
        {
        public:
            /*!
             *    constructor
             */
            VtkOutputProcessor();
            
            virtual ~VtkOutputProcessor();
            
            /*!
             *   Write \p mesh to the output stream \p output in a Vtk readable format
             */
            virtual void writeMesh(std::ostream& output, const FESystem::Mesh::MeshBase& mesh, const FESystem::DegreeOfFreedom::DegreeOfFreedomMap& dof_map);
            
            /*!
             *   Writes the variables specified in \p variables_to_write to the output stream. The variables should be 
             *   corresponding to the \p mesh and associated \p dof_map. The values of which are stored in the vector \p vec. 
             */
            virtual void writeSolution(std::ostream& output, 
                                       const std::string& data_name,
                                       const FESystem::Mesh::MeshBase& mesh, 
                                       const FESystem::DegreeOfFreedom::DegreeOfFreedomMap& dof_map,
                                       const std::vector<FESystemUInt>& variables_to_write,
                                       const FESystem::Numerics::VectorBase<FESystemDouble>& vec);
            
            
        protected:
            
            /*!
             *   Returns the VTK cell type number for the given FESystem element \p type. Since VTK does not have a quadratic element, 
             *   all elements are written as their linear counterparts, and only the first \p n_nodes_to_write nodes will be written for the cell
             */
            void getVtkElemTypeNum(FESystem::Mesh::ElementType type, FESystemUInt& elem_type_num, FESystemUInt& n_nodes_to_write);
            
            /*!
             *   flag of whether the previous data was a vector(=3) or scalar(=1) or mesh(=0)
             */
            FESystemUInt n_previous_dofs;
        };
    }
    
}

#endif // __fesystem_vtk_output_processor_h__

