//
//  DegreeOfFreedomMap.h
//  FESystem
//
//  Created by Manav Bhatia on 4/6/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_degree_of_freedom_map_h__
#define __fesystem_degree_of_freedom_map_h__

// C++ includes
#include <vector>

// FESystem includes
#include "Base/FESystemTypes.h"
#include "Base/FESystemExceptions.h"


namespace FESystem
{
    // Forward declerations
    namespace Mesh {class MeshBase;}
    namespace Mesh {class ElemBase;}
    namespace Numerics {template <typename ValType> class MatrixBase;}
    namespace Numerics {template <typename ValType> class VectorBase;}
    namespace Numerics {class SparsityPattern;}
    
    namespace Base
    {
        // Forward declerations
        class DegreeOfFreedomObject;
        
        struct Variable
        {
            FESystemUInt system_id;
            FESystemUInt variable_id;
            FESystemUInt order;
            std::string name;
        };
        
        
        class DegreeOfFreedomMap
        {
        public:
            /*!
             *    constructor with \p m as the reference to the Mesh object for which this object handles the dofs
             */
            DegreeOfFreedomMap(FESystem::Mesh::MeshBase& m);
            
            ~DegreeOfFreedomMap();
            
            /*!
             *    Adds a variable to this DegreeOfFreedomMap, with the name \p var_name and with polynomial \p order and returns 
             *    the internal ID Of the variable. Note that the order defines the interpolation order for this variable. This is used
             *    in different ways depending on the type of finite element in use. For Lagrange elements, the order of the polynomial 
             *    is defined by the number of nodes in an element. Thus, the order should be less than or equal to the Lagrange polynomial 
             *    order in use. For a hierarchich basis like the Legendre polynomials, the \p order is used as the initial order of the
             *    polynomial. 
             */
            FESystemUInt addVariable(std::string& var_name, FESystemUInt order);

            /*!
             *    Returns the number of variables in this degree of freedom map
             */ 
            FESystemUInt getNVariables() const;

            
            /*!
             *    Returns the Variable referred to by \p var_id
             */ 
            const Variable& getVariable(FESystemUInt var_id) const;
            

            /*!
             *    Returns the number of degrees of freedom in this DegreeOfFreedomMap
             */
            FESystemUInt getNDofs() const;

            /*!
             *    Returns the number of degrees of freedom associated with the nodes and element 
             */
            FESystemUInt getNDofsForElem(const FESystem::Mesh::ElemBase& elem) const;

            /*!
             *    Returns the SparsityPattern associated with this map
             */
            const FESystem::Numerics::SparsityPattern& getSparsityPattern() const;
            
            /*!
             *    Initialize the degree of freedom distribution, create the conneectivity map, and the sparsity data structure
             *    for the mesh. The user should have added the system and variables, and finished adding nodes and elements to
             *    the mesh data structure. 
             */
            void reinit();
            
            
            /*!
             *    This method assembles the element matrix into the global matrix. The sequence of the local matrix entries should be
             *    the following: node1_var1, node2_var1, ... , nodeN_var1, node1_var2, ..., nodeN_varM. 
             *    
             */ 
            template <typename ValType>
            void addToGlobalMatrix(const FESystem::Mesh::ElemBase& elem, const FESystem::Numerics::MatrixBase<ValType>& elem_mat, FESystem::Numerics::MatrixBase<ValType>& global_mat) const;

            
            /*!
             *    This method assembles the element vector into the global vector. The sequence of the local matrix entries should be
             *    the following: node1_var1, node2_var1, ... , nodeN_var1, node1_var2, ..., nodeN_varM. 
             *    
             */ 
            template <typename ValType>
            void addToGlobalVector(const FESystem::Mesh::ElemBase& elem, const FESystem::Numerics::VectorBase<ValType>& elem_vec, FESystem::Numerics::VectorBase<ValType>& global_vec) const;
            
            /*!
             *    This method extracts the values of dofs for the element from the global vector. The sequence of the local matrix entries should be
             *    the following: node1_var1, node2_var1, ... , nodeN_var1, node1_var2, ..., nodeN_varM. 
             *    
             */ 
            template <typename ValType>
            void getFromGlobalVector(const FESystem::Mesh::ElemBase& elem, const FESystem::Numerics::VectorBase<ValType>& global_vec, FESystem::Numerics::VectorBase<ValType>& elem_vec) const;

        protected:
            
            /*!
             *   Initializes the sparsity pattern
             */
            void initializeSparsityPattern();
            
            /*!
             *   Sparsity pattern correspding to the for which this DOF map is set up
             */
            FESystem::Numerics::SparsityPattern* sparsity_pattern;
                         
            /*!
             *   Mesh object for which this degree of freeom map handles the dofs
             */
            FESystem::Mesh::MeshBase& mesh;
            
            /*!
             *   If the analysis is h-adaptable or not
             */  
            FESystemBoolean if_h_adaptable;

            
            /*!
             *   If the analysis is p-adaptable or not
             */  
            FESystemBoolean if_p_adaptable;

            /*!
             *   Number of degrees of freedom
             */  
            FESystemUInt n_dofs;
            
            /*!
             *   Vector storage of the variables
             */
            std::vector<FESystem::Base::Variable> variable_vec;
            
            /*!
             *   vector of degree of freedom objects for each global ID.
             */
            std::vector<FESystem::Base::DegreeOfFreedomObject*> dof_object_vector;
            
        };
        
        /*!
         *    Exception is thrown if the new variable name already exists in the system
         */
        DeclareException1(VariableNameAlreadyUsed, 
						  std::string, 
						  << "Variable with name already exists : " << Arg1);

    }
}



#endif //__fesystem_degree_of_freedom_map_h__
