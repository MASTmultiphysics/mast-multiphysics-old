//
//  UGFlutterSolver.h
//  FESystem
//
//  Created by Manav Bhatia on 7/24/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_ug_flutter_solver_h__
#define __fesystem_ug_flutter_solver_h__

// C++ includes
#include <vector>


// FESystem includes
#include "Disciplines/Aeroelasticity/FrequencyDomainFlutterSolverBase.h"


namespace FESystem
{
    namespace Aeroelasticity
    {
        template <typename ValType>
        class UGFlutterSolver: public FESystem::Aeroelasticity::FrequencyDomainFlutterSolverBase<ValType>
        {
        public:
            /*!
             *    Constructor
             */
            UGFlutterSolver();
            
            virtual ~UGFlutterSolver();
            
            /*!
             *   method to specify the frequency range for flutter solution. Use the other method to specify the range so that the flutter solver will automatically 
             *   select the k_values
             */
            void setReducedFrequency(std::vector<typename RealOperationType(ValType)>& k_vals);
            
            /*!
             *   method to specify the range of k_values for performing and analysis
             */
            void setReducedFrequency(typename RealOperationType(ValType) k_low, typename RealOperationType(ValType) k_up, FESystemUInt n_divs = 100);
            
            /*!
             *    sets the aerodynamic matrix, which should be for the current value of the reduced frequency
             */
            void setAerodynamicMatrix(const FESystem::Numerics::MatrixBase<ValType>& q_mat);
            
            /*!
             *    solves for the flutter root
             */
            virtual FESystem::Aeroelasticity::FlutterSolutionCallBack solve();
            
        protected:
            
            /*!
             *    performs the eigensolution at the latest reduced frequency and performs other operations
             */
            void eigenSolution();
            
            /*!
             *    Boolean to check if the k_ref data is available
             */
            FESystemBoolean if_k_ref_initialized;
            
            /*!
             *    Boolean to check if the k_ref values to use for flutter solution have been provided, or a range is to be used
             */
            FESystemBoolean if_given_k_ref_vals;

            /*!
             *   range of k_ref values for analysis
             */
            std::pair<typename RealOperationType(ValType), typename RealOperationType(ValType)> k_ref_range;

            /*!
             *    number of steps to use for k_refs between the specified range
             */
            FESystemUInt n_k_ref_steps;
            
            /*!
             *   vector of k_ref value for analysis
             */
            std::vector<typename RealOperationType(ValType)> k_ref_vals;
            
            /*!
             *    Current value of the reduced frequency for which the user should update the aerodynamic matrix
             */
            typename RealOperationType(ValType)  current_k_ref;

            /*!
             *    current iteration number for flutter solution
             */
            FESystemUInt current_iter_num;
            
            /*!
             *   pointer to the aerodynamic matrix
             */
            const FESystem::Numerics::MatrixBase<ValType>* aerodynamic_matrix;

            /*!
             *   temporary matrix
             */
            const FESystem::Numerics::MatrixBase<ValType>* tmp_mat;
            
        };
    }
}

#endif // __fesystem_ug_flutter_solver_h__
