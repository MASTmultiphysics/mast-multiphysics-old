////
////  RadiationCavityAnalysis.h
////  FESystem
////
////  Created by Manav Bhatia on 1/25/13.
////
////
//
//#ifndef __FESystem__RadiationCavityAnalysis__
//#define __FESystem__RadiationCavityAnalysis__
//
////C++ includes
//#include <memory>
//#include <vector>
//
//// FESystem includes
//
//
//
//namespace FESystem
//{
//    namespace Radiation
//    {
//        
//        class RadiationCavityAnalysis
//        {
//        public:
//            /// enum for the matrix names for this analysis
//            enum RadiationQty
//            {
//                F_MATRIX,
//                F_MATRIX_SENSITIVITY,
//                A_MATRIX,
//                A_INVERSE_MATRIX,
//                A_MATRIX_SENSITIVITY,
//                B_MATRIX,
//                B_MATRIX_SENSITIVITY,
//                B_INVERSE_MATRIX,
//                G_MATRIX,
//                G_MATRIX_SENSITIVITY,
//                FACTOR_MATRIX,
//                AB_INV_MATRIX
//            };
//            
//            /*!
//             *   constructor
//             */
//            RadiationCavityAnalysis(FESystem::FESystemController& fesys_controller,
//                                    Discipline::ThermalAnalysis& thermal_anal,
//                                    RadiationCavity& rad_cavity);
//            
//            /*!
//             *    destructor
//             */
//            ~RadiationCavityAnalysis();
//            
//            /*!
//             *    method returns a reference to the quantity for this cavity analysis
//             *    @param qty_name name of the quantity that is to be returned
//             *    @param DV pointer to the design variable if sensitivity is seeked. Should
//             *    NULL if no sensitivity is required
//             */
//            void getMatrixQuantity(FESystem::Numerics::MatrixBase<FESystemDouble>& matrix,
//                                   const RadiationQty qty_name);
//            
//            
//            /*!   method returns the FE nodal dof IDs for the nodes of the FE
//             *    that are participating in radiation. Values in the nodal load
//             *    vector will be arranged in the same way as the dofs in this
//             *    vector
//             *    @param vector the IDs will be returned in this vector
//             */
//            void getFENodalDofIDs(std::vector<unsigned int>& vector);
//            
//            /*!
//             *    method returns the FE nodal heat load vector
//             *    @param nodal_temp_val value of the nodal temperature arranged in the same way
//             *    as the dof vector obtrained from the method getFENodalDofIDs()
//             *    @param nodal_load_vector the nodal heat load vector is returned in this
//             */
//            void getFENodalLoadVector(FESystem::Numerics::VectorBase<FESystemDouble>& nodal_temp_val,
//                                      FESystem::Numerics::VectorBase<FESystemDouble>& nodal_load_vector);
//            
//            
//            /*!
//             *    method returns the FE nodal heat load vector
//             *    @param nodal_temp_val value of the nodal temperature arranged in the same way
//             *    as the dof vector obtrained from the method getFENodalDofIDs()
//             *    @param nodal_load_vector the nodal heat load vector is returned in this
//             *    @param DV pointer to design variable for which sensitivity is sought
//             */
//            void getFENodalLoadVectorSensitivity(FESystem::Numerics::VectorBase<FESystemDouble>& nodal_temp_val,
//                                                 FESystem::Numerics::VectorBase<FESystemDouble>& nodal_load_vector);
//            
//            
//            /*!
//             *    method returns the FE jacobian matrix
//             *    @param nodal_temp_val value of the nodal temperature arranged in the same way
//             *    as the dof vector obtrained from the method getFENodalDofIDs()
//             *    @param jacobian the jacobian is returned in this
//             */
//            void getFEJacobianMatrix(FESystem::Numerics::VectorBase<FESystemDouble>& nodal_temp_val,
//                                     FESystem::Numerics::MatrixBase<FESystemDouble>& jacobian);
//            
//        protected:
//            
//            /*!   
//             *   calculates the vector of emissivities for the current temperature of the cavity
//             */
//            void calculateEmissivityVector(FESystem::Numerics::VectorBase<FESystemDouble>* vector_ptr);
//            
//            /*!
//             *    calculates the matrix product \f$ [A][B]^{-1} \f$, and approximates it based on the
//             *    settings provided by the user about whether to approximate or to always recalculate the
//             *    factor matrix
//             */
//            void calculateABInvFactorMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>* matrix_ptr);
//            
//            /*!
//             *    only retrieves the base epsilon vector stored in the database. This method
//             *    does not make sure that the calculated quantity is current. Hence, it should
//             *    be ensured in the appropriate context
//             */
//            void getBaseEpsilonVector(FESystem::Numerics::VectorBase<FESystemDouble>& vec);
//            
//            /*!
//             *    calculates the shape factor matrix
//             *    @param matrix_ptr pointer to matrix where the calculated matrix is stored
//             *    @param DV pointer to the design variable for sensitivity. (should NULL if
//             *    no sensitivity is needed)
//             */
//            void calculateShapeFactorMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>* matrix_ptr,
//                                            const bool sensitivity);
//            
//            /*!
//             *    calculates the A matrix
//             *    @param matrix_ptr pointer to matrix where the calculated matrix is stored
//             *    @param DV pointer to the design variable for sensitivity. (should NULL if
//             *    no sensitivity is needed)
//             */
//            void calculateA_Matrix(FESystem::Numerics::MatrixBase<FESystemDouble>* matrix_ptr,
//                                   const bool sensitivity);
//            
//            /*!
//             *    calculates inverse of A matrix
//             *    @param matrix_ptr pointer to matrix where the calculated matrix is stored
//             */
//            void calculateA_InverseMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>* matrix_ptr);
//            
//            /*!
//             *    calculates B matrix
//             *    @param matrix_ptr pointer to matrix where the calculated matrix is stored
//             *    @param DV pointer to the design variable for sensitivity. (should NULL if
//             *    no sensitivity is needed)
//             */
//            void calculateB_Matrix(FESystem::Numerics::MatrixBase<FESystemDouble>* matrix_ptr,
//                                   const bool sensitivity);
//            
//            
//            /*!
//             *    calculates inverse of B matrix
//             *    @param matrix_ptr pointer to matrix where the calculated matrix is stored
//             */
//            void calculateB_InverseMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>* matrix_ptr);
//            
//            
//            /*!
//             *   calculate the matrix \f$ \frac{\partial B_{ij}}{\partial T_{rad_k}} \lambda_j \f$
//             */
//            void calculateB_TemperatureDependentJacobianContribution
//            (FESystem::Numerics::VectorBase<FESystemDouble>& vector,
//             const FESystem::Numerics::VectorBase<FESystemDouble>& rad_temp_vector,
//             FESystem::Numerics::MatrixBase<FESystemDouble>& A_inv_mat,
//             FESystem::Numerics::MatrixBase<FESystemDouble>& AB_inv_mat);
//            
//            
//            /*!
//             *   calculates the vector \f$ \sum_{j} [((AB^{-1})_{ii})^{j-1} (\gamma_i)^j ]_{diag} \f$ .
//             *   This is needed for the case of temperature dependent material property when the
//             *    \f$ [A][B]^{-1} \f$ is being approximated
//             */
//            void getTemperatureDependentFluxCorrectionVector
//            (FESystem::Numerics::VectorBase<FESystemDouble>& correction_vec,
//             FESystem::Numerics::MatrixBase<FESystemDouble>& AB_inv_mat);
//            
//            /*!
//             *    calculated G matrix
//             *    @param matrix_ptr pointer to matrix where the calculated matrix is stored
//             */
//            void calculateG_Matrix(FESystem::Numerics::MatrixBase<FESystemDouble>* matrix_ptr,
//                                   const bool sensitivity);
//            
//            /*! 
//             *    calculate factor matrix
//             *    @param matrix_ptr pointer to matrix where the calculated matrix is stored
//             */
//            void calculateFactorMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>* matrix_ptr);
//            
//            
//            /*!
//             *    calculates the radiation temperature vector from the FE nodal load vector
//             *    @param rad_temp_vector_ptr radiation vector
//             */
//            void calculateRadElemTempVector(FESystem::Numerics::VectorBase<FESystemDouble>* rad_temp_vector_ptr,
//                                            FESystem::Numerics::VectorBase<FESystemDouble>& nodal_temp_vector,
//                                            const bool reduce_power_by_four,
//                                            const bool sensitivity);
//            
//            
//            /*! 
//             *    initializes the property card at the specified temperature
//             */
//            void initializePropertyCardAtElementTemperature(ElemDataCard& elem_data_card,
//                                                            const unsigned int rad_elem_ID);
//            
//            
//            /*!
//             *    calculates and returns the quantity in the matrix that is passed as the
//             *    arguement. This method is merely an interface that calls the respective
//             *    method based on the quantity that is being seeked
//             *    @param matrix pointer to matrix where the quantity will be stored
//             *    @param qty name of the quantity that is to be calculated
//             *    @param DV pointer to a design variable for which sensitivity is being seeked.
//             *    If no sensitivity is required, DV should be NULL
//             */
//            void calculateQty(FESystem::Numerics::MatrixBase<FESystemDouble>* matrix,
//                              const RadiationQty qty);
//            
//            
//            /*!
//             *   @returns true if the quantity is temperature dependent in case of temperature dependent
//             *    properties.
//             */
//            bool checkIfTemperatureDependent(RadiationQty quantity);
//            
//            /*!
//             *   clears the property card initialization
//             */
//            void clearPropertyCardInitialization(ElemDataCard& elem_data_card);
//            
//            
//            /*!
//             *    @returns a DataInfoBase object for the given quantity, for the current analysis in progress.
//             *    The object will point to the current time iteration for a transient analysis, and the
//             *    current DV for a sensitivity quantity
//             */
//            std::auto_ptr<FESystemDatabase::DataInfoBase> getDataInfoForQty(RadiationQty qty);
//            
//            /*!
//             *   fesystem controller for this analysis
//             */
//            FESystem::FESystemController& fesystem_controller;
//            
//            /*!
//             *   reference to the thermal analysis discipline that owns this object
//             */
//            Discipline::ThermalAnalysis& thermal_discipline;
//            
//            /*!
//             *   cavity for this analysis
//             */
//            RadiationCavity& radiation_cavity;
//            
//            /*!
//             *   number of times the ABInv matrix has been calculated
//             */
//            unsigned int n_ABinv_calculations;
//            
//            /*!
//             *   pointer to the FE nodal temperature vector for the analysis
//             */
//            FESystem::Numerics::VectorBase<FESystemDouble>* fe_nodal_temp_vector;
//            
//            /*!
//             *   performance logging header
//             */
//            const std::string performance_logging_header;
//        };
//        
//    }
//}
//
//
//
//#endif /* defined(__FESystem__RadiationCavityAnalysis__) */
