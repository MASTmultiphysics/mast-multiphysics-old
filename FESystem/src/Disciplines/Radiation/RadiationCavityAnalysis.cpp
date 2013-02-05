////
////  RadiationCavityAnalysis.cpp
////  FESystem
////
////  Created by Manav Bhatia on 1/25/13.
////
////
//
////C++ includes
//#include <sstream>
//#include <fstream>
//
//// FESystem includes
//#include "Disciplines/Radiation/RadiationCavityAnalysis.h"
//
//
//
//FESystem::Raditaion::RadiationCavityAnalysis::RadiationCavityAnalysis(FESystem::FESystemController& fesys_controller,
//                                                                      Discipline::ThermalAnalysis& thermal_anal,
//                                                                      RadiationCavity& rad_cavity):
//fesystem_controller(fesys_controller),
//thermal_discipline(thermal_anal),
//radiation_cavity(rad_cavity),
//n_ABinv_calculations(0),
//fe_nodal_temp_vector(NULL),
//performance_logging_header("RadiationCavityAnalysis")
//{
//    this->fesystem_controller.performance_logging->setEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::RadiationCavityAnalysis()",
//     this->performance_logging_header);
//    
//    // set the mesh data for the radiation cavity and tell it to initialize its data
//    // structures
//    unsigned int mesh_ID = this->radiation_cavity.getFiniteElementMeshID();
//    
//    MeshDS::FEMeshData* mesh_data =
//    this->fesystem_controller.mesh_list->getMeshDataFromID(mesh_ID);
//    
//    this->radiation_cavity.setMeshData(mesh_data);
//    
//    // now, iterate over all the dsign variables, and set the mesh data for each
//    std::auto_ptr<std::vector<DesignData::DesignParameter*> > dv_vector =
//    this->fesystem_controller.design_database->getParameters();
//    
//    std::vector<DesignData::DesignParameter*>::const_iterator dv_it, dv_end;
//    dv_it = dv_vector->begin();
//    dv_end = dv_vector->end();
//    
//    
//    for (; dv_it != dv_end; dv_it++)
//    {
//        if ((*dv_it)->getParameterTypeEnumID() == DesignData::SHAPE_PARAMETER::num())
//        {
//            DesignData::ShapeParameter* shape_dv = dynamic_cast<DesignData::ShapeParameter*>(*dv_it);
//            mesh_ID = shape_dv->getPerturbedMeshID(Discipline::THERMAL_DISCIPLINE::num());
//            mesh_data =
//            this->fesystem_controller.mesh_list->getMeshDataFromID(mesh_ID);
//            
//            this->radiation_cavity.setMeshData(mesh_data, shape_dv->getID());
//        }
//    }
//    
//    this->fesystem_controller.performance_logging->unsetEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::RadiationCavityAnalysis()",
//     this->performance_logging_header);
//}
//
//
//
//
//
//
//
//FESystem::Raditaion::RadiationCavityAnalysis::~RadiationCavityAnalysis()
//{
//    
//}
//
//
//
//
//
//void
//FESystem::Raditaion::RadiationCavityAnalysis::getMatrixQuantity(FESystem::Numerics::MatrixBase<FESystemDouble>& matrix,
//                                                                RadiationQty qty_name)
//{
//    // check if the matrix exists. If not, then calculate it, store it and
//    // return it.
//    
//    // create a data info object based on the kind of analysis
//    std::auto_ptr<FESystemDatabase::DataInfoBase> data_info(this->getDataInfoForQty(qty_name).release());
//    
//    // this bool will check if this quantity needs to be recalculated and resaved
//    bool recalculate_and_store = false;
//    if (this->thermal_discipline.checkPropertyDependenceOnTemperature() &&
//        this->checkIfTemperatureDependent(qty_name))
//        recalculate_and_store = true;
//    
//    
//    bool qty_exists =
//    this->fesystem_controller.global_data_storage->checkIfQtyExists(*data_info);
//    
//    
//    // if the quantity is not present or if the qty is dependent on temperature
//    // and
//    if ( !qty_exists || recalculate_and_store)
//    {
//        this->calculateQty(&matrix, qty_name);
//        this->fesystem_controller.global_data_storage->storeMatrix(*data_info, matrix);
//        
//        //**************** this is only temporary *********************
//        /*      if (qty_name == FESystem::Raditaion::RadiationCavityAnalysis::F_MATRIX)
//         {
//         FESystem::Numerics::MatrixBase<FESystemDouble>& matrix = *(insert_return_pair.first->second);
//         
//         Mat petsc_matrix;
//         MatCreateSeqAIJ(PETSC_COMM_SELF,
//         matrix.m(),
//         matrix.n(),
//         0,
//         PETSC_NULL,
//         &petsc_matrix);
//         
//         for (unsigned int i=0; i < matrix.m(); i++)
//         for (unsigned int j=0; j < matrix.n(); j++)
//         MatSetValue(petsc_matrix, i,j, matrix(i,j), INSERT_VALUES);
//         
//         MatAssemblyBegin(petsc_matrix, MAT_FINAL_ASSEMBLY);
//         MatAssemblyEnd(petsc_matrix, MAT_FINAL_ASSEMBLY);
//         
//         PetscViewer binary_viewer;
//         
//         std::string file_name = name;
//         std::ostringstream osstream;
//         osstream << this->radiation_cavity.getCavityID();
//         file_name += osstream.str();
//         
//         PetscViewerBinaryOpen( FESystem::COMM_WORLD,
//         file_name.c_str(),
//         PETSC_FILE_CREATE,
//         &binary_viewer);
//         
//         
//         MatView(petsc_matrix, binary_viewer);
//         
//         PetscViewerDestroy(binary_viewer);
//         
//         MatDestroy(petsc_matrix);
//         }*/
//        //**************** this is only temporary *********************
//        //**************** this is only temporary *********************
//        /*          if (qty_name == FESystem::Raditaion::RadiationCavityAnalysis::F_MATRIX )
//         {
//         FESystem::Numerics::MatrixBase<FESystemDouble>& matrix = *(insert_return_pair.first->second);
//         
//         Mat Fmatrix, dFmatrix;
//         
//         PetscViewer binary_viewer;
//         
//         std::string base_name, sens_name;
//         base_name = this->getStringNameForQty(qty_name);
//         sens_name = this->getStringNameForQtySensitivity(qty_name, 1);
//         
//         std::ostringstream osstream;
//         osstream << this->radiation_cavity.getCavityID();
//         
//         base_name += osstream.str();
//         sens_name += osstream.str();
//         
//         // Read in F
//         PetscViewerBinaryOpen( FESystem::COMM_WORLD,
//         base_name.c_str(),
//         PETSC_FILE_RDONLY,
//         &binary_viewer);
//         MatLoad(binary_viewer, MATSEQAIJ, &Fmatrix);
//         PetscViewerDestroy(binary_viewer);
//         MatAssemblyBegin(Fmatrix, MAT_FINAL_ASSEMBLY);
//         MatAssemblyEnd(Fmatrix, MAT_FINAL_ASSEMBLY);
//         
//         // now read in dF/dalpha
//         PetscViewerBinaryOpen( FESystem::COMM_WORLD,
//         sens_name.c_str(),
//         PETSC_FILE_RDONLY,
//         &binary_viewer);
//         MatLoad(binary_viewer, MATSEQAIJ, &dFmatrix);
//         PetscViewerDestroy(binary_viewer);
//         MatAssemblyBegin(dFmatrix, MAT_FINAL_ASSEMBLY);
//         MatAssemblyEnd(dFmatrix, MAT_FINAL_ASSEMBLY);
//         
//         
//         double ddv = this->fesystem_controller.analysis_case->getRealParameter("DDV");
//         
//         MatAXPY(Fmatrix , ddv, dFmatrix, SAME_NONZERO_PATTERN);
//         MatAssemblyBegin(Fmatrix, MAT_FINAL_ASSEMBLY);
//         MatAssemblyEnd(Fmatrix, MAT_FINAL_ASSEMBLY);
//         
//         
//         double value;
//         for (unsigned int i=0; i < matrix.m(); i++)
//         for (unsigned int j=0; j < matrix.n(); j++)
//         {
//         const int  ii  = i;
//         const int jj = j;
//         if (DV == NULL)
//         MatGetValues(Fmatrix, 1, &ii, 1, &jj, &value);
//         else
//         MatGetValues(dFmatrix, 1, &ii, 1, &jj, &value);
//         
//         matrix(i,j) = value;
//         }
//         
//         MatDestroy(Fmatrix);
//         MatDestroy(dFmatrix);
//         
//         }*/
//        //**************** this is only temporary *********************
//        
//    }
//    else
//        this->fesystem_controller.global_data_storage->fillMatrix(*data_info, matrix);
//}
//
//
//
//
//
//void FESystem::Raditaion::RadiationCavityAnalysis::calculateQty(FESystem::Numerics::MatrixBase<FESystemDouble>* matrix,
//                                                                RadiationQty qty)
//{
//    switch (qty)
//    {
//        case FESystem::Raditaion::RadiationCavityAnalysis::F_MATRIX:
//            this->calculateShapeFactorMatrix(matrix, false);
//            break;
//            
//        case FESystem::Raditaion::RadiationCavityAnalysis::F_MATRIX_SENSITIVITY:
//            this->calculateShapeFactorMatrix(matrix, true);
//            break;
//            
//        case FESystem::Raditaion::RadiationCavityAnalysis::A_MATRIX:
//            this->calculateA_Matrix(matrix, false);
//            break;
//            
//        case FESystem::Raditaion::RadiationCavityAnalysis::A_INVERSE_MATRIX:
//            this->calculateA_InverseMatrix(matrix);
//            break;
//            
//        case FESystem::Raditaion::RadiationCavityAnalysis::A_MATRIX_SENSITIVITY:
//            this->calculateA_Matrix(matrix, true);
//            break;
//            
//        case FESystem::Raditaion::RadiationCavityAnalysis::B_MATRIX:
//            this->calculateB_Matrix(matrix, false);
//            break;
//            
//        case FESystem::Raditaion::RadiationCavityAnalysis::B_MATRIX_SENSITIVITY:
//            this->calculateB_Matrix(matrix, true);
//            break;
//            
//        case FESystem::Raditaion::RadiationCavityAnalysis::B_INVERSE_MATRIX:
//            this->calculateB_InverseMatrix(matrix);
//            break;
//            
//        case FESystem::Raditaion::RadiationCavityAnalysis::G_MATRIX:
//            this->calculateG_Matrix(matrix, false);
//            break;
//            
//        case FESystem::Raditaion::RadiationCavityAnalysis::G_MATRIX_SENSITIVITY:
//            this->calculateG_Matrix(matrix, true);
//            break;
//            
//        case FESystem::Raditaion::RadiationCavityAnalysis::FACTOR_MATRIX:
//            this->calculateFactorMatrix(matrix);
//            break;
//            
//        case FESystem::Raditaion::RadiationCavityAnalysis::AB_INV_MATRIX:
//            this->calculateABInvFactorMatrix(matrix);
//            break;
//            
//        default:
//            abort();
//            break;
//    }
//}
//
//
//
//void FESystem::Raditaion::RadiationCavityAnalysis::calculateShapeFactorMatrix
//(FESystem::Numerics::MatrixBase<FESystemDouble>* matrix_ptr,
// bool sensitivity)
//{
//    this->fesystem_controller.performance_logging->setEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::calculateShapeFactorMatrix()",
//     this->performance_logging_header);
//    
//    FESystem::Numerics::MatrixBase<FESystemDouble>& matrix = *(matrix_ptr);
//    
//    unsigned int n_elems = this->radiation_cavity.nRadiationElems();
//    if (matrix.m() != n_elems || matrix.n() != n_elems)
//        matrix.resize(n_elems, n_elems);
//    matrix.zero();
//    
//    // for shape design variable, the shape factor matrix will be calculated for
//    // the perturbed mesh. The base shape factor matrix, is then retrived and the
//    // sensitivity is then calculated
//    unsigned int DV_ID = 0;
//    double perturbation = 0.0;
//    if (sensitivity)
//    {
//        // for a Property DV, the sensitivity of this quantity is zero
//        const DesignData::DesignParameter& dparam =
//        this->thermal_discipline.getAnalysisDriver().getCurrentDesignParameter();
//        if (dparam.getParameterTypeEnumID() == DesignData::PROPERTY_PARAMETER::num())
//            return;
//        else
//        {
//            DV_ID = dparam.getID();
//            perturbation = dparam.getPerturbationStepSize();
//        }
//    }
//    
//    // get the iterators to the elements of the radiation mesh
//    RadiationElemConstIterator elem_it, elem_end, nested_elem_it;
//    elem_it = this->radiation_cavity.getElemBeginIterator(DV_ID);
//    elem_end = this->radiation_cavity.getElemEndIterator(DV_ID);
//    
//    unsigned int dof_elem1 =0, dof_elem2 = 0;
//    RadiationElement *rad_elem1 = NULL, *rad_elem2 = NULL;
//    RadiationElementPair elem_pair;
//    
//    // do a nested iteration on the elements
//    for ( ; elem_it != elem_end; elem_it++)
//    {
//        rad_elem1 = *elem_it;
//        
//        dof_elem1 = rad_elem1->ID();
//        
//        // inside the loop, create a pair of the two elements for which shape
//        // factor calculations will have to be performed
//        nested_elem_it = elem_it;
//        nested_elem_it++;
//        
//        for (; nested_elem_it != elem_end; nested_elem_it++)
//        {
//            // create an element pair for the two elements, and calculate the shape factor
//            rad_elem2 = *nested_elem_it;
//            
//            dof_elem2 = rad_elem2->ID();
//            
//            elem_pair.reinit(rad_elem1, rad_elem2);
//#ifdef WITH_MATHEMATICA
//            elem_pair.setMathLink(FESystem::mathlink_link);
//#endif
//            
//            // get the shape factors from this element pair, and set the appropriate
//            // matrix elements in the shape factor matrix. note that both f12 and f21
//            // need to be obtained.
//            std::pair<double, double> shape_factors =
//            elem_pair.getShapeFactor(RadiationElementPair::MITALAS_CONTOUR_INT);
//            
//            // next, insert the two shape factors into the shape factor matrix
//            matrix.set(dof_elem1, dof_elem2, shape_factors.first);
//            matrix.set(dof_elem2, dof_elem1, shape_factors.second);
//            
//            // clear the radiation element pair initialization
//            elem_pair.clear();
//        }
//    }
//    
//    if (sensitivity) // implies shape sensitivity is being sought
//    {
//        FESystem::Numerics::MatrixBase<FESystemDouble> matrix_base(n_elems, n_elems);
//        this->getMatrixQuantity(matrix_base,
//                                FESystem::Raditaion::RadiationCavityAnalysis::F_MATRIX);
//        matrix.add(-1.0, matrix_base, true);
//        matrix.scale(1.0/perturbation);
//    }
//    
//    this->fesystem_controller.performance_logging->unsetEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::calculateShapeFactorMatrix()",
//     this->performance_logging_header);
//}
//
//
//
//
//
//void FESystem::Raditaion::RadiationCavityAnalysis::calculateA_Matrix(FESystem::Numerics::MatrixBase<FESystemDouble>* matrix_ptr,
//                                                                     bool sensitivity)
//{
//    this->fesystem_controller.performance_logging->setEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::calculateA_Matrix()",
//     this->performance_logging_header);
//    
//    // get a pointer to the matrices
//    FESystem::Numerics::MatrixBase<FESystemDouble>& matrix = *matrix_ptr;
//    
//    unsigned int n_elems = this->radiation_cavity.nRadiationElems();
//    if (matrix.m() != n_elems || matrix.n() != n_elems)
//        matrix.resize(n_elems, n_elems);
//    matrix.zero();
//    
//    // for shape design variable, the shape factor matrix will be calculated for
//    // the perturbed mesh. The base shape factor matrix, is then retrived and the
//    // sensitivity is then calculated
//    if (sensitivity)
//    {
//        // for a Property DV, the sensitivity of this quantity is zero
//        const DesignData::DesignParameter& dparam =
//        this->thermal_discipline.getAnalysisDriver().getCurrentDesignParameter();
//        if (dparam.getParameterTypeEnumID() == DesignData::PROPERTY_PARAMETER::num())
//            return;
//    }
//    
//    
//    // in case the DV is NULL, the F_matrix will contain the actual shape factor matrix,
//    // otherwise, it will contain the sensitivity of the shape factor matrix
//    FESystem::Numerics::MatrixBase<FESystemDouble> F_matrix(n_elems, n_elems);
//    if (sensitivity)
//        this->getMatrixQuantity(F_matrix, FESystem::Raditaion::RadiationCavityAnalysis::F_MATRIX_SENSITIVITY);
//    else
//        this->getMatrixQuantity(F_matrix, FESystem::Raditaion::RadiationCavityAnalysis::F_MATRIX);
//    
//    // get the value of sigma
//    double sigma =
//    this->fesystem_controller.analysis_case->getRealParameter("S_B_CONSTANT");
//    
//    // calculate [A] = sigma * ([I] - [F]) or its sensitivity
//    matrix.add(-1.0, F_matrix, true);
//    
//    if (!sensitivity)
//        matrix.shift(1.0);
//    
//    // finally, scale the parameter with sigma
//    matrix.scale(sigma);
//    
//    
//    this->fesystem_controller.performance_logging->unsetEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::calculateA_Matrix()",
//     this->performance_logging_header);
//}
//
//
//
//void FESystem::Raditaion::RadiationCavityAnalysis::calculateB_Matrix(FESystem::Numerics::MatrixBase<FESystemDouble>* matrix_ptr,
//                                                                     bool sensitivity)
//{
//    this->fesystem_controller.performance_logging->setEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::calculateB_Matrix()",
//     this->performance_logging_header);
//    
//    // get a pointer to the matrices
//    FESystem::Numerics::MatrixBase<FESystemDouble>& matrix = *matrix_ptr;
//    
//    unsigned int n_elems = this->radiation_cavity.nRadiationElems();
//    if (matrix.m() != n_elems || matrix.n() != n_elems)
//        matrix.resize(n_elems, n_elems);
//    matrix.zero();
//    
//    const DesignData::DesignParameter* dparam = NULL;
//    unsigned int DV_ID = 0;
//    if (sensitivity)
//    {
//        // for a Property DV, the sensitivity of this quantity is zero
//        dparam =
//        &(this->thermal_discipline.getAnalysisDriver().getCurrentDesignParameter());
//        DV_ID = dparam->getID();
//        
//        if (dparam->getParameterTypeEnumID() == DesignData::PROPERTY_PARAMETER::num())
//            this->getMatrixQuantity(matrix, FESystem::Raditaion::RadiationCavityAnalysis::F_MATRIX);
//        else
//            this->getMatrixQuantity(matrix, FESystem::Raditaion::RadiationCavityAnalysis::F_MATRIX_SENSITIVITY);
//    }
//    else
//        this->getMatrixQuantity(matrix, FESystem::Raditaion::RadiationCavityAnalysis::F_MATRIX);
//    
//    
//    Property::PropertyDatabase& property_database =
//    *(this->fesystem_controller.property_database.get());
//    MeshDS::FEMeshData& mesh_data = *(this->radiation_cavity.getMeshData());
//    
//    // calculate [B] = [1/epsilon] - [(1 - epsilon)/epsilon] [F]
//    unsigned int n_rad_elems = this->radiation_cavity.nRadiationElems(),
//    fe_elem_ID = 0,  property_ID = 0;
//    double  factor1 = 0.0, factor2 = 0.0;
//    
//    FESystem::Numerics::VectorBase<FESystemDouble> factor_1_vec(n_elems),
//    factor_2_vec(n_elems);
//    
//    for (unsigned int i=0; i < n_rad_elems; i++)
//    {
//        // get the ID of the finite elem, get the property card for this elem
//        fe_elem_ID = this->radiation_cavity.getElemFromID(i).FeElementID();
//        const Elem* elem = mesh_data.getElemFromForeignID(fe_elem_ID);
//        property_ID = mesh_data.getElemPropertyID(elem);
//        ElemDataCard& elem_data_card = property_database.getElemDataCardFromID(property_ID);
//        
//        // if the property is temperature dependent, then initialize the property card
//        // at the average temperature of the finite element to which this radiation
//        // element belongs
//        this->initializePropertyCardAtElementTemperature(elem_data_card, i);
//        
//        if (sensitivity)
//        {
//            if (dparam->getParameterTypeEnumID() == DesignData::PROPERTY_PARAMETER::num())
//            {
//                // sensitivity wrt global property parameter
//                elem_data_card.getFactorSensitivityForGlobalParameter(factor1,
//                                                                      RADIATION_EPSILON_FACTOR_1::num(),
//                                                                      DV_ID);
//                elem_data_card.getFactorSensitivityForGlobalParameter(factor2,
//                                                                      RADIATION_EPSILON_FACTOR_2::num(),
//                                                                      DV_ID);
//            }
//            else
//            {
//                // shape sensitivity
//                elem_data_card.getFactor(factor1, RADIATION_EPSILON_FACTOR_1::num());
//                elem_data_card.getFactor(factor2, RADIATION_EPSILON_FACTOR_2::num());
//            }
//        }
//        else
//        {
//            elem_data_card.getFactor(factor1, RADIATION_EPSILON_FACTOR_1::num());
//            elem_data_card.getFactor(factor2, RADIATION_EPSILON_FACTOR_2::num());
//            if (factor1 <= 1.0)
//                std::cout << "emissivity greater than 1, factor : " << factor1 << std::endl;
//        }
//        
//        
//        
//        factor_2_vec.set(i, -factor2);
//        factor_1_vec.set(i, factor1);
//        
//        this->clearPropertyCardInitialization(elem_data_card);
//    }
//    
//    matrix.leftDiagonalScale(factor_2_vec);
//    
//    if (!(sensitivity &&
//          dparam->getParameterTypeEnumID() == DesignData::SHAPE_PARAMETER::num()))
//        matrix.diagonalAdd(1.0, factor_1_vec);
//    
//    this->fesystem_controller.performance_logging->unsetEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::calculateB_Matrix()",
//     this->performance_logging_header);
//}
//
//
//
//void FESystem::Raditaion::RadiationCavityAnalysis::calculateG_Matrix
//(FESystem::Numerics::MatrixBase<FESystemDouble>* matrix_ptr,
// bool sensitivity)
//{
//    this->fesystem_controller.performance_logging->setEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::calculateG_Matrix()",
//     this->performance_logging_header);
//    
//    // get a pointer to the matrices
//    FESystem::Numerics::MatrixBase<FESystemDouble>& matrix = *matrix_ptr;
//    
//    unsigned int n_elems = this->radiation_cavity.nRadiationElems();
//    unsigned int n_nodes = this->radiation_cavity.getFENodes().size();
//    if (matrix.m() != n_elems || matrix.n() != n_nodes)
//        matrix.resize(n_elems, n_nodes);
//    matrix.zero();
//    
//    // for shape design variable, the shape factor matrix will be calculated for
//    // the perturbed mesh. The base shape factor matrix, is then retrived and the
//    // sensitivity is then calculated
//    unsigned int DV_ID = 0;
//    double perturbation = 0.0;
//    if (sensitivity)
//    {
//        // for a Property DV, the sensitivity of this quantity is zero
//        const DesignData::DesignParameter& dparam =
//        this->thermal_discipline.getAnalysisDriver().getCurrentDesignParameter();
//        if (dparam.getParameterTypeEnumID() == DesignData::PROPERTY_PARAMETER::num())
//            return;
//        else
//        {
//            DV_ID = dparam.getID();
//            perturbation = dparam.getPerturbationStepSize();
//        }
//    }
//    
//    
//    // iterate over all the radiation elements, get G matrix from each,
//    // and add them to the global G matrix. The rows where the G matrix rows
//    // will be added will be the dof of the radiation element
//    // and the column will be the node number for the FE to which this
//    // radiation element belongs
//    unsigned int elem_dof = 0, n_elem_nodes = 0;
//    std::vector<unsigned int> fe_node_ID_vec;
//    FESystem::Numerics::VectorBase<FESystemDouble> elem_G_vec;
//    
//    // if sensitivity is requested, then the DV_ID will be nonzero, due to which
//    // the elem iterators will point to the perturbed geom elements
//    RadiationElemConstIterator elem_it, elem_end;
//    elem_it = this->radiation_cavity.getElemBeginIterator(DV_ID);
//    elem_end = this->radiation_cavity.getElemEndIterator(DV_ID);
//    
//    for (; elem_it != elem_end; elem_it++)
//    {
//        elem_dof = (*elem_it)->ID();
//        this->radiation_cavity.getInternalNodeIDForRadElem
//        (elem_dof, fe_node_ID_vec);
//        
//        n_elem_nodes = fe_node_ID_vec.size();
//        
//        elem_G_vec.resize(n_elem_nodes);
//        
//        (*elem_it)->getRadElemG_Vector(elem_G_vec);
//        
//        matrix.setRow(elem_dof, fe_node_ID_vec, elem_G_vec);
//    }
//    
//    if (sensitivity) // implies shape sensitivity is being sought
//    {
//        FESystem::Numerics::MatrixBase<FESystemDouble> matrix_base(n_elems, n_nodes);
//        this->getMatrixQuantity(matrix_base, FESystem::Raditaion::RadiationCavityAnalysis::G_MATRIX);
//        matrix.add(-1.0, matrix_base, true);
//        matrix.scale(1.0/perturbation);
//    }
//    
//    this->fesystem_controller.performance_logging->unsetEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::calculateG_Matrix()",
//     this->performance_logging_header);
//    
//}
//
//
//void FESystem::Raditaion::RadiationCavityAnalysis::calculateA_InverseMatrix
//(FESystem::Numerics::MatrixBase<FESystemDouble>* matrix_ptr)
//{
//    this->fesystem_controller.performance_logging->setEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::calculateA_InverseMatrix()",
//     this->performance_logging_header);
//    
//    //  std::cout<< "asked for A_inv" << std::endl;
//    
//    FESystem::Numerics::MatrixBase<FESystemDouble>& A_inverse_matrix = *matrix_ptr;
//    
//    unsigned int n_elems = this->radiation_cavity.nRadiationElems();
//    if (A_inverse_matrix.m() != n_elems || A_inverse_matrix.n() != n_elems)
//        A_inverse_matrix.resize(n_elems, n_elems);
//    A_inverse_matrix.zero();
//    
//    // get the A_matrix
//    FESystem::Numerics::MatrixBase<FESystemDouble> A_matrix(n_elems, n_elems);
//    FESystem::Numerics::VectorBase<FESystemDouble> tmp_vector(n_elems), sol_vector(n_elems);
//    this->getMatrixQuantity(A_matrix, FESystem::Raditaion::RadiationCavityAnalysis::A_MATRIX);
//    A_matrix.finalAssemble();
//    tmp_vector.zero();
//    
//    // create a petsc mat and set its values to the B matrix.
//    Mat matrix = A_matrix.getMat();
//    
//    // create a PC context and create the LU factorization
//    PC petsc_pc;
//    
//    PCCreate(PETSC_COMM_SELF, &petsc_pc);
//    PCSetType(petsc_pc , PCLU);
//    PCFactorSetUseInPlace(petsc_pc);
//    PCSetOperators(petsc_pc, matrix, matrix, SAME_NONZERO_PATTERN);
//    PCSetUp(petsc_pc);
//    
//    Mat lu_mat;
//	MatStructure str;
//    PCGetOperators(petsc_pc, PETSC_NULL, &lu_mat, &str);
//    
//    for (unsigned int i=0; i < n_elems; i++)
//    {
//        tmp_vector.set(i,1.0);
//        tmp_vector.assemble();
//        sol_vector.zero();
//        MatSolve(lu_mat, tmp_vector.getVec(), sol_vector.getVec());
//        tmp_vector.set(i,0.0);
//        
//        for (unsigned int j=0; j<n_elems; j++)
//            A_inverse_matrix.set(j,i,sol_vector.el(j));
//    }
//    
//    PCDestroy(petsc_pc);
//    
//    this->fesystem_controller.performance_logging->unsetEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::calculateA_InverseMatrix()",
//     this->performance_logging_header);
//}
//
//
//
//void FESystem::Raditaion::RadiationCavityAnalysis::calculateB_InverseMatrix
//(FESystem::Numerics::MatrixBase<FESystemDouble>* matrix_ptr)
//{
//    this->fesystem_controller.performance_logging->setEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::calculateB_InverseMatrix()",
//     this->performance_logging_header);
//    
//    //  std::cout<< "asked for B_inv" << std::endl;
//    
//    FESystem::Numerics::MatrixBase<FESystemDouble>& B_inverse_matrix = *matrix_ptr;
//    
//    unsigned int n_elems = this->radiation_cavity.nRadiationElems();
//    if (B_inverse_matrix.m() != n_elems || B_inverse_matrix.n() != n_elems)
//        B_inverse_matrix.resize(n_elems, n_elems);
//    B_inverse_matrix.zero();
//    
//    // get the B_matrix and the B_inverse matrices.
//    FESystem::Numerics::MatrixBase<FESystemDouble> B_matrix(n_elems, n_elems);
//    FESystem::Numerics::VectorBase<FESystemDouble> tmp_vector(n_elems), sol_vector(n_elems);
//    this->getMatrixQuantity(B_matrix, FESystem::Raditaion::RadiationCavityAnalysis::B_MATRIX);
//    B_matrix.finalAssemble();
//    tmp_vector.zero();
//    
//    // create a petsc mat and set its values to the B matrix.
//    Mat matrix = B_matrix.getMat();
//    
//    // create a PC context and create the LU factorization
//    PC petsc_pc;
//    
//    PCCreate(PETSC_COMM_SELF, &petsc_pc);
//    PCSetType(petsc_pc , PCLU);
//    PCFactorSetUseInPlace(petsc_pc);
//    PCSetOperators(petsc_pc, matrix, matrix, SAME_NONZERO_PATTERN);
//    PCSetUp(petsc_pc);
//    
//    Mat lu_mat;
//    MatStructure str;
//    PCGetOperators(petsc_pc, PETSC_NULL, &lu_mat, &str);
//    
//    for (unsigned int i=0; i < n_elems; i++)
//    {
//        tmp_vector.set(i,1.0);
//        tmp_vector.assemble();
//        sol_vector.zero();
//        MatSolve(lu_mat, tmp_vector.getVec(), sol_vector.getVec());
//        tmp_vector.set(i,0.0);
//        
//        for (unsigned int j=0; j<n_elems; j++)
//            B_inverse_matrix.set(j,i,sol_vector.el(j));
//    }
//    
//    PCDestroy(petsc_pc);
//    
//    this->fesystem_controller.performance_logging->unsetEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::calculateB_InverseMatrix()",
//     this->performance_logging_header);
//}
//
//
//
//
//
//void FESystem::Raditaion::RadiationCavityAnalysis::calculateEmissivityVector
//(FESystem::Numerics::VectorBase<FESystemDouble>* vector_ptr)
//{
//    this->fesystem_controller.performance_logging->setEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::calculateEmissivityVector()",
//     this->performance_logging_header);
//    
//    FESystem::Numerics::VectorBase<FESystemDouble>& eps_vec = *vector_ptr;
//    
//    unsigned int n_elems = this->radiation_cavity.nRadiationElems();
//    if (eps_vec.size() != n_elems)
//        eps_vec.resize(n_elems);
//    eps_vec.zero();
//    
//    
//    Property::PropertyDatabase& property_database =
//    *(this->fesystem_controller.property_database.get());
//    MeshDS::FEMeshData& mesh_data = *(this->radiation_cavity.getMeshData());
//    
//    unsigned int fe_elem_ID = 0,  property_ID = 0;
//    double  factor1 = 0.0;
//    
//    for (unsigned int i=0; i < n_elems; i++)
//    {
//        // get the ID of the finite elem, get the property card for this elem
//        fe_elem_ID = this->radiation_cavity.getElemFromID(i).FeElementID();
//        const Elem* elem = mesh_data.getElemFromForeignID(fe_elem_ID);
//        property_ID = mesh_data.getElemPropertyID(elem);
//        ElemDataCard& elem_data_card = property_database.getElemDataCardFromID(property_ID);
//        
//        // if the property is temperature dependent, then initialize the property card
//        // at the average temperature of the finite element to which this radiation
//        // element belongs
//        this->initializePropertyCardAtElementTemperature(elem_data_card, i);
//        
//        elem_data_card.getPropertyValueFromMaterialCard(EMISSIVITY::num(), factor1);
//        if (factor1 > 1.0)
//            std::cout << "emissivity greater than 1, value : " << factor1 << std::endl;
//        
//        
//        eps_vec.set(i, factor1);
//        
//        this->clearPropertyCardInitialization(elem_data_card);
//    }
//    
//    
//    this->fesystem_controller.performance_logging->unsetEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::calculateEmissivityVector()",
//     this->performance_logging_header);
//    
//}
//
//
//
//void
//FESystem::Raditaion::RadiationCavityAnalysis::calculateABInvFactorMatrix
//(FESystem::Numerics::MatrixBase<FESystemDouble>* matrix_ptr)
//{
//    this->fesystem_controller.performance_logging->setEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::calculateABInvFactorMatrix()",
//     this->performance_logging_header);
//    
//    // this is now being limited only to cases with temperature "dependent"
//    // material properties
//    Assert(this->thermal_discipline.checkPropertyDependenceOnTemperature(),
//           ExcInternalError());
//    
//    FESystem::Numerics::MatrixBase<FESystemDouble>& factor_matrix = *matrix_ptr;
//    
//    unsigned int n_elems = this->radiation_cavity.nRadiationElems();
//    factor_matrix.zero();
//    
//    bool recalculate = true;
//    bool store_base = false;
//    
//    if (!this->thermal_discipline.checkPropertyDependenceOnTemperature())
//    {
//        recalculate = true;
//        store_base = false;
//    }
//    else
//    {
//        if(!this->radiation_cavity.ifApproximateABInvFactor())
//        {
//            recalculate = true;
//            store_base = false;
//        }
//        else
//        {
//            if (fmod(1.0*this->n_ABinv_calculations,
//                     1.0*this->radiation_cavity.getNItertionsToRecalculateABInv()) == 0.0)
//            {
//                recalculate = true;
//                store_base = true;
//            }
//            else
//            {
//                
//                recalculate = false;
//                store_base = false;
//            }
//        }
//    }
//    
//    //  std::cout << "if approx  " << this->radiation_cavity.ifApproximateABInvFactor() << std::endl;
//    //  std::cout << "n_AB_inv,n_iter_to_recalc = " << this->n_ABinv_calculations << " "
//    //    << this->radiation_cavity.getNItertionsToRecalculateABInv() << std::endl;
//    //  std::cout << "fmod =" << fmod(1.0*this->n_ABinv_calculations,
//    //                    1.0*this->radiation_cavity.getNItertionsToRecalculateABInv()) << std::endl;
//    
//    if (recalculate)
//    {
//        // simply calculate
//        this->getMatrixQuantity(factor_matrix,
//                                FESystem::Raditaion::RadiationCavityAnalysis::B_INVERSE_MATRIX);
//        FESystem::Numerics::MatrixBase<FESystemDouble> matrix1(n_elems, n_elems);
//        
//        this->getMatrixQuantity(matrix1, FESystem::Raditaion::RadiationCavityAnalysis::A_MATRIX);
//        
//        factor_matrix.leftMultiply(matrix1);
//        
//        std::cout<< "recalculating AB_inv" << std::endl;
//        
//        if (store_base)
//        {
//            std::string name;
//            std::ostringstream cavity_id;
//            cavity_id << this->radiation_cavity.getCavityID();
//            
//            name.clear();
//            name = "BaseEpsilonVector";
//            name += "_Cavity_";
//            name += cavity_id.str();
//            
//            // calculate the epsilon vector and store it
//            FESystem::Numerics::VectorBase<FESystemDouble> vec1(n_elems);
//            this->calculateEmissivityVector(&vec1);
//            FESystemDatabase::TimeIndependentDataInfo time_independent_data_info;
//            time_independent_data_info.setDisciplineEnumID(Discipline::RADIATION_DISCIPLINE::num());
//            time_independent_data_info.setLoadCase(this->thermal_discipline.getAnalysisDriver().
//                                                   getCurrentLoadCase());
//            time_independent_data_info.setName(name);
//            
//            this->fesystem_controller.global_data_storage->storeVector(time_independent_data_info,
//                                                                       vec1);
//            
//            // create a name for storing the matrix
//            std::ostringstream num;
//            name.clear();
//            name = "BaseABInvMatrix";
//            name += "_Cavity_";
//            name += cavity_id.str();
//            
//            time_independent_data_info.clear();
//            time_independent_data_info.setDisciplineEnumID(Discipline::RADIATION_DISCIPLINE::num());
//            time_independent_data_info.setLoadCase(this->thermal_discipline.getAnalysisDriver().
//                                                   getCurrentLoadCase());
//            time_independent_data_info.setName(name);
//            
//            this->fesystem_controller.global_data_storage->storeMatrix(time_independent_data_info,
//                                                                       factor_matrix);
//        }
//    }
//    else
//    {
//        std::cout<< "returning base AB_inv" << std::endl;
//        
//        // create a name for storing the matrix
//        std::string name;
//        std::ostringstream cavity_id;
//        cavity_id << this->radiation_cavity.getCavityID();
//        std::ostringstream num;
//        name.clear();
//        name = "BaseABInvMatrix";
//        name += "_Cavity_";
//        name += cavity_id.str();
//        
//        FESystemDatabase::TimeIndependentDataInfo time_independent_data_info;
//        time_independent_data_info.clear();
//        time_independent_data_info.setDisciplineEnumID(Discipline::RADIATION_DISCIPLINE::num());
//        time_independent_data_info.setLoadCase(this->thermal_discipline.getAnalysisDriver().
//                                               getCurrentLoadCase());
//        time_independent_data_info.setName(name);
//        
//        this->fesystem_controller.global_data_storage->fillMatrix(time_independent_data_info,
//                                                                  factor_matrix);
//    }
//    
//    // now increment the counter
//    this->n_ABinv_calculations++;
//    
//    this->fesystem_controller.performance_logging->unsetEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::calculateABInvFactorMatrix()",
//     this->performance_logging_header);
//    
//}
//
//
//void FESystem::Raditaion::RadiationCavityAnalysis::getBaseEpsilonVector
//(FESystem::Numerics::VectorBase<FESystemDouble>& vec)
//{
//    std::string name;
//    std::ostringstream cavity_id;
//    cavity_id << this->radiation_cavity.getCavityID();
//    
//    name.clear();
//    name = "BaseEpsilonVector";
//    name += "_Cavity_";
//    name += cavity_id.str();
//    
//    FESystemDatabase::TimeIndependentDataInfo time_independent_data_info;
//    time_independent_data_info.setDisciplineEnumID(Discipline::RADIATION_DISCIPLINE::num());
//    time_independent_data_info.setLoadCase(this->thermal_discipline.getAnalysisDriver().
//                                           getCurrentLoadCase());
//    time_independent_data_info.setName(name);
//    
//    this->fesystem_controller.global_data_storage->fillVector(time_independent_data_info,
//                                                              vec);
//}
//
//
//
//void FESystem::Raditaion::RadiationCavityAnalysis::calculateFactorMatrix
//(FESystem::Numerics::MatrixBase<FESystemDouble>* matrix_ptr)
//{
//    this->fesystem_controller.performance_logging->setEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::calculateFactorMatrix()",
//     this->performance_logging_header);
//    
//    // this is now being limited only to cases with temperature "independent"
//    // material properties
//    Assert(!this->thermal_discipline.checkPropertyDependenceOnTemperature(),
//           ExcInternalError());
//    
//    FESystem::Numerics::MatrixBase<FESystemDouble>& factor_matrix = *matrix_ptr;
//    
//    unsigned int n_elems = this->radiation_cavity.nRadiationElems();
//    factor_matrix.zero();
//    
//    // get references to the A and G matrices
//    // the factor matrix will be calculated using the expression
//    // factor = [G]^T  [A] [B]^(-1)
//    this->getMatrixQuantity(factor_matrix,
//                            FESystem::Raditaion::RadiationCavityAnalysis::B_INVERSE_MATRIX);
//    FESystem::Numerics::MatrixBase<FESystemDouble> matrix1(n_elems, n_elems);
//    
//    this->getMatrixQuantity(matrix1, FESystem::Raditaion::RadiationCavityAnalysis::A_MATRIX);
//    
//    factor_matrix.leftMultiply(matrix1);
//    
//    this->getMatrixQuantity(matrix1, FESystem::Raditaion::RadiationCavityAnalysis::G_MATRIX);
//    
//    factor_matrix.leftMultiplyTranspose(matrix1);
//    
//    this->fesystem_controller.performance_logging->unsetEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::calculateFactorMatrix()",
//     this->performance_logging_header);
//    
//}
//
//
//
//
//void FESystem::Raditaion::RadiationCavityAnalysis::getFENodalDofIDs(std::vector<unsigned int>& vector)
//{
//    // get the vector of FE nodes from the cavity
//    const std::vector<Node*>& node_vector =
//    this->radiation_cavity.getFENodes();
//    
//    std::vector<Node*>::const_iterator it, end;
//    it = node_vector.begin();
//    end = node_vector.end();
//    
//    unsigned int dof_ID = 0;
//    
//    // clear the vector, and then insert the nodes in it
//    vector.clear();
//    for (; it != end; it++)
//    {
//        dof_ID = (*it)->dof_number(0,0,0);
//        vector.push_back(dof_ID);
//    }
//}
//
//
//
//void FESystem::Raditaion::RadiationCavityAnalysis::calculateRadElemTempVector
//(FESystem::Numerics::VectorBase<FESystemDouble>* rad_temp_vector_ptr,
// FESystem::Numerics::VectorBase<FESystemDouble>& nodal_temp_vector,
// const bool reduce_power_by_four,
// const bool sensitivity)
//{
//    FESystem::Numerics::VectorBase<FESystemDouble>& rad_temp_vector = *rad_temp_vector_ptr;
//    
//    unsigned int n_elems = this->radiation_cavity.nRadiationElems();
//    if (rad_temp_vector.size() != n_elems)
//        rad_temp_vector.resize(n_elems);
//    rad_temp_vector.zero();
//    
//    // get the value of absolute temperature
//    const double temp_abs =
//    this->fesystem_controller.analysis_case->getRealParameter("ABSOLUTE_TEMP");
//    
//    unsigned int DV_ID = 0;
//    double perturbation = 0.0;
//    if (sensitivity)
//    {
//        // for a Property DV, the sensitivity of this quantity is zero
//        const DesignData::DesignParameter& dparam =
//        this->thermal_discipline.getAnalysisDriver().getCurrentDesignParameter();
//        if (dparam.getParameterTypeEnumID() == DesignData::PROPERTY_PARAMETER::num())
//            return;
//        else
//        {
//            DV_ID = dparam.getID();
//            perturbation = dparam.getPerturbationStepSize();
//        }
//    }
//    
//    // first calculate the temperatures for the radiation elements
//    // this is done using the expression
//    // T_rad^4 = 1/A * int_A (N^T Temp)^4 dA
//    
//    unsigned int elem_dof = 0;
//    
//    std::vector<Point> points;
//    double qp_value = 0.0, value = 0.0;
//    std::vector<unsigned int> fe_node_ID_vec;
//    
//    std::auto_ptr<QBase> quadrature(new QGauss(2,FIFTH));
//    std::auto_ptr<FEBase> fe_rad(new FE<2, LAGRANGE>(FEType()));
//    std::auto_ptr<FEBase> global_fe(new FE<2,LAGRANGE>(FEType()));
//    
//    fe_rad->attach_quadrature_rule(quadrature.get());
//    
//    RadiationElemConstIterator elem_it, elem_end;
//    elem_it = this->radiation_cavity.getElemBeginIterator(DV_ID);
//    elem_end = this->radiation_cavity.getElemEndIterator(DV_ID);
//    
//    
//    for (; elem_it != elem_end; elem_it++)
//    {
//        elem_dof = (*elem_it)->ID();
//        
//        fe_node_ID_vec.clear();
//        
//        // get the vector of internal node IDs
//        this->radiation_cavity.getInternalNodeIDForRadElem
//        (elem_dof, fe_node_ID_vec);
//        
//        // get the local elem and the nodes in the local coordinates
//        Elem* local_rad_elem =
//        (*elem_it)->getLocalRadiationGeometricElem();
//        const std::vector<Node*> local_nodes =
//        (*elem_it)->getLocalCoordinateNodes();
//        
//        fe_rad->reinit(local_rad_elem);
//        
//        // now obtain the value of the shape function for
//        // the associated FE
//        const std::vector<std::vector<Real> >& phi_rad = fe_rad->get_phi();
//        const std::vector<Real>& JxW_rad = fe_rad->get_JxW();
//        
//        points.clear();
//        
//        for (unsigned int qp = 0; qp < phi_rad[0].size(); qp++)
//        {
//            Point pt;
//            for (unsigned int i=0; i < phi_rad.size(); i++)
//            {
//                Node& node = *(local_nodes[i]);
//                pt(0) += phi_rad[i][qp] * node(0);
//                pt(1) += phi_rad[i][qp] * node(1);
//            }
//            points.push_back(pt);
//        }
//        
//        // next, get the values of the shape function for the global
//        // FE at the quad points for the local rad elem
//        global_fe->reinit(local_rad_elem, &points);
//        const std::vector<std::vector<Real> >& phi_global = global_fe->get_phi();
//        
//        // now, calculate the G_vector
//        value = 0.0;
//        for (unsigned int qp=0; qp < JxW_rad.size(); qp++)
//        {
//            qp_value = 0.0;
//            for (unsigned int i=0; i < phi_global.size(); i++)
//                qp_value += phi_global[i][qp] *
//                (nodal_temp_vector.el(fe_node_ID_vec[i]) + temp_abs);
//            
//            value += pow(qp_value,4) * JxW_rad[qp];
//        }
//        
//        rad_temp_vector.set(elem_dof, value / (*elem_it)->getArea());
//    }
//    
//    if (reduce_power_by_four)
//        rad_temp_vector.power(0.25);
//    
//    if (sensitivity) // shape sensitivity is calculated using finite differences
//    {
//        FESystem::Numerics::VectorBase<FESystemDouble> rad_temp_vector_base(n_elems);
//        this->calculateRadElemTempVector(&rad_temp_vector_base,
//                                         nodal_temp_vector,
//                                         reduce_power_by_four,
//                                         false);
//        rad_temp_vector.add(-1.0, rad_temp_vector_base);
//        rad_temp_vector.scale(1.0/perturbation);
//    }
//}
//
//
//
//
//
//void FESystem::Raditaion::RadiationCavityAnalysis::getFENodalLoadVector
//(FESystem::Numerics::VectorBase<FESystemDouble>& nodal_temp_vector,
// FESystem::Numerics::VectorBase<FESystemDouble>& nodal_load_vector)
//{
//    this->fesystem_controller.performance_logging->setEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::getFENodalLoadVector()",
//     this->performance_logging_header);
//    
//    unsigned int n_elems = this->radiation_cavity.nRadiationElems();
//    unsigned int n_nodes = this->radiation_cavity.getFENodes().size();
//    
//    Assert(nodal_temp_vector.size() == n_nodes, ExcInternalError());
//    Assert(nodal_load_vector.size() == n_nodes, ExcInternalError());
//    
//    // first set the pointer to the nodal load vector
//    this->fe_nodal_temp_vector = &(nodal_temp_vector);
//    
//    FESystem::Numerics::VectorBase<FESystemDouble> rad_temp_vector(n_elems);
//    rad_temp_vector.zero();
//    this->calculateRadElemTempVector(&rad_temp_vector, nodal_temp_vector,
//                                     false, false);
//    
//    // get the factor_matrix, multiply to get the nodal load vector and return
//    FESystem::Numerics::MatrixBase<FESystemDouble> factor_matrix(n_elems, n_nodes);
//    nodal_load_vector.zero();
//    
//    if (!this->thermal_discipline.checkPropertyDependenceOnTemperature())
//    {
//        this->getMatrixQuantity(factor_matrix, FESystem::Raditaion::RadiationCavityAnalysis::FACTOR_MATRIX);
//        factor_matrix.rightMultiplyVector(rad_temp_vector,
//                                          nodal_load_vector);
//    }
//    else
//    {
//        this->getMatrixQuantity(factor_matrix, FESystem::Raditaion::RadiationCavityAnalysis::G_MATRIX);
//        FESystem::Numerics::MatrixBase<FESystemDouble> mat1(rad_temp_vector.size(),
//                                                           rad_temp_vector.size());
//        
//        this->getMatrixQuantity(mat1, FESystem::Raditaion::RadiationCavityAnalysis::AB_INV_MATRIX);
//        FESystem::Numerics::VectorBase<FESystemDouble>
//        vector1(n_elems), vector2(n_elems), vector3(n_elems);
//        mat1.rightMultiplyVector(rad_temp_vector, vector1);
//        
//        if(this->radiation_cavity.ifApproximateABInvFactor())
//        {
//            this->getTemperatureDependentFluxCorrectionVector(vector2, mat1);
//            vector2.pointwiseMultiply(vector1);
//            mat1.rightMultiplyVector(vector2, vector3);
//            vector1.add(1.0, vector3);
//        }
//        
//        factor_matrix.leftMultiplyVector(vector1, nodal_load_vector);
//    }
//    
//    
//    // now clear the pointer
//    this->fe_nodal_temp_vector = NULL;
//    
//    this->fesystem_controller.performance_logging->unsetEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::getFENodalLoadVector()",
//     this->performance_logging_header);
//}
//
//
//
//
//void FESystem::Raditaion::RadiationCavityAnalysis::getFENodalLoadVectorSensitivity
//(FESystem::Numerics::VectorBase<FESystemDouble>& nodal_temp_vector,
// FESystem::Numerics::VectorBase<FESystemDouble>& nodal_load_vector_sens)
//{
//    this->fesystem_controller.performance_logging->setEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::getFENodalLoadVectorSensitivity()",
//     this->performance_logging_header);
//    
//    
//    unsigned int n_elems = this->radiation_cavity.nRadiationElems();
//    unsigned int n_nodes = this->radiation_cavity.getFENodes().size();
//    
//    Assert(nodal_temp_vector.size() == n_nodes, ExcInternalError());
//    Assert(nodal_load_vector_sens.size() == n_nodes, ExcInternalError());
//    
//    // first set the pointer to the nodal load vector
//    this->fe_nodal_temp_vector = &(nodal_temp_vector);
//    nodal_load_vector_sens.zero();
//    
//    // calculate T_rad and d(T_rad)/d(Alpha)
//    FESystem::Numerics::VectorBase<FESystemDouble> rad_temp_vector(n_elems);
//    FESystem::Numerics::VectorBase<FESystemDouble> rad_temp_vector_sens(n_elems);
//    
//    rad_temp_vector.zero();
//    rad_temp_vector_sens.zero();
//    
//    this->calculateRadElemTempVector(&rad_temp_vector,
//                                     nodal_temp_vector,
//                                     false,
//                                     false);
//    
//    const DesignData::DesignParameter& dparam =
//    this->thermal_discipline.getAnalysisDriver().getCurrentDesignParameter();
//    if (dparam.getParameterTypeEnumID() == DesignData::SHAPE_PARAMETER::num())
//        this->calculateRadElemTempVector(&rad_temp_vector_sens,
//                                         nodal_temp_vector,
//                                         false,
//                                         true);
//    
//    
//    FESystem::Numerics::MatrixBase<FESystemDouble> G_matrix(n_elems, n_nodes);
//    this->getMatrixQuantity(G_matrix, FESystem::Raditaion::RadiationCavityAnalysis::G_MATRIX);
//    FESystem::Numerics::MatrixBase<FESystemDouble> G_matrix_sens(n_elems, n_nodes);
//    this->getMatrixQuantity(G_matrix_sens, FESystem::Raditaion::RadiationCavityAnalysis::G_MATRIX_SENSITIVITY);
//    
//    FESystem::Numerics::MatrixBase<FESystemDouble> A_matrix(n_elems, n_elems);
//    this->getMatrixQuantity(A_matrix, FESystem::Raditaion::RadiationCavityAnalysis::A_MATRIX);
//    FESystem::Numerics::MatrixBase<FESystemDouble> A_matrix_sens(n_elems, n_elems);
//    this->getMatrixQuantity(A_matrix_sens, FESystem::Raditaion::RadiationCavityAnalysis::A_MATRIX_SENSITIVITY);
//    
//    FESystem::Numerics::MatrixBase<FESystemDouble> B_matrix_sens(n_elems, n_elems);
//    this->getMatrixQuantity(B_matrix_sens, FESystem::Raditaion::RadiationCavityAnalysis::B_MATRIX_SENSITIVITY);
//    
//    // out of these vectors, vector4 has size n_nodes, hence, result of any vector multiplied by
//    // G matrix will go into this vector
//    FESystem::Numerics::VectorBase<FESystemDouble> vector1(n_elems), vector2(n_elems), vector3(n_elems),
//    vector4(n_nodes);
//    vector1.zero(); vector2.zero(); vector3.zero(); vector4.zero();
//    
//    std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemDouble> > AB_inv, A_inv, factor_matrix;
//    
//    if (!this->thermal_discipline.checkPropertyDependenceOnTemperature())
//    {
//        factor_matrix.reset(new FESystem::Numerics::MatrixBase<FESystemDouble>(n_nodes, n_elems));
//        this->getMatrixQuantity(*factor_matrix, FESystem::Raditaion::RadiationCavityAnalysis::FACTOR_MATRIX);
//        FESystem::Numerics::MatrixBase<FESystemDouble> B_inverse_matrix(n_elems, n_elems);
//        this->getMatrixQuantity(B_inverse_matrix,
//                                FESystem::Raditaion::RadiationCavityAnalysis::B_INVERSE_MATRIX);
//        B_inverse_matrix.rightMultiplyVector(rad_temp_vector,
//                                             vector1);
//    }
//    else
//    {
//        AB_inv.reset(new FESystem::Numerics::MatrixBase<FESystemDouble>(n_elems, n_elems));
//        A_inv.reset(new FESystem::Numerics::MatrixBase<FESystemDouble>(n_elems, n_elems));
//        this->getMatrixQuantity(*AB_inv, FESystem::Raditaion::RadiationCavityAnalysis::AB_INV_MATRIX);
//        this->getMatrixQuantity(*A_inv, FESystem::Raditaion::RadiationCavityAnalysis::A_INVERSE_MATRIX);
//        
//        AB_inv->rightMultiplyVector(rad_temp_vector, vector1);
//        A_inv->rightMultiplyVector(vector1, vector2);
//    }
//    
//    // calculate the first term
//    if (!this->thermal_discipline.checkPropertyDependenceOnTemperature())
//    {
//        vector3.zero();
//        A_matrix.rightMultiplyVector(vector1, vector3);
//        G_matrix_sens.leftMultiplyVector(vector3, nodal_load_vector_sens);
//    }
//    else
//        G_matrix_sens.leftMultiplyVector(vector1, nodal_load_vector_sens);
//    
//    // second term
//    vector3.zero();
//    if (!this->thermal_discipline.checkPropertyDependenceOnTemperature())
//        A_matrix_sens.rightMultiplyVector(vector1, vector3);
//    else
//        A_matrix_sens.rightMultiplyVector(vector2, vector3);
//    G_matrix.leftMultiplyVector(vector3, vector4);
//    nodal_load_vector_sens.add(1.0, vector4);
//    
//    // third term
//    if (!this->thermal_discipline.checkPropertyDependenceOnTemperature())
//    {
//        B_matrix_sens.rightMultiplyVector(vector1, vector3);
//        rad_temp_vector_sens.add(-1.0, vector3);
//        factor_matrix->rightMultiplyVector(rad_temp_vector_sens, vector4);
//    }
//    else
//    {
//        B_matrix_sens.rightMultiplyVector(vector2, vector3);
//        rad_temp_vector_sens.add(-1.0, vector3);
//        AB_inv->rightMultiplyVector(rad_temp_vector_sens, vector3);
//        G_matrix.leftMultiplyVector(vector3, vector4);
//    }
//    
//    nodal_load_vector_sens.add(1.0, vector4);
//    
//    // clear the nodal load vector pointer
//    this->fe_nodal_temp_vector = NULL;
//    
//    this->fesystem_controller.performance_logging->unsetEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::getFENodalLoadVectorSensitivity()",
//     this->performance_logging_header);
//}
//
//
//
//
//
//void FESystem::Raditaion::RadiationCavityAnalysis::getFEJacobianMatrix
//(FESystem::Numerics::VectorBase<FESystemDouble>& nodal_temp_val,
// FESystem::Numerics::MatrixBase<FESystemDouble>& jacobian)
//{
//    this->fesystem_controller.performance_logging->setEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::getFEJacobianMatrix()",
//     this->performance_logging_header);
//    
//    // first set the pointer to the nodal load vector
//    this->fe_nodal_temp_vector = &(nodal_temp_val);
//    
//    unsigned int elem_dof = 0,
//    n_elems = this->radiation_cavity.nRadiationElems(),
//    n_nodes = this->radiation_cavity.getFENodes().size();
//    
//    // get the factor_matrix, multiply it with the nodal load vector and return
//    FESystem::Numerics::MatrixBase<FESystemDouble> factor_matrix(n_nodes, n_elems);
//    
//    assert (jacobian.m() == jacobian.n());
//    assert (jacobian.m() == factor_matrix.m());
//    assert (jacobian.m() == nodal_temp_val.size());
//    jacobian.zero();
//    
//    
//    // get the value of absolute temperature
//    const double temp_abs =
//    this->fesystem_controller.analysis_case->getRealParameter("ABSOLUTE_TEMP");
//    
//    // create some temporary set of values for the procedure below
//    FESystem::Numerics::MatrixBase<FESystemDouble> elem_temp_matrix(n_elems, n_nodes);
//    elem_temp_matrix.zero();
//    
//    std::vector<Point> points;
//    double qp_value = 0.0;
//    FESystem::Numerics::VectorBase<FESystemDouble> value;
//    std::vector<unsigned int> fe_node_ID_vec;
//    
//    std::auto_ptr<QBase> quadrature(new QGauss(2,FIFTH));
//    std::auto_ptr<FEBase> fe_rad(new FE<2, LAGRANGE>(FEType()));
//    std::auto_ptr<FEBase> global_fe(new FE<2,LAGRANGE>(FEType()));
//    
//    fe_rad->attach_quadrature_rule(quadrature.get());
//    
//    RadiationElemConstIterator elem_it, elem_end;
//    elem_it = this->radiation_cavity.getElemBeginIterator();
//    elem_end = this->radiation_cavity.getElemEndIterator();
//    
//    
//    // this for loop calculates and stores in elem_temp_matrix the following quantity
//    // [d T_R^4 / d T_FE] = [int_Omega 4 ({N}^T {T_FE}+T_abs)^3 {N}^T d Omega]
//    for (; elem_it != elem_end; elem_it++)
//    {
//        elem_dof = (*elem_it)->ID();
//        
//        fe_node_ID_vec.clear();
//        
//        // get the vector of internal node IDs
//        this->radiation_cavity.getInternalNodeIDForRadElem
//        (elem_dof, fe_node_ID_vec);
//        
//        // get the local elem and the nodes in the local coordinates
//        Elem* local_rad_elem =
//        (*elem_it)->getLocalRadiationGeometricElem();
//        const std::vector<Node*> local_nodes =
//        (*elem_it)->getLocalCoordinateNodes();
//        
//        fe_rad->reinit(local_rad_elem);
//        
//        // now obtain the value of the shape function for
//        // the associated FE
//        const std::vector<std::vector<Real> >& phi_rad = fe_rad->get_phi();
//        const std::vector<Real>& JxW_rad = fe_rad->get_JxW();
//        
//        points.clear();
//        
//        for (unsigned int qp = 0; qp < phi_rad[0].size(); qp++)
//        {
//            Point pt;
//            for (unsigned int i=0; i < phi_rad.size(); i++)
//            {
//                Node& node = *(local_nodes[i]);
//                pt(0) += phi_rad[i][qp] * node(0);
//                pt(1) += phi_rad[i][qp] * node(1);
//            }
//            points.push_back(pt);
//        }
//        
//        // next, get the values of the shape function for the global
//        // FE at the quad points for the local rad elem
//        global_fe->reinit(local_rad_elem, &points);
//        const std::vector<std::vector<Real> >& phi_global = global_fe->get_phi();
//        
//        if (value.size() != phi_global.size())
//            value.resize(phi_global.size());
//        
//        // now, calculate the G_vector
//        value.zero();
//        for (unsigned int qp=0; qp < JxW_rad.size(); qp++)
//        {
//            // calculate the temperaure at the quadrature point first
//            qp_value = 0.0;
//            for (unsigned int i=0; i < phi_global.size(); i++)
//                qp_value += phi_global[i][qp] *
//                (nodal_temp_val.el(fe_node_ID_vec[i]) + temp_abs);
//            
//            qp_value = pow(qp_value, 3);
//            // next, add this to the vector
//            for (unsigned int i=0; i < phi_global.size(); i++)
//                value.add(i, 4 * phi_global[i][qp] * qp_value * JxW_rad[qp]);
//        }
//        
//        value.scale(1.0 / (*elem_it)->getArea());
//        // now insert the value in the elem matrix
//        elem_temp_matrix.setRow(elem_dof, fe_node_ID_vec, value);
//    }
//    
//    
//    
//    // now multiply the factor matrix and the elem_temp_matrix
//    // together to get the jacobian
//    if (!this->thermal_discipline.checkPropertyDependenceOnTemperature())
//    {
//        // for temp independent case, we have [d T_R^4 / d T_FE] in the elem temp matrix
//        // this only needs to be multiplied with the factor matrix G^T A B^{-1} on the left
//        this->getMatrixQuantity(factor_matrix, FESystem::Raditaion::RadiationCavityAnalysis::FACTOR_MATRIX);
//        elem_temp_matrix.leftMultiply(factor_matrix);
//        jacobian.copy(elem_temp_matrix);
//    }
//    else
//    {
//        FESystem::Numerics::MatrixBase<FESystemDouble>
//        A_inv(n_elems, n_elems), AB_inv(n_elems, n_elems), G_mat(n_elems, n_nodes);
//        this->getMatrixQuantity(G_mat, FESystem::Raditaion::RadiationCavityAnalysis::G_MATRIX);
//        this->getMatrixQuantity(A_inv, FESystem::Raditaion::RadiationCavityAnalysis::A_INVERSE_MATRIX);
//        this->getMatrixQuantity(AB_inv, FESystem::Raditaion::RadiationCavityAnalysis::AB_INV_MATRIX);
//        factor_matrix.copy(AB_inv);
//        factor_matrix.leftMultiplyTranspose(G_mat);
//        
//        FESystem::Numerics::VectorBase<FESystemDouble> temp_dep_vector(n_elems), rad_temp_vector;
//        temp_dep_vector.zero();
//        rad_temp_vector.zero();
//        this->calculateRadElemTempVector(&rad_temp_vector, nodal_temp_val,
//                                         false, // false for reduce power by 4
//                                         false); // false for sensitivity
//        
//        // get the contribution matrix due to temperature dependent properties
//        // this willl return the term [\frac{\partial B_{ij}}{\partial T_{R_k}}\psi_j]_diag
//        this->calculateB_TemperatureDependentJacobianContribution(temp_dep_vector, rad_temp_vector,
//                                                                  A_inv, AB_inv);
//        
//        // multiply it with the sensitivity of the radiation temperatures, with respect to the
//        // FE temperatures
//        // [d {T_rad} / d {T_fe}] = (1/4) * diag([(T_rad + T_a)^3])^(-1) * [elem_temp_matrix]
//        
//        // to perform this calculation, first obtain the temperature vector raised to 3rd power
//        rad_temp_vector.power(0.75);
//        rad_temp_vector.scale(-4.0); // -4.0 to take care of the -1.0 of this term
//        rad_temp_vector.invertValues();
//        
//        // this pointwise multiplication operation gives the vector
//        // -1.0 * [\frac{\partial B_{ij}}{\partial T_{R_k}}\psi_j]_diag * [(1/4) * (T_rad + T_a)^-3]_diag
//        temp_dep_vector.pointwiseMultiply(rad_temp_vector);
//        // this shift operation takes care of the [d((T_{rad}+T_{abs})^4)/d T_{rad}]
//        temp_dep_vector.shift(1.0);
//        
//        // now use this to calculate the following term using the approximation
//        // [G]^T [A] [B]^{-1} [lambda] [dT_R^4/dT_{FE}]
//        if(!this->radiation_cavity.ifApproximateABInvFactor())
//        {
//            elem_temp_matrix.leftDiagonalScale(temp_dep_vector);
//            elem_temp_matrix.leftMultiply(factor_matrix);
//        }
//        else
//        {
//            FESystem::Numerics::VectorBase<FESystemDouble> vector1(n_elems), vector2(n_elems);
//            this->getTemperatureDependentFluxCorrectionVector(vector2, AB_inv);
//            AB_inv.getDiagonal(vector1);
//            vector2.pointwiseMultiply(vector1);
//            vector2.pointwiseMultiply(temp_dep_vector);
//            vector2.add(1.0, temp_dep_vector);
//            AB_inv.rightDiagonalScale(temp_dep_vector);
//            
//            AB_inv.leftMultiplyTranspose(G_mat);
//            elem_temp_matrix.leftMultiply(AB_inv);
//        }
//        
//        jacobian.copy(elem_temp_matrix);
//    }
//    
//    this->fesystem_controller.performance_logging->unsetEvent
//    ("FESystem::Raditaion::RadiationCavityAnalysis::getFEJacobianMatrix()",
//     this->performance_logging_header);
//}
//
//
//
//void
//FESystem::Raditaion::RadiationCavityAnalysis::getTemperatureDependentFluxCorrectionVector
//(FESystem::Numerics::VectorBase<FESystemDouble>& correction_vec,
// FESystem::Numerics::MatrixBase<FESystemDouble>& AB_inv_mat)
//{
//    unsigned int n_elems = this->radiation_cavity.nRadiationElems();
//    
//    Assert(correction_vec.size() == n_elems, ExcInternalError());
//    Assert(AB_inv_mat.m() == n_elems, ExcInternalError());
//    Assert(AB_inv_mat.n() == n_elems, ExcInternalError());
//    
//    correction_vec.zero();
//    
//    // calculate the epsilon vector and store it
//    FESystem::Numerics::VectorBase<FESystemDouble> vec1(n_elems), vec2(n_elems);
//    
//    std::string name;
//    std::ostringstream cavity_id;
//    cavity_id << this->radiation_cavity.getCavityID();
//    
//    name.clear();
//    name = "BaseEpsilonVector";
//    name += "_Cavity_";
//    name += cavity_id.str();
//    
//    FESystemDatabase::TimeIndependentDataInfo time_independent_data_info;
//    time_independent_data_info.setDisciplineEnumID(Discipline::RADIATION_DISCIPLINE::num());
//    time_independent_data_info.setLoadCase(this->thermal_discipline.getAnalysisDriver().
//                                           getCurrentLoadCase());
//    time_independent_data_info.setName(name);
//    
//    this->fesystem_controller.global_data_storage->fillVector(time_independent_data_info,
//                                                              vec1);
//    this->calculateEmissivityVector(&vec2);
//    
//    // next, calculate the gamma vector
//    vec2.add(-1.0, vec1);
//    double base_eps = 0.0, delta_eps = 0.0, val = 0.0;
//    for (unsigned int i=0; i<n_elems; i++)
//    {
//        base_eps = vec1.el(i);
//        delta_eps = vec2.el(i);
//        
//        val = -delta_eps/base_eps/(base_eps + delta_eps);
//        vec2.set(i, val);
//    }
//    double sigma =
//    this->fesystem_controller.analysis_case->getRealParameter("S_B_CONSTANT");
//    // store a copy of vec2 for future use
//    vec2.scale(1.0/sigma);
//    
//    
//    AB_inv_mat.getDiagonal(vec1);
//    
//    // get the number of terms to use in the taylors series for approximation
//    const unsigned int n_taylor_terms =
//    this->radiation_cavity.getNTaylorTermsForABInvApproximation();
//    
//    // now iterate over the terms, and calculate the correction vector
//    for (unsigned int i=0; i < n_taylor_terms; i++)
//        for (unsigned int j=0; j<n_elems; j++)
//        {
//            val = pow(-vec2.el(j)*vec1.el(j),i+1) / vec1.el(j);
//            correction_vec.add(j,val);
//        }
//}
//
//
//
//void
//FESystem::Raditaion::RadiationCavityAnalysis::calculateB_TemperatureDependentJacobianContribution
//(FESystem::Numerics::VectorBase<FESystemDouble>& vector,
// const FESystem::Numerics::VectorBase<FESystemDouble>& rad_temp_vector,
// FESystem::Numerics::MatrixBase<FESystemDouble>& A_inv_mat,
// FESystem::Numerics::MatrixBase<FESystemDouble>& AB_inv_mat)
//{
//    unsigned int  n_elems = this->radiation_cavity.nRadiationElems();
//    
//    // create some temporary set of values for the procedure below
//    if (vector.size() != n_elems)
//        vector.resize(n_elems);
//    Assert(A_inv_mat.m() == n_elems &&
//           A_inv_mat.n() == n_elems, ExcInternalError());
//    Assert(AB_inv_mat.m() == n_elems &&
//           AB_inv_mat.n() == n_elems, ExcInternalError());
//    
//    vector.zero();
//    
//    FESystem::Numerics::VectorBase<FESystemDouble> vector1(n_elems), vector2(n_elems), vector3(n_elems);
//    vector1.zero(); vector2.zero(); vector3.zero();
//    
//    FESystem::Numerics::MatrixBase<FESystemDouble> F_matrix(n_elems, n_elems);
//    this->getMatrixQuantity(F_matrix, FESystem::Raditaion::RadiationCavityAnalysis::F_MATRIX);
//    
//    // first, get a vector for [B]^(-1) * {(T_rad + T_a)^4}
//    AB_inv_mat.rightMultiplyVector(rad_temp_vector, vector2);
//    
//    // if the cavity AB_inv matrix is being approximated, then calculated
//    // vector needs to be corrected
//    if(this->radiation_cavity.ifApproximateABInvFactor())
//    {
//        this->getTemperatureDependentFluxCorrectionVector(vector1, AB_inv_mat);
//        vector1.pointwiseMultiply(vector2);
//        AB_inv_mat.rightMultiplyVector(vector1, vector3);
//        vector2.add(1.0, vector3);
//    }
//    
//    A_inv_mat.rightMultiplyVector(vector2, vector1);
//    
//    
//    Property::PropertyDatabase& property_database =
//    *(this->fesystem_controller.property_database.get());
//    MeshDS::FEMeshData& mesh_data = *(this->radiation_cavity.getMeshData());
//    
//    unsigned int fe_elem_ID, property_ID;
//    double factor;
//    
//    // iterate over each radiation element, and calculate the contribution
//    // this loop calculates the term [\frac{\partial B_{ij}}{\partial T_{R_k}}\psi_j]_diag
//    for (unsigned int i=0; i < n_elems; i++)
//    {
//        // get the ID of the finite elem, get the property card for this elem
//        fe_elem_ID = this->radiation_cavity.getElemFromID(i).FeElementID();
//        const Elem* elem = mesh_data.getElemFromForeignID(fe_elem_ID);
//        property_ID = mesh_data.getElemPropertyID(elem);
//        ElemDataCard& elem_data_card = property_database.getElemDataCardFromID(property_ID);
//        
//        // if the property is temperature dependent, then initialize the property card
//        // at the average temperature of the finite element to which this radiation
//        // element belongs
//        this->initializePropertyCardAtElementTemperature(elem_data_card, i);
//        
//        // this is the term (\frac{-1}{epsilon^2_{R_i}}
//        //  \frac{\partial \epsilon_{R_i}(T_{R_i})}{\partial T_{R_k}} \delta_{ij} \delta_{ik})
//        // sensitivities for both factor1 and factor2 are the same, hence, only one is
//        // obtained
//        elem_data_card.getFactorSensitivityForLocalParameter(factor,
//                                                             RADIATION_EPSILON_FACTOR_1::num(),
//                                                             Property::TEMPERATURE::num());
//        
//        F_matrix.getRow(i, vector2);
//        vector2.scale(-factor);
//        vector2.add(i, factor);
//        vector.set(i, vector1.dot(vector2));
//        
//        
//        // once this card has been used, clear the initialization
//        this->clearPropertyCardInitialization(elem_data_card);
//    }
//    
//}
//
//
//
//
//bool
//FESystem::Raditaion::RadiationCavityAnalysis::checkIfTemperatureDependent(RadiationQty quantity)
//{
//    switch (quantity)
//    {
//        case  FESystem::Raditaion::RadiationCavityAnalysis::F_MATRIX:
//        case  FESystem::Raditaion::RadiationCavityAnalysis::F_MATRIX_SENSITIVITY:
//            return false;
//            break;
//			
//        case  FESystem::Raditaion::RadiationCavityAnalysis::A_MATRIX:
//        case  FESystem::Raditaion::RadiationCavityAnalysis::A_MATRIX_SENSITIVITY:
//            return false;
//            break;
//			
//        case  FESystem::Raditaion::RadiationCavityAnalysis::A_INVERSE_MATRIX:
//            return false;
//            break;
//            
//        case  FESystem::Raditaion::RadiationCavityAnalysis::G_MATRIX:
//        case  FESystem::Raditaion::RadiationCavityAnalysis::G_MATRIX_SENSITIVITY:
//            return false;
//            break;
//            
//        case  FESystem::Raditaion::RadiationCavityAnalysis::B_MATRIX:
//        case  FESystem::Raditaion::RadiationCavityAnalysis::B_MATRIX_SENSITIVITY:
//            return true;
//            break;
//			
//        case  FESystem::Raditaion::RadiationCavityAnalysis::B_INVERSE_MATRIX:
//            return true;
//            break;
//			
//        case  FESystem::Raditaion::RadiationCavityAnalysis::FACTOR_MATRIX:
//            return false; // false, since it is now defined only for temperature independent cases
//            break;
//            
//        case  FESystem::Raditaion::RadiationCavityAnalysis::AB_INV_MATRIX:
//            return true; // this is defined only for temperature dependent cases
//            break;
//            
//        default:
//            abort();
//            break;
//    }
//    
//    // should never get here.
//    Assert(false, ExcInternalError());
//    return false;
//}
//
//
//
//std::auto_ptr<FESystemDatabase::DataInfoBase>
//FESystem::Raditaion::RadiationCavityAnalysis::getDataInfoForQty(RadiationQty quantity)
//{
//    std::auto_ptr<FESystemDatabase::DataInfoBase> data_info(NULL);
//    std::ostringstream rad_num;
//    rad_num << this->radiation_cavity.getCavityID();
//    std::string name;
//    
//    switch (quantity)
//    {
//        case  FESystem::Raditaion::RadiationCavityAnalysis::F_MATRIX:
//        {
//            name = "F_matrix";
//            name += "_Cavity_";
//            name += rad_num.str();
//            name += "_";
//            data_info.reset(new FESystemDatabase::TimeIndependentDataInfo);
//            data_info->setDisciplineEnumID(Discipline::RADIATION_DISCIPLINE::num());
//            data_info->setName(name);
//        }
//            break;
//			
//        case  FESystem::Raditaion::RadiationCavityAnalysis::F_MATRIX_SENSITIVITY:
//        {
//            data_info.reset(this->getDataInfoForQty(FESystem::Raditaion::RadiationCavityAnalysis::F_MATRIX).release());
//            data_info->setDVID(this->thermal_discipline.getAnalysisDriver().getCurrentDesignParameter().getID());
//        }
//            break;
//            
//        case  FESystem::Raditaion::RadiationCavityAnalysis::A_MATRIX:
//        {
//            name = "A_matrix";
//            name += "_Cavity_";
//            name += rad_num.str();
//            name += "_";
//            data_info.reset(new FESystemDatabase::TimeIndependentDataInfo);
//            data_info->setDisciplineEnumID(Discipline::RADIATION_DISCIPLINE::num());
//            data_info->setName(name);
//        }
//            break;
//            
//        case  FESystem::Raditaion::RadiationCavityAnalysis::A_INVERSE_MATRIX:
//        {
//            name = "A_Inverse_matrix";
//            name += "_Cavity_";
//            name += rad_num.str();
//            name += "_";
//            data_info.reset(new FESystemDatabase::TimeIndependentDataInfo);
//            data_info->setDisciplineEnumID(Discipline::RADIATION_DISCIPLINE::num());
//            data_info->setName(name);
//        }
//            break;
//            
//        case  FESystem::Raditaion::RadiationCavityAnalysis::A_MATRIX_SENSITIVITY:
//        {
//            data_info.reset(this->getDataInfoForQty(FESystem::Raditaion::RadiationCavityAnalysis::A_MATRIX).release());
//            data_info->setDVID(this->thermal_discipline.getAnalysisDriver().getCurrentDesignParameter().getID());
//        }
//            break;
//            
//        case  FESystem::Raditaion::RadiationCavityAnalysis::G_MATRIX:
//        {
//            name = "G_matrix";
//            name += "_Cavity_";
//            name += rad_num.str();
//            name += "_";
//            data_info.reset(new FESystemDatabase::TimeIndependentDataInfo);
//            data_info->setDisciplineEnumID(Discipline::RADIATION_DISCIPLINE::num());
//            data_info->setName(name);
//        }
//            break;
//            
//        case  FESystem::Raditaion::RadiationCavityAnalysis::G_MATRIX_SENSITIVITY:
//        {
//            data_info.reset(this->getDataInfoForQty(FESystem::Raditaion::RadiationCavityAnalysis::G_MATRIX).release());
//            data_info->setDVID(this->thermal_discipline.getAnalysisDriver().getCurrentDesignParameter().getID());
//        }
//            break;
//            
//        case  FESystem::Raditaion::RadiationCavityAnalysis::B_MATRIX:
//        {
//            switch(this->thermal_discipline.getSolution().getDisciplineTransientNatureEnumID
//                   (Discipline::THERMAL_DISCIPLINE::num()))
//            {
//                case STEADY_STATE_SOLUTION_ENUM_ID:
//                {
//                    name = "B_matrix";
//                    name += "_Cavity_";
//                    name += rad_num.str();
//                    name += "_";
//                    data_info.reset(new FESystemDatabase::TimeIndependentDataInfo);
//                    data_info->setDisciplineEnumID(Discipline::RADIATION_DISCIPLINE::num());
//                    data_info->setLoadCase(this->thermal_discipline.getAnalysisDriver().
//                                           getCurrentLoadCase());
//                    data_info->setName(name);
//                }
//                    break;
//                    
//                case TRANSIENT_SOLUTION_ENUM_ID:
//                {
//                    name = "B_matrix";
//                    name += "_Cavity_";
//                    name += rad_num.str();
//                    name += "_";
//                    Driver::TransientAnalysisDriver& transient_driver =
//                    dynamic_cast<Driver::TransientAnalysisDriver&>
//                    (this->thermal_discipline.getAnalysisDriver());
//                    
//                    data_info.reset(new FESystemDatabase::TimeDependentDataInfo);
//                    data_info->setDisciplineEnumID(Discipline::RADIATION_DISCIPLINE::num());
//                    data_info->setLoadCase(transient_driver.getCurrentLoadCase());
//                    dynamic_cast<FESystemDatabase::TimeDependentDataInfo*>
//                    (data_info.get())->setTransientIterationInfo(transient_driver.getCurrentIterationNumber(),
//                                                                 transient_driver.getCurrentAnalysisTime());
//                    dynamic_cast<FESystemDatabase::TimeDependentDataInfo*>(data_info.get())->setOrder(0);
//                    data_info->setName(name);
//                }
//                    break;
//                    
//                default:
//                    Assert(false, ExcInternalError());
//                    break;
//            }
//        }
//            break;
//            
//        case  FESystem::Raditaion::RadiationCavityAnalysis::B_MATRIX_SENSITIVITY:
//        {
//            data_info.reset(this->getDataInfoForQty(FESystem::Raditaion::RadiationCavityAnalysis::B_MATRIX).release());
//            data_info->setDVID(this->thermal_discipline.getAnalysisDriver().getCurrentDesignParameter().getID());
//        }
//            break;
//            
//        case  FESystem::Raditaion::RadiationCavityAnalysis::B_INVERSE_MATRIX:
//        {
//            switch(this->thermal_discipline.getSolution().getDisciplineTransientNatureEnumID
//                   (Discipline::THERMAL_DISCIPLINE::num()))
//            {
//                case STEADY_STATE_SOLUTION_ENUM_ID:
//                {
//                    name = "B_Inverse_Matrix";
//                    name += "_Cavity_";
//                    name += rad_num.str();
//                    name += "_";
//                    data_info.reset(new FESystemDatabase::TimeIndependentDataInfo);
//                    data_info->setDisciplineEnumID(Discipline::RADIATION_DISCIPLINE::num());
//                    data_info->setLoadCase(this->thermal_discipline.getAnalysisDriver().
//                                           getCurrentLoadCase());
//                    data_info->setName(name);
//                }
//                    break;
//                    
//                case TRANSIENT_SOLUTION_ENUM_ID:
//                {
//                    name = "B_Inverse_Matrix";
//                    name += "_Cavity_";
//                    name += rad_num.str();
//                    name += "_";
//                    Driver::TransientAnalysisDriver& transient_driver =
//                    dynamic_cast<Driver::TransientAnalysisDriver&>
//                    (this->thermal_discipline.getAnalysisDriver());
//                    
//                    data_info.reset(new FESystemDatabase::TimeDependentDataInfo);
//                    data_info->setDisciplineEnumID(Discipline::RADIATION_DISCIPLINE::num());
//                    data_info->setLoadCase(transient_driver.getCurrentLoadCase());
//                    dynamic_cast<FESystemDatabase::TimeDependentDataInfo*>
//                    (data_info.get())->setTransientIterationInfo(transient_driver.getCurrentIterationNumber(),
//                                                                 transient_driver.getCurrentAnalysisTime());
//                    dynamic_cast<FESystemDatabase::TimeDependentDataInfo*>(data_info.get())->setOrder(0);
//                    data_info->setName(name);
//                }
//                    break;
//                    
//                default:
//                    Assert(false, ExcInternalError());
//                    break;
//            }
//        }
//            break;
//			
//        case  FESystem::Raditaion::RadiationCavityAnalysis::FACTOR_MATRIX:
//        {
//            switch(this->thermal_discipline.getSolution().getDisciplineTransientNatureEnumID
//                   (Discipline::THERMAL_DISCIPLINE::num()))
//            {
//                case STEADY_STATE_SOLUTION_ENUM_ID:
//                {
//                    name = "Factor_matrix";
//                    name += "_Cavity_";
//                    name += rad_num.str();
//                    name += "_";
//                    data_info.reset(new FESystemDatabase::TimeIndependentDataInfo);
//                    data_info->setDisciplineEnumID(Discipline::RADIATION_DISCIPLINE::num());
//                    data_info->setLoadCase(this->thermal_discipline.getAnalysisDriver().
//                                           getCurrentLoadCase());
//                    data_info->setName(name);
//                }
//                    break;
//                    
//                case TRANSIENT_SOLUTION_ENUM_ID:
//                {
//                    Driver::TransientAnalysisDriver& transient_driver =
//                    dynamic_cast<Driver::TransientAnalysisDriver&>
//                    (this->thermal_discipline.getAnalysisDriver());
//                    
//                    name = "Factor_matrix";
//                    name += "_Cavity_";
//                    name += rad_num.str();
//                    name += "_";
//                    data_info.reset(new FESystemDatabase::TimeDependentDataInfo);
//                    data_info->setDisciplineEnumID(Discipline::RADIATION_DISCIPLINE::num());
//                    data_info->setLoadCase(transient_driver.getCurrentLoadCase());
//                    dynamic_cast<FESystemDatabase::TimeDependentDataInfo*>
//                    (data_info.get())->setTransientIterationInfo(transient_driver.getCurrentIterationNumber(),
//                                                                 transient_driver.getCurrentAnalysisTime());
//                    dynamic_cast<FESystemDatabase::TimeDependentDataInfo*>(data_info.get())->setOrder(0);
//                    data_info->setName(name);
//                }
//                    break;
//                    
//                default:
//                    Assert(false, ExcInternalError());
//                    break;
//            }
//        }
//            break;
//            
//        case  FESystem::Raditaion::RadiationCavityAnalysis::AB_INV_MATRIX:
//        {
//            switch(this->thermal_discipline.getSolution().getDisciplineTransientNatureEnumID
//                   (Discipline::THERMAL_DISCIPLINE::num()))
//            {
//                case STEADY_STATE_SOLUTION_ENUM_ID:
//                {
//                    name = "AB_Inv_matrix";
//                    name += "_Cavity_";
//                    name += rad_num.str();
//                    name += "_";
//                    data_info.reset(new FESystemDatabase::TimeIndependentDataInfo);
//                    data_info->setDisciplineEnumID(Discipline::RADIATION_DISCIPLINE::num());
//                    data_info->setLoadCase(this->thermal_discipline.getAnalysisDriver().
//                                           getCurrentLoadCase());
//                    data_info->setName(name);
//                }
//                    break;
//                    
//                case TRANSIENT_SOLUTION_ENUM_ID:
//                {
//                    Driver::TransientAnalysisDriver& transient_driver =
//                    dynamic_cast<Driver::TransientAnalysisDriver&>
//                    (this->thermal_discipline.getAnalysisDriver());
//                    
//                    name = "AB_Inv_matrix";
//                    name += "_Cavity_";
//                    name += rad_num.str();
//                    name += "_";
//                    data_info.reset(new FESystemDatabase::TimeDependentDataInfo);
//                    data_info->setDisciplineEnumID(Discipline::RADIATION_DISCIPLINE::num());
//                    data_info->setLoadCase(transient_driver.getCurrentLoadCase());
//                    dynamic_cast<FESystemDatabase::TimeDependentDataInfo*>
//                    (data_info.get())->setTransientIterationInfo(transient_driver.getCurrentIterationNumber(),
//                                                                 transient_driver.getCurrentAnalysisTime());
//                    dynamic_cast<FESystemDatabase::TimeDependentDataInfo*>(data_info.get())->setOrder(0);
//                    data_info->setName(name);
//                }
//                    break;
//                    
//                default:
//                    Assert(false, ExcInternalError());
//                    break;
//            }
//        }
//            break;
//            
//            
//        default:
//            Assert(false, ExcInternalError());
//            break;
//    }
//    
//    return data_info;
//}
//
//
//
//
//
//void
//FESystem::Raditaion::RadiationCavityAnalysis::initializePropertyCardAtElementTemperature
//(ElemDataCard& elem_data_card,
// const unsigned int rad_elem_ID)
//{
//    // now initialize the property card at the specified value.
//    // if the properties are temperature dependent, then initialize the
//    // element property cards at the element average temperature
//    if (this->thermal_discipline.checkPropertyDependenceOnTemperature())
//    {
//        Assert(this->fe_nodal_temp_vector != NULL, ExcEmptyObject());
//        
//        // calculate the average of the nodal temperature of the finite element
//        // to which the specified radiation element belongs.
//        static std::vector<unsigned int> fe_node_ID_vec;
//        fe_node_ID_vec.clear();
//        
//        // get the vector of internal node IDs
//        this->radiation_cavity.getInternalNodeIDForRadElem
//        (rad_elem_ID, fe_node_ID_vec);
//        
//        double value = 0.0;
//        for (unsigned int i=0; i < fe_node_ID_vec.size(); i++)
//            value += this->fe_nodal_temp_vector->el(fe_node_ID_vec[i]);
//        
//        value /= (1.0 * fe_node_ID_vec.size());
//        
//        static std::map<unsigned int, double> local_param_map;
//        
//        local_param_map[Property::TEMPERATURE::num()] = value;
//        elem_data_card.reinitElemAndMaterialCardForLocalParameters(&local_param_map);
//    }
//    else
//    {
//        elem_data_card.reinitElemAndMaterialCardForLocalParameters();
//    }
//    
//}
//
//
//void
//FESystem::Raditaion::RadiationCavityAnalysis::clearPropertyCardInitialization(ElemDataCard& elem_data_card)
//{
//    static std::vector<unsigned int> param_vec;
//    static std::vector<unsigned int> no_param_vec;
//    
//    if (param_vec.size() == 0)
//        param_vec.push_back(Property::TEMPERATURE::num());
//    
//    // if the properties are temperature dependent, then initialize the
//    // element property cards at the element average temperature
//    if(this->thermal_discipline.checkPropertyDependenceOnTemperature())
//        elem_data_card.partialClearElemAndMaterialCardLocalParameterInitialization
//        (param_vec);
//    else
//        elem_data_card.partialClearElemAndMaterialCardLocalParameterInitialization
//        (no_param_vec);
//}
//
