////
////  RadiationElementPair.h
////  FESystem
////
////  Created by Manav Bhatia on 1/25/13.
////
////
//
//#ifndef __FESystem__RadiationElementPair__
//#define __FESystem__RadiationElementPair__
//
//// FESystem includes
//#include "Base/FESystemTypes.h"
//
//
//// forward declerations
//
//
//namespace FESystem
//{
//    
//    namespace Mesh {class FaceElemBase;}
//    
//    namespace Radiation
//    {
//
//        /*!
//         *   enum describing the method to be used
//         */
//        enum ShapeFactorMethod
//        {
//            DOUBLE_AREA_INT,
//            CONTOUR_INT,
//            MITALAS_CONTOUR_INT
//        };
//
//        
//        class RadiationElementPair
//        {
//        public:
//            
//            /*!
//             *   constructor
//             */
//            RadiationElementPair();
//            
//            /*!
//             *  destructor
//             */
//            ~RadiationElementPair();
//            
//            /*!
//             *  initializes the pair with two elements
//             */
//            void reinit(const FESystem::Mesh::FaceElemBase& el1, const FESystemInt orient1, const FESystem::Mesh::FaceElemBase& el2, const FESystemInt orient2);
//            
//            /*!
//             *  clears the initialization
//             */
//            void clear();
//            
//            /*!
//             *   returns the shape factors between the elements as a pair
//             */
//            void getShapeFactor(ShapeFactorMethod method, FESystemDouble& f1, FESystemDouble& f2);
//            
//            
//        protected:
//            
//            /*!
//             *   calculates and returns the double area integration shape factor
//             */
//            void doubleAreaIntShapeFactor(FESystemDouble& f1, FESystemDouble& f2);
//            
//            /*!
//             *   calculates and returns the contour integration shape factor
//             */
//            void contourIntShapeFactor(FESystemDouble& f1, FESystemDouble& f2);
//            
//            /*!
//             *   calculates and returns the shape factors using mitalas' method of intgration
//             */
//            void mitalasContourIntShapeFactor(FESystemDouble& f1, FESystemDouble& f2);
//                        
//            /*!
//             *   boolean to keep track of whether the pair has been initialized
//             */
//            FESystemBoolean initialized;
//            
//            /*!
//             *   Radiation elems
//             */
//            const FESystem::Mesh::FaceElemBase *rad_elem1, *rad_elem2;
//            
//            
//            FESystemInt rad_elem1_orientation, rad_elem2_orientation;
//        };
//    }
//}
//
//
//
//
//#endif /* defined(__FESystem__RadiationElementPair__) */
