//
//  analysis_driver.h
//  MAST
//
//  Created by Manav Bhatia on 11/25/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_analysis_driver_h__
#define __MAST_analysis_driver_h__

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"


namespace MAST {
    class AnalysisDriver {
    public:
        AnalysisDriver();
        
        virtual ~AnalysisDriver();
        
        
        void init(LibMeshInit& init, GetPot& infile,
                  int argc, char* const argv[]);
        
        
        void solve();
        
    protected:

        
        void init_geometry();
        
        
        void init_solver();
        
        
        void perform_amr();
        
    };
}


#endif // __MAST_analysis_driver_h__
