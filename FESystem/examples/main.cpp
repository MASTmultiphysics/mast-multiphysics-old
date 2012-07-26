
// C++ includes
#include <iostream>
#include <memory>
#include <string>

// FESystem includes
#include "Base/FESystemBase.h"
#include "Base/FESystemTypes.h"
#include "Utils/CommandLineArguments.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "Field/LocalField.h"
#include "Surrogates/LeastSquareSurrogate.h"
#include "Surrogates/LeastSquareSurrogate.h"

// Plplot includes
//#include "plstream.h"

int old_main (int argc, char * const argv[]) {
    
	std::auto_ptr<FESystem::Base::InitializerBase> 
	fes_library(FESystem::Base::createInitializer(FESystem::Base::FESystem)); 
	
	
	FESystem::Utility::CommandLineArguments cargs;
	std::string str;
	str = "-np 1 ";
	cargs.appendArgument(str);
    
	fes_library->setCommandArguments(cargs);
	fes_library->initialize();
    
	FESystem::Numerics::LocalVector<FESystemDouble> vec1, vec2;
    FESystem::Field::LocalField<FESystemDouble> field;
    FESystem::Surrogates::LeastSquareSurrogate<FESystemDouble> ls_surr;
    
    FESystemUInt n_params = 2, n_resps = 1, n_points = 3;
    
    field.setDimensions(n_points*n_points, n_params, n_resps);
    
    vec1.resize(n_params);
    vec2.resize(n_resps);
	
    FESystemUInt npt = 0;
    for (FESystemUInt i=0; i<n_points; i++)
        for (FESystemUInt j=0; j<n_points; j++)
        {
            vec1.setVal(0,1.0*i);
            vec1.setVal(1,1.0*j);
            
            vec2.setVal(0,1.0*(i*i+j*j));
            field.addParameterAndResponses(npt++, vec1, vec2);
        }
    
    ls_surr.setField(field);
    ls_surr.setPolynomialType(FESystem::Surrogates::SECOND_ORDER_COMPLETE);
    ls_surr.initialize();
    
    for (FESystemUInt i=0; i<n_points*n_points; i++)
    {
        field.getParametersAtPoint(i, vec1);
        field.getResponsesAtPoint(i, vec2);
        
        std::cout << "Point: ";
        vec1.write(std::cout);
        std::cout << "\n";
        
        std::cout << "Actual Response: " << vec2.getVal(0) << "\n";
        std::cout << "Predicted Response: " << ls_surr.getEstimate(0, vec1) << "\n";    
    }
    
    //
    // plotting
    //
    //  std::auto_ptr<plstream> plot(new plstream());
    //  
    //  plot->sdev("xwin");
    //  plot->init();
    //  plot->adv(0);
    //  plot->vsta();
    //  plot->prec(1,1);
    //  plot->wind(t0,t_final,-40.0,40.0);
    //  plot->col0(1);
    //  plot->box("bcgnst",2,2,"bcgnst",0,0);
    //  //  plot->col0(3);
    //  //  plot->line(n_steps,&(x[0]),&(y[0]));
    //  //  plot->col0(4);
    //  //  plot->line(n_steps,&(x[0]),&(y[n_steps+1]));
    //  //  plot->col0(5);
    //  //  plot->line(n_steps,&(x[0]),&(p[0]));
    //  //  plot->col0(6);
    //  //  plot->line(n_steps,&(x[0]),&(p[n_steps+1]));
    //  //  plot->col0(7);
    //  //  plot->line(n_steps,&(x[0]),&(p[2*n_steps+2]));
    //  plot->col0(8);
    //  plot->line(n_steps,&(x[0]),&(p[3*n_steps+3]));
    //  plot->lab("x values","y values","title of the plot");
    
    
    return 0;
}
