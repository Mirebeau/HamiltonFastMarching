//
//  main2D.cpp
//  
//
//  Created by Jean-Marie Mirebeau on 14/09/2016.
//
//


#include "Output/ExportMacros.h"
#include "Output/MexIO.h"
typedef IO_<Mex::BaseIO> IO;
typedef typename IO::Msg Msg;
typedef typename IO::WarnMsg WarnMsg;
#include "DispatchAndRun.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[] ){
    if(nrhs!=1 || nlhs!=1){
        IO::WarnMsg() << "Exactly one input and one output are expected (structures).";
        return;
    }
    
    try {
        try {
            IO io(prhs[0],plhs);
            io.arrayOrdering = ArrayOrdering::Transposed;
            Run(io);
        } catch (const std::exception & e) {
            IO::WarnMsg() << "Lifted Fast Marching exception.\n " << e.what() << "\n";
        }
    } catch(...){
        IO::WarnMsg() << "MatlabLFM error : unspecified exception caught.\n";
    }
}
