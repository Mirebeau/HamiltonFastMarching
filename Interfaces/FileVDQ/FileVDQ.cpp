// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#define Voronoi6 // Uncomment for 6 dimensional reduction

#include "JMM_CPPLibs/Macros/ExportArrow.h"
#include "JMM_CPPLibs/Output/FileIO.h"
typedef IO_<FileIO> IO;
typedef typename IO::Msg Msg;
typedef typename IO::WarnMsg WarnMsg;
#include "DispatchAndRun.h"

/* // Set directory under windows
#include <direct.h>
char path[] = "C:/Users/Jean-Marie Mirebeau/Documents/Cpp/FileVDE/Debug";
_chdir(path);
*/


int main(int argc, const char * argv[]) { 
    std::string inputPrefix, outputPrefix;
    
    inputPrefix  = argc > 1 ? argv[1] : "input";
    outputPrefix = argc > 2 ? argv[2] : "output";
    
    try {
        IO io(inputPrefix,outputPrefix);
        Run(io);
    } catch(const std::logic_error & e) {
        IO::WarnMsg() << "Exception caught. " << e.what();
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
