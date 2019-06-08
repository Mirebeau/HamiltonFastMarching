// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.
//
// Wolfram library link interface by Erik Bekers, Eindhoven University of Technology, Eindhoven, the Netherlands
// 2017-04-05








/* ***************************************************************************** */
/* **************************** Set up Mathematica IO IO *********************** */
/* ***************************************************************************** */

// Include main header 
#include "JMM_CPPLibs/Output/MathematicaIO.h"
// Define IO type
typedef IO_<MathematicaIO> IO;
typedef typename IO::Msg Msg;
typedef typename IO::WarnMsg WarnMsg;
// Make IO object (initialize with libData=NULL)
IO io(NULL);
//io.ArrayOrdering = ArrayOrdering::Reversed;
// Define the extern "C" functions that rely on io
#include "JMM_CPPLibs/Output/MathematicaIOExternC.h"


#include "DispatchAndRun.h"

extern "C"
{//Begin of Extern "C"
    
    
    EXTERN_C DLLEXPORT int RunModel(WolframLibraryData libData,
                                    mint Argc, MArgument *Args, MArgument Res) {
        MathematicaTryCatch(0,
                            io.SetWolframLibraryData(libData);
                            Run(io);
                            io.UsageReport();
                            )
    }
}






