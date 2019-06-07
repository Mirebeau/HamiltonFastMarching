// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Copyright 2017 Erik Bekkers, Eindhoven University of Technology, the Netherlands
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef MathematicaIOExternC_h
#define MathematicaIOExternC_h


#define MathematicaTryCatch(nArgs,code) \
if(nArgs!=Argc) { \
IO::WarnMsg() << "Library call error : number of arguments " << Argc << " differs from expected " << nArgs << ".\n"; \
return LIBRARY_TYPE_ERROR;} \
try{code; \
return LIBRARY_NO_ERROR; \
} catch (const std::logic_error & e) { \
IO::WarnMsg() << "Exception caught. " << e.what(); \
return LIBRARY_FUNCTION_ERROR; \
}




// All functions that interact with Wolfram Mathematica are included in -extern "C"
extern "C"
{//Begin of Extern "C"

    DLLEXPORT bool HasField(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
        MathematicaTryCatch(1,
                            io.SetWolframLibraryData(libData);
                            const std::string key = MArgument_getUTF8String(Args[0]);
                            MArgument_setBoolean(Res,io.HasField(key));
                            )
        
    }
    
    DLLEXPORT bool EraseField(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
        MathematicaTryCatch(1,
                            io.SetWolframLibraryData(libData);
                            const std::string key = MArgument_getUTF8String(Args[0]);
                            MArgument_setBoolean(Res,io.EraseField(key));
                            )
    }
    
	  /* ***************************************************************************** */
	  /* ********************* All the set/get functions ***************************** */
	  /* ***************************************************************************** */

	  // GetScalar/SetScalar
	DLLEXPORT int GetScalar(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
        MathematicaTryCatch(1,
                            io.SetWolframLibraryData(libData);
                            const std::string key = MArgument_getUTF8String(Args[0]);
                            MArgument_setReal(Res, io.Get<double>(key));
                            )
    }
    
	DLLEXPORT int SetScalar(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
        MathematicaTryCatch(2,
                            io.SetWolframLibraryData(libData);
                            const std::string key = MArgument_getUTF8String(Args[0]);
                            const double val = MArgument_getReal(Args[1]);
                            io.Set<double>(key, val);
                            )
    }

	// GetString/SetString
	DLLEXPORT int GetString(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)
	{
        MathematicaTryCatch(1,
                            io.SetWolframLibraryData(libData);
                            const std::string key = MArgument_getUTF8String(Args[0]);
                            const std::string stringSTD = io.GetString(key);
                            MArgument_setUTF8String(Res, (char*)stringSTD.c_str());
        )
//			char * string = new char[stringSTD.size() + 1];
//			std::copy(stringSTD.begin(), stringSTD.end(), string);
	}
	DLLEXPORT int SetString(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
        MathematicaTryCatch(2,
                            io.SetWolframLibraryData(libData);
                            const std::string key = MArgument_getUTF8String(Args[0]);
                            const std::string val = MArgument_getUTF8String(Args[1]);
                            io.SetString(key, val);
                            )
    }

	// GetVector/SetVector
	DLLEXPORT int GetVector(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
        MathematicaTryCatch(1,
                            io.SetWolframLibraryData(libData);
                            const std::string key = MArgument_getUTF8String(Args[0]);
                            MArgument_setMTensor(Res, io.MathGetArray(key,1));
                            )
	}
    
	DLLEXPORT int SetVector(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)
	{
        MathematicaTryCatch(2,
                            io.SetWolframLibraryData(libData);
                            const std::string key = MArgument_getUTF8String(Args[0]);
                            MTensor val = MArgument_getMTensor(Args[1]);
                            io.MathSetArray(key, val);
                            )
	}

	// GetArray/SetArray
	DLLEXPORT int GetArray(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
        MathematicaTryCatch(2,
                            io.SetWolframLibraryData(libData);
                            const std::string key = MArgument_getUTF8String(Args[0]);
                            const int ndims = MArgument_getInteger(Args[1]); // Array dimension
                            MArgument_setMTensor(Res,io.MathGetArray(key,ndims));
                            )
	}
    DLLEXPORT int SetArray(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
        MathematicaTryCatch(2,
                            io.SetWolframLibraryData(libData);
                            const std::string key = MArgument_getUTF8String(Args[0]);
                            const MTensor val = MArgument_getMTensor(Args[1]);
                            io.MathSetArray(key, val);
                            )
    }
}//End of Extern "C"


#endif /* MathematicaIOExternC_h */
