// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Copyright 2017 Erik Bekkers, Eindhoven University of Technology, the Netherlands
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef MathematicaIO_h
#define MathematicaIO_h


//#include "MathLink.h"
#include "WolframLibrary.h"
#include <map>
#include <fstream>
#include <sstream>
#include "BaseIO.h"

/*
A simple interface for making Mathematica data available to the c++ code. 
*/

struct MathematicaIO : BaseIO {
    typedef double ScalarType;
    
    // An adaptation of the file IO:
    MathematicaIO(WolframLibraryData libDataIn){SetWolframLibraryData(libDataIn);arrayOrdering=TraitsIO::ArrayOrdering::RowMajor;}
    
    void MathSetArray(KeyCRef,const MTensor);
    MTensor MathGetArray(KeyCRef,int d) const;
    std::string GetComputedKeys() const;
    
    void SetWolframLibraryData(WolframLibraryData libDataIn) { libData=libDataIn; latestLibData=libDataIn;};
protected:
    WolframLibraryData libData = nullptr;// libData is necessary in the type conversions
    static WolframLibraryData latestLibData;
    
    static void StaticSendMsg(bool warn, const std::string & msg, WolframLibraryData libData);
    static void StaticSendMsg(bool warn, const std::string & msg) {
        return StaticSendMsg(warn,msg,latestLibData);}
    virtual void SendMsg(bool warn, const std::string & msg) const override {
        return StaticSendMsg(warn,msg,libData);}
    template<bool,typename> friend struct _Msg;
};
WolframLibraryData MathematicaIO::latestLibData = nullptr;


#include "MathematicaIO.hxx"

#endif /* MathematicaIO_h */
