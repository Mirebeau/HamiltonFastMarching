// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef FileIO_h
#define FileIO_h

#include <fstream>
#include <sstream>
#include "BaseIO.h"
/* 
 A simple interface for importing/exporting strings and numeric arrays (with double type) to files.
 
 Mostly used as an interface with Mathematica(R) and Python. (Direct, library based interfaces are also available.)
 */

struct FileIO : BaseIO {
    FileIO(std::string inputPrefix, std::string outputPrefix_);
    ~FileIO();
protected:
    const std::string outputPrefix;
};
    
#include "FileIO.hxx"

#endif /* FileIO_h */
