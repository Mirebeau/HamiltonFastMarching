// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef PythonIO_h
#define PythonIO_h

#include "BaseIO.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;


/*
A simple interface for making Python data available to the c++ code.
*/

struct PythonIO : BaseIO {
    // An adaptation of the file IO:
    PythonIO(){arrayOrdering = ArrayOrdering::YXZ_RowMajor;};
    typedef py::array_t<ScalarType, py::array::c_style | py::array::forcecast > ndarray;
    ndarray PyGetArray(KeyCRef) const;
    void PySetArray(KeyCRef, ndarray);
    std::string GetComputedKeys() const;
protected:
    static void StaticSendMsg(bool warn, std::string msg){
		py::print(warn ? ("Warning : "+msg) : msg, py::arg("end")="" );}
    void SendMsg(bool warn, const std::string & msg) const {return StaticSendMsg(warn,msg);}
    template<bool,typename> friend struct _Msg;
};


#include "PythonIO.hxx"

#endif /* PythonIO_h */
