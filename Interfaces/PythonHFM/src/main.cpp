#include <iostream>
#include <vector>
#include "JMM_CPPLibs/Macros/PPCat.h"
#include "JMM_CPPLibs/Output/PythonIO.h"
typedef IO_<PythonIO> IO;
Redeclare2Types(IO,Msg,WarnMsg)

#define HFMIO PPCAT(ModelName,_PY11) // Name of py11 exported class must be unique to model
void Run(IO &); // Defined in DispatchAndRun.h below, which also redefines ModelName

struct HFMIO {
	Redeclare3Types(IO,ScalarType,KeyCRef,ndarray);
	IO io;
	void Run(){::Run(io); io.UsageReport();}
	
	// Input output
	std::string GetString(KeyCRef key) const {return io.GetString(key);}
	void SetString(KeyCRef key, const std::string & val){io.SetString(key,val);}
	ScalarType GetScalar(KeyCRef key) const {return io.Get<ScalarType>(key);}
	void SetScalar(KeyCRef key, ScalarType val){io.Set<ScalarType>(key,val);}
	ndarray GetArray(KeyCRef key) const {return io.PyGetArray(key);}
	void SetArray(KeyCRef key, ndarray val) {io.PySetArray(key,val);}
	
	// Field inspection
	std::string GetComputedKeys() {return io.GetComputedKeys();}
	bool HasField(KeyCRef key) const {return io.HasField(key);}
	bool EraseField(KeyCRef key) {return io.EraseField(key);}
};

PYBIND11_MODULE(PythonModuleName, m){
	py::register_exception<JMMCppException>(m, "PyExp");
	
	m.doc() = "HFM python plugin"; // optional module docstring
	
	py::class_<HFMIO>(m, "HFMIO")
	.def(py::init<>())
	.def("run",&HFMIO::Run)
	.def("get_string",&HFMIO::GetString)
	.def("set_string",&HFMIO::SetString)
	.def("get_scalar",&HFMIO::GetScalar)
	.def("set_scalar",&HFMIO::SetScalar)
	.def("get_array", &HFMIO::GetArray)
	.def("set_array", &HFMIO::SetArray)
	.def("computed_keys", &HFMIO::GetComputedKeys)
	.def("has_field", &HFMIO::HasField)
	.def("erase_field",&HFMIO::EraseField)
	;
	    /*
    py::class_<HFMIO>(m, "HFMIO")
    .def(py::init<>())
    .def("Run",&HFMIO::Run)
    .def("GetString",&HFMIO::GetString)
    .def("SetString",&HFMIO::SetString)
    .def("GetScalar",&HFMIO::GetScalar)
    .def("SetScalar",&HFMIO::SetScalar)
    .def("GetArray", &HFMIO::GetArray)
    .def("SetArray", &HFMIO::SetArray)
    .def("GetComputedKeys", &HFMIO::GetComputedKeys)
    .def("HasField", &HFMIO::HasField)
    .def("EraseField",&HFMIO::EraseField)
    ;*/
}

#include "DispatchAndRun.h"
