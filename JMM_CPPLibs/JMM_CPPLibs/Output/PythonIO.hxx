// Copyright 2018 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef PythonIO_hxx
#define PythonIO_hxx

void PythonIO::PySetArray(KeyCRef key, ndarray arr) {
    
    py::buffer_info info = arr.request();
    RawElement & raw = CreateElement(key);
    raw.Clear(currentSetter);
    const int ndims = arr.ndim();
    raw.dims.resize(ndims);

    for(int i=0; i<ndims; ++i)
        raw.dims[i]=info.shape[i];
    
    const DiscreteType flattenedLength = raw.FlattenedLength(true); // Empty field
    raw.data.resize(flattenedLength);
    
    ScalarType * beg = reinterpret_cast<ScalarType*>(info.ptr);
    std::copy(beg,beg+flattenedLength, raw.data.begin());
}

auto PythonIO::PyGetArray(KeyCRef key) const -> ndarray {
    const RawElement & raw = GetRaw(key);
    if(raw.IsString()) {ExceptionMacro("PythonIO error : field " << key << " is a string, not an array");}
    if(raw.dims.empty()){ExceptionMacro("PythonIO error : field " << key << " is a scalar, not an array");}
    ndarray result(raw.dims);
    
    assert(raw.FlattenedLength()==raw.data.size());
    ScalarType * beg = static_cast<ScalarType*>(result.request().ptr);
    std::copy(raw.data.begin(),raw.data.end(),beg);
    return result;
}

std::string PythonIO::GetComputedKeys() const {
    std::ostringstream oss;
    oss << "[";
    for(const auto & a : rawElems){
        if(a.second.setter == SetterTag::Compute) {
            oss << "['" << a.first << "','" <<
            (a.second.IsString() ? "string" : a.second.IsScalar() ? "float" : "array")
            << "'],";}
    }
    oss << "]";
    return oss.str();
}
#endif /* PythonIO_hxx */
