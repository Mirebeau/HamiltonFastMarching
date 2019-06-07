// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Copyright 2017 Erik Bekkers, Eindhoven University of Technology, the Netherlands
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef MathematicaIO_hxx
#define MathematicaIO_hxx

void MathematicaIO::StaticSendMsg(bool warn, const std::string & msg, WolframLibraryData libData){
    if(!libData) return;
    std::ostringstream oss;
    oss << "WriteString[\"stdout\",\"";
    if(warn) oss << "***** Warning ! *****\n" << msg << "********************\n";
    else oss << msg;
    if(!libData) return;
    oss << "\"]";
    libData->evaluateExpression(libData,(char*)oss.str().c_str(),6,0,0);
}

void MathematicaIO::MathSetArray(KeyCRef key,const MTensor val){
    // Get the MTensor data
    double* arrayMath = libData->MTensor_getRealData(val);// The data
    int arrayMathLen = libData->MTensor_getFlattenedLength(val);// Flattened length
    const int ndims = libData->MTensor_getRank(val); // Array dimension
    mint const* dimsMath = libData->MTensor_getDimensions(val);// Dimensions
    
    RawElement & raw = CreateElement(key);
    raw.dims.resize(ndims);
    for(int i=0; i<ndims; ++i) {
        raw.dims[i]=dimsMath[i];} // Row major
    
    // Copy data
    raw.data.resize(arrayMathLen);
    std::copy(arrayMath, arrayMath+arrayMathLen,raw.data.begin());
}

MTensor MathematicaIO::MathGetArray(KeyCRef key,int ndims) const{
    const RawElement & raw = GetRaw(key);
    if(ndims!=raw.dims.size()){
        ExceptionMacro("MathGetArray error: tensor " << key << " has rank " << raw.dims.size() << " and not " << ndims);}
    
    std::vector<mint> dims(ndims);
    for(int i=0; i<ndims; ++i) dims[i] = raw.dims[i]; // Row major
    
    MTensor vectorMath_T;
    libData->MTensor_new(MType_Real, ndims, &dims[0], &vectorMath_T);
    
    // Copy data and return
    std::copy(raw.data.begin(),raw.data.begin()+raw.FlattenedLength(),libData->MTensor_getRealData(vectorMath_T));
    return vectorMath_T;
}

std::string MathematicaIO::GetComputedKeys() const {
    std::ostringstream oss;
    oss << "{";
    for(const auto & a : rawElems){
        if(a.second.setter == SetterTag::Compute) {
            oss << "{" << a.first << ",";
            if(a.second.IsString()) {oss << "String";}
            else if(a.second.IsScalar()) {oss << "Real";}
            else oss << "{Real," << a.second.dims.size() << "}";
            oss << "},";}
    }
    if(!rawElems.empty()){oss.seekp(-1,oss.cur);}
    oss << "}";
    return oss.str();
}

#endif /* MathematicaIO_hxx */
