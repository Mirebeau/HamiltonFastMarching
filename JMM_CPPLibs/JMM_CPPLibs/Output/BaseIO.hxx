// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0


#ifndef BaseIO_hxx
#define BaseIO_hxx

// ---- Data container ----

struct BaseIO::RawElement {
    // Case of a string field
    std::string str;
    
    // case of a numerical field
    std::vector<DiscreteType> dims;
    std::vector<ScalarType> data;
    
    SetterTag setter = SetterTag::Unknown;
    
    bool IsString() const {return data.empty() && dims.size()==0;}
    bool IsScalar() const {return !data.empty() && dims.size()==0;}
    void Clear(SetterTag tag) {str.clear(); dims.clear(); data.clear();setter=tag;}
    RawElement(SetterTag tag):setter(tag){};
    DiscreteType FlattenedLength(const bool empty=false) const {
        /* When the data is empty, and the size is empty, 
        (a.k.a the field will hold a yet unspecified scalar) 
        IsString returns a false positive.*/
        if(!empty && IsString()) {ExceptionMacro("BaseIO::RawElement::FlattenedLength error, field is a string");}
        DiscreteType res=1;
        for(auto dim : dims) res*=dim;
        return res;
    }
};


// ----- Public fields -----

bool BaseIO::HasField(KeyCRef key) const {
    if (rawElems.find(key) != rawElems.end()) return true;
    visitedUnset.insert(key);
    return false;
}

bool BaseIO::EraseField(KeyCRef key) {
    unused.erase(key);
    return rawElems.erase(key)!=0;
}

std::string BaseIO::GetString(KeyCRef key) const {
    const auto & val = GetRaw(key);
    if (!val.data.empty())
        ExceptionMacro("BaseIO import error : field " << key << " is not a string.");
    unused.erase(key);
    return val.str;
}

void BaseIO::SetString(KeyCRef key, const std::string & val) {
    RawElement & raw = CreateElement(key);
    raw.str = val;
}

// ----- Protected fields -----

void BaseIO::UsageReport(){
    currentSetter = SetterTag::Compute;
    if(!unused.empty()){
        std::ostringstream ossUser, ossCompute;
        for(KeyCRef key : unused) {
            switch(GetSetter(key)){
                case SetterTag::User: ossUser << key << " "; break;
                case SetterTag::Compute: ossCompute << key << " "; break;
                case SetterTag::Unknown: ExceptionMacro("BaseIO Error : unknown setter for key " << key);
            }
        }
        const std::string & fromUser = ossUser.str(), fromCompute = ossCompute.str();
        if(!fromUser.empty()){
            SetString("unusedFromUser", fromUser);
            if(verbosity>=1) _Msg<true,BaseIO>(this) << "Unused fields from user: " << fromUser << "\n";
        }
        if(!fromCompute.empty()){
            SetString("unusedFromCompute", fromCompute);
            if(verbosity>=2) _Msg<false,BaseIO>(this) << "Unused fields from compute : " << fromCompute << "\n";
        }
    }
    if(!defaulted.empty()){
        std::ostringstream oss;
        for(KeyCRef key : defaulted) oss << key << " ";
        SetString("defaulted", oss.str());
        if(verbosity>=2) _Msg<false,BaseIO>(this) << "Defaulted fields : " << oss.str() << "\n";
    }
    if(!visitedUnset.empty()){
        std::ostringstream oss;
        for(KeyCRef key : visitedUnset) oss << key << " ";
        SetString("visitedUnset", oss.str());
        if(verbosity>=3) _Msg<false,BaseIO>(this) << "Visited but unset fields : " << oss.str() << "\n";
    }
/*	{
		std::ostringstream oss;
		for(KeyCRef key : unusedHelp) oss << key << " ";
		SetString("unusedHelp",oss.str());
	}*/
    currentSetter = SetterTag::User;
	
}

auto BaseIO::GetSetter(KeyCRef key) const -> SetterTag {
    return GetRaw(key,false).setter;
}

const BaseIO::RawElement & BaseIO::GetRaw(KeyCRef key, bool forUse) const {
    const auto it = rawElems.find(key);
    if (it == rawElems.end()) ExceptionMacro("BaseIO import error : field " << key << " not found.");
    if(forUse) unused.erase(key);
    return it->second;
}


auto BaseIO::CreateElement(KeyCRef key) -> RawElement & {
    auto it = rawElems.find(key);
    if (it != rawElems.end() && verbosity>=1){
        _Msg<true,BaseIO>(this) << "BaseIO: redefining field " << key << ".\n";
        it->second.Clear(currentSetter);
    } else {
        it= rawElems.insert(std::pair<std::string,RawElement>(key,RawElement(currentSetter))).first;
    }
    unused.insert(key);
    return it->second;
}

template<typename T> auto BaseIO::GetDimsPtr(KeyCRef key) const
-> std::pair<std::vector<BaseIO::DiscreteType>, FCI<T,const ScalarType> > {
    const auto & raw = GetRaw(key);
    if(raw.IsString()) {ExceptionMacro("IO error : key " << key
									   << " is a string, not a numerical field.");}
    auto dims = raw.dims;
	typedef FCI<T,const ScalarType> FCIT;
    if (!std::is_same<T, ScalarType>::value) {
		if (dims.empty() || dims.back() != FCIT::nComp()) // Row major
            ExceptionMacro("BaseIO input error first dimension "
						   << (dims.empty() ? 0 : dims.back()) << " of field "
                           << key << " does not match expected value " << FCIT::nComp() << ".");
        dims.pop_back();
    }
	return {dims, FCIT(&raw.data[0]) };
}

template<typename T, size_t d, typename F> void BaseIO::Set(KeyCRef key, DimType<d> dims, const F & vals) {
    const DiscreteType size = dims.Product();
	typedef FCI<T, ScalarType> FCIT;
    
    RawElement & raw = CreateElement(key);
	raw.data.resize(FCIT::nComp()*size);
	
	FCIT input{&raw.data[0]};
    for (DiscreteType i = 0; i<size; ++i, ++input)
        input.Set(vals(i));
    
    for(auto dim : dims) raw.dims.push_back(dim);
    if(!std::is_same<T, ScalarType>::value) raw.dims.push_back((DiscreteType)FCIT::nComp()); // Row major
}




#endif /* BaseIO_h */
