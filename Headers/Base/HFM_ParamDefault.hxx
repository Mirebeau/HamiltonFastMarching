//
//  HFM_ParamDefault.hxx
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 04/12/2017.
//

#ifndef HFM_ParamDefault_hxx
#define HFM_ParamDefault_hxx

// --- Default domain parametrization ----

template<typename Traits> template<typename Dummy> struct
HamiltonFastMarching<Traits>::_ParamDefault<true,Dummy> : HFM::ParamInterface {
    PointType origin = PointType::Constant(0);
    // Two distinct scales, for "physical" and "bundle" parts respectively.
    ScalarType gridScale=1, dependScale=1;
    virtual PointType ADim(const PointType & p) const override;
    virtual PointType ReDim(const PointType & p) const override;
    virtual VectorType ReDim(const VectorType & u) const override;
    void Setup(IO & io, ScalarType _dependScale){
        // Setting the scale and origin. This is really specialized for R^n x S^d structures.
        gridScale = io.template Get<ScalarType>("gridScale");
        dependScale = _dependScale;
        origin.back() = -dependScale/2.;
        
        if(io.HasField("origin")){ // Only importing the physical dimensions of origin
            const DiscreteType DimIndep = Dimension-Traits::nStencilDependencies;
            const auto indepOrigin = io.template Get<std::array<ScalarType,DimIndep> >("origin");
            auto indIt = indepOrigin.begin(); const auto dep=Traits::stencilDependencies; auto depIt = dep.begin();
            for(int i=0; i<Dimension;++i){
                if(depIt!=dep.end() && i==*depIt){++depIt;}
                else {origin[i]=*indIt; ++indIt;}
            }
        }
    }
    void Setup(HFMI * that, ScalarType _dependScale) {Setup(that->io,_dependScale);}
};

template<typename Traits> template<typename Dummy> struct
HamiltonFastMarching<Traits>::_ParamDefault<false,Dummy> : HFM::ParamInterface {
    PointType origin = PointType::Constant(0);
    ScalarType gridScale=1;
    virtual PointType ADim(const PointType & p) const override {
        return PointType::FromOrigin((p-origin)/gridScale);}
    virtual PointType ReDim(const PointType & p) const override {
        return VectorType::FromOrigin(p)*gridScale+origin;}
    virtual VectorType ReDim(const VectorType & v) const override {
        return v*gridScale;}
    void Setup(HFMI * that) {
        origin=that->io.template Get<PointType>("origin",origin);
        gridScale = that->io.template Get<ScalarType>("gridScale");
    }
};


template<typename T> template<typename Dummy> auto
HamiltonFastMarching<T>::_ParamDefault<true,Dummy>::
ADim(const PointType & p) const -> PointType {
    PointType result;
    for(int i=0, j=0; i<Dimension; ++i){
        if(j<Traits::nStencilDependencies && i==Traits::stencilDependencies[j]){
            result[i] = (p[i]-origin[i])/dependScale;
            ++j;
        } else {
            result[i] = (p[i]-origin[i])/gridScale;
        }
    }
    return result;
}

template<typename T> template<typename Dummy> auto
HamiltonFastMarching<T>::_ParamDefault<true,Dummy>::
ReDim(const PointType & p) const -> PointType {
    PointType result;
    for(int i=0,j=0; i<Dimension; ++i){
        if(j<T::nStencilDependencies && i==Traits::stencilDependencies[j]){
            result[i] = p[i]*dependScale+origin[i];
            ++j;
        } else {
            result[i] = p[i]*gridScale+origin[i];
        }
    }
    return result;
}

template<typename T> template<typename Dummy> auto
HamiltonFastMarching<T>::_ParamDefault<true,Dummy>::
ReDim(const VectorType & p) const -> VectorType {
    VectorType result;
    for(int i=0,j=0; i<Dimension; ++i){
        if(j<T::nStencilDependencies && i==Traits::stencilDependencies[j]){
            result[i] = p[i]*dependScale;
            ++j;
        } else {
            result[i] = p[i]*gridScale;
        }
    }
    return result;
}


/*// --- Default domain parametrization ----
 
 template<typename Traits> template<typename Dummy> struct
 HamiltonFastMarching<Traits>::_StencilDataType<true,Dummy>::ParamType : HFM::ParamInterface {
 PointType origin = PointType::Constant(0);
 // Two distinct scales, for "physical" and "bundle" parts respectively.
 ScalarType gridScale=1, dependScale=1;
 virtual PointType ADim(const PointType & p) const override;
 virtual PointType ReDim(const PointType & p) const override;
 virtual VectorType ReDim(const VectorType & u) const override;
 void Setup(IO & io, ScalarType _dependScale){
 // Setting the scale and origin. This is really specialized for R^n x S^d structures.
 gridScale = io.template Get<ScalarType>("gridScale");
 dependScale = _dependScale;
 origin.back() = -dependScale/2.;
 
 if(io.HasField("origin")){ // Only importing the physical dimensions of origin
 const DiscreteType DimIndep = Dimension-Traits::nStencilDependencies;
 const auto indepOrigin = io.template Get<std::array<ScalarType,DimIndep> >("origin");
 auto indIt = indepOrigin.begin(); const auto dep=Traits::stencilDependencies; auto depIt = dep.begin();
 for(int i=0; i<Dimension;++i){
 if(depIt!=dep.end() && i==*depIt){++depIt;}
 else {origin[i]=*indIt; ++indIt;}
 }
 }
 }
 void Setup(HFMI * that, ScalarType _dependScale) {Setup(that->io,_dependScale);}
 };
 
 template<typename Traits> template<typename Dummy> struct
 HamiltonFastMarching<Traits>::_StencilDataType<false,Dummy>::ParamType : HFM::ParamInterface {
 PointType origin = PointType::Constant(0);
 ScalarType gridScale=1;
 virtual PointType ADim(const PointType & p) const override {
 return PointType::FromOrigin((p-origin)/gridScale);}
 virtual PointType ReDim(const PointType & p) const override {
 return VectorType::FromOrigin(p)*gridScale+origin;}
 virtual VectorType ReDim(const VectorType & v) const override {
 return v*gridScale;}
 void Setup(HFMI * that) {
 origin=that->io.template Get<PointType>("origin",origin);
 gridScale = that->io.template Get<ScalarType>("gridScale");
 }
 };
 
 
 template<typename T> template<typename Dummy> auto
 HamiltonFastMarching<T>::_StencilDataType<true,Dummy>::ParamType::
 ADim(const PointType & p) const -> PointType {
 PointType result;
 for(int i=0, j=0; i<Dimension; ++i){
 if(j<Traits::nStencilDependencies && i==Traits::stencilDependencies[j]){
 result[i] = (p[i]-origin[i])/dependScale;
 ++j;
 } else {
 result[i] = (p[i]-origin[i])/gridScale;
 }
 }
 return result;
 }
 
 template<typename T> template<typename Dummy> auto
 HamiltonFastMarching<T>::_StencilDataType<true,Dummy>::ParamType::
 ReDim(const PointType & p) const -> PointType {
 PointType result;
 for(int i=0,j=0; i<Dimension; ++i){
 if(j<T::nStencilDependencies && i==Traits::stencilDependencies[j]){
 result[i] = p[i]*dependScale+origin[i];
 ++j;
 } else {
 result[i] = p[i]*gridScale+origin[i];
 }
 }
 return result;
 }
 
 template<typename T> template<typename Dummy> auto
 HamiltonFastMarching<T>::_StencilDataType<true,Dummy>::ParamType::
 ReDim(const VectorType & p) const -> VectorType {
 VectorType result;
 for(int i=0,j=0; i<Dimension; ++i){
 if(j<T::nStencilDependencies && i==Traits::stencilDependencies[j]){
 result[i] = p[i]*dependScale;
 ++j;
 } else {
 result[i] = p[i]*gridScale;
 }
 }
 return result;
 }
 */


#endif /* HFM_ParamDefault_hxx */
