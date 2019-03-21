//
//  HFM_ParamDefault.hxx
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 04/12/2017.
//

#ifndef HFM_ParamDefault_hxx
#define HFM_ParamDefault_hxx

// --- Default domain parametrization ----

template<typename TTraits> template<typename Dummy> struct
HamiltonFastMarching<TTraits>::_ParamDefault<1,Dummy> : HFM::ParamInterface {
    typedef HamiltonFastMarching<TTraits> HFM;
    Redeclare6Types(HFM,PointType,VectorType,ScalarType,DiscreteType,HFMI,Traits);
    Redeclare1Constant(HFM,Dimension)
    PointType origin = PointType::Constant(0);
    // Two distinct scales, for "physical" and "bundle" parts respectively.
    ScalarType gridScale=1, dependScale=1;
    virtual PointType ADim(const PointType & p) const override;
    virtual VectorType ADim(const VectorType & u) const override;
    virtual PointType ReDim(const PointType & p) const override;
    virtual VectorType ReDim(const VectorType & u) const override;
    void Setup(IO & io, ScalarType _dependScale){
        // Setting the scale and origin. This is really specialized for R^n x S^d structures.
        gridScale = io.template Get<ScalarType>("gridScale");
		if(gridScale<=0){ExceptionMacro("Parametrization error : gridscale "
										<< gridScale << " should be positive");}
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

template<typename TTraits> template<typename Dummy> struct
HamiltonFastMarching<TTraits>::_ParamDefault<0,Dummy> : HFM::ParamInterface {
    typedef HamiltonFastMarching<TTraits> HFM;
    Redeclare5Types(HFM,PointType,VectorType,ScalarType,HFMI,Traits);
    PointType origin = PointType::Constant(0);
    ScalarType gridScale=1;
    virtual PointType ADim(const PointType & p) const override {
        return PointType::FromOrigin((p-origin)/gridScale);}
    virtual VectorType ADim(const VectorType & v) const override {
        return v/gridScale;}
    virtual PointType ReDim(const PointType & p) const override {
        return VectorType::FromOrigin(p)*gridScale+origin;}
    virtual VectorType ReDim(const VectorType & v) const override {
        return v*gridScale;}
    void Setup(HFMI * that) {
        origin=that->io.template Get<PointType>("origin",origin);
        gridScale = that->io.template Get<ScalarType>("gridScale");
		if(gridScale<=0){ExceptionMacro("Parametrization error : gridscale "
										<< gridScale << " should be positive");}
    }
};

template<typename TTraits> template<typename Dummy> struct
HamiltonFastMarching<TTraits>::_ParamDefault<2,Dummy> : HFM::ParamInterface {
    typedef HamiltonFastMarching<TTraits> HFM;
    Redeclare5Types(HFM,PointType,VectorType,ScalarType,HFMI,Traits);
    Redeclare1Constant(HFM,Dimension)

    PointType origin = PointType::Constant(0);
    PointType gridScales = PointType::Constant(1);
    
    virtual PointType ADim(const PointType & p) const override {
        PointType result;
        for(int i=0; i<Dimension; ++i) result[i]=(p[i]-origin[i])/gridScales[i];
        return result;}
    virtual VectorType ADim(const VectorType & v) const override {
        VectorType result;
        for(int i=0; i<Dimension; ++i) result[i]=v[i]/gridScales[i];
        return result;}
    
    virtual PointType ReDim(const PointType & p) const override {
        PointType result;
        for(int i=0; i<Dimension; ++i) result[i]=p[i]*gridScales[i]+origin[i];
        return result;}
    virtual VectorType ReDim(const VectorType & v) const override {
        VectorType result;
        for(int i=0; i<Dimension; ++i) result[i]=v[i]*gridScales[i];
        return result;}
    void Setup(HFMI * that) {
		auto & io = that->io;
        origin=io.template Get<PointType>("origin",origin);
		if(io.HasField("gridScales")){gridScales = io.template Get<PointType>("gridScales");}
		else {gridScales = PointType::Constant(io.template Get<ScalarType>("gridScale"));}
		if(!gridScales.IsPositive()){ExceptionMacro("Parametrization error : gridscales "
													<< gridScales << " should be positive");}
    }
};

// ---------- Implementation ----------

template<typename T> template<typename Dummy> auto
HamiltonFastMarching<T>::_ParamDefault<1,Dummy>::
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
HamiltonFastMarching<T>::_ParamDefault<1,Dummy>::
ADim(const VectorType & p) const -> VectorType {
    VectorType result;
    for(int i=0,j=0; i<Dimension; ++i){
        if(j<T::nStencilDependencies && i==Traits::stencilDependencies[j]){
            result[i] = p[i]/dependScale;
            ++j;
        } else {
            result[i] = p[i]/gridScale;
        }
    }
    return result;
}


template<typename T> template<typename Dummy> auto
HamiltonFastMarching<T>::_ParamDefault<1,Dummy>::
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
HamiltonFastMarching<T>::_ParamDefault<1,Dummy>::
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


#endif /* HFM_ParamDefault_hxx */
