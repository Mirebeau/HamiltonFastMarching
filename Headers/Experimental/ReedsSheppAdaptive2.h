// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef ReedsSheppAdaptive2_h
#define ReedsSheppAdaptive2_h

/*
 
 The model ReedsSheppForwardAdaptive2 uses an adaptive angular discretization,
 not uniform contrary to the models in curvature2 and curvature3.
 The angles are chosen such that (cos t, sin t) is proportional to a vector with small integer coordinates.
 - Pros : this avoids the artefact where the vehicle advances sideways, i.e. the constraint on the direction is not sufficiently enforced.
 - Cons : unfortunately, there are too few such angles, hence the discretization is not very accurate.
 
 TODO : find another way to avoid sideways motions, which has more directions. (i.e. the related tensor decomposition is mostly aligned).
 e.g. build the tensors manually, so as to guarantee that their obtuse superbase only contains well aligned vectors.
 */

/*
 For convenience, I also add a ReedsShepp model where one can independently set the forward, reverse, and angular speed.
 */

#include "Specializations/Curvature2.h"

// -----------

struct TraitsReedsSheppThreeSpeeds2 : TraitsR2S1 {
    typedef Difference<3> DifferenceType;
    static const int nForward = 3+3, nSymmetric=1;
};

struct StencilReedsSheppThreeSpeeds2 : HamiltonFastMarching<TraitsReedsSheppThreeSpeeds2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsReedsSheppThreeSpeeds2> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::ParamDefault param;
    ScalarType eps=0.1, xi=1; // xi is the typical curvature radius
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        const ScalarType theta = param.dependScale*index[2];
        typedef Traits::BasisReduction<2> ReductionType;
        const typename ReductionType::VectorType v{cos(theta),sin(theta)};
        Voronoi1Vec<ReductionType>(&stencil.forward[0], v/param.gridScale,eps);
        Voronoi1Vec<ReductionType>(&stencil.forward[3],-v/param.gridScale,eps);
        
        for(int i=0; i<6; ++i) stencil.forward[i].multIndex = (i<3 ? 0 : 1);
        stencil.symmetric[0] = {1./square(xi*param.dependScale),OffsetType{0,0,1},2};
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {Superclass::Setup(that);
        param.Setup(that,2.*Traits::mathPi/dims[2]);
        eps=that->io.template Get<ScalarType>("eps",eps);
        xi=that->io.template Get<ScalarType>("xi",xi);
    }
};

// -----------
struct TraitsReedsSheppAdaptive2 : TraitsR2S1 {
    typedef Difference<3> DifferenceType; // Physical and angular speed functions
    static const int nForward = 4;
    
    // Interpolate last coordinate w.r.t. provided weights
    template<typename E> struct InterpolatedDataSource : DataSource<E> {
        std::vector<std::pair<std::pair<DiscreteType,DiscreteType>,ScalarType> > weights;
        std::unique_ptr<DataSource<E> > data;
        virtual bool CheckDims(const IndexType & index) const {
            if(data==nullptr) ExceptionMacro("InterpolatedDataSource error : missing data");
            return data->CheckDims(index);}
        
        typedef typename DataSource<E>::ReturnType ReturnType;
        virtual ReturnType operator()(const IndexType & index) const {
            assert(data);
            assert(index.back()<weights.size());
            const auto w = weights[index.back()];
            IndexType ind0=index, ind1=index;
            ind0.back()=w.first.first;
            ind1.back()=w.first.second;
            const auto r0 = (*data)(ind0), r1=(*data)(ind1);
            return r0+(r1-r0)*w.second;
        }
    };
};

struct StencilReedsSheppAdaptive2 : HamiltonFastMarching<TraitsReedsSheppAdaptive2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsReedsSheppAdaptive2> HFM;
    typedef HFM::StencilDataType Superclass;
    struct ParamType : HFM::ParamInterface {
        std::vector<ScalarType> angles;
        ScalarType ADim(ScalarType) const;
        ScalarType ReDim(ScalarType) const;
        

        ScalarType gridScale=1;
        typedef LinearAlgebra::Point<ScalarType, 2> OriginType; // Only used for operator << instead of std::array
        OriginType origin{0,0};
        virtual PointType ADim(const PointType & p) const override {
            return {(p[0]-origin[0])/gridScale,(p[1]-origin[1])/gridScale,ADim(p[2])};}
        virtual PointType ReDim(const PointType & p) const override {
            return {p[0]*gridScale+origin[0],p[1]*gridScale+origin[0],ReDim(p[2])};}
        virtual VectorType ReDim(const VectorType & v) const override {
            return {v[0]*gridScale,v[1]*gridScale,ReDim(v[2])};}
    } param;
    
    ScalarType xi=1;
    typedef std::array<DiscreteType,2> PhysicalOffsetType;
    std::vector<PhysicalOffsetType> angularOffsets;
    static constexpr std::array<PhysicalOffsetType,11> defaultOffsets =
    {{ {{2,1}}, {{3, 1}}, {{3, 2}}, {{5, 1}}, {{5, 4}}, {{8, 1}},
        {{5,2}}, {{4, 1}}, {{5, 3}}, {{12, 1}}, {{8, 7}} }};
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        const DiscreteType iTheta = index[2], nTheta=dims[2];
        assert(0<=iTheta && iTheta<nTheta);
        assert(angularOffsets.size()==nTheta && param.angles.size()==nTheta);
        
        // Physical offset
        const auto offset = angularOffsets[iTheta];
        stencil.forward[0]=
        {1./(square(param.gridScale)*(square((ScalarType)offset[0])+square((ScalarType)offset[1]))),
            {(Traits::ShortType)offset[0],(Traits::ShortType)offset[1],0},0};
        stencil.forward[1]=
        {1./(square(param.gridScale)*(square((ScalarType)offset[0])+square((ScalarType)offset[1]))),
            {(Traits::ShortType)-offset[0],(Traits::ShortType)-offset[1],0},1};

        
        // Angular offsets
        auto stencilIt = stencil.forward.begin()+2;
        const ScalarType theta = param.angles[iTheta];
        for(int eps=-1; eps<=1; eps+=2, ++stencilIt){
            ScalarType dTheta = param.angles[PosMod(iTheta+eps,nTheta)]-theta;
            if(dTheta>=mathPi) dTheta-=2*mathPi;
            if(dTheta<=-mathPi) dTheta+=2*mathPi;
            *stencilIt = {1./square(xi*dTheta),{0,0,(Traits::ShortType)eps},2};
        }
        assert(stencilIt==stencil.forward.end());
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI*) override;
};

constexpr decltype(StencilReedsSheppAdaptive2::defaultOffsets) StencilReedsSheppAdaptive2::defaultOffsets;

auto StencilReedsSheppAdaptive2::ParamType::ADim(ScalarType t) const -> ScalarType {
    const ScalarType tPer = fPosMod(t, 2.*mathPi);
    auto it = std::upper_bound(angles.begin(), angles.end(), t); // *it>=t
    const ScalarType theta1 = (it==angles.end()) ? (angles.front()+2.*mathPi) : *it;
    assert(it!=angles.begin());
    --it;
    const DiscreteType iTheta = it-angles.begin();
    const ScalarType theta0 = *it;
    
    return (t-tPer)*angles.size()/(2.*mathPi) + iTheta + (tPer-theta0)/(theta1-theta0) + 0.5;
}

auto StencilReedsSheppAdaptive2::ParamType::ReDim(ScalarType t) const -> ScalarType {
    t-=0.5;
    const DiscreteType nTheta = angles.size();
    const ScalarType tPer = fPosMod(t, (ScalarType)nTheta);
    const DiscreteType iTheta = std::floor(tPer);
    const ScalarType theta0 = angles[iTheta], theta1=angles[PosMod(iTheta+1, nTheta)];
    return (t-tPer)*(2.*mathPi)/nTheta+theta0+(theta1-theta0)*(tPer-iTheta);
}

void StencilReedsSheppAdaptive2::Setup(HFMI*that) {
    Superclass::Setup(that);
    auto & io = that->io;
    param.gridScale = io.template Get<ScalarType>("gridScale");
    typedef StencilReedsSheppAdaptive2::ParamType::OriginType OriginType;
    param.origin = io.template Get<OriginType>("origin",param.origin);
    xi = io.template Get<ScalarType>("xi",xi);
    
    // Setup angular offsets
    auto & angles  = param.angles;
    auto & offsets = angularOffsets;
    
    typedef PhysicalOffsetType V;
    typedef std::array<ScalarType,2> W;
    auto Angle = [](const V & v)->ScalarType {ScalarType t = atan2((ScalarType)v[1],(ScalarType)v[0]);
        return t<0 ? t+2.*mathPi : t;};
    
    if(io.HasField("angularOffsets")){
        auto angOff = io.template GetVector<W>("angularOffsets");
        angularOffsets.reserve(angOff.size());
        for(const W & o : angOff) angularOffsets.push_back({{(DiscreteType)o[0],(DiscreteType)o[1]}});
    } else {
        const int n=dims[2];
        if(n<8 || n%8!=0) ExceptionMacro("Error : number of angles (" << n << ") should be a multiple of 8.\n");
        const int maxN = 8+8*defaultOffsets.size();
        if(n>maxN) ExceptionMacro("Angular offsets defaulted up to " << maxN << "directions only.\n");
        offsets.reserve(n);
        offsets.insert(offsets.end(),{ {{1,0}}, {{0,1}}, {{-1,0}}, {{0,-1}} });
        offsets.insert(offsets.end(),{ {{1,1}}, {{-1,1}}, {{-1,-1}}, {{1,-1}} });
        for(int i=0; i<n/8-1; ++i) {const V v = defaultOffsets[i];
            offsets.insert(offsets.end(),{
                {{v[0],v[1]}}, {{v[1],v[0]}}, {{-v[0],-v[1]}}, {{-v[1],-v[0]}},
                {{v[0],-v[1]}}, {{-v[1],v[0]}}, {{-v[0],v[1]}}, {{v[1],-v[0]}}
            });}
        std::sort(offsets.begin(),offsets.end(),[&](V u,V v)->bool{return Angle(u)<Angle(v);});
        std::vector<W> angOff; angOff.reserve(offsets.size());
        for(const V v : offsets) angOff.push_back({{(ScalarType)v[0],(ScalarType)v[1]}});
        io.template SetVector<W>("angularOffsets",angOff);
    }
    angles.reserve(offsets.size());
    ScalarType lastAngle = -Traits::Infinity();
    for(const auto & v : offsets){
        const ScalarType theta = Angle(v);
        if(theta<=lastAngle) ExceptionMacro("Error : invalid angular offsets");
        angles.push_back(theta);
        lastAngle=theta;
    }
    
    
    if(!io.HasField("uniformlySampledSpeed")){
        pMultSource = that->template GetField<MultiplierType>("speed");
    } else {
        DiscreteType nTheta = io.template Get<ScalarType>("uniformlySampledSpeed");
        if(nTheta<=0) ExceptionMacro("Error : uniformlySampledSpeed should be positive (=number of angular orientations)");
        auto * interp = new typename Traits::template InterpolatedDataSource<MultiplierType>;
        pMultSource = std::unique_ptr<MultSourceType>(interp);
        std::swap(nTheta,dims[2]); // For CheckDims
        interp->data = that->template GetField<MultiplierType>("speed");
        std::swap(nTheta,dims[2]);
        
        interp->weights.reserve(angles.size());
        for(const ScalarType theta : angles){
            const ScalarType adAngle = theta*nTheta/(2.*mathPi);
            const DiscreteType i = std::floor(adAngle);
            assert(i<nTheta);
            interp->weights.push_back({{i,PosMod(i+1, nTheta)},adAngle-i});
        }
    }
}
#endif /* ReedsSheppAdaptive2_h */
