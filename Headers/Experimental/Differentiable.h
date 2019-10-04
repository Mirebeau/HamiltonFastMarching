// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2018 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef Differentiable_h
#define Differentiable_h

/* In this file, we define differentiable implementations of some common models.
 More precisely we replace symmetric upwind finite differences
    max(0,U(x)-U(x-e),U(x)-U(x+e))^2
 with two upwind finite differences
    max(0,U(x)-U(x-e))^2 + max(0,U(x)-U(x+e))^2
 
 This new implementation is expected to be slightly less accurate near the cut locus, and marginally slower.
 However the value function should be continuously differentiable w.r.t. the parameters (cost, boundary conditions),
 which may prove useful in optimization tasks.
 */

#include "Specializations/CommonTraits.h"

// ------------- Isotropic metric, everywhere differentiable implementation ------------
template<size_t VDimension>
struct TraitsIsotropicDiff : TraitsBase<VDimension> {
    typedef TraitsBase<VDimension> Superclass;
    Redeclare3Types(Superclass,DiscreteType,ScalarType,OffsetType)
    Redeclare1Constant(Superclass,Dimension)

    typedef EulerianDifference<OffsetType,ScalarType,1> DifferenceType;
    typedef EulerianStencil<DifferenceType,0,2*Dimension> StencilType;
    
    typedef PeriodicGrid<TraitsIsotropicDiff> DomainType;
};

template<size_t VDimension>
struct StencilIsotropicDiff : HamiltonFastMarching<TraitsIsotropicDiff<VDimension> >::StencilDataType {
    typedef HamiltonFastMarching<TraitsIsotropicDiff<VDimension> > HFM;
    typedef typename HFM::StencilDataType Superclass;
    Redeclare5Types(HFM,ParamDefault,IndexType,StencilType,ParamInterface,HFMI)
    Redeclare1Constant(HFM,Dimension)
    ParamDefault param;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        int n=0;
        for(int i=0; i<Dimension; ++i){
            for(int s=-1; s<=1; s+=2){
                auto & diff = stencil.forward[0][n]; ++n;
                diff.baseWeight=1./square(param.gridScale);
                diff.offset.fill(0);
                diff.offset[i]=s;
            }
        }
        assert(n==stencil.forward[0].size());
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {Superclass::Setup(that); param.Setup(that);}
};


// --------- Differentiation of riemannian distances ------

// ------------- 2D Riemannian metrics ------------
struct TraitsRiemannDiff2 : TraitsBase<2> {
    typedef EulerianDifference<OffsetType,ScalarType,3> DifferenceType;
    typedef EulerianStencil<DifferenceType,3> StencilType;

    static const DiscreteType nStencilDependencies=2;
    constexpr static const std::array<Boundary, Dimension>  boundaryConditions =
    {{Boundary::Closed, Boundary::Closed}};
    
    constexpr static std::array<DiscreteType, nStencilDependencies>
    stencilDependencies = {{0,1}};
    
    typedef PeriodicGrid<TraitsRiemannDiff2> DomainType;
};
// Linker wants the following line for some obscure reason.
constexpr const decltype(TraitsRiemannDiff2::boundaryConditions) TraitsRiemannDiff2::boundaryConditions;
constexpr const decltype(TraitsRiemannDiff2::stencilDependencies) TraitsRiemannDiff2::stencilDependencies;


struct StencilRiemannDiff2 : HamiltonFastMarching<TraitsRiemannDiff2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsRiemannDiff2> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::ParamDefault param;
    ScalarType minWeightRatio = 1e-4;
    
    typedef Traits::BasisReduction<2> ReductionType;
    typedef ReductionType::SymmetricMatrixType SymmetricMatrixType;
    typedef SymmetricMatrixType MetricElementType;
    typedef Traits::DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pDualMetric;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        typedef Traits::ShortType ShortType;
        const auto diff = (*pDualMetric)(index);
        auto sb = ReductionType::CanonicalSuperBase();
        ReductionType::ObtuseSuperbase(diff, sb);
        auto & st=stencil.symmetric[0];
        for(int i=0; i<3; ++i)
            st[i]={OffsetType{(ShortType)-sb[i][1],(ShortType)sb[i][0]},
                -diff.ScalarProduct(sb[PosMod(i+1,3)], sb[PosMod(i+2,3)]),
                (ShortType)i};
        
        ScalarType maxWeight=0;
        for(auto & diff : st){maxWeight=std::max(maxWeight,diff.baseWeight);}
        for(auto & diff : st){
            diff.baseWeight = std::max(diff.baseWeight, maxWeight*minWeightRatio);
            diff.baseWeight/=square(param.gridScale);
        }
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI*) override;
};

void StencilRiemannDiff2::Setup(HFMI*that){
    auto & io = that->io;
    dims = IndexType::CastCoordinates(io.template Get<PointType>("dims"));
    param.Setup(that);
    pMultSource = std::unique_ptr<MultSourceType>(new HFMI::template DataSource_Value<MultiplierType>(MultiplierType::Constant(1)));
    pDualMetric = that->template GetField<SymmetricMatrixType>("dualMetric");
    minWeightRatio = io.template Get<ScalarType>("minWeightRatio",minWeightRatio);
    
    if(io.HasField("dualMetricVariation")){
        const auto & dualVariation = io.template GetArray<SymmetricMatrixType,Dimension+1>("dualMetricVariation");
        for(int i=0; i<Dimension; ++i)
            if(dualVariation.dims[i]!=dims[i])
                ExceptionMacro("dualMetricVariation has inconsistent dims");
        typedef Traits::template Array<MultiplierType, Dimension+1> DeepArrayType;
        DeepArrayType variation;
        variation.dims = dualVariation.dims; variation.resize(dualVariation.size());
        typename Traits::template Array<ScalarType,Dimension> a; a.dims = dims;
        const DiscreteType nPts = a.dims.Product();
        for(int linearIndex=0; linearIndex<nPts; ++linearIndex){
            const IndexType index = a.Convert(linearIndex);
            StencilType stencil;
            SetStencil(index,stencil);
            auto & symmetric = stencil.symmetric[0];
            typename ReductionType::SuperbaseType sb;
            for(int k=0; k<3; ++k) {
                const OffsetType o=symmetric[k].offset;
                sb[k] = {(DiscreteType)o[1],(DiscreteType)-o[0]};
            }
            
            typename DeepArrayType::IndexType deepIndex;
            for(int i=0; i<Dimension; ++i) deepIndex[i]=index[i];
            for(int j=0; j<dualVariation.dims[Dimension]; ++j){
                deepIndex[Dimension]=j;
                const SymmetricMatrixType dDiff = dualVariation(deepIndex);
                for(int k=0; k<3; ++k){
                    variation(deepIndex)[k] =
                    dDiff.ScalarProduct(sb[PosMod(k+1,3)],sb[PosMod(k+2,3)])
                    / (2.*symmetric[k].baseWeight*square(param.gridScale));
                } // for k
            } // for j
        } // for linearIndex
        
        io.SetArray("speedVariation",variation);
    } // if dualForwardVariations
}
/*
 template<typename TS,typename TH> template<typename Dummy> struct
 HFMInterface<TS,TH>::SpecializedStencil<StencilRiemannDiff2,Dummy> {
 typedef HFMInterface<TS,TH> HFMI;
 typedef TH HFM; typedef TS StencilDataType;
 Redeclare5Types(HFM,StencilType,ScalarType,PointType,MultiplierType,OffsetType)
 typedef HFMI::template DataSource<MultiplierType> MultSourceType;
 typedef typename StencilDataType::ReductionType ReductionType;
 typedef typename StencilDataType::SymmetricMatrixType SymmetricMatrixType;
 
 static void Setup(HFMI * that){
 auto & stencil = *that->pStencil; auto & io = that->io;
 stencil.param.gridScale = io.template Get<ScalarType>("gridScale");
 stencil.param.origin = io.template Get<PointType>("origin",PointType::Constant(0));
 stencil.pMultSource = std::unique_ptr<MultSourceType>(new HFMI::template DataSource_Value<MultiplierType>(MultiplierType::Constant(1)));
 stencil.pDualMetric = that->template GetField<SymmetricMatrixType>("dualMetric");
 
 //        stencil.pMultSource = that->template GetField<MultiplierType>("speed");
 stencil.minWeightRatio = io.template Get<ScalarType>("minWeightRatio",stencil.minWeightRatio);
 
 if(io.HasField("dualForwardVariation")){
 const auto & dualVariation = io.template GetArray<SymmetricMatrixType,Dimension+1>("dualForwardVariation");
 for(int i=0; i<Dimension; ++i)
 if(dualVariation.dims[i]!=stencil.dims[i])
 ExceptionMacro("dualForwardVariation has inconsistent dims");
 typedef Array<MultiplierType, Dimension+1> DeepArrayType;
 DeepArrayType variation;
 variation.dims = dualVariation.dims; variation.resize(dualVariation.size());
 Array<ScalarType,Dimension> a; a.dims = stencil.dims;
 const DiscreteType nPts = a.dims.ProductOfCoordinates();
 for(int linearIndex=0; linearIndex<nPts; ++linearIndex){
 const IndexType index = a.Convert(linearIndex);
 StencilType st;
 stencil.SetStencil(index,st);
 typename ReductionType::SuperbaseType sb;
 for(int k=0; k<3; ++k) {
 const OffsetType o=st.symmetric[k].offset;
 sb[k] = {(DiscreteType)o[1],(DiscreteType)-o[0]};
 }
 
 typename DeepArrayType::IndexType deepIndex;
 for(int i=0; i<Dimension; ++i) deepIndex[i]=index[i];
 for(int j=0; j<dualVariation.dims[Dimension]; ++j){
 deepIndex[Dimension]=j;
 const SymmetricMatrixType dDiff = dualVariation(deepIndex);
 for(int k=0; k<3; ++k){
 variation(deepIndex)[k] =
 dDiff.ScalarProduct(sb[PosMod(k+1,3)],sb[PosMod(k+2,3)])
 / (2.*st.symmetric[k].baseWeight*square(stencil.param.gridScale));
 } // for k
 } // for j
 } // for linearIndex
 
 io.SetArray("forwardVariation",variation);
 } // if dualForwardVariations
 }
 };
 */

#endif /* Differentiable_h */
