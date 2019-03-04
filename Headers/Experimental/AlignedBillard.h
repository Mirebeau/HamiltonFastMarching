//
//  AlignedBillard.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 18/07/2018.
//

#ifndef AlignedBillard_h
#define AlignedBillard_h

#include "Specializations/CommonTraits.h"


// ----- Non square domain, periodic (in a weird way)
// TODO : similar class, with transformation matrices (i.e. not aligned edges)
// In that case, IndexFromPoint would likely better return an index corresponding to the same orientation.

// TODO : find a better way to initialize, current approach is ugly (and get rid of mutable).
// Probably, the initialization of HFM and its StencilData needs to be modified, in a single function call, with reference to the domain. Or better, the domain should be able to initialize itself, but then some modifications are needed.

// TODO : check if works with the other stencil cache policy.

template<typename TTraits>
struct AlignedBillardGrid : TTraits::BaseDomain {
    typedef TTraits Traits;
    Redeclare6Types(Traits,ScalarType,PointType,VectorType,DiscreteType,IndexType,IndexCRef)
    Redeclare1Constant(Traits,Dimension)
    static_assert(Dimension==2, "AlignedBillardGrid error: two dimensional only");
    
    typedef typename Traits::BaseDomain Superclass;
    using Superclass::Superclass;
    
    struct Transform;
    
    static const bool periodizeUsesBase = true;
    Transform Periodize(IndexType & target, IndexCRef base) const;
    Transform Periodize(PointType & target, PointType base) const;
    
    Transform PeriodizeNoBase(IndexType & target) const;
    Transform PeriodizeNoBase(PointType & target) const;

    struct EdgeType {
        PointType p;
        VectorType v, trans;
    };

    mutable std::vector<EdgeType> edges;
    void SetPeriodizeLazy() const;
protected:
    mutable typename Traits::template Array<bool,Dimension> inside; // for lazy periodization

    ScalarType LoopIndex(const PointType &) const;
    ScalarType Cut(PointType &, const VectorType &, ScalarType) const; //modifies point
    ScalarType Cut(const PointType &, const VectorType &, const EdgeType &) const;
};

// ------------- Riemannian metric on such a domain ------------
struct TraitsAlignedBillard : TraitsBase<2> {
    typedef TraitsBase<2> Superclass;
    typedef EulerianDifference<OffsetType,ScalarType,0> DifferenceType;
    typedef EulerianStencil<DifferenceType,(Dimension*(Dimension+1))/2> StencilType;
    typedef AlignedBillardGrid<TraitsAlignedBillard> DomainType;
};

// Transform does essentially nothing
template<typename T> struct
AlignedBillardGrid<T>::Transform {
	bool IsValid() const {return valid;}
	template<typename TVec> void PullVector(TVec & v) const {};
	void Invalidate(){valid=false;}
	//    static Transform Error() {Transform result; result.Invalidate(); return result;}
protected:
	bool valid = true;
};


struct StencilAlignedBillard final :
HamiltonFastMarching<TraitsAlignedBillard>::StencilDataType {
	
    typedef HamiltonFastMarching<TraitsAlignedBillard> HFM;
    typedef typename HFM::StencilDataType Superclass;
    Redeclare5Types(HFM,ParamDefault,IndexType,StencilType,ParamInterface,HFMI)
    Redeclare1Constant(HFM,Dimension)
    ParamDefault param;
    
    typedef typename Traits::template BasisReduction<Dimension> ReductionType;
    typedef typename ReductionType::SymmetricMatrixType SymmetricMatrixType;
    typedef SymmetricMatrixType MetricElementType;
    typedef typename Traits::template DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pDualMetric;

    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0][0],(*pDualMetric)(index));
        const ScalarType hm2 = 1/square(param.gridScale);
        for(auto & diff : stencil.symmetric[0]) diff.baseWeight*=hm2;
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        param.Setup(that);
        pDualMetric = that->template GetField<MetricElementType>("dualMetric");
        pHFMI=that;
    }
    
protected:
    HFMI * pHFMI=nullptr;
    virtual void Initialize(const HFM * pHFM) override {
        auto & io = pHFMI->io;
        const auto p = io.GetVector<PointType>("edge_p");
        const auto v = io.GetVector<VectorType>("edge_v");
        const auto trans = io.GetVector<VectorType>("edge_trans");
        
        const size_t n = p.size();
        if(v.size()!=n || trans.size()!=n){
            ExceptionMacro("Aligned billard error : inconsistent edge data lengths");}
        
        auto & edges = pHFM->dom.edges;
        edges.reserve(n);
        for(size_t i=0; i<n; ++i){
            edges.push_back({param.ADim(p[i]),param.ADim(v[i]),param.ADim(trans[i])});
        }
        pHFM->dom.SetPeriodizeLazy();
        
        Superclass::Initialize(pHFM);
    }
	 
};

#include "Implementation/AlignedBillard.hxx"
#endif /* AlignedBillard_h */
