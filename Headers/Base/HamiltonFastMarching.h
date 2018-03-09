// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef HamiltonFastMarching_h
#define HamiltonFastMarching_h

#include <array>
#include <bitset>
#include <queue>
#include <map>
#include <cmath>
#include <memory>

#include "DataStructures/ShallowMap.h"
#include "DataStructures/RangeAccessor.h"
#include "DataStructures/CappedVector.h"
#include "PeriodicGrid.h"
#include "Output/ExportMacros.h"

/*
// TODO : remove unnecessary virtual tags.
// Note : walls must not be taken into account at initialization with mult based.
 */

#define FromHFM(x) HFM:: x

// Silly replacement for a constexpr log. TODO : hide somewhere.

template<size_t i> struct FloorLog2 {static const size_t value = 1+FloorLog2<i/2>::value;};
template<> struct FloorLog2<1> {static const size_t value = 0;};
template<size_t i> class CeilLog2 {
    static const size_t floor = FloorLog2<i>::value;
public:
    static const size_t value = (i== (1<<floor) ? floor : (floor+1));
};

template<typename T> struct HFMInterface;

template<typename TTraits>
struct HamiltonFastMarching {
    typedef TTraits Traits;
    Redeclare1Constant(FromTraits,Dimension)
    Redeclare3Types(FromTraits,ScalarType,DiscreteType,ShortType)
    Redeclare5Types(FromTraits,PointType,VectorType,IndexType,OffsetType,DifferenceType)
    
    typedef const IndexType & IndexCRef;
    typedef const OffsetType & OffsetCRef;
    struct FullIndexType {IndexType index; DiscreteType linear;};
    typedef const FullIndexType & FullIndexCRef;
    typedef HamiltonFastMarching HFM;
    typedef HFMInterface<Traits> HFMI;
    
    template<typename E, size_t n> using Array = typename Traits::template Array<E,n>;
    template<typename E> using DataSource = typename Traits::template DataSource<E>;
    template<size_t n> using BasisReduction=typename Traits::template BasisReduction<n>;

    // Domain and Running
    typedef PeriodicGrid<Traits> DomainType;
    typedef typename DomainType::ReverseFlag ReverseFlag;
    const DomainType dom;
    virtual ReverseFlag VisibleOffset(IndexCRef, OffsetCRef, IndexType &) const;

    std::map<IndexType,ScalarType,typename IndexType::LexicographicCompare> seeds;
    void Run();
    bool sndOrder=false;
    Array<ScalarType,Dimension> values;
    
    // Geodesic related functions. Require values and activeNeighs to be correctly set,
    // either by running the algorithm or restoring them from a previous run
    struct GeodesicSolverInterface;
    struct FlowDataType;
    FlowDataType GeodesicFlow(IndexCRef) const;
    static const DiscreteType nNeigh =
    Traits::nForward + 2*Traits::nSymmetric + Traits::nMax*(Traits::nMaxForward+2*Traits::nMaxSymmetric);
    static const DiscreteType nMaxBits = CeilLog2<Traits::nMax>::value;
    typedef std::bitset<nNeigh+nMaxBits> ActiveNeighFlagType;
    Array<ActiveNeighFlagType,Dimension> activeNeighs;

    static const DiscreteType nActiveNeigh =
    Traits::nForward + Traits::nSymmetric + Traits::nMaxForward + Traits::nMaxSymmetric;
    struct DiscreteFlowElement {OffsetType offset; ScalarType weight;};
    typedef CappedVector<DiscreteFlowElement, nActiveNeigh> DiscreteFlowType;
    struct RecomputeType {ScalarType value,width;};
    RecomputeType Recompute(IndexCRef, DiscreteFlowType &) const;

    // StencilDataType must be subclassed.
    struct StencilType;
    typedef ParamInterface_<PointType,VectorType> ParamInterface;
    typedef typename DifferenceType::MultiplierType MultiplierType;
    static const bool hasMultiplier = DifferenceType::multSize>0;
    template<bool b=hasMultiplier, typename Dummy=void> struct _StencilDataType;
    typedef _StencilDataType<> StencilDataType;
    DiscreteType MaxStencilWidth() const;
    std::unique_ptr<StencilDataType> pStencilData;
    
    // Default domain parametrization (Currently used in all instantiations)
    static const bool hasBundle = (0<Traits::nStencilDependencies) && (Traits::nStencilDependencies<Dimension);
    template<int b=int(hasBundle), typename Dummy=void> struct _ParamDefault;
    typedef _ParamDefault<> ParamDefault;
    
    // Extra algorithms may be inserted at different points
    enum Decision {kAccept=0, kContinue=1, kTerminate=2, kRecompute=1<<10};
    struct ExtraAlgorithmInterface;
    struct ExtraAlgorithmPtrs;
    ExtraAlgorithmPtrs extras;
    
    HamiltonFastMarching(std::unique_ptr<StencilDataType>);
protected:
    template<size_t n=Traits::nMax> struct _QuadType;
    typedef _QuadType<> QuadType;
    
    struct QueueElement;
    std::priority_queue<QueueElement> queue;
    
    virtual void RunInit();
    bool RunOnce(); // returns true if should stop
    
    int PostProcess(IndexCRef); // returns a combination of decitions. Priority: kTerminate > kContinue > kAccept
    int PostProcessWithRecompute(IndexCRef, const RecomputeType &, const DiscreteFlowType &);
    
    Array<bool,Dimension> acceptedFlags;
    void ConditionalUpdate(IndexCRef,OffsetType,ScalarType);
    void Update(FullIndexCRef, OffsetCRef, ScalarType);
    friend struct _StencilDataType<>;
};

// ------- Some sub-structures of interest -------

template<typename T> struct HamiltonFastMarching<T>::ExtraAlgorithmInterface {
    virtual ~ExtraAlgorithmInterface(){};
    virtual void Setup(HFMI *){};
    virtual void Finally(HFMI *){};
    virtual bool ImplementIn(HFM *)=0; // Insert in adequate field of extras. Returns false if algorithm is vacuous.
protected:
    virtual int PostProcess(IndexCRef) {return 0;};
    virtual int PostProcessWithRecompute(IndexCRef, const RecomputeType &, const DiscreteFlowType &) {return 0;};
    virtual bool Visible(IndexCRef, OffsetCRef, IndexCRef) const {return true;}
    virtual void BeforeRecompute(IndexCRef) const {};
    friend struct HamiltonFastMarching<T>;
};

template<typename T> struct HamiltonFastMarching<T>::ExtraAlgorithmPtrs {
    std::vector<ExtraAlgorithmInterface*> postProcess;
    std::vector<ExtraAlgorithmInterface*> postProcessWithRecompute;
    std::vector<ExtraAlgorithmInterface*> visible;
    std::vector<ExtraAlgorithmInterface*> beforeRecompute;
};

template<typename T> struct HamiltonFastMarching<T>::GeodesicSolverInterface {
    const HFM & fm;
    virtual bool Run(std::vector<PointType> &) const = 0; // true if failed. Subclass inspection may yield more info on failure.
    virtual void Setup(HFMI *) = 0;
    virtual std::vector<std::vector<PointType> > Run(HFMI *, const std::vector<PointType> &) = 0;
    GeodesicSolverInterface(const HFM & _fm):fm(_fm){};
    virtual ~GeodesicSolverInterface(){};
};


template<typename Traits> struct
HamiltonFastMarching<Traits>::FlowDataType{
    VectorType flow;
    ScalarType value; // map value at point
    ScalarType width; // stencil width
    PrintSelfMacro(FlowDataType)
};

template<typename Traits> struct HamiltonFastMarching<Traits>::
StencilType {
    std::array<DifferenceType, Traits::nForward> forward;
    std::array<DifferenceType, Traits::nSymmetric> symmetric;
    std::array<std::array<DifferenceType, Traits::nMaxForward>, Traits::nMax> maxForward;
    std::array<std::array<DifferenceType, Traits::nMaxSymmetric>, Traits::nMax> maxSymmetric;
    PrintSelfMacro(StencilType);
};

// ******** Stencil data - multiplier based *********

template<typename T> template<typename Dummy>
struct HamiltonFastMarching<T>::_StencilDataType<true,Dummy>{
    typedef HamiltonFastMarching<T> HFM;
    Redeclare8Types(FromHFM,IndexCRef,OffsetCRef,StencilType,QuadType,DifferenceType,MultiplierType,Traits,FullIndexCRef)
    Redeclare3Types(FromHFM,ParamInterface,HFMI,DomainType)
    Redeclare6Types(FromTraits,DiscreteType,ScalarType,PointType,VectorType,IndexType,OffsetType)
    Redeclare2Constants(FromTraits,Dimension,mathPi)
    
    typedef HFM::DataSource<MultiplierType> MultSourceType;
    IndexType dims; // Needs value
    virtual void SetStencil(IndexCRef, StencilType &) = 0; // Needs specialization
    std::unique_ptr<MultSourceType> pMultSource = nullptr; // Needs assignment
    
    struct RecomputeDataType {const StencilType & stencil;MultiplierType mult;};
    RecomputeDataType RecomputeData(IndexCRef); // Needed in some subclasses of HFM...
    virtual void Setup(HFMI *);
    virtual const ParamInterface & Param() const = 0;
protected:
    friend struct HamiltonFastMarching<Traits>;
    void EraseCache(DiscreteType index) {shallowMultQuads.erase(index);}
    struct UpdateDataType {
        const StencilType & stencil;
        const MultiplierType & mult;
        QuadType & quad;
    };
    // Basic data. Index, and linearIndex given to avoid multiple conversions
    UpdateDataType UpdateData(FullIndexCRef);
    RangeAccessor<OffsetType*> ReversedOffsets(FullIndexCRef);
    virtual void Initialize(const HFM *);
private:
    typedef HFM::Array<StencilType, T::nStencilDependencies> StencilArrayType;
    StencilArrayType stencils;
    ShallowMap<DiscreteType, std::pair<MultiplierType, QuadType> > shallowMultQuads;
    typedef typename StencilArrayType::IndexType ShortIndexType;
    IndexType IndexFromShortIndex(const ShortIndexType &) const;
    ShortIndexType ShortIndexFromIndex(IndexCRef) const;
    std::vector<OffsetType> reversedOffsets;
    std::vector<DiscreteType> reversedOffsetsSplits;
};

// ******** Stencil data - standard *********

template<typename T> template<typename Dummy>
struct HamiltonFastMarching<T>::_StencilDataType<false,Dummy>{
    typedef HamiltonFastMarching<T> HFM;
    Redeclare8Types(FromHFM,IndexCRef,OffsetCRef,StencilType,QuadType,DifferenceType,Traits,FullIndexCRef,DomainType)
    Redeclare3Types(FromHFM,ParamInterface,HFMI,MultiplierType)
    Redeclare6Types(FromTraits,DiscreteType,ScalarType,PointType,VectorType,IndexType,OffsetType)
    Redeclare2Constants(FromTraits,Dimension,mathPi)

    IndexType dims; // Needs value
    virtual void SetStencil(IndexCRef, StencilType &) = 0; // Needs specialization
    virtual void Setup(HFMI *);
    virtual const ParamInterface & Param() const = 0;
    struct RecomputeDataType {StencilType stencil; MultiplierType mult;}; // mult is dummy here
    RecomputeDataType RecomputeData(IndexCRef);
protected:
    friend struct HamiltonFastMarching<Traits>;
    void EraseCache(DiscreteType index) {shallowStencilQuads.erase(index);}
    struct UpdateDataType {
        const StencilType & stencil; MultiplierType mult; // mult is dummy here
        QuadType & quad;
    };
    UpdateDataType UpdateData(FullIndexCRef);
    RangeAccessor<OffsetType*> ReversedOffsets(FullIndexCRef);

    virtual void Initialize(const HFM *);
private:
    ShallowMap<DiscreteType, std::pair<StencilType, QuadType> > shallowStencilQuads;
    std::vector<OffsetType> reversedOffsets;
    std::vector<DiscreteType> reversedOffsetsSplits;
};

#include "HamiltonFastMarching.hxx"
#include "HFM_StencilDataType.hxx"
#include "HFM_ParamDefault.hxx"

#endif /* HamiltonFastMarching_h */
