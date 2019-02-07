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

#include "JMM_CPPLibs/DataStructures/ShallowMap.h"
#include "JMM_CPPLibs/DataStructures/RangeAccessor.h"
#include "JMM_CPPLibs/DataStructures/CappedVector.h"
#include "BaseGrid.h"
#include "JMM_CPPLibs/Macros/ExportArrow.h"
#include "JMM_CPPLibs/Macros/TemplateLog2.h"
/*
// Note : walls must not be taken into account at initialization with mult based.
 */

template<typename T> struct HFMInterface;
template<typename T> struct Factoring;

enum class StencilStoragePolicy {Share,Recomp,Lag2,Lag3};
enum class FactoringMethod {None, Static, Dynamic};

template<typename TTraits>
struct HamiltonFastMarching {
    typedef TTraits Traits;
    Redeclare1Constant(Traits,Dimension)
    Redeclare6Types(Traits,ScalarType,DiscreteType,ShortType,DomainType,
					StencilType,DistanceGuess)
    Redeclare6Types(Traits,PointType,VectorType,IndexType,OffsetType,
					DifferenceType,IndexDiff)
    
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
    typedef typename DomainType::Transform DomainTransformType;
    const DomainType dom;
    DomainTransformType VisibleOffset(IndexCRef, OffsetCRef, IndexType &) const;

    std::map<IndexType,ScalarType,typename IndexType::LexicographicCompare> seeds;
    void Run();
	int order = 1;
    Array<ScalarType,Dimension> values;
    
    // Geodesic related functions. Require values and activeNeighs to be correctly set,
    // either by running the algorithm or restoring them from a previous run
    struct GeodesicSolverInterface;
    struct FlowDataType;
    FlowDataType GeodesicFlow(IndexCRef) const;
    
    Redeclare4Types(StencilType,ActiveNeighFlagType,DiscreteFlowElement,
					DiscreteFlowType,RecomputeType)
    Array<ActiveNeighFlagType,Dimension> activeNeighs;
    
//    struct DiscreteFlowElement {OffsetType offset; ScalarType weight;};
//    typedef CappedVector<DiscreteFlowElement, nActiveNeigh> DiscreteFlowType;
//    struct RecomputeType {ScalarType value,width;};
    RecomputeType Recompute(IndexCRef, DiscreteFlowType &) const;

    typedef ParamInterface_<PointType,VectorType> ParamInterface;
    typedef typename DifferenceType::MultiplierType MultiplierType;

    static const StencilStoragePolicy policy =
    DifferenceType::multSize>0 ? StencilStoragePolicy::Share :
    DifferenceType::multSize==0 ? StencilStoragePolicy::Recomp :
    StencilStoragePolicy::Lag2;
    template<StencilStoragePolicy, typename Dummy> struct _StencilDataType;
    typedef _StencilDataType<policy,void> StencilDataType;
    DiscreteType MaxStencilWidth() const;
    StencilDataType & stencilData;
    
    // Default domain parametrization (Currently used in all instantiations)
    static const bool hasBundle = (0<Traits::nStencilDependencies) && (Traits::nStencilDependencies<Dimension);
    template<int hasBundle_, typename Dummy> struct _ParamDefault;
    typedef _ParamDefault<int(hasBundle),void> ParamDefault;
    
    // Extra algorithms may be inserted at different points
    enum Decision {kAccept=0, kContinue=1, kTerminate=2, kRecompute=1<<10};
    struct ExtraAlgorithmInterface;
    struct ExtraAlgorithmPtrs;
    ExtraAlgorithmPtrs extras;
	typedef Factoring<Traits> FactoringType;
	mutable FactoringType factoring;
    
    HamiltonFastMarching(StencilDataType &);
protected:
    struct QueueElement;
    std::priority_queue<QueueElement> queue;
    
    void RunInit();
    bool RunOnce(); // returns true if should stop
    
    int PostProcess(IndexCRef); // returns a combination of decitions. Priority: kTerminate > kContinue > kAccept
    int PostProcessWithRecompute(IndexCRef, const RecomputeType &, const DiscreteFlowType &);
    
    Array<bool,Dimension> acceptedFlags;
    void ConditionalUpdate(IndexCRef,OffsetType,ScalarType);
    void Update(FullIndexCRef, OffsetCRef, ScalarType);
    friend struct _StencilDataType<policy, void>;
};

// ------- Some sub-structures of interest -------

template<typename T> struct HamiltonFastMarching<T>::ExtraAlgorithmInterface {
    virtual ~ExtraAlgorithmInterface(){};
    virtual void Setup(HFMI *){};
    virtual void Finally(HFMI *){};
    virtual bool ImplementIn(HFM *)=0; // Insert in adequate field of extras. Returns false if algorithm is vacuous.
protected:
	// Returns a Decision enum
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

typedef StencilStoragePolicy SSP;

// ******** Stencil data - multiplier based *********

template<typename T> template<typename Dummy>
struct HamiltonFastMarching<T>::_StencilDataType<SSP::Share,Dummy>{
    typedef HamiltonFastMarching<T> HFM;
    Redeclare6Types(HFM,IndexCRef,OffsetCRef,StencilType,DifferenceType,MultiplierType,Traits)
    Redeclare5Types(HFM,ParamInterface,HFMI,DomainType,FullIndexCRef,DistanceGuess)
    Redeclare6Types(Traits,DiscreteType,ScalarType,PointType,VectorType,IndexType,OffsetType)
    Redeclare1Type(StencilType,QuadType)
    Redeclare2Constants(Traits,Dimension,mathPi)
    
    typedef typename HFM::template DataSource<MultiplierType> MultSourceType;
    IndexType dims; // Needs value
    virtual void SetStencil(IndexCRef, StencilType &) = 0; // Needs specialization
    std::unique_ptr<MultSourceType> pMultSource = nullptr; // Needs assignment
    
    struct RecomputeDataType {const StencilType & stencil;MultiplierType mult;};
    RecomputeDataType RecomputeData(IndexCRef); // Needed in some subclasses of HFM...
    
    virtual void Setup(HFMI *);
    virtual const ParamInterface & Param() const = 0;
    virtual DistanceGuess GetGuess(IndexCRef) const {ExceptionMacro("Equation factoring error : no guess");};
	virtual ~_StencilDataType(){};
protected:
    friend struct HamiltonFastMarching<Traits>;
    void EraseCache(DiscreteType index) {shallowMultQuads.erase(index);}
    
    ScalarType HopfLaxUpdate(FullIndexCRef, OffsetCRef, ScalarType, ActiveNeighFlagType &);
    template<typename F> RecomputeType HopfLaxRecompute(const F &, IndexCRef, ActiveNeighFlagType, DiscreteFlowType &);
    RangeAccessor<OffsetType*> ReversedOffsets(FullIndexCRef);
    
    virtual void Initialize(const HFM *);
private:
    typedef typename HFM::template Array<StencilType, T::nStencilDependencies> StencilArrayType;
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
struct HamiltonFastMarching<T>::_StencilDataType<SSP::Recomp,Dummy>{
    typedef HamiltonFastMarching<T> HFM;
    Redeclare6Types(HFM,IndexCRef,OffsetCRef,StencilType,DifferenceType,Traits,FullIndexCRef)
    Redeclare5Types(HFM,ParamInterface,HFMI,MultiplierType,DomainType,DistanceGuess)
    Redeclare6Types(Traits,DiscreteType,ScalarType,PointType,VectorType,IndexType,OffsetType)
    Redeclare1Type(StencilType,QuadType)
    Redeclare2Constants(Traits,Dimension,mathPi)

    IndexType dims; // Needs value
    virtual void SetStencil(IndexCRef, StencilType &) = 0; // Needs specialization
    virtual void Setup(HFMI *);
    virtual const ParamInterface & Param() const = 0;
    struct RecomputeDataType {StencilType stencil; MultiplierType mult;}; // mult is dummy here
    RecomputeDataType RecomputeData(IndexCRef);
    virtual DistanceGuess GetGuess(IndexCRef) const {ExceptionMacro("Equation factoring error : no guess");};
	virtual ~_StencilDataType(){};
protected:
    friend struct HamiltonFastMarching<Traits>;
    void EraseCache(DiscreteType index) {shallowStencilQuads.erase(index);}
    
    ScalarType HopfLaxUpdate(FullIndexCRef, OffsetCRef, ScalarType, ActiveNeighFlagType &);
    template<typename F> RecomputeType HopfLaxRecompute(const F &, IndexCRef, ActiveNeighFlagType, DiscreteFlowType &);
    RangeAccessor<OffsetType*> ReversedOffsets(FullIndexCRef);

    virtual void Initialize(const HFM *);
private:
    ShallowMap<DiscreteType, std::pair<StencilType, QuadType> > shallowStencilQuads;
    std::vector<OffsetType> reversedOffsets;
    std::vector<DiscreteType> reversedOffsetsSplits;
};

// ********** Stencil data - Semi-Lagrangian2 *******

template<typename T> template<typename Dummy>
struct HamiltonFastMarching<T>::_StencilDataType<SSP::Lag2,Dummy>{
    typedef HamiltonFastMarching<T> HFM;
    Redeclare6Types(HFM,IndexCRef,OffsetCRef,StencilType,DifferenceType,Traits,FullIndexCRef)
    Redeclare5Types(HFM,ParamInterface,HFMI,MultiplierType,DomainType,DistanceGuess)
    Redeclare6Types(Traits,DiscreteType,ScalarType,PointType,VectorType,IndexType,OffsetType)
    Redeclare1Constant(Traits,Dimension)
    
    IndexType dims; // Needs value
    void SetStencil(IndexCRef, StencilType &);
    virtual void Setup(HFMI *);
    virtual const ParamInterface & Param() const = 0;
    virtual void Initialize(const HFM *);

    virtual void SetNeighbors(IndexCRef, std::vector<OffsetType> &) = 0;
	virtual DistanceGuess GetGuess(IndexCRef) const {ExceptionMacro("Equation factoring error : no guess");};
	virtual ~_StencilDataType(){};
protected:
    friend struct HamiltonFastMarching<Traits>;
    ScalarType HopfLaxUpdate(FullIndexCRef, OffsetCRef, ScalarType, ActiveNeighFlagType &);
    template<typename F> RecomputeType HopfLaxRecompute(const F &, IndexCRef, ActiveNeighFlagType, DiscreteFlowType &);
    RangeAccessor<OffsetType*> ReversedOffsets(FullIndexCRef);
    void EraseCache(DiscreteType index){};
    
    typedef CappedVector<std::pair<OffsetType,ScalarType>, 3> OffsetVal3;
    virtual std::pair<ScalarType,int> HopfLaxUpdate(IndexCRef,const OffsetVal3 &) = 0;
    virtual RecomputeType HopfLaxRecompute(IndexCRef,DiscreteFlowType &) = 0;
private:
    const HFM * pFM = nullptr;
    HFM::Array<ScalarType,Dimension> indexConverter;
    std::vector<OffsetType> reversedOffsets;
    std::vector<DiscreteType> reversedOffsetsSplits;

    std::vector<OffsetType> directOffsets;
    std::vector<DiscreteType> directOffsetsSplits;
};

// ********** Stencil data - Semi-Lagrangian3 *********

template<typename T> template<typename Dummy>
struct HamiltonFastMarching<T>::_StencilDataType<SSP::Lag3,Dummy>{
	typedef HamiltonFastMarching<T> HFM;
	
	IndexType dims; // Needs value
	void SetStencil(IndexCRef, StencilType &);
	virtual void Setup(HFMI *);
	virtual const ParamInterface & Param() const = 0;
	virtual void Initialize(const HFM *);
	
	virtual void SetNeighbors(IndexCRef, std::vector<OffsetType> &) = 0;
	virtual DistanceGuess GetGuess(IndexCRef) const {ExceptionMacro("Equation factoring error : no guess");};
	virtual ~_StencilDataType(){};
protected:
	// TODO
};

#include "HamiltonFastMarching.hxx"
#include "HFM_StencilDataType.hxx"
#include "HFM_ParamDefault.hxx"

#endif /* HamiltonFastMarching_h */
