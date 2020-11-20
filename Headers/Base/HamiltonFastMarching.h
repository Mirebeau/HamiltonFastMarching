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
#include "JMM_CPPLibs/DataStructures/CappedVector.h"
#include "JMM_CPPLibs/DataStructures/MapWeightedSum.h"
#include "JMM_CPPLibs/DataStructures/MultiVector.h"
#include "BaseGrid.h"
#include "JMM_CPPLibs/Macros/ExportArrow.h"
#include "JMM_CPPLibs/Macros/TemplateLog2.h"

template<typename T> struct HFMInterface;
template<typename T> struct StaticFactoring;

enum class FactoringMethod {None, Static, Dynamic};
enum class StencilStoragePolicy {Share,Recomp,Lag2,Lag3};

/// This structure provides default definitions for several traits of the models. Can be specialized.
template<typename Traits, typename Dummy=void> struct ExtraTraits {
    Redeclare1Constant(Traits,Dimension)
	Redeclare1Type(Traits,DifferenceType)
	
	static constexpr bool hasBundle =
	(0<Traits::nStencilDependencies) && (Traits::nStencilDependencies<Dimension);
		
	/* multSize<0 usually identifies semi-Lagrangian schemes, which have
	specific storage policies, and enjoy strict causality.*/
    static constexpr StencilStoragePolicy policy =
    DifferenceType::multSize>0 ? StencilStoragePolicy::Share :
    DifferenceType::multSize==0 ? StencilStoragePolicy::Recomp :
	Dimension==2 ? StencilStoragePolicy::Lag2 :
	StencilStoragePolicy::Lag3;
	
	/// Wether the scheme is strictly causal. Used when considering high order schemes.
	static constexpr bool strictlyCausal = DifferenceType::multSize<0;
	/// Wether the scheme uses singularity factoring at the first pass, before recomputation.
	static constexpr bool factorFirstPass =
	(policy==StencilStoragePolicy::Lag2 ||
	 policy==StencilStoragePolicy::Lag3);
};


template<typename TTraits>
struct HamiltonFastMarching {
    typedef TTraits Traits;
    Redeclare1Constant(Traits,Dimension)
    Redeclare12Types(Traits,ScalarType,DiscreteType,ShortType,DomainType,
					 StencilType,DistanceGuess,PointType,VectorType,IndexType,OffsetType,
					 DifferenceType,IndexDiff)
	Redeclare4Constants(ExtraTraits<Traits>,hasBundle,policy,strictlyCausal,
						factorFirstPass)
	
    using IndexCRef = const IndexType &;
    using OffsetCRef = const OffsetType &;
    struct FullIndexType {IndexType index; DiscreteType linear; PrintSelfMacro(FullIndexType)};
    using FullIndexCRef = const FullIndexType &;
    using HFM = HamiltonFastMarching;
    using HFMI = HFMInterface<Traits>;
    
    template<typename E, size_t n> using Array = typename Traits::template Array<E,n>;
    template<typename E> using DataSource = typename Traits::template DataSource<E>;
    template<size_t n> using BasisReduction=typename Traits::template BasisReduction<n>;

    // Domain and Running
    using DomainTransformType = typename DomainType::Transform;
    const DomainType dom;
    DomainTransformType VisibleOffset(IndexCRef, OffsetCRef, IndexType &) const;

    std::map<IndexType,ScalarType,typename IndexType::LexicographicCompare> seeds;
    void Run();
	Array<ScalarType,Dimension> values;

	// High order scheme extension.
	int order = 1;
//	static constexpr bool strictlyCausal = DifferenceType::multSize<0; // A.k.a. semi-Lagrangian.
	ScalarType maxRatioOrder2 = 0.5, maxRatioOrder3 = 0.25; // Fallback to low order. Req above
	
    // Geodesic related functions. Require values and activeNeighs to be correctly set,
    // either by running the algorithm or restoring them from a previous run
    struct GeodesicSolverInterface;
    struct FlowDataType;
    FlowDataType GeodesicFlow(IndexCRef) const;
    
    Redeclare4Types(StencilType,ActiveNeighFlagType,DiscreteFlowElement,
					DiscreteFlowType,RecomputeType)
    Array<ActiveNeighFlagType,Dimension> activeNeighs;
    
    RecomputeType Recompute(IndexCRef, DiscreteFlowType &) const;

    using ParamInterface = ParamInterface_<PointType,VectorType>;
    using MultiplierType = typename DifferenceType::MultiplierType;
	
	struct _StencilDataTypeBase;
    template<StencilStoragePolicy, typename Dummy> struct _StencilDataType;
    using StencilDataType = _StencilDataType<policy,void>;
    DiscreteType MaxStencilWidth() const;
    StencilDataType & stencilData;
    
    // Default domain parametrization (Currently used in all instantiations)
    template<int hasBundle_, typename Dummy> struct _ParamDefault;
    using ParamDefault = _ParamDefault<int(hasBundle),void>;
    
    // Extra algorithms may be inserted at different points
    enum Decision {kAccept=0, kContinue=1, kTerminate=2, kRecompute=1<<10};
    struct ExtraAlgorithmInterface;
    struct ExtraAlgorithmPtrs;
    ExtraAlgorithmPtrs extras;
	mutable StaticFactoring<Traits> factoring;
    
    HamiltonFastMarching(StencilDataType &);
	
	// Access to values around a point
	template<bool useFactoring, bool smallCorrection>
	void SetIndex(IndexCRef) const;
	template<bool useFactoring, bool smallCorrection, int maxOrder>
	ScalarType GetNeighborValue(OffsetType,int&) const;
//protected:
    struct QueueElement;
    std::priority_queue<QueueElement> queue;
    
    void RunInit();
    bool RunOnce(); // returns true if should stop
	
	// returns a combination of decisions. Priority: kTerminate > kContinue > kAccept
    int PostProcess(IndexCRef);
    int PostProcessWithRecompute(IndexCRef, const RecomputeType &, const DiscreteFlowType &);
    
    Array<bool,Dimension> acceptedFlags;
    void ConditionalUpdate(IndexCRef,OffsetType,ScalarType);
    void Update(FullIndexCRef, OffsetCRef, ScalarType);
    friend struct _StencilDataType<policy, void>;
	
private:
	struct {
		IndexType index;
		ScalarType value;
	} mutable getNeighborValue_tmp;
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

// ******** Stencil data, base class *********

template<typename T> struct HamiltonFastMarching<T>::_StencilDataTypeBase {
	IndexType dims;
	virtual ~_StencilDataTypeBase(){};

	virtual void SetStencil(IndexCRef, StencilType &) = 0; 
	virtual const ParamInterface & Param() const = 0;

	virtual DistanceGuess GetGuess(const PointType &) const {
		ExceptionMacro("Equation factoring error : no guess");};
	virtual DistanceGuess GetGuess(const IndexType & index) const {
		return GetGuess(pFM->dom.PointFromIndex(index));}
	
	virtual void Setup(HFMI * that){ // Import data from user interface
		dims = IndexType::CastCoordinates( that->io.template Get<PointType>("dims") );}
	
	using MultiOffsetVector = MultiVector<OffsetType,DiscreteType>;
	using ConstOffsetRange = typename MultiOffsetVector::ConstRangeType;
	virtual ConstOffsetRange ReversedOffsets(FullIndexCRef full) const {
		return reversedOffsets[full.linear];}
protected:
	friend struct HamiltonFastMarching<Traits>;
	virtual ScalarType HopfLaxUpdate(FullIndexCRef, OffsetCRef, ScalarType, ActiveNeighFlagType &) = 0;
	// Also, a templated HopfLaxRecompute
	virtual void EraseCache(DiscreteType index) {}
	virtual void Initialize(const HFM * _pFM){pFM = _pFM;} // Prepare for fast marching
	
	const HFM* pFM;
	MultiOffsetVector reversedOffsets;
	
};








#include "Implementation/HamiltonFastMarching.hxx"
#include "Implementation/HFM_StencilShare.hxx"
#include "Implementation/HFM_StencilRecompute.hxx"
#include "Implementation/HFM_StencilLag2.hxx"
#include "Implementation/HFM_StencilLag3.hxx"
#include "Implementation/HFM_ParamDefault.hxx"

#endif /* HamiltonFastMarching_h */
