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
/*
// Note : walls must not be taken into account at initialization with mult based.
 
 Failsafe mechanism enableMaxRatioOrder, for high order schemes, requires an estimate of the offset primal norm. For now, we use a rough estimate, requiring the strict causality of the scheme.
 */

template<typename T> struct HFMInterface;
template<typename T> struct Factoring;

enum class StencilStoragePolicy {Share,Recomp,Lag2,Lag3};
enum class FactoringMethod {None, Static, Dynamic};

template<typename TTraits>
struct HamiltonFastMarching {
    typedef TTraits Traits;
    Redeclare1Constant(Traits,Dimension)
    Redeclare12Types(Traits,ScalarType,DiscreteType,ShortType,DomainType,
					 StencilType,DistanceGuess,PointType,VectorType,IndexType,OffsetType,
					 DifferenceType,IndexDiff)
	
    typedef const IndexType & IndexCRef;
    typedef const OffsetType & OffsetCRef;
    struct FullIndexType {IndexType index; DiscreteType linear; PrintSelfMacro(FullIndexType)};
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
	Array<ScalarType,Dimension> values;

	// High order scheme extension.
	int order = 1;
	static const bool strictlyCausal = DifferenceType::multSize<0; // A.k.a. semi-Lagrangian.
	ScalarType maxRatioOrder2 = 0.5, maxRatioOrder3 = 0.25; // Fallback to low order. Req above
	
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
	Dimension==2 ? StencilStoragePolicy::Lag2 :
	StencilStoragePolicy::Lag3;
	
	struct _StencilDataTypeBase;
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
	
	// Access to values around a point
	template<bool useFactoring, bool smallCorrection>
	void SetIndex(IndexCRef) const;
	template<bool useFactoring, bool smallCorrection, int maxOrder>
	ScalarType GetNeighborValue(OffsetType,int&) const;
protected:
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
