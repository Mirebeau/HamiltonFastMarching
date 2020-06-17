//
//  HFM_StencilLag2.hxx
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 18/02/2019.
//

#ifndef HFM_StencilLag2_hxx
#define HFM_StencilLag2_hxx

#include "Base/Lagrangian2Stencil.h"

/*
 This class implements a stencil data structure for the Hamilton-Fast-Marching class,
 with the following features :
 - Stencils have a variable size, and involve only offsets. (Semi-Lagrangian type.)
 - Each point has an independent stencil.
 - They are computed only once.
 */

// ********** Stencil data - Semi-Lagrangian2 *******

template<typename T> template<typename Dummy>
struct HamiltonFastMarching<T>::_StencilDataType<SSP::Lag2,Dummy> :
HamiltonFastMarching<T>::_StencilDataTypeBase {
	using HFM=HamiltonFastMarching<T>;
	using Superclass = HFM::_StencilDataTypeBase;
	Redeclare17Types(HFM,IndexCRef,OffsetCRef,StencilType,DifferenceType,Traits,FullIndexCRef,
					ParamInterface,HFMI,MultiplierType,DomainType,DistanceGuess,
					DiscreteType,ScalarType,PointType,VectorType,IndexType,OffsetType)
	Redeclare1Type(Superclass,  ConstOffsetRange)
	Redeclare1Constant(Traits,Dimension)
	
	virtual void Setup(HFMI *) override;
	virtual void Initialize(const HFM *) override final;
	
	virtual void SetStencil(IndexCRef, StencilType &) override final;
	virtual void SetNeighbors(IndexCRef, std::vector<OffsetType> &) = 0;
	virtual ConstOffsetRange ReversedOffsets(FullIndexCRef) const override final;
protected:
	friend struct HamiltonFastMarching<Traits>;
	virtual ScalarType HopfLaxUpdate(FullIndexCRef, OffsetCRef, ScalarType, ActiveNeighFlagType &) override;
	template<typename F> RecomputeType
	HopfLaxRecompute(const F &, IndexCRef, ActiveNeighFlagType, DiscreteFlowType &);
	virtual RecomputeType // Only in TTI
	_HopfLaxRecompute(IndexCRef,ActiveNeighFlagType,DiscreteFlowType &) {
		assert(false); return RecomputeType{};}
	
	Lagrangian2StencilGeometry geom = Lagrangian2StencilGeometry::None;
	using OffsetVal3 = CappedVector<std::pair<OffsetType,ScalarType>, 3>;
	virtual std::pair<ScalarType,int> HopfLaxUpdate(IndexCRef,const OffsetVal3 &) = 0;
	virtual RecomputeType HopfLaxRecompute(IndexCRef,DiscreteFlowType &) = 0;
	HFM::Array<ScalarType,Dimension> indexConverter;
	Redeclare1Type(ActiveNeighFlagType, SectorIndexType)
private:
	MultiVector<OffsetType,DiscreteType> directOffsets;
};

// --------------- Semi-Lagrangian FM-ASR scheme --------------

template<typename Traits> template<typename Dummy> void
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag2, Dummy>::
SetStencil(IndexCRef index, StencilType & stencil){
	const auto rg =
	directOffsets[ (geom==Lagrangian2StencilGeometry::None) ? indexConverter.Convert(index) : 0];
	stencil.pOffsets = &rg.front();
	stencil.nOffsets = rg.size();
}

template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag2, Dummy>::
ReversedOffsets(FullIndexCRef full) const -> ConstOffsetRange {
	return this->reversedOffsets[ (geom==Lagrangian2StencilGeometry::None) ? full.linear : 0];
}


// Setup and initialization
template<typename Traits> template<typename Dummy> void
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag2, Dummy>::
Setup(HFMI * that){
	Superclass::Setup(that);
	indexConverter.dims = this->dims;

	if(that->io.HasField("stencilGeometry")){
		geom = enumFromString<Lagrangian2StencilGeometry>(that->io.GetString("stencilGeometry"));
		if( (int)geom==-1 ) ExceptionMacro("Stencil Lagrangian 2 error : unrecognized stencil Geometry");
	}
}

// Compute the direct and reversed offsets
template<typename Traits> template<typename Dummy> void
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag2, Dummy>::
Initialize(const HFM * _pFM){
	Superclass::Initialize(_pFM);
	
	if(geom!=Lagrangian2StencilGeometry::None){
		const auto stencil = StencilType::MakeStencil(geom);
		auto & dir = directOffsets, & rev = this->reversedOffsets;
		const size_t size = stencil.size();
		dir.values.reserve(size);
		rev.values.reserve(size);
		dir.values.insert(dir.values.end(),stencil.begin(),stencil.end());
		rev.values.insert(rev.values.end(),stencil.begin(),stencil.end());
		dir.splits.push_back(size);
		rev.splits.push_back(size);
		return;
	}
	
	const DiscreteType size = this->dims.Product();
	
	assert(directOffsets.empty());
	directOffsets.splits.reserve(size+1);
	std::vector<std::pair<DiscreteType,OffsetType> > targets;
	
	for(DiscreteType linearIndex=0; linearIndex<size; ++linearIndex){
		const IndexType index = indexConverter.Convert(linearIndex);
		SetNeighbors(index,directOffsets.values);
		directOffsets.splits.push_back((DiscreteType)directOffsets.values.size());
		
		for(const OffsetType & offset : directOffsets.back()){
			IndexType neighbor = index+IndexDiff::CastCoordinates(offset);
			if(this->pFM->dom.Periodize(neighbor,index).IsValid()){
				targets.push_back({indexConverter.Convert(neighbor), offset});}
		}
	}
	
	std::sort(targets.begin(),targets.end());
	assert(this->reversedOffsets.empty());
	this->reversedOffsets.insert(targets,size);
}

// HofLaxUpdate
template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag2, Dummy>::
HopfLaxUpdate(FullIndexCRef full, OffsetCRef acceptedOffset,
			  ScalarType __acceptedValue__, ActiveNeighFlagType & active) -> ScalarType {
	const IndexType updatedIndex = full.index;
	
	// Get the stencil
	StencilType stencil;
	SetStencil(updatedIndex, stencil);
	
	// Get the relevant sectors
	std::array<int,2> sectors;
	{
		int i = 0;
		for(;stencil.Sector(i,0)!=acceptedOffset; ++i){
			assert(i<stencil.NSectors());}
		sectors[0] = i;
		sectors[1] = i==0 ? stencil.NSectors()-1 : i-1;
	}
	
	// Get the relevant offsets, in addition to acceptedOffset
	const std::array<OffsetType,2> neighOffsets {
		stencil.Sector(sectors[0],1),
		stencil.Sector(sectors[1],0)
	};
	
	/*
	// A failed experiment. TODO : remove ?
	constexpr bool useSecondOrderOnFirstPass = false;
	if constexpr(useSecondOrderOnFirstPass){ // -----------------------
		auto & fm = *(this->pFM);
		fm.template SetIndex<true,true>(updatedIndex); // useFactoring, smallCorrection
		int order=2;
		// Corrected value (with second order and factoring), hides function parameter
		const ScalarType acceptedValue =
		fm.template GetNeighborValue<true,true,2>(acceptedOffset,order);
		assert(order>=1);
		
		ScalarType oldValue = fm.values[full.linear];
		
		for(int i=0; i<2; ++i){
			int ord = order;
			const ScalarType neighVal =
			fm.template GetNeighborValue<true,true,2>(neighOffsets[i],ord);
			if(ord==0) continue;
			const ScalarType acceptedVal =
			ord==order ? acceptedValue :
			fm.template GetNeighborValue<true,true,1>(acceptedOffset,ord);
			OffsetVal3 offsetVal{ {acceptedOffset,acceptedVal}, {neighOffsets[i], neighVal} };
			
			const auto [newValue,newActive] = HopfLaxUpdate(updatedIndex,offsetVal);
			if(newValue<oldValue){
				oldValue = newValue;
				active.sectorIndex = sectors[i];
			}
		}
		return oldValue;
		
	} else { 
	*/
	// ------- do not use second order on first pass... -----------
	auto & fm = *(this->pFM);
	fm.template SetIndex<true,false>(updatedIndex); // useFactoring, smallCorrection
	int order=1;

	// Get accepted value, with factoring correction
	OffsetVal3 offsetVal;
	offsetVal.push_back({acceptedOffset,
		fm.template GetNeighborValue<true,false,1>(acceptedOffset,order)});
	assert(order==1);
	
	std::array<ShortType, 3> act;
	act[0] = sectors[0];
	
	// Get neighbor values
	for(int i=0; i<2; ++i){
		int ord=order;
		const ScalarType neighVal =
		fm.template GetNeighborValue<true,false,1>(neighOffsets[i],ord);
		if(ord==1){
			act[offsetVal.size()]=sectors[i];
			offsetVal.push_back({neighOffsets[i],neighVal});
		}
	}
	const auto [newValue,newActive] = HopfLaxUpdate(updatedIndex,offsetVal);
	const ScalarType oldValue = fm.values[full.linear];
	if(newValue>=oldValue) return oldValue;
	
	if constexpr(std::is_same_v<void,SectorIndexType>){assert(false);}
	else {active.sectorIndex = act[newActive];}
	return newValue;
}

// HopfLaxRecompute
template<typename Traits> template<typename Dummy> template<typename F> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag2, Dummy>::
HopfLaxRecompute(const F & f, IndexCRef index, ActiveNeighFlagType active,
				 DiscreteFlowType & discreteFlow) -> RecomputeType {
	
	if constexpr(std::is_same_v<void,SectorIndexType>){ // -> TTINorm
		return _HopfLaxRecompute(index, active, discreteFlow);}
	
	assert(!active.none());
	StencilType stencil;
	SetStencil(index, stencil);
	
	int sectorIndex;
	if constexpr(std::is_same_v<void,SectorIndexType>){assert(false); sectorIndex=0;}
	else {sectorIndex = active.sectorIndex;}
	
	// Get data of neighbor 0
	int ord0 = 3;
	const OffsetType offset0 = stencil.Sector(sectorIndex,0);
	ScalarType val0 = f(offset0,ord0);
	
	
	// Get data of neighbor 1
	int ord1 = ord0==0 ? 3 : ord0;
	const OffsetType offset1 = stencil.Sector(sectorIndex,1);
	const ScalarType val1 = f(offset1,ord1);
	
	
	// Semi-Lagrangian scheme requires same scheme order for all neighbor values
	if(0<ord1 && ord1<ord0){
		ord0=ord1;
		val0 = f(offset0,ord0);
		assert(ord0==ord1);
	}
	
	assert(ord0>=0 && ord1>=0);
	assert(ord0+ord1>0);
	const int ord = std::max(ord0,ord1); // snd order used ?
	
	assert(discreteFlow.empty());
	const std::array<ScalarType,4> mult = {0.,1.,3./2.,11./6.};
	if(ord0>0) {discreteFlow.push_back({offset0,mult[ord]*val0});}
	if(ord1>0) {discreteFlow.push_back({offset1,mult[ord]*val1});}
	
	RecomputeType result = HopfLaxRecompute(index,discreteFlow);
	const std::array<ScalarType,4> div = {0.,1.,2./3.,6./11.};
	result.width*=div[ord]; result.value*=div[ord];
	
	return result;
};

#endif /* HFM_StencilLag2_hxx */
