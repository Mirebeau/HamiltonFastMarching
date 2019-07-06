//
//  HFM_StencilLag3.hxx
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 18/02/2019.
//

#ifndef HFM_StencilLag3_hxx
#define HFM_StencilLag3_hxx

#include "Base/Lagrangian3Stencil.h"

// ********** Stencil data - Semi-Lagrangian3 *********

template<typename T> template<typename Dummy>
struct HamiltonFastMarching<T>::_StencilDataType<SSP::Lag3,Dummy>
: HamiltonFastMarching<T>::_StencilDataTypeBase {
	using HFM = HamiltonFastMarching<T>;
	using Superclass = HFM::_StencilDataTypeBase;
	Redeclare1Type(Superclass,  ConstOffsetRange)

	virtual void Setup(HFMI *) override;
	virtual void Initialize(const HFM *) override final;
	
	virtual ConstOffsetRange ReversedOffsets(FullIndexCRef) const override final;
protected:
	friend struct HamiltonFastMarching<Traits>;
	virtual ScalarType HopfLaxUpdate(FullIndexCRef, OffsetCRef, ScalarType, ActiveNeighFlagType &) override final;
	template<typename F> RecomputeType
	HopfLaxRecompute(const F &, IndexCRef, ActiveNeighFlagType, DiscreteFlowType &);
	
	typedef CappedVector<std::pair<OffsetType,ScalarType>, 9> OffsetVals;
	virtual std::pair<ScalarType,int> HopfLaxUpdate(IndexCRef,const OffsetVals &) = 0;
	virtual RecomputeType HopfLaxRecompute(IndexCRef,DiscreteFlowType &) = 0;
	Lagrangian3StencilGeometry geom = Lagrangian3StencilGeometry::None;
private:
	HFM::Array<StencilType,Dimension> stencils;
	std::vector<OffsetType> tmpHLU;
};

// -------------- Implementation --------


// Setup and initialization
template<typename Traits> template<typename Dummy> void
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag3, Dummy>::
Setup(HFMI * that){
	Superclass::Setup(that);
	stencils.dims = this->dims;
	
	if(that->io.HasField("stencilGeometry")){
		geom = enumFromString<Lagrangian3StencilGeometry>(that->io.GetString("stencilGeometry"));
		if( (int)geom==-1 ) ExceptionMacro("Stencil Lagrangian 3 error : unrecognized stencil Geometry");
	}
}

template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag3, Dummy>::
ReversedOffsets(FullIndexCRef full) const -> ConstOffsetRange {
	return this->reversedOffsets[ (geom==Lagrangian3StencilGeometry::None) ? full.linear : 0];
}

// Compute the direct and reversed offsets
template<typename Traits> template<typename Dummy> void
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag3, Dummy>::
Initialize(const HFM * pFM_){
	Superclass::Initialize(pFM_);
	
	if(geom!=Lagrangian3StencilGeometry::None){
		const auto stencil = StencilType::MakeStencil(geom);
		stencils.push_back(stencil);
		
		auto & rev = this->reversedOffsets;
		stencil.Neighbors(rev.values);
		rev.splits.push_back(rev.values.size());
		return;
	}
	
	const DiscreteType size = this->dims.Product();
	stencils.resize(size);
	std::vector<std::pair<DiscreteType,OffsetType> > targets;
	std::vector<OffsetType> directOffsets;
	
	for(DiscreteType linearIndex=0; linearIndex<size; ++linearIndex){
		const IndexType index = stencils.Convert(linearIndex);
		StencilType & stencil = stencils[linearIndex];
		this->SetStencil(index,stencil);
		
		directOffsets.clear();
		stencil.Neighbors(directOffsets);
		for(const OffsetType & offset : directOffsets){
			IndexType neighbor = index+IndexDiff::CastCoordinates(offset);
			if(this->pFM->dom.Periodize(neighbor,index).IsValid()){
				targets.push_back({stencils.Convert(neighbor), offset});}
		}
	}
	std::sort(targets.begin(),targets.end());
	assert(this->reversedOffsets.empty());
	this->reversedOffsets.insert(targets,this->dims.Product());
}


// HofLaxUpdate
template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag3, Dummy>::
HopfLaxUpdate(FullIndexCRef full, OffsetCRef acceptedOffset, ScalarType acceptedValue, ActiveNeighFlagType & active) -> ScalarType {
	const IndexType updatedIndex = full.index;
	
	// Get the index of the accepted offset
	const StencilType & stencil =
	geom==Lagrangian3StencilGeometry::None ? stencils[full.linear] : stencils.front();
	const ShortType neighborIndex = stencil.NeighborIndex(acceptedOffset);
	
	// Iterate over neighbors, and fill relative data
	OffsetVals offsetVal;
	const HFM & fm = *this->pFM;
	fm.template SetIndex<true,false>(updatedIndex); // useFactoring, smallCorrection
	
	tmpHLU.clear();
	stencil.NeighborsAround(neighborIndex, tmpHLU);
	for(auto & offset : tmpHLU){
		IndexType neigh = updatedIndex+IndexDiff::CastCoordinates(offset);
		const auto transform = fm.dom.Periodize(neigh,updatedIndex);
		if(transform.IsValid()){
			const int neighLin = fm.values.Convert(neigh);
			if(fm.acceptedFlags[neighLin]){
				
				//offsetVal.push_back({offset,fm.values[neighLin]}); // No factorization
				int ord=1;
				offsetVal.push_back({offset,
					fm.template GetNeighborValue<true,false,1>(offset,ord)});
				assert(ord==1);
				continue;
			}
		}
		offsetVal.push_back({offset,Traits::Infinity()});
	}
	
	// Put the data of the accepted value
	int ord=1;
	offsetVal.push_back({acceptedOffset,
		fm.template GetNeighborValue<true,false,1>(acceptedOffset,ord)});
	assert(ord==1);
//	offsetVal.push_back({acceptedOffset, acceptedValue});
		
	// Evaluate the Hopf-Lax update, and compare
	const auto [value,sectorIndex] = HopfLaxUpdate(updatedIndex,offsetVal);
	const ScalarType oldValue = fm.values[full.linear];
	
	/*
	if(full.index==IndexType{0,0,0}){
	std::cout << "In HLU "
	ExportVarArrow(full.index)
	ExportVarArrow(acceptedOffset)
	ExportVarArrow(value)
	ExportVarArrow(oldValue)
	ExportArrayArrow(tmpHLU)
	ExportArrayArrow(offsetVal)
		ExportVarArrow(this->pFM->values(updatedIndex+IndexDiff::CastCoordinates(acceptedOffset)))
	<< std::endl;
		}*/
		
	
	if(value>=oldValue) return oldValue;
	
	active.neighborIndex = neighborIndex;
	active.sectorIndex = sectorIndex;
	return value;
}

// HopfLaxRecompute
template<typename Traits> template<typename Dummy> template<typename F> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag3, Dummy>::
HopfLaxRecompute(const F & f, IndexCRef updatedIndex, ActiveNeighFlagType active, DiscreteFlowType & discreteFlow) -> RecomputeType {
	assert(!active.none());
	const StencilType & stencil =
	geom==Lagrangian3StencilGeometry::None ? stencils(updatedIndex) : stencils.front();
	
	// Get the offsets toward the active neighbors
	// TODO : check if these array allocations are perf limiting, and act if necessary.
	tmpHLU.clear();
	stencil.Neighbors(tmpHLU);
	const OffsetType offset0 = tmpHLU[active.neighborIndex];
	
	tmpHLU.clear();
	stencil.NeighborsAround(active.neighborIndex,tmpHLU);
	const OffsetType offset1 = tmpHLU[active.sectorIndex];
	const OffsetType offset2 = tmpHLU[(active.sectorIndex+1)%tmpHLU.size()];
	
	// Get the values at the active neighbors, and their approximation order
	int ord0 = 3;
	ScalarType val0 = f(offset0,ord0);
	assert(ord0!=0); // This value must be valid, hence with positive approx order
	int ord1 = ord0;
	ScalarType val1 = f(offset1,ord1);
	int ord2 = ord1!=0 ? ord1 : ord0;
	ScalarType val2 = f(offset2,ord2);
	
	// All neighbor values must be at the same approximation order
	const int ord = ord2!=0 ? ord2 : ord1!=0 ? ord1 : ord0;
	if(ord<ord0) {ord0=ord; val0 = f(offset0,ord0); assert(ord0==ord);}
	if(ord<ord1) {ord1=ord; val1 = f(offset1,ord1); assert(ord1==ord);}
	assert(ord==ord2 || ord2==0);
	
	assert(discreteFlow.empty());
	const std::array<ScalarType,4> mult = {0.,1.,3./2.,11./6.};
	discreteFlow.push_back({offset0,mult[ord]*val0});
	if(ord1>0) {discreteFlow.push_back({offset1,mult[ord]*val1});}
	if(ord2>0) {discreteFlow.push_back({offset2,mult[ord]*val2});}
	
	RecomputeType result = HopfLaxRecompute(updatedIndex,discreteFlow);
	const std::array<ScalarType,4> div = {0.,1.,2./3.,6./11.};
	result.width*=div[ord]; result.value*=div[ord];
	
	/*
	if(updatedIndex==IndexType{0,0,0}){
	std::cout << "In HL recompute, stencil 3 "
	ExportVarArrow(updatedIndex)
	ExportVarArrow(val0)
	ExportVarArrow(val1)
	ExportVarArrow(val2)
	ExportVarArrow(offset0)
	ExportVarArrow(offset1)
	ExportVarArrow(offset2)
	ExportVarArrow(result.value)
	ExportVarArrow(active)
	ExportArrayArrow(tmpHLU)
	ExportVarArrow(this->pFM->values(updatedIndex+IndexDiff::CastCoordinates(offset0)))
	ExportVarArrow(this->pFM->values(updatedIndex+IndexDiff::CastCoordinates(offset1)))
	ExportVarArrow(this->pFM->values(updatedIndex+IndexDiff::CastCoordinates(offset2)))
	<< std::endl;
		}*/
	
	return result;
};

#endif /* HFM_StencilLag3_hxx */
