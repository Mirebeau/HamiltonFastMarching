//
//  CommonStencil.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 18/02/2019.
//

#ifndef CommonStencil_h
#define CommonStencil_h
#include "JMM_CPPLibs/DataStructures/CappedVector.h"
#include "JMM_CPPLibs/Macros/ExportArrow.h"


template<typename TOffset, typename TScalar, int nActiveNeigh>
struct CommonStencil {
	using OffsetType = TOffset;
	using ScalarType = TScalar;
	static_assert(nActiveNeigh>0);

	struct DiscreteFlowElement {OffsetType offset; ScalarType weight;
		PrintSelfMacro(DiscreteFlowElement);};
	using DiscreteFlowType = CappedVector<DiscreteFlowElement, nActiveNeigh>;
	struct RecomputeType {ScalarType value,width; PrintSelfMacro(RecomputeType);};
};


template<typename TO, typename TS, int nA> void
CommonStencil<TO,TS,nA>::DiscreteFlowElement::PrintSelf(std::ostream & os) const {
	os << "{" << offset << "," << weight << "}";
}

template<typename TO, typename TS, int nA> void
CommonStencil<TO,TS,nA>::RecomputeType::PrintSelf(std::ostream & os) const {
	os << "{" << value << "," << width << "}";
}


#endif /* CommonStencil_h */
