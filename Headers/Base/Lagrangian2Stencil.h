//
//  Lagrangian2Stencil.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 16/08/2018.
//

#ifndef Lagrangian2Stencil_h
#define Lagrangian2Stencil_h

#include "CommonStencil.h"
#include "JMM_CPPLibs/LinearAlgebra/VectorType.h"
#include "JMM_CPPLibs/Output/EnumToString.h"
#include "JMM_CPPLibs/Macros/DependentFalse.h"

// ------ Default stencils ------
enum class Lagrangian2StencilGeometry {Diamond,Square,SpikySquare,Voronoi,None};
template<> char const * enumStrings<Lagrangian2StencilGeometry>::data[] =
{"Diamond","Square","SpikySquare","Voronoi","None"};


// ----------- Semi-Lagrangian scheme ------------

template<typename TOffset, typename TScalar, typename TDiscrete>
struct Lagrangian2Stencil {
    using OffsetType = TOffset;
	using DiscreteType = TDiscrete;
	using ScalarType = TScalar;
    using ShortType = typename OffsetType::ComponentType;
    static constexpr DiscreteType Dimension = OffsetType::Dimension;
//	static_assert(Dimension==2,"Two dimensional stencil class");
    
	DiscreteType NSectors() const {return nOffsets;}
    const OffsetType & Sector(DiscreteType n, DiscreteType k) const {
        assert(0<=k && k<Dimension);
        assert(0<=n && n<NSectors()); // slightly abusive use case
        return pOffsets[(n+k)%nOffsets];
    }
    
    PrintSelfMacro(Lagrangian2Stencil);
    struct ActiveNeighFlagType {
		using SectorIndexType = ShortType;
        SectorIndexType sectorIndex;
        ActiveNeighFlagType():sectorIndex(-1){};
        bool none() const {return sectorIndex==-1;}
		explicit ActiveNeighFlagType(ScalarType a):sectorIndex(a){};
		explicit operator ScalarType() const {return sectorIndex;}		
    };
    
    static constexpr int nActiveNeigh = Dimension;
	using CommonStencilType = CommonStencil<OffsetType,ScalarType,nActiveNeigh>;
	Redeclare3Types(CommonStencilType,DiscreteFlowElement,DiscreteFlowType,RecomputeType);
	
    OffsetType * pOffsets;
    DiscreteType nOffsets;
	
	static std::vector<OffsetType>
	MakeStencil(Lagrangian2StencilGeometry,
				OffsetType=OffsetType::CanonicalBasis(0),
				OffsetType=OffsetType::CanonicalBasis(1));
};


template<typename TO, typename TS, typename TD> void
Lagrangian2Stencil<TO,TS,TD>::PrintSelf(std::ostream & os) const {
    os << "{";
    for(int i=0; i<nOffsets; ++i) os << pOffsets[i] << ",";
    os << "}";
}

template<typename TO, typename TS, typename TD> auto
Lagrangian2Stencil<TO,TS,TD>::MakeStencil(Lagrangian2StencilGeometry geom, OffsetType u, OffsetType v)
-> std::vector<OffsetType> {
	if constexpr(Dimension==2){
		assert(std::abs(LinearAlgebra::Determinant(u,v) )==1);
	} else {assert(false);}
	
	switch(geom){
		case Lagrangian2StencilGeometry::Diamond: return {u,v,-u,-v};
		case Lagrangian2StencilGeometry::Square: return {
			u,   u+v,  v,  v-u,
			-u,-(u+v),-v,-(v-u)};
		case Lagrangian2StencilGeometry::SpikySquare: return {
			u, 2*u+v, u+v, u+2*v,
			v, 2*v-u, v-u, v-2*u,
			-u,-2*u-v,-u-v,-u-2*v,
			-v,-2*v+u,-v+u,-v+2*u
		};
		case Lagrangian2StencilGeometry::Voronoi: return {
			u,   u+v,  v,
			-u,-(u+v),-v
		};
		default:
			ExceptionMacro("Lagrangian2Stencil::GetStencil error : unsupported geometry"
						   << enumToRealString(geom) << "\n")
	}
}

// Refines the orig stencil, according to the stop predicate, and stores the result in refined
// (with order reversed).
template<typename Predicate, typename VectorType>
void SternBrocotRefine(const Predicate & stop,
					   std::vector<VectorType> & refined, std::vector<VectorType> & orig){
	assert(!orig.empty());
	refined.push_back(orig.front());
	int size=orig.size();
	while(!orig.empty()){
		const VectorType & u = refined.back(), & v = orig.back();
		if(stop(u,v)){
			refined.push_back(v);
			orig.pop_back();
		} else {
			orig.push_back(u+v);
			++size;
			if(size>=127){ExceptionMacro("Stern-Brocot refine error : excessive stencil size. "
										 "Metric is either non-definite or too anisotropic.");}
		}
	}
	refined.pop_back();
}

// An optimized implementation, avoiding costly gradient recomputations,
// when the predicate is known to be cosAngleMin

template<typename NormType,
typename ScalarType,
typename OffsetType,
typename VectorType
> void
SternBrocotRefine_AcuteBound(const NormType & norm,
							 const ScalarType & cosAngleMin,
							 std::vector<OffsetType> & refined,
							 std::vector<OffsetType> & orig,
							 std::vector<VectorType> & grads,
							 std::vector<ScalarType> & norms) {
	assert(!orig.empty());
	assert(grads.empty());
	assert(norms.empty());
		
	auto grad = [&norm](const OffsetType & u) -> VectorType {
		return norm.Gradient(VectorType::CastCoordinates(u));};
	auto scal = [](const VectorType & u, const OffsetType & v) -> ScalarType {
		return u.ScalarProduct(VectorType::CastCoordinates(v));};
	
	grads.reserve(orig.size());
	for(const OffsetType & offset : orig){
		const VectorType g = grad(offset);
		grads.push_back(g);
		norms.push_back(scal(g,offset));
	}
	
	refined.push_back(orig.front());
	VectorType grad_u = grads.front();
	ScalarType norm_u = norms.front();
	
	// Refinement
	int size=orig.size();
	while(!orig.empty()){
		const OffsetType & u = refined.back(), & v = orig.back();
		const VectorType & grad_v = grads.back();
		const ScalarType & norm_v = norms.back();
		
		// CosAngleMin predicate
		const bool stop =
		scal(grad_u,v) >= cosAngleMin*norm_v &&
		scal(grad_v,u) >= cosAngleMin*norm_u;
		
		if(stop){
			refined.push_back(v);
			grad_u = grad_v;
			norm_u = norm_v;
			
			orig.pop_back();
			grads.pop_back();
			norms.pop_back();
		} else {
			const OffsetType w = u+v; // Inserted offset
			const VectorType grad_w = grad(w);
			const ScalarType norm_w = scal(grad_w,w);

			orig.push_back(w);
			grads.push_back(grad_w);
			norms.push_back(norm_w);
			
			++size;
			if(size>=127){ExceptionMacro("Stern-Brocot refine error : excessive stencil size. "
										 "Metric is either non-definite or too anisotropic.");}
		}
	}
	refined.pop_back();

	grads.clear();
	norms.clear();
}

// -----
/*
 // Alternative version using a forward list, slightly less efficient.
 template<typename TPred, typename TVec>
 void SternBrocotRefine(const TPred & stop, std::forward_list<TVec> & l){
 typedef TVec VectorType;
 assert(!l.empty());
 int size = 8; // There is no l.size(), but this is a typical upper bound.
 auto it0 = l.begin(), it1=it0; ++it1;
 for(; it1!=l.end(); ){
 const VectorType & u=*it0, & v=*it1;
 if(stop(u,v)){
 ++it0;
 ++it1;
 } else {
 it1 = l.insert_after(it0,u+v);
 ++size;
 if(size>=127){ExceptionMacro("Stern-Brocot refine error : excessive stencil size. "
 "Metric is either non-definite or too anisotropic.");}
 }
 }
 }
 */
#endif /* Lagrangian2Stencil_h */
