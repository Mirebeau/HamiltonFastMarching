//
//  Lagrangian3Stencil.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 20/11/2018.
//

#ifndef Lagrangian3Stencil_h
#define Lagrangian3Stencil_h

#include "Specializations/CommonTraits.h"
#include "JMM_CPPLibs/Output/EnumToString.h"


// Possible additions : VoronoiRefined, or CubeRefined (several ways)
enum class Lagrangian3StencilGeometry {Cube,Voronoi};
template<> char const * enumStrings<Lagrangian3StencilGeometry>::data[] =
{"Cube","Voronoi"};

template<typename TOff, typename TScalar, typename TDiscrete>
struct Lagrangian3Stencil {
	using OffsetType = TOff;
	using ScalarType = TScalar;
	using DiscreteType = TDiscrete;
	using ShortType = typename OffsetType::ComponentType ;
	static const DiscreteType Dimension = OffsetType::Dimension;
	static_assert(Dimension==3,"Three dimensional stencil class");
	
	// Lists all the direct neighbors appearing in the scheme
	void Neighbors(std::vector<OffsetType> &) const;
	ShortType NeighborIndex(const OffsetType &) const;
	// Lists all the neighbors around a given one, given by its index
	void NeighborsAround(ShortType, std::vector<OffsetType> &) const;
	
	
	PrintSelfMacro(Lagrangian3Stencil);
	// Three main orientations of the stencil
	std::array<OffsetType,3> offsets;
	Lagrangian3StencilGeometry geom;
	
	struct ActiveNeighFlagType {
		// Indices to be used with NeighborsAround
		ShortType neighborIndex, sectorIndex;
		ActiveNeighFlagType():sectorIndex(-1){};
		bool none() const {return sectorIndex==-1;}
		ActiveNeighFlagType(long a)
		:neighborIndex(a & 255), sectorIndex((a>>8) & 255) {}
		unsigned long to_ulong() const {
			return (unsigned long)(neighborIndex) | ((unsigned long)(sectorIndex)<<8);}
	};

	static const int nActiveNeigh = Dimension;
	using CommonStencilType = CommonStencil<OffsetType,ScalarType,nActiveNeigh>;
	Redeclare3Types(CommonStencilType,DiscreteFlowElement,DiscreteFlowType,RecomputeType);
};

template<typename TO, typename TS, typename TD> void
Lagrangian3Stencil<TO,TS,TD>::PrintSelf(std::ostream & os) const {
	os << "{" <<
	offsets[0] << "," <<
	offsets[1] << "," <<
	offsets[2] << "," <<
	enumToString(geom) << "}";
}

template<typename TO, typename TS, typename TD> void
Lagrangian3Stencil<TO,TS,TD>::Neighbors(
	std::vector<OffsetType> & neigh) const {
	
	switch (geom) {
		case Lagrangian3StencilGeometry::Cube:{
			const OffsetType &
			u=offsets[0],
			v=offsets[1],
			w=offsets[2];
			neigh.reserve(neigh.size()+26);
			neigh.insert(neigh.end(),{
				u,-u,v,-v,w,-w,
				u+v,u-v,-u+v,-u-v,
				v+w,v-w,-v+w,-v-w,
				w+u,w-u,-w+u,-w-u,
				u+v+w,u+v-w,u-v+w,u-v-w,
				-u+v+w,-u+v-w,-u-v+w,-u-v-w
			});
			return;
		}
			
//		case Lagrangian3StencilGeometry::VoronoiRefined:
//			neigh.reserve(neigh.size()+18);
		case Lagrangian3StencilGeometry::Voronoi:{
//			const bool refined = geom==Lagrangian3StencilGeometry::VoronoiRefined;
			const OffsetType
			& a=offsets[0],
			& b=offsets[1],
			& c=offsets[2];
			const OffsetType d= -(a+b+c);
			neigh.insert(neigh.end(),{
				a,b,c,d,
				-a,-b,-c,-d,
				a+b,a+c,a+d,b+c,b+d,c+d
			});
//			if(refined){
//				const OffsetType x = a-b, y=c-d;
//				neigh.insert(neigh.end(),{x,-x,y,-y});
//			}
			return;
		}
			
		default:
			assert(false);
			return;
	}
}

template<typename TO, typename TS, typename TD> auto
Lagrangian3Stencil<TO,TS,TD>::NeighborIndex(const OffsetType & offset)
const -> ShortType {
	const OffsetType &
	u=offsets[0],
	v=offsets[1],
	w=offsets[2];
	assert(LinearAlgebra::Determinant(u,v,w) == 1);
	
	const auto
	a = LinearAlgebra::Determinant(offset,v,w),
	b = LinearAlgebra::Determinant(u,offset,w),
	c = LinearAlgebra::Determinant(u,v,offset);

	ShortType index=-1;
	
	switch(geom){
		case Lagrangian3StencilGeometry::Cube:{
			static const ShortType ind[3][3][3] = {
				{{25, 9, 24}, {17, 1, 15}, {23, 8, 22}},
				{{13, 3, 12}, {5, -1, 4}, {11, 2, 10}},
				{{21, 7, 20}, {16, 0, 14}, {19, 6, 18}}
			};
			index = ind[a+1][b+1][c+1];
			break;
		}
			
		case Lagrangian3StencilGeometry::Voronoi:{
			static const ShortType ind[3][3][3] = {
				{{3, 13, -1}, {12, 4, -1}, {-1, -1, -1}},
				{{10, 5, -1}, {6, -1, 2}, {-1, 1, 11}},
				{{-1, -1, -1}, {-1, 0, 9}, {-1, 8, 7}}
			};
			index = ind[a+1][b+1][c+1];
			break;
		}
			
		default:
			assert(false);
	}
	assert(index!=-1);
#ifdef DEBUG
	std::vector<OffsetType> offsets;
	Neighbors(offsets);
	assert(offsets[index]==offset);
#endif
	return index;
}


template<typename TO, typename TS, typename TD> void
Lagrangian3Stencil<TO,TS,TD>::NeighborsAround(
	ShortType index, std::vector<OffsetType> & neigh) const {
	switch (geom) {
		case Lagrangian3StencilGeometry::Cube:
			// i,j,k permutation of 0,1,2
			if(index<6){
				// Around a face center
				const int eps = index%2, i=index/2;
				const OffsetType
				u = (eps ? -1 : 1)*offsets[i];
				const OffsetType &
				v = offsets[(i+1)%3],
				w = offsets[(i+2)%3];
				neigh.insert(neigh.end(),
							 {u+v,u+v+w,u+w,u-v+w,u-v,u-v-w,u-w,u+v-w});
				return;
			}
			index-=6;
			
			if(index<12){
				// Around an edge center
				const int eps1=index%2; index/=2;
				const int eps0=index%2; index/=2;
				const int i=index;
				const OffsetType u = (eps0 ? -1 : 1)*offsets[i];
				const OffsetType v = (eps1 ? -1 : 1)*offsets[(i+1)%3];
				const OffsetType & w = offsets[(i+2)%3];
				neigh.insert(neigh.end(),
						{u,u+v+w,v,u+v-w});
				return;
			}
			index-=12;
	
			assert(index<8);
			if(true){
				// Around a cube vertex
				const int eps2=index%2; index/=2;
				const int eps1=index%2; index/=2;
				const int eps0=index;
				const OffsetType u = (eps0 ? -1 : 1)*offsets[0];
				const OffsetType v = (eps1 ? -1 : 1)*offsets[1];
				const OffsetType w = (eps2 ? -1 : 1)*offsets[2];
				neigh.insert(neigh.end(),
							 {u,u+v,v,v+w,w,w+u});
				return;
			}
			assert(false);
			
		case Lagrangian3StencilGeometry::Voronoi:{
			const OffsetType opp = -(offsets[0]+offsets[1]+offsets[2]);
			auto sb = [opp,this](int i)->OffsetType {
				return i<Dimension ? offsets[i]: opp;};
			
			if(index<=Dimension){
				//superbase a,u,v,w
				// around a
				const int i=index;
				const OffsetType
				a=sb(i),
				u=sb((i+1)%4),
				v=sb((i+2)%4),
				w=sb((i+3)%4);
				neigh.insert(neigh.end(),
//							 {a+u,a+u+v,a+v,a+v+w,a+w,a+w+u}
							 {a+u,-w,a+v,-u,a+w,-v});
				return;
			}
			index-=(Dimension+1);
			
			if(index<=Dimension){
				//superbase : a,u,v,w
				// around -a
				const int i=index;
				const OffsetType
//				a=sb[i],
				u=sb((i+1)%4),
				v=sb((i+2)%4),
				w=sb((i+3)%4);

				neigh.insert(neigh.end(),
							 {u,u+v,v,v+w,w,w+u});
				return;
			}
			index-=(Dimension+1);

			assert(index<6);
			if(true){
				// superbase a,b,c,d
				// around a+b
				typedef std::array<int,4> PermType;
				const std::array<PermType,6> perms = {{
					{{0,1,2,3}},{{0,2,1,3}},{{0,3,1,2}},
					{{1,2,0,3}},{{1,3,0,2}},{{2,3,0,1}}
				}};
				const PermType & p = perms[index];
				const OffsetType
				a=sb(p[0]),b=sb(p[1]),
				c=sb(p[2]),d=sb(p[3]);
				neigh.insert(neigh.end(),
							 //{a,a+b+c,b,a+b+d}
							 {a,-d,b,-c}
							 );
				return;
			}
			assert(false);
		}
			
//		case Lagrangian3StencilGeometry::VoronoiRefined:
//			assert(false);
			// TODO
			
		default:
			assert(false);
			return;
	}
}
#endif /* Lagrangian3Stencil_h */
