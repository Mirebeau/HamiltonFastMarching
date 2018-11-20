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


enum class Lagrangian3StencilGeometry {
	Cube,Voronoi,VoronoiRefined};
template<> char const * enumStrings<Lagrangian3StencilGeometry>::data[] =
{"Cube","Voronoi","VoronoiRefined"};

template<typename TOff, typename TScalar, typename TDiscrete>
struct Lagrangian3Stencil {
	typedef TOff OffsetType;
	typedef TScalar ScalarType;
	typedef TDiscrete DiscreteType;
	typedef typename OffsetType::ComponentType ShortType;
	static const DiscreteType Dimension = OffsetType::Dimension;
	static_assert(Dimension==3,"Three dimensional stencil class");
	
	void Neighbors(std::vector<OffsetType> &) const;
	void NeighborsAround(ShortType, std::vector<OffsetType> &) const;
	
	PrintSelfMacro(Lagrangian3Stencil);
	OffsetType * pOffsets = nullptr;
	Lagrangian3StencilGeometry geom;
};

template<typename TO, typename TS, typename TD> void
Lagrangian3Stencil<TO,TS,TD>::PrintSelf(std::ostream & os) const {
	assert(pOffsets!=nullptr);
	os << "{" <<
	pOffsets[0] << "," <<
	pOffsets[1] << "," <<
	pOffsets[2] << "," <<
	geom << "}";
}

template<typename TO, typename TS, typename TD> void
Lagrangian3Stencil<TO,TS,TD>::Neighbors(
	std::vector<OffsetType> & neigh) const {
	
	switch (geom) {
		case Lagrangian3StencilGeometry::Cube:{
			const OffsetType &
			u=pOffsets[0],
			v=pOffsets[0],
			w=pOffsets[0];
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
			
		case Lagrangian3StencilGeometry::VoronoiRefined:
			neigh.reserve(neigh.size()+18);
		case Lagrangian3StencilGeometry::Voronoi:{
			const bool refined = geom==Lagrangian3StencilGeometry::VoronoiRefined;
			const OffsetType
			& a=pOffsets[0],
			& b=pOffsets[1],
			& c=pOffsets[2];
			const OffsetType d= -(a+b+c);
			neigh.insert(neigh.end(),{
				a,b,c,d,
				-a,-b,-c,-d,
				a+b,a+c,a+d,b+c,b+d,c+d
			});
			if(refined){
				const OffsetType x = a-b, y=c-d;
				neigh.insert(neigh.end(),{x,-x,y,-y});
			}
			return;
		}
			
		default:
			assert(false);
			return;
	}
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
				const OffsetType u = (eps ? -1 : 1)*pOffsets[i];
				const OffsetType &
				v = pOffsets[(i+1)%3],
				w = pOffsets[(i+2)%3];
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
				const OffsetType u = (eps0 ? -1 : 1)*pOffsets[i];
				const OffsetType v = (eps1 ? -1 : 1)*pOffsets[(i+1)%3];
				const OffsetType & w = pOffsets[(i+2)%3];
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
				const OffsetType u = (eps0 ? -1 : 1)*pOffsets[0];
				const OffsetType v = (eps1 ? -1 : 1)*pOffsets[1];
				const OffsetType w = (eps2 ? -1 : 1)*pOffsets[2];
				neigh.insert(neigh.end(),
							 {u,u+v,v,v+w,w,w+u});
				return;
			}
			assert(false);
			
		case Lagrangian3StencilGeometry::Voronoi:{
			const OffsetType opp = -(pOffsets[0]+pOffsets[1]+pOffsets[2]);
			auto sb = [opp,this](int i)->OffsetType {
				return i<Dimension ? pOffsets[i]: opp;};
			
			if(index<=Dimension){
				//superbase a,u,v,w
				// around a
				const int i=index;
				const OffsetType
				a=sb[i],
				u=sb[(i+1)%4],
				v=sb[(i+2)%4],
				w=sb[(i+3)%4];
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
				u=sb[(i+1)%4],
				v=sb[(i+2)%4],
				w=sb[(i+3)%4];

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
				a=sb[p[0]],b=sb[p[1]],
				c=sb[p[2]],d=sb[p[3]];
				neigh.insert(neigh.end(),
							 //{a,a+b+c,b,a+b+d}
							 {a,-d,b,-c}
							 );
				return;
			}
			assert(false);
		}
			
		case Lagrangian3StencilGeometry::VoronoiRefined:
			assert(false);
			// TODO
			
		default:
			assert(false);
			return;
	}
}
#endif /* Lagrangian3Stencil_h */
