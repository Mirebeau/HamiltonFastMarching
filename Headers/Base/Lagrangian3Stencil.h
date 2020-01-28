//
//  Lagrangian3Stencil.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 20/11/2018.
//

#ifndef Lagrangian3Stencil_h
#define Lagrangian3Stencil_h

#include "JMM_CPPLibs/LinearAlgebra/VectorType.h"
#include "JMM_CPPLibs/Output/EnumToString.h"

#include "CommonStencil.h"

// Possible additions : VoronoiRefined, or CubeRefined (several ways)
enum class Lagrangian3StencilGeometry {Diamond,CutCube,Cube,SpikyCube,Voronoi,None};
template<> char const * enumStrings<Lagrangian3StencilGeometry>::data[] =
{"Diamond","CutCube","Cube","SpikyCube","Voronoi","None"};

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
	Lagrangian3StencilGeometry geom;
	std::array<OffsetType,3> offsets;
	
	struct ActiveNeighFlagType {
		// Indices to be used with NeighborsAround
		ShortType neighborIndex, sectorIndex;
		ActiveNeighFlagType():sectorIndex(-1){};
		bool none() const {return sectorIndex==-1;}
		explicit ActiveNeighFlagType(ScalarType a){
			unsigned long u(a);
			neighborIndex= u & 255;
			sectorIndex = (u>>8) & 255;
		}
		explicit operator ScalarType() const {
			using uShort = std::make_unsigned_t<ShortType>;
			using uLong = unsigned long;
			return uLong(uShort(neighborIndex)) | uLong(uShort(sectorIndex)<<8);}
		PrintSelfMacro(ActiveNeighFlagType);
	};

	static const int nActiveNeigh = Dimension;
	using CommonStencilType = CommonStencil<OffsetType,ScalarType,nActiveNeigh>;
	Redeclare3Types(CommonStencilType,DiscreteFlowElement,DiscreteFlowType,RecomputeType);
	
	bool IsValid() const;
	
	static Lagrangian3Stencil MakeStencil(Lagrangian3StencilGeometry geom,
										  OffsetType u={1,0,0}, OffsetType v={0,1,0}, OffsetType w={0,0,1}){
		return Lagrangian3Stencil{geom,{u,v,w}};}
};

template<typename TO, typename TS, typename TD> void
Lagrangian3Stencil<TO,TS,TD>::PrintSelf(std::ostream & os) const {
	os << "{" <<
	offsets[0] << "," <<
	offsets[1] << "," <<
	offsets[2] << "," <<
	enumToString(geom) << "}";
}

template<typename TO, typename TS, typename TD> bool
Lagrangian3Stencil<TO,TS,TD>::IsValid() const {
	return geom!=Lagrangian3StencilGeometry::None
	&& LinearAlgebra::Determinant(offsets[0],offsets[1],offsets[2])==1;
}


template<typename TO, typename TS, typename TD> void
Lagrangian3Stencil<TO,TS,TD>::ActiveNeighFlagType::PrintSelf(std::ostream & os) const {
	os << "{" << (int)neighborIndex << "," << (int)sectorIndex << "}";
}

template<typename TO, typename TS, typename TD> void
Lagrangian3Stencil<TO,TS,TD>::Neighbors(
	std::vector<OffsetType> & neigh) const {
	assert(IsValid());
	
	// Voronoi class stencil
	// TODO : Voronoi refined
	if(geom==Lagrangian3StencilGeometry::Voronoi){
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
		
	} else {
		// ------ Cube class stencil ------
		const int stencilSize =
		geom == Lagrangian3StencilGeometry::Diamond ? 6 :
		geom == Lagrangian3StencilGeometry::CutCube ? 18 :
		geom == Lagrangian3StencilGeometry::Cube ? 26 :
		50; // geom == Lagrangian3StencilGeometry::SpikyCube
		neigh.reserve(neigh.size()+stencilSize);
		
		for(int i=0; i<3; ++i){
			for(int s : {-1,1}){
				neigh.push_back(s*offsets[i]);
			}
		}
		if(geom == Lagrangian3StencilGeometry::Diamond) return;

		for(int i=0; i<3; ++i){
			const int j=(i+1)%3;
			for(int si : {-1,1}){
				for(int sj : {-1,1}){
					neigh.push_back(si*offsets[i]+sj*offsets[j]);
				}
			}
		}
		if(geom == Lagrangian3StencilGeometry::CutCube) return;

		for(int s0 : {-1,1}){
			for(int s1 : {-1,1}){
				for(int s2 : {-1,1}){
					neigh.push_back(s0*offsets[0]+s1*offsets[1]+s2*offsets[2]);
				}
			}
		}
		if(geom == Lagrangian3StencilGeometry::Cube) return;

		for(int i=0; i<3; ++i){
			const int j=(i+1)%3, k=(i+2)%3;
			for(int si : {-1,1}){
				for(int sj : {-1,1}){
					for(int sk : {-1,1}){
						neigh.push_back(2*si*offsets[i]+sj*offsets[j]+sk*offsets[k]);
					}
				}
			}
		}
		assert(geom == Lagrangian3StencilGeometry::SpikyCube); return;

		
	}
	/*
			
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
	 */
}

template<typename TO, typename TS, typename TD> auto
Lagrangian3Stencil<TO,TS,TD>::NeighborIndex(const OffsetType & offset)
const -> ShortType {
	assert(IsValid());
	
	const OffsetType &
	u=offsets[0],
	v=offsets[1],
	w=offsets[2];
	
	const auto
	a = LinearAlgebra::Determinant(offset,v,w),
	b = LinearAlgebra::Determinant(u,offset,w),
	c = LinearAlgebra::Determinant(u,v,offset);
	
	assert(std::abs(a)<=2 && std::abs(b)<=2 && std::abs(c)<=2);

	ShortType index=-1;
	
	if(geom==Lagrangian3StencilGeometry::Voronoi){
		static const ShortType ind[3][3][3] = {
			{{3, 13, -1}, {12, 4, -1}, {-1, -1, -1}},
			{{10, 5, -1}, {6, -1, 2}, {-1, 1, 11}},
			{{-1, -1, -1}, {-1, 0, 9}, {-1, 8, 7}}
		};
		index = ind[a+1][b+1][c+1];
	} else {
		static const ShortType ind[5][5][5]={{{-1,-1,-1,-1,-1},{-1,26,-1,27,-1},{-1,-1,-1,-1,-1},{-1,28,-1,29,-1},{-1,-1,-1,-1,-1}},{{-1,34,-1,36,-1},{42,18,6,19,46},{-1,14,0,16,-1},{43,20,7,21,47},{-1,38,-1,40,-1}},{{-1,-1,-1,-1,-1},{-1,10,2,11,-1},{-1,4,-1,5,-1},{-1,12,3,13,-1},{-1,-1,-1,-1,-1}},{{-1,35,-1,37,-1},{44,22,8,23,48},{-1,15,1,17,-1},{45,24,9,25,49},{-1,39,-1,41,-1}},{{-1,-1,-1,-1,-1},{-1,30,-1,31,-1},{-1,-1,-1,-1,-1},{-1,32,-1,33,-1},{-1,-1,-1,-1,-1}}};
		index = ind[a+2][b+2][c+2];
	}
	
	
	assert(index!=-1);
#ifndef NDEBUG
	std::vector<OffsetType> neighs;
	Neighbors(neighs);
	assert(neighs[index]==offset);
#endif
	return index;
}


template<typename TO, typename TS, typename TD> void
Lagrangian3Stencil<TO,TS,TD>::NeighborsAround(
	ShortType index, std::vector<OffsetType> & neigh) const {
	assert(index>=0);
	assert(IsValid());
	
	
	
#ifndef NDEBUG
	std::vector<OffsetType> neighs;
	Neighbors(neighs);
	const OffsetType center = neighs[index];
#endif

	
	if(geom==Lagrangian3StencilGeometry::Voronoi){
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
	
	// ----- Cube class ----
	if(index<6){
		const int i = index/2;
		const int j = (i+1)%3, k=(i+2)%3;
		const OffsetType & v = offsets[j], w = offsets[k];
		
		if(geom==Lagrangian3StencilGeometry::Diamond){
			neigh.insert(neigh.end(),{v,w,-v,-w}); return;}

		const int s=2*(index%2)-1;
		const OffsetType u = s*offsets[i];
		//  Around u
		assert(center==u);

		if(geom==Lagrangian3StencilGeometry::CutCube){
			neigh.insert(neigh.end(),{u+v,u+w,u-v,u-w}); return;}
		
		if(geom==Lagrangian3StencilGeometry::Cube){
			neigh.insert(neigh.end(),{
				u+v,u+v+w,u+w,u+w-v,
				u-v,u-v-w,u-w,u-w+v}); return;}
		
		assert(geom==Lagrangian3StencilGeometry::SpikyCube);{
			neigh.insert(neigh.end(),{
				u+v,2*u+v+w,u+w,2*u+w-v,
				u-v,2*u-v-w,u-w,2*u-w+v}); return;}
	}
	
	index-=6;
	if(index<12){
		const int sj = 2*(index%2)-1; index/=2;
		const int si = 2*(index%2)-1; index/=2;
		const int i = index;
		const int j = (i+1)%3, k=(i+2)%3;
		
		const OffsetType u = si*offsets[i], v=sj*offsets[j], w =offsets[k];
		// Around u+v
		assert(center==u+v);

		if(geom==Lagrangian3StencilGeometry::CutCube){
			neigh.insert(neigh.end(),{u,u+w,v+w,v,v-w,u-w}); return;}
		
		if(geom==Lagrangian3StencilGeometry::Cube){
			neigh.insert(neigh.end(), {u,u+v+w,v,u+v-w}); return;}
		
		assert(geom==Lagrangian3StencilGeometry::SpikyCube);{
			neigh.insert(neigh.end(), {
				u,2*u+v+w,u+v+w,u+2*v+w,
				v,u+2*v-w,u+v-w,2*u+v-w}); return;}
	}
	
	index-=12;
	if(index<8){
		const int s2 = 2*(index%2)-1; index/=2;
		const int s1 = 2*(index%2)-1; index/=2;
		const int s0 = 2*(index%2)-1; index/=2;
		
		const OffsetType u = s0*offsets[0], v=s1*offsets[1], w = s2*offsets[2];
		// Around u+v+w
		assert(center==u+v+w);

		
		if(geom==Lagrangian3StencilGeometry::Cube){
			neigh.insert(neigh.end(), {u,u+v,v,v+w,w,w+u}); return;}
		
		assert(geom==Lagrangian3StencilGeometry::SpikyCube);{
			neigh.insert(neigh.end(), {2*u+v+w,u+v,u+2*v+w,v+w,u+v+2*w,w+u}); return;}
	}
	
	index-=8;
	assert(index<24);{
		const int sk = 2*(index%2)-1; index/=2;
		const int sj = 2*(index%2)-1; index/=2;
		const int si = 2*(index%2)-1; index/=2;
		const int i = index;
		const int j = (i+1)%3, k = (i+2)%3;
		
		const OffsetType u = si*offsets[i], v = sj*offsets[j], w = sk*offsets[k];
		// Around 2*u + v + w
		assert(center==2*u+v+w);

		assert(geom==Lagrangian3StencilGeometry::SpikyCube);{
			neigh.insert(neigh.end(), {u,u+v,u+v+w,u+w}); return;}
	}
/*
			
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
	 */
}
#endif /* Lagrangian3Stencil_h */
