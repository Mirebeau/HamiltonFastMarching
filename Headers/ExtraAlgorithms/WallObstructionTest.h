// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef WallObstructionTest_h
#define WallObstructionTest_h

#include "Base/HFMInterface.h"

template<typename T> struct WallObstructionTest :
HamiltonFastMarching<T>::ExtraAlgorithmInterface {
    typedef HamiltonFastMarching<T> HFM;
    typedef typename HFM::ExtraAlgorithmInterface Superclass;
    Redeclare11Types(HFM,IndexCRef,IndexType,ScalarType,Traits,HFMI,
					 PointType,ShortType,DiscreteType,OffsetCRef,OffsetType,DomainTransformType)
    Redeclare1Constant(HFM,Dimension)
    
    virtual void Setup(HFMI*) override;
    virtual bool ImplementIn(HFM*) override;
    typedef typename HFM::template DataSource<bool> BoolField;
    std::unique_ptr<BoolField> pWalls;
protected:
    virtual bool Visible(IndexCRef,OffsetCRef,IndexCRef) const override;
    const HFM * pFM=nullptr;
    typedef typename HFM::template Array<ShortType,Dimension> ShortArray;
    ShortArray wallDist;
    bool refineStencilAtWallBoundary = false;
};

template<typename T> void
WallObstructionTest<T>::Setup(HFMI*that){
    if(that->io.HasField("walls")){
        pWalls = that->template GetIntegralField<bool>("walls");}
    refineStencilAtWallBoundary = (bool)
    that->io.template Get<ScalarType>("refineStencilAtWallBoundary",false,3);
}


template<typename T> bool
WallObstructionTest<T>::ImplementIn(HFM*_pFM) {
    if(pWalls==nullptr) return false;
    assert(_pFM!=nullptr);
    pFM=_pFM;
    _pFM->extras.visible.push_back(this);
    const IndexType dims = pFM->values.dims;
    if(!pWalls->CheckDims(dims)) ExceptionMacro("Error : walls have inconsistent dims");
    
    // ------------ Running a Dikjstra to compute the distance to walls ---------

    wallDist.dims = dims;
    wallDist.resize(dims.Product(),std::numeric_limits<ShortType>::max());
    const ShortType maxDist = pFM->MaxStencilWidth();

    struct QueueElementType {
        ShortType val;
        DiscreteType linearIndex;
        bool operator < (const QueueElementType & other) const {
            return val>other.val || (val==other.val && linearIndex> other.linearIndex);}
    };
    
    std::priority_queue<QueueElementType> queue;
    
    for(DiscreteType linearIndex=0; linearIndex<wallDist.size(); ++linearIndex)
        if((*pWalls)(wallDist.Convert(linearIndex))){
            queue.push({0,linearIndex});
            wallDist[linearIndex]=0;
        }
    
    while(!queue.empty()){
        const auto top = queue.top();
        queue.pop();
        const IndexType index = wallDist.Convert(top.linearIndex);
        if(wallDist[top.linearIndex]!=top.val) continue;
        const ShortType neighBound = top.val+1;
        if(neighBound>maxDist) continue;
        for(int i=0; i<Dimension; ++i)
            for(int eps=-1; eps<=1; eps+=2){
                IndexType neighbor = index;
                neighbor[i]+=eps;
                const auto transform = pFM->dom.Periodize(neighbor,index);
                if(!transform.IsValid()) continue;
                const DiscreteType neighLinearIndex = wallDist.Convert(neighbor);
                ShortType & neighVal = wallDist[neighLinearIndex];
                if(neighVal<=neighBound) continue;
                neighVal = neighBound;
                queue.push({neighBound, neighLinearIndex});
            }
    }
    
    return true;
}


// ----------- Walking on a line to see if walls are met ---------
template<typename T> bool
WallObstructionTest<T>::Visible(IndexCRef base, OffsetCRef v, IndexCRef sum) const  {
    assert(!wallDist.empty()); assert(pFM!=nullptr);
    
    // Are we far from the wall ?
    DiscreteType wd = (DiscreteType)wallDist(sum);
    if(wd==std::numeric_limits<ShortType>::max()) return true;
    else if(wd==0) return false;
    
    // If stencil is refined at wall boundary, then these cases pass though
    if(refineStencilAtWallBoundary && wallDist(base)==1) return true;
    
    // Begin walking on the line from source to target
    OffsetType vSign = OffsetType::Constant(0);
    DiscreteType vSum=0, vProd=1, vMax=0;
    for(int i=0; i<Dimension; ++i){
        const DiscreteType vi=v[i];
        if(vi!=0){
            vSign[i] = vi>0 ? 1 : -1;
            const DiscreteType avi = std::abs(vi);
            vSum += avi;
            vProd *= avi;
            vMax = std::max(vMax,avi);
        }
    }
    
    // Are we far enough ?
    if(vSum<=wd) return true;
    // We next construct a straight discrete path from source to seed, and see if it intersects a wall.
    
    IndexType pos=base, vNext, vInc;
    for(int i=0; i<Dimension; ++i){
        vNext[i] = (v[i]==0 ?
                    std::numeric_limits<DiscreteType>::max() :
                    vProd/std::abs(v[i]));
        vInc[i]=2*vNext[i];
    }
    
    for(int i=1; i<vMax; ++i){
        const DiscreteType threshold = (2*i*vProd)/vMax;
        for(int j=0; j<Dimension; ++j){
            if(vNext[j]<=threshold){
                vNext[j]+=vInc[j];
                --vSum;
                pos[j]+=vSign[j];
            }
        }
        IndexType posPer=pos;
        const DomainTransformType transform = pFM->dom.Periodize(posPer,base);
        assert(transform.IsValid()); (void)transform;
        wd=(DiscreteType)wallDist(posPer);
        if(vSum<=wd) return true;
        else if(wd==0) return false;
    }
    return true;
};

#endif /* WallObstructionTest_h */
