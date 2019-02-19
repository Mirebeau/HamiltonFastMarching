// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef GeodesicODESolver_h
#define GeodesicODESolver_h

// ------------------------------------------
// -------- Approach based on ODEs ----------
// ------------------------------------------


#include "HamiltonFastMarching.h"
#include "JMM_CPPLibs/LinearAlgebra/SquareCube.h"


template <typename TTraits> struct GeodesicODESolver
: HamiltonFastMarching<TTraits>::GeodesicSolverInterface {
    typedef TTraits Traits; typedef HamiltonFastMarching<Traits> HFM;
    Redeclare7Types(Traits,DiscreteType,ScalarType,ShortType,IndexType,OffsetType,VectorType,PointType)
    Redeclare4Types(HFM,HFMI,GeodesicSolverInterface,DomainTransformType,FlowDataType);
    template<typename E, size_t n> using Array = typename Traits::template Array<E,n>;
    Redeclare1Constant(Traits,Dimension);
    
    ScalarType geodesicStep = 0.25;
    DiscreteType targetTolerance = 6;
    ScalarType causalityTolerance = 4;
    ScalarType stationarityThreshold = 6*sqrt(Dimension);
    
    virtual bool Run(std::vector<PointType> &) const override;
    virtual void Setup(HFMI *) override;
    virtual std::vector<std::vector<PointType> > Run(HFMI *, const std::vector<PointType> &) override;
    GeodesicODESolver(const HFM & hfm);
protected:
    typedef std::map<DiscreteType,std::pair<FlowDataType,int> > FlowCacheType;
    struct FlowAvgType;
    FlowAvgType GeodesicFlow(const PointType &, const Array<ShortType,Dimension> &, FlowCacheType &) const;
    Array<ShortType,Dimension> LInfDistance(const std::vector<IndexType> &, ShortType) const;
    Array<ShortType,Dimension> targetDistances;
};

#include "Implementation/GeodesicODESolver.hxx"

#endif /* GeodesicODESolver_h */
