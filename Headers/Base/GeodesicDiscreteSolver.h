// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef GeodesicDiscreteSolver_h
#define GeodesicDiscreteSolver_h

#include "HamiltonFastMarching.h"

template<typename TTraits> struct GeodesicDiscreteSolver
: HamiltonFastMarching<TTraits>::GeodesicSolverInterface {
    typedef TTraits Traits; typedef HamiltonFastMarching<Traits> HFM;
    Redeclare8Types(Traits,DiscreteType,ScalarType,ShortType,IndexType,OffsetType,VectorType,PointType,IndexDiff)
    Redeclare6Types(HFM,DiscreteFlowType,MultiplierType,IndexCRef,DifferenceType,HFMI,GeodesicSolverInterface);
    template<typename E, size_t n> using Array = typename Traits::template Array<E,n>;
    Redeclare1Constant(Traits,Dimension);
    
    ScalarType geodesicStep = 0.25;
    ScalarType weightThreshold = 0.001; // smaller contributions eliminated
    ScalarType volumeBound= 5.*pow(1.3, Dimension);  // restart if too much spread
    DiscreteType nRestartsBeforeVolumeIncrease = 5;
    ScalarType volumeIncreaseRatio = 1.3;
    DiscreteType nMaxRestarts = 50;

    mutable DiscreteType nRestarts; mutable ScalarType effectiveVolumeBound; // Setup by Run
    virtual bool Run(std::vector<PointType> &) const override;
    virtual void Setup(HFMI *) override;
    virtual std::vector<std::vector<PointType> > Run(HFMI *, const std::vector<PointType> &) override;
    using GeodesicSolverInterface::GeodesicSolverInterface;
protected:
    bool GeodesicDiscrete(std::vector<PointType> &) const;
    // TODO separate virtual method for pruning.
};


#include "Implementation/GeodesicDiscreteSolver.hxx"

#endif /* GeodesicDiscreteSolver_h */
