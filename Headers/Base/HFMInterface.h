// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef HFMInterface_h
#define HFMInterface_h

#include <memory>
#include <time.h>
#include "GeodesicODESolver.h"
#include "GeodesicDiscreteSolver.h"
#include "ExtraAlgorithms/TimeDependentFields.h"
#include "ExtraAlgorithms/StaticFactoring.h"

template<typename TTraits> struct HFMInterface {
    typedef TTraits Traits;
    typedef HamiltonFastMarching<Traits> HFM;
    Redeclare7Types(HFM,ActiveNeighFlagType,StencilDataType,DiscreteFlowType,
					ExtraAlgorithmInterface,GeodesicSolverInterface,IndexCRef,OffsetCRef);
    Redeclare6Types(Traits,DiscreteType,ScalarType,PointType,VectorType,IndexType,
					IndexDiff);
	Redeclare1Type(IO,KeyCRef)
    Redeclare2Constants(Traits,Dimension,mathPi)

    template<typename E, size_t n> using Array = typename Traits::template Array<E,n>;
    template<typename E> using DataSource = typename Traits::template DataSource<E>;
    template<typename E> struct DataSource_Value;
    
    IO & io;
    StencilDataType & stencil;
    std::unique_ptr<HFM> pFM;
    std::vector<std::unique_ptr<ExtraAlgorithmInterface> > extras;
    std::unique_ptr<TimeDependentFields<Traits> > pTime;
    std::unique_ptr<GeodesicSolverInterface> pGeodesicSolver;
	
    HFMInterface(IO & _io, StencilDataType & _stencil) :io(_io), stencil(_stencil) {};
    virtual void Run();

    template<typename E> std::unique_ptr<DataSource<E> > GetField(KeyCRef, bool=true); // bool field: allow time dependency
	int FieldElementSize(KeyCRef) const;
    template<typename E> std::unique_ptr<DataSource<E> > GetIntegralField(KeyCRef);
    void ExportGeodesics(KeyCRef, const std::vector<PointType> &);
    template<typename E> struct DataSource_Inverse; // A data source which inverses the values of another data source.
    template<bool b, typename Dummy> struct SpecializationsDefault_; // Mostly internal methods
	typedef SpecializationsDefault_<HFM::hasBundle,void> SpecializationsDefault;
protected:
    template<typename E> struct DataSource_Array;
    template<typename E, bool=HFM::hasBundle> struct DataSource_Indep; // A data source independent of the "bundle" coords.
    template<typename E, bool=HFM::hasBundle> struct DataSource_Dep;   // A data source depending only on the "bundle" coords.
    template<typename E> struct TimeDependentSource;

//    std::unique_ptr<StencilDataType> _pStencil;
    template<typename Alg> Alg * SetupSingleAlgorithm();
    
    virtual void Run_SetupIO();
    virtual void Run_SetupStencil(); // Also creates HFM
    virtual void Run_SetupSolver();
    virtual void Run_SetupExtraAlgorithms();
    virtual bool Run_RunSolver();
    virtual void Run_ExtractGeodesics();
    virtual void Run_ExportData();	
};

// Please add more extra algorithms if needed
#include "ExtraAlgorithms/CommonStoppingCriteria.h"
#include "ExtraAlgorithms/WallObstructionTest.h"
#include "ExtraAlgorithms/EuclideanPathLength.h"
#include "ExtraAlgorithms/VoronoiDiagram.h"
#include "ExtraAlgorithms/FirstVariation.h"

template<typename T> void HFMInterface<T>::
Run_SetupExtraAlgorithms(){
    SetupSingleAlgorithm<CommonStoppingCriteria<T> >();
    SetupSingleAlgorithm<WallObstructionTest<T> >();
    SetupSingleAlgorithm<EuclideanPathLength<T> >();
    SetupSingleAlgorithm<VoronoiDiagram<T> >();
    auto pVar = SetupSingleAlgorithm<FirstVariation<T> >();
    if(!pTime->ImplementIn(pFM.get())) pTime.reset();
    if(pVar && pTime) pVar->pCurrentTime = &(pTime->currentTime);	
}

#include "Implementation/HFMInterface.hxx"

#endif /* HFMInterface_h */
