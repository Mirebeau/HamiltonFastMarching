// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef DispatchAndRun_h
#define DispatchAndRun_h

#define PPCAT_NX(A, B) A ## B
#define PPCAT(A, B) PPCAT_NX(A, B)
#define STRING_NX(s) #s
#define STRING(s) STRING_NX(s)

// **** Do we need high dimensional Voronoi reduction ***

#define HighVoronoi_Riemann4                1
#define HighVoronoi_Riemann5                1
#define HighVoronoi_AsymmetricQuadratic4    1

#if PPCAT(HighVoronoi_,ModelName)
#define HighVoronoi
#endif

// **** Which header contains which model ****
// Standard specializations

#define Isotropic_Isotropic1            1
#define Isotropic_Isotropic2            1
#define Isotropic_Isotropic3            1
#define Isotropic_Isotropic4            1
#define Isotropic_Isotropic5            1
#define Isotropic_Diagonal2             1
#define Isotropic_Diagonal3             1
#define Isotropic_Diagonal4             1
#define Isotropic_Diagonal5             1

#define Riemannian_Riemann2             1
#define Riemannian_Riemann3             1
#define Riemannian_Riemann4        1
#define Riemannian_Riemann5        1

#define Curvature2_ReedsShepp2          1
#define Curvature2_ReedsSheppForward2   1
#define Curvature2_Dubins2              1
#define Curvature2_Elastica2            1
#define Curvature2_Elastica2_9          1

#define Curvature3_ReedsShepp3          1
#define Curvature3_ReedsSheppForward3   1

// Experimental specializations

#define PrescribedCurvature2_ReedsSheppExt2         1
#define PrescribedCurvature2_ReedsSheppForwardExt2  1
#define PrescribedCurvature2_DubinsExt2             1
#define PrescribedCurvature2_ElasticaExt2_5         1

#define Differentiable_IsotropicDiff2   1
#define Differentiable_IsotropicDiff3   1
#define Differentiable_IsotropicDiff4   1
#define Differentiable_RiemannDiff2     1

#define AsymmetricQuadratic_AsymmetricQuadratic2    1
#define AsymmetricQuadratic_AsymmetricQuadratic3    1
#define AsymmetricQuadratic_AsymmetricQuadratic4    1
#define AsymmetricQuadratic_AsymmetricQuadratic3p1  1

#define RiemannLifted_RiemannLifted2_Closed         1
#define RiemannLifted_RiemannLifted2_Periodic       1
#define RiemannLifted_RiemannLifted3                1

#define HalfDisk_HalfDisk2              1
#define HalfDisk_HalfDisk3              1
#define HalfDisk_HalfDisk3p1            1



// **** Include the correct header ****
// Standard specializations
#if PPCAT(Isotropic_,ModelName)
#include "Specializations/Isotropic.h"
#elif PPCAT(Riemannian_,ModelName)
#include "Specializations/Riemannian.h"
#elif PPCAT(Curvature2_,ModelName)
#include "Specializations/Curvature2.h"
#elif PPCAT(Curvature3_,ModelName)
#include "Specializations/Curvature3.h"

// Experimental specializations
#elif PPCAT(PrescribedCurvature2_,ModelName)
#include "Experimental/PrescribedCurvature2.h"
#elif PPCAT(Differentiable_,ModelName)
#include "Experimental/Differentiable.h"
#elif PPCAT(AsymmetricQuadratic_,ModelName)
#include "Experimental/AsymmetricQuadratic.h"
#elif PPCAT(RiemannLifted_,ModelName)
#include "Experimental/RiemannLifted.h"
#elif PPCAT(RollingBall_,ModelName)
#include "Experimental/RollingBall.h"

// Very experimental specializations
#else
#include "Experimental/ReedsSheppAdaptive2.h"
#include "Experimental/Quaternionic.h"
#include "Experimental/RollingBall.h"
#endif

// Saving the model name as a string for future reference
const std::string ModelNameString=STRING(ModelName);

// Key redefinitions for templated classes
#define Isotropic1 Isotropic<1>
#define Isotropic2 Isotropic<2>
#define Isotropic3 Isotropic<3>
#define Isotropic4 Isotropic<4>
#define Isotropic5 Isotropic<5>
#define Diagonal2 Diagonal<2>
#define Diagonal3 Diagonal<3>
#define Diagonal4 Diagonal<4>
#define Diagonal5 Diagonal<5>

#define Riemann2 Riemann<2>
#define Riemann3 Riemann<3>
#define Riemann4 Riemann<4>
#define Riemann5 Riemann<5>

#define IsotropicDiff2 IsotropicDiff<2>
#define IsotropicDiff3 IsotropicDiff<3>
#define IsotropicDiff4 IsotropicDiff<4>

#define AsymmetricQuadratic2 AsymmetricQuadratic<2>
#define AsymmetricQuadratic3 AsymmetricQuadratic<3>
#define AsymmetricQuadratic4 AsymmetricQuadratic<4>

#define Elastica2 Elastica2<5>
#define Elastica2_9 Elastica2<9>
#define RiemannLifted2_Closed RiemannLifted2<Boundary::Closed>
#define RiemannLifted2_Periodic RiemannLifted2<Boundary::Periodic>


// ------- Custom invocation, with multiple models.  ---------

#define HFMSpecializationMacro(modelName) \
{typedef HFMInterface<Traits ## modelName> HFMI; typedef Stencil ## modelName StencilDataType;\
if(model== #modelName){ \
    io.currentSetter=IO::SetterTag::Compute;\
    HFMI(io, std::unique_ptr<StencilDataType>(new StencilDataType)).Run();\
    io.currentSetter=IO::SetterTag::User; return;} }

#ifdef Custom
#include "Experimental/Differentiable.h"
#include "Experimental/RiemannLifted.h"
#endif
/*
#ifdef AllBaseModels
#define Riemann
#define Curvature2
#define Curvature3
#define Isotropic
#endif
*/
//#define HFMSpecializationMacro(hfmName,modelName) \
//{ if(model== #modelName) {HFMInterface<Stencil ## modelName,hfmName<Traits ## modelName> >(io).Run(); return;} }

template<typename TTraits> struct HFMInterface2 {
    typedef TTraits Traits;
    typedef HamiltonFastMarching<Traits> HFM;
    Redeclare4Types(FromHFM,ActiveNeighFlagType,StencilDataType,ExtraAlgorithmInterface,GeodesicSolverInterface);
    
    template<typename E, size_t n> using Array = typename Traits::template Array<E,n>;
    template<typename E> using DataSource = typename Traits::template DataSource<E>;
    template<typename E> struct DataSource_Value;
    
    IO & io;
    StencilDataType & stencil;
    HFMInterface2(IO & _io, StencilDataType & _stencil) : io(_io), stencil(_stencil) {}
    //_pStencil(std::move(pStencil)) {};
    virtual void Run(){};
};


void Run(IO & io){

    // ------- Run a single model --------
#ifdef ModelName
    
    if(io.HasField("model")){
        if(io.GetString("model")!=ModelNameString)
            ExceptionMacro("Executable applies to " << ModelNameString
                           << ", not to " << io.GetString("model") << ".");
    }
    typedef HFMInterface<PPCAT(Traits,ModelName)> HFMI;
    typedef PPCAT(Stencil,ModelName) StencilDataType;
    
    io.currentSetter=IO::SetterTag::Compute;
    StencilDataType stencil;
    HFMI(io, stencil).Run();
    io.currentSetter=IO::SetterTag::User;
    return;
#endif
    
    // ------- Run one of several models (mostly debug purposes) ------
    const std::string model = io.GetString("model");
#ifdef Custom
// This custom executable is here to let the user choose the adequate combination of (FastMarchingClass, Model) for his/her application.
    HFMSpecializationMacro(IsotropicDiff<2>)
    HFMSpecializationMacro(RiemannLifted2<Boundary::Closed>)
#endif
    
// ------------- Model options ----------
    /*
#ifdef Riemann
    HFMSpecializationMacro(Riemann2);
    HFMSpecializationMacro(Riemann3);
#endif
    
#ifdef Curvature2
    HFMSpecializationMacro(ReedsShepp2);
    HFMSpecializationMacro(ReedsSheppForward2);
    HFMSpecializationMacro(Dubins2);
    HFMSpecializationMacro(Elastica2<5>);
//    HFMSpecializationMacro(Elastica2<9>);
#endif
    
#ifdef Curvature3
    HFMSpecializationMacro(ReedsShepp3)
    HFMSpecializationMacro(ReedsSheppForward3)
#endif
    
#ifdef Isotropic
    HFMSpecializationMacro(Isotropic2)
    HFMSpecializationMacro(Diagonal2)
    HFMSpecializationMacro(Isotropic3)
    HFMSpecializationMacro(Diagonal3)
#endif
    


// ----------- Experimental -----------
 
#ifdef PrescribedCurvature2
    HFMSpecializationMacro(ReedsSheppExt2)
    HFMSpecializationMacro(ReedsSheppForwardExt2)
    HFMSpecializationMacro(DubinsExt2)
    HFMSpecializationMacro(ElasticaExt2<5>)
#endif
    
#ifdef RiemannExtra
    // HalfDisk models
    HFMSpecializationMacro(HalfDisk2)
    HFMSpecializationMacro(HalfDisk3)
    HFMSpecializationMacro(HalfDisk3p1)
    
    // Differentiation with riemannian metrics
    HFMSpecializationMacro(RiemannDiff2)
    
    // Continuous differentiation of isotropic metrics
    HFMSpecializationMacro(IsotropicDiff2)
    
    // Lifted riemannian metrics
    HFMSpecializationMacro(RiemannLifted2<Boundary::Closed>)
    HFMSpecializationMacro(RiemannLifted2<Boundary::Periodic>)
    HFMSpecializationMacro(RiemannLifted3)
#endif
    
#ifdef Experimental0
    // Quaternions
    HFMSpecializationMacro(ReedsSheppSO3)
    HFMSpecializationMacro(ReedsSheppForwardSO3)
    
    // Variants of Reeds-Shepp, with adaptive angular resolution, distinct foward/backward/angular speeds
    HFMSpecializationMacro(ReedsSheppAdaptive2)
    HFMSpecializationMacro(ReedsSheppThreeSpeeds2)
#endif


#ifdef Experimental1
    HFMSpecializationMacro(Riemann4)
    HFMSpecializationMacro(Riemann5)
#endif
*/

    ExceptionMacro("Unrecognized model : " << model);
}


#endif /* DispatchAndRun_h */
