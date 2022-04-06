// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef DispatchAndRun_h
#define DispatchAndRun_h

#include "JMM_CPPLibs/Macros/String.h"
#include "JMM_CPPLibs/Macros/PPCat.h"

// Saving the model name as a string for future reference
const std::string ModelNameString=STRING(ModelName);

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
#define Riemannian_Riemann4        		1
#define Riemannian_Riemann5        		1
#define Riemannian_Riemann3_Periodic	1

#define Curvature2_ReedsShepp2          1
#define Curvature2_ReedsSheppForward2   1
#define Curvature2_Dubins2              1
#define Curvature2_Elastica2            1
#define Curvature2_Elastica2_9          1

#define ConvexCurvature2_ConvexReedsSheppForward2  1
#define ConvexCurvature2_ConvexDubins2             1
#define ConvexCurvature2_ConvexElastica2           1

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

//#define AsymmetricQuadratic_AsymmetricQuadratic2    1 // Now use semi-lag scheme QuadLinLag2
#define AsymmetricQuadratic_AsymmetricQuadratic3    1
#define AsymmetricQuadratic_AsymmetricQuadratic4    1
#define AsymmetricQuadratic_AsymmetricQuadratic3p1  1

#define RiemannLifted_RiemannLifted2_Closed         1
#define RiemannLifted_RiemannLifted2_Periodic       1
#define RiemannLifted_RiemannLifted3                1

#define QuadLinLag2_Rander2 1
#define QuadLinLag2_AsymmetricQuadratic2 1

#define Seismic2_Seismic2 1
#define Seismic2_SeismicTopographic2 1
#define Seismic3_Seismic3 1
#define Seismic3_SeismicTopographic3 1
#define Seismic3SIMD_Seismic3SIMD 1

#define TTI_TTI2 1
#define TTI_TTI3 1

#define AlignedBillard_AlignedBillard 1

#define AsymRander2_AsymRander2 1

// **** Include the correct header ****
// Standard specializations
#if PPCAT(Isotropic_,ModelName)
#include "Specializations/Isotropic.h"
using StencilIsotropic1 = StencilIsotropic<1>;
using StencilIsotropic2 = StencilIsotropic<2>;
using StencilIsotropic3 = StencilIsotropic<3>;
using StencilIsotropic4 = StencilIsotropic<4>;
using StencilIsotropic5 = StencilIsotropic<5>;
using StencilDiagonal2 = StencilDiagonal<2>;
using StencilDiagonal3 = StencilDiagonal<3>;
using StencilDiagonal4 = StencilDiagonal<4>;
using StencilDiagonal5 = StencilDiagonal<5>;
#elif PPCAT(Riemannian_,ModelName)
#include "Specializations/Riemannian.h"
using StencilRiemann2 = StencilRiemann<2>;
using StencilRiemann3 = StencilRiemann<3>;
using StencilRiemann4 = StencilRiemann<4>;
using StencilRiemann5 = StencilRiemann<5>;
using StencilRiemann3_Periodic = StencilRiemann<3,Boundary::Periodic>;
#elif PPCAT(Curvature2_,ModelName)
#include "Specializations/Curvature2.h"
using StencilElastica2_9 = StencilElastica2<9>;
using StencilElastica2_5 = StencilElastica2<5>;
#define Elastica2 Elastica2_5 // Default discretization
#elif PPCAT(Curvature3_,ModelName)
#include "Specializations/Curvature3.h"
#elif PPCAT(QuadLinLag2_,ModelName)
#include "Specializations/QuadLinLag2.h"

// Experimental specializations
#elif PPCAT(PrescribedCurvature2_,ModelName)
#include "Experimental/PrescribedCurvature2.h"
using StencilElasticaExt2_5 = StencilElasticaExt2<5>;
#elif PPCAT(ConvexCurvature2_,ModelName)
#define convex_curvature_macro 1
#include "Specializations/Curvature2.h"
using StencilConvexReedsSheppForward2 = StencilReedsSheppForward2;
using StencilConvexDubins2 = StencilDubins2;
using StencilConvexElastica2 = StencilElastica2<5>;
#elif PPCAT(Differentiable_,ModelName)
#include "Experimental/Differentiable.h"
using StencilIsotropicDiff2 = StencilIsotropicDiff<2>;
using StencilIsotropicDiff3 = StencilIsotropicDiff<3>;
using StencilIsotropicDiff4 = StencilIsotropicDiff<4>;
#elif PPCAT(AsymmetricQuadratic_,ModelName)
#include "Experimental/AsymmetricQuadratic.h"
using StencilAsymmetricQuadratic2 = StencilAsymmetricQuadratic<2>;
using StencilAsymmetricQuadratic3 = StencilAsymmetricQuadratic<3>;
using StencilAsymmetricQuadratic4 = StencilAsymmetricQuadratic<4>;
#elif PPCAT(RiemannLifted_,ModelName)
#include "Experimental/RiemannLifted.h"
using StencilRiemannLifted2_Closed = StencilRiemannLifted2<Boundary::Closed>;
using StencilRiemannLifted2_Periodic = StencilRiemannLifted2<Boundary::Periodic>;
#elif PPCAT(RollingBall_,ModelName)
#include "Experimental/RollingBall.h"
#elif PPCAT(Seismic2_,ModelName)
#include "Experimental/Seismic2.h"
#elif PPCAT(Seismic3_,ModelName)
#include "Experimental/Seismic3.h"
#elif PPCAT(Seismic3SIMD_,ModelName)
#include "xsimd/xsimd.hpp"
#include "Experimental/Seismic3.h"
using StencilSeismic3SIMD = StencilSeismic3;
#elif PPCAT(TTI_,ModelName)
#include "Experimental/TTI.h"
using StencilTTI2 = StencilTTI<2>;
using StencilTTI3 = StencilTTI<3>;
#elif PPCAT(AsymRander2_,ModelName)
#include "Experimental/AsymRander.h"
#elif PPCAT(AlignedBillard_,AlignedBillard)
#include "Experimental/AlignedBillard.h"
// Very experimental specializations
#else
/*
#include "Experimental/ReedsSheppAdaptive2.h"
#include "Experimental/Quaternionic.h"
#include "Experimental/RollingBall.h"
 */
#endif

// ------- Custom invocation, with multiple models.  ---------

#define HFMSpecializationMacro(modelName) \
{ \
using StencilDataType = Stencil ## modelName ;\
using HFMI = StencilDataType::HFMI; \
if(model== #modelName){ \
    io.currentSetter=IO::SetterTag::Compute;\
    StencilDataType stencil; \
    HFMI(io, stencil).Run();\
    io.currentSetter=IO::SetterTag::User; return;} \
}

#ifdef Custom
//#include "Experimental/Differentiable.h"
//#include "Experimental/RiemannLifted.h"

#include "Experimental/AlignedBillard.h"
#include "Experimental/TTI.h"
#include "Experimental/AsymRander.h"

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
    Redeclare4Types(HFM,ActiveNeighFlagType,StencilDataType,ExtraAlgorithmInterface,GeodesicSolverInterface);
    
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
//    typedef HFMInterface<PPCAT(Traits,ModelName)> HFMI;
    typedef PPCAT(Stencil,ModelName) StencilDataType;
	using HFMI = StencilDataType::HFMI;
    
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

	HFMSpecializationMacro(AsymRander2)
//	HFMSpecializationMacro(TTI2)
//	HFMSpecializationMacro(TTI3)
//    HFMSpecializationMacro(IsotropicDiff<2>)
//    HFMSpecializationMacro(RiemannLifted2<Boundary::Closed>)
	
	/*
	 HFMSpecializationMacro(AlignedBillard)
	 */
/*
    {
        typedef HFMInterface<TraitsLagrangian2> HFMI;
        typedef StencilRanderLag2 StencilDataType;
        if(model== "RanderLag2"){
            io.currentSetter=IO::SetterTag::Compute;
            StencilDataType stencil;
            HFMI(io, stencil).Run();
            io.currentSetter=IO::SetterTag::User; return;}
    }
    
    {
        typedef HFMInterface<TraitsLagrangian2> HFMI;
        typedef StencilAsymmetricQuadraticLag2 StencilDataType;
        if(model== "AsymmetricQuadraticLag2"){
            io.currentSetter=IO::SetterTag::Compute;
            StencilDataType stencil;
            HFMI(io, stencil).Run();
            io.currentSetter=IO::SetterTag::User; return;}
    }*/

    
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
