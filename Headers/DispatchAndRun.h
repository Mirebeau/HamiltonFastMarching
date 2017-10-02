// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef DispatchAndRun_h
#define DispatchAndRun_h

/*#include "Base/HFMUtilities.h"
#include "Experimental/HFMTimeVar.h"
#include "Experimental/HFMState.h"*/
#include "Base/HFMInterface.h"

#include "Specializations/Isotropic.h"
#include "Specializations/Riemannian.h"
#include "Specializations/Curvature2.h"
#include "Specializations/Curvature3.h"

#include "Experimental/HalfDisk.h"
#include "Experimental/PrescribedCurvature2.h"
#include "Experimental/RiemannLifted.h"
#include "Experimental/ReedsSheppAdaptive2.h"
#include "Experimental/Quaternionic.h"


#define HFMSpecializationMacro(modelName) \
{typedef HFMInterface<Traits ## modelName> HFMI; typedef Stencil ## modelName StencilDataType;\
if(model== #modelName){ \
    io.currentSetter=IO::SetterTag::Compute;\
    HFMI(io, std::unique_ptr<StencilDataType>(new StencilDataType)).Run();\
    io.currentSetter=IO::SetterTag::User; return;} }

#ifdef AllBaseModels
#define Riemann
#define Curvature2
#define Curvature3
#define Isotropic
#endif

//#define HFMSpecializationMacro(hfmName,modelName) \
//{ if(model== #modelName) {HFMInterface<Stencil ## modelName,hfmName<Traits ## modelName> >(io).Run(); return;} }


void Run(IO & io){
    const std::string model = io.GetString("model");
    
#ifdef Custom
// This custom executable is here to let the user choose the adequate combination of (FastMarchingClass, Model) for his/her application.
#endif
    
// ------------- Model options ----------
    
#ifdef Riemann
    HFMSpecializationMacro(Riemann2);
    HFMSpecializationMacro(Riemann3);
#endif
    
#ifdef Curvature2
    HFMSpecializationMacro(ReedsShepp2);
    HFMSpecializationMacro(ReedsSheppForward2);
    HFMSpecializationMacro(Dubins2);
    HFMSpecializationMacro(Elastica2<5>);
    HFMSpecializationMacro(Elastica2<9>);
#endif
    
#ifdef Curvature3
    HFMSpecializationMacro(ReedsShepp3)
    HFMSpecializationMacro(ReedsSheppForward3)
#endif
    
#ifdef Isotropic
    HFMSpecializationMacro(IsotropicBox2<Boundary::Closed>);
    HFMSpecializationMacro(IsotropicBox3<Boundary::Closed>);
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
    
    // Differentiation with riemannian metrics
    HFMSpecializationMacro(RiemannDiff2)
    
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
    
    // Checking boundary conditions
    HFMSpecializationMacro(IsotropicBox2<Boundary::Periodic>);
    HFMSpecializationMacro(IsotropicBox2<Boundary::Sphere2_0>);
    HFMSpecializationMacro(Sphere2);
#endif
    
#ifdef Experimental1
    
#endif


    ExceptionMacro("Unrecognized model : " << model);
}


#endif /* DispatchAndRun_h */
