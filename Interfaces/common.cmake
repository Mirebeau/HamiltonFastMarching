# ---- Choose the models to be compiled ---- 
SET(TestCompilationModelNames "Isotropic2;Riemann2" CACHE STRING "TestCompilationModelNames")

SET(StandardModelNames "Isotropic2;Isotropic3;Diagonal2;Diagonal3;Riemann2;Riemann3;ReedsShepp2;ReedsSheppForward2;Elastica2;Dubins2;ReedsShepp3;ReedsSheppForward3;IsotropicDiff2;DubinsExt2;TTI2;TTI3" CACHE STRING "ModelNames")

# Also, Seismic3SIMD, which requires additional libraries

Set(ExperimentalModelNames "Seismic2;SeismicTopographic2;Seismic3;SeismicTopographic3;AlignedBillard;Riemann4;Riemann5;Elastica2_9;ElasticaExt2_5;ReedsSheppExt2;ReedsSheppForwardExt2;RiemannDiff2;RiemannLifted2_Periodic;Rander2;AsymmetricQuadratic2;AsymmetricQuadratic3;AsymmetricQuadratic3p1;Riemann3_Periodic;AsymRander2;ConvexReedsSheppForward2;ConvexDubins2;ConvexElastica2" CACHE STRING "ExperimentalModelNames")
Set(CustomModelNames "" CACHE STRING "CustomModelNames")

option(IncludeStandardModels "IncludeStandardModels" TRUE)
option(IncludeExperimentalModels "IncludeExperimentalModels" TRUE)

option(CustomExecutable "CustomExecutable")
option(TestCompilation "TestCompilation" FALSE)

Set(ModelNames ${CustomModelNames})

# Setting the model names
if(TestCompilation)
Set(ModelNames ${ModelNames} ${TestCompilationModelNames})
else()

	if(IncludeStandardModels)
	Set(ModelNames ${ModelNames} ${StandardModelNames})
	endif()

	if(IncludeExperimentalModels)
	Set(ModelNames ${ModelNames} ${ExperimentalModelNames})
	endif()

endif()


list(REMOVE_DUPLICATES ModelNames)
message(STATUS "Compiled models : ${ModelNames}")

# ---- More options -----

option(OPTIMIZE_FOR_NATIVE "Build with -march=native" OFF)
if(OPTIMIZE_FOR_NATIVE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()