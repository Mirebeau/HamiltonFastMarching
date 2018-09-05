# ---- Choose the models to be compiled ---- 
SET(StandardModelNames "Isotropic2;Isotropic3;Diagonal2;Diagonal3;Riemann2;Riemann3;ReedsShepp2;ReedsSheppForward2;Elastica2;Dubins2;ReedsShepp3;ReedsSheppForward3;IsotropicDiff2;DubinsExt2" CACHE STRING "ModelNames")
Set(ExperimentalModelNames "RiemannLifted2_Periodic" CACHE STRING "ExperimentalModelNames")
Set(CustomModelNames "" CACHE STRING "CustomModelNames")

option(IncludeStandardModels "IncludeStandardModels" TRUE)
option(IncludeExperimentalModels "IncludeExperimentalNames" FALSE)

option(CustomExecutable "CustomExecutable")

Set(ModelNames ${CustomModelNames})

# Setting the model names
if(IncludeStandardModels)
Set(ModelNames ${ModelNames} ${StandardModelNames})
endif()

if(IncludeExperimentalNames)
Set(ModelNames ${ModelNames} ${ExperimentalModelNames})
endif()

list(REMOVE_DUPLICATES ModelNames)
message(STATUS "Compiled models : ${ModelNames}")

# # ---- set a few machine dependent paths ----
# set(JMM_CPPLibs_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../../JMM_CPPLibs-master" CACHE FILEPATH "JMM_CPPLibs directory")

# #paths required for compilation of Lifted Fast Marching
# set(LinearAlgebra_dir 		"${JMM_CPPLibs_dir}/LinearAlgebra")
# set(Output_dir 				"${JMM_CPPLibs_dir}/Output")
# set(DataStructures_dir 		"${JMM_CPPLibs_dir}/DataStructures")


# set(DummyBinDir "${CMAKE_CURRENT_BINARY_DIR}/Dummy")


# # ---- Include headers -----
# add_subdirectory(${LinearAlgebra_dir} "${DummyBinDir}/LinearAlgebra")
# add_subdirectory(${Output_dir} "${DummyBinDir}/Output")
# add_subdirectory(${DataStructures_dir} 	"${DummyBinDir}/DataStructures")
# add_subdirectory("../../Headers" "${DummyBinDir}/HFMHeaders")

# set(DataStructures_Headers
# 	${DataStructures_Headers}
# 	${DataStructures_dir}/DataStructures/CappedVector.h
# 	${DataStructures_dir}/DataStructures/ShallowMap.h
# 	${DataStructures_dir}/DataStructures/RedeclareTypesMacro.h
# )

# set(LinearAlgebra_Headers 
# 	${LinearAlgebra_Headers} 
# 	${LinearAlgebra_dir}/LinearAlgebra/ArrayType.h
# 	${LinearAlgebra_dir}/LinearAlgebra/BasisReduction.h	
# )

# set(Output_Headers 
# 	${Output_Headers} 
# 	"${Output_dir}/Output/BaseIO.h"
# 	"${Output_dir}/Output/BaseIO.hxx"
# 	"${Output_dir}/Output/FileIO.h"
# 	"${Output_dir}/Output/FileIO.hxx"
# 	"${Output_dir}/Output/IO.h"
# 	"${Output_dir}/Output/IO.hxx"
# 	"${Output_dir}/Output/ExceptionMacro.h"
# 	"${Output_dir}/Output/MexIO.h"
# 	"${Output_dir}/Output/MexIO.hxx"
# 	"${Output_dir}/Output/PythonIO.h"
# 	"${Output_dir}/Output/PythonIO.hxx"
# )

# set(Project_Headers 
# 	${Base_Headers}
# 	${Base_Implementation_Headers}
# 	${Specializations_Headers}
# 	${Experimental_Headers}
# 	${ExtraAlgorithms_Headers}
# 	${Root_Headers}

# 	${DataStructures_Headers}
# 	${LinearAlgebra_Headers}
# 	${Output_Headers}
# )


# # ---- Headers IDE layout -----
# source_group("DataStructures" FILES ${DataStructures_Headers})
# source_group("LinearAlgebra" FILES ${LinearAlgebra_Headers})
# source_group("Output" FILES ${Output_Headers})

# source_group("Base" FILES ${Base_Headers})
# source_group("Base\\Implementation" FILES ${Base_Implementation_Headers})
# source_group("Specializations" FILES ${Specializations_Headers})
# source_group("Experimental" FILES ${Experimental_Headers})
# source_group("ExtraAlgorithms" FILES ${ExtraAlgorithms_Headers})