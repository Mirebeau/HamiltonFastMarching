# ------- THIS DOES NOT BUILD MATLAB HFM ------

# This file is only used for compilation tests with xcode
# Please use the provided Matlab(R) script to compile the project

cmake_minimum_required(VERSION 3.1)
set (CMAKE_CXX_STANDARD 14)

project(NotMatlabHFM)


#set a few machine dependent paths
set(Matlab_dir "/Applications/MATLAB_R2012b.app/extern/include" CACHE FILEPATH "Matlab directory")
set(JMM_CPPLibs_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../../../JMM_CPPLibs-master" CACHE FILEPATH "JMM_CPPLibs directory")

#paths required for compilation of Lifted Fast Marching
set(LinearAlgebra_dir 		"${JMM_CPPLibs_dir}/LinearAlgebra")
set(Output_dir 				"${JMM_CPPLibs_dir}/Output")
set(DataStructures_dir 		"${JMM_CPPLibs_dir}/DataStructures")


set(DummyBinDir "${CMAKE_CURRENT_BINARY_DIR}/Dummy")

#LiftedFastMarching headers
add_subdirectory(${LinearAlgebra_dir} "${DummyBinDir}/LinearAlgebra")
add_subdirectory(${Output_dir} "${DummyBinDir}/Output")
add_subdirectory(${DataStructures_dir} 	"${DummyBinDir}/DataStructures")
add_subdirectory("../../../Headers" "${DummyBinDir}/HFMHeaders")
#add_subdirectory("Base")
#add_subdirectory("Specializations")
#add_subdirectory("Experimental")
#add_subdirectory("ExtraAlgorithms")

set(DataStructures_Headers
	${DataStructures_Headers}
	${DataStructures_dir}/DataStructures/CappedVector.h
	${DataStructures_dir}/DataStructures/ShallowMap.h
	${DataStructures_dir}/DataStructures/RedeclareTypesMacro.h
)

set(LinearAlgebra_Headers 
	${LinearAlgebra_Headers} 
	${LinearAlgebra_dir}/LinearAlgebra/ArrayType.h
	${LinearAlgebra_dir}/LinearAlgebra/BasisReduction.h	
)

set(Output_Headers 
	${Output_Headers} 
	"${Output_dir}/Output/BaseIO.h"
	"${Output_dir}/Output/BaseIO.hxx"
	"${Output_dir}/Output/MexIO.h"
	"${Output_dir}/Output/MexIO.hxx"
	"${Output_dir}/Output/IO.h"
	"${Output_dir}/Output/IO.hxx"
	"${Output_dir}/Output/ExceptionMacro.h"
)

set(Project_Headers 
	${Base_Headers}
	${Base_Implementation_Headers}
	${Specializations_Headers}
	${Experimental_Headers}
	${ExtraAlgorithms_Headers}
	${Root_Headers}

	${DataStructures_Headers}
	${LinearAlgebra_Headers}
	${Output_Headers}
)

include_directories(${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${Matlab_dir})

# List of executables, as specified in DispachAndRun

option(AllBaseExecutable "AllBaseExecutable" TRUE)
if(AllBaseExecutable)
add_executable(NotMatlabHFM_AllBase "../MatlabHFM.cpp" ${Project_Headers})
target_compile_definitions(NotMatlabHFM_AllBase PRIVATE AllBaseModels)
endif()

option(RiemannExtraExecutable "RiemannExtraExecutable")
if(RiemannExtraExecutable)
add_executable(NotMatlabHFM_RiemannExtra "MatlabHFM.cpp")
target_compile_definitions(NotMatlabHFM_RiemannExtra PRIVATE RiemannExtra)
endif()

# Executables that are part of NotMatlabHFM_AllBase
# (On some systems, NotMatlabHFM_AllBase does not compile because it is too big (?))


option(IsotropicExecutable "IsotropicExecutable")
if(IsotropicExecutable)
add_executable(NotMatlabHFM_Isotropic "MatlabHFM.cpp")
target_compile_definitions(NotMatlabHFM_Isotropic PRIVATE Isotropic)
endif()

option(RiemannExecutable "RiemannExecutable")
if(RiemannExecutable)
add_executable(NotMatlabHFM_Riemann "MatlabHFM.cpp")
target_compile_definitions(NotMatlabHFM_Riemann PRIVATE Riemann)
endif()


# Other optional executables

option(CustomExecutable "CustomExecutable")
if(CustomExecutable)
add_executable(NotMatlabHFM_Custom "MatlabHFM.cpp")
target_compile_definitions(NotMatlabHFM_Custom PRIVATE Custom)
endif()


option(Experimental0Executable "Experimental0Executable")
if(Experimental0Executable)
add_executable(NotMatlabHFM_Experimental0 "MatlabHFM.cpp")
target_compile_definitions(NotMatlabHFM_Experimental0 PRIVATE Experimental0)
endif()

# Headers IDE layout
source_group("DataStructures" FILES ${DataStructures_Headers})
source_group("LinearAlgebra" FILES ${LinearAlgebra_Headers})
source_group("Output" FILES ${Output_Headers})

source_group("Base" FILES ${Base_Headers})
source_group("Base\\Implementation" FILES ${Base_Implementation_Headers})
source_group("Specializations" FILES ${Specializations_Headers})
source_group("Experimental" FILES ${Experimental_Headers})
source_group("ExtraAlgorithms" FILES ${ExtraAlgorithms_Headers})