cmake_minimum_required(VERSION 3.1)
set (CMAKE_CXX_STANDARD 14)

project(MathematicaHFM)

#set a few machine dependent paths
#add_subdirectory("../LocalSettings" "${CMAKE_CURRENT_BINARY_DIR}/LocalSettings")
set(JMM_CPPLibs_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../../JMM_CPPLibs-master" CACHE FILEPATH "JMM_CPPLibs directory")

#paths required for compilation of Lifted Fast Marching
set(LinearAlgebra_dir 		"${JMM_CPPLibs_dir}/LinearAlgebra")
set(Output_dir 				"${JMM_CPPLibs_dir}/Output")
set(DataStructures_dir 		"${JMM_CPPLibs_dir}/DataStructures")


set(DummyBinDir "${CMAKE_CURRENT_BINARY_DIR}/Dummy")

#LiftedFastMarching headers
add_subdirectory(${LinearAlgebra_dir} "${DummyBinDir}/LinearAlgebra")
add_subdirectory(${Output_dir} "${DummyBinDir}/Output")
add_subdirectory(${DataStructures_dir} 	"${DummyBinDir}/DataStructures")
add_subdirectory("../../Headers" "${DummyBinDir}/HFMHeaders")
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
	"${Output_dir}/Output/FileIO.h"
	"${Output_dir}/Output/FileIO.hxx"
	"${Output_dir}/Output/IO.h"
	"${Output_dir}/Output/IO.hxx"
	"${Output_dir}/Output/ExceptionMacro.h"
	"${Output_dir}/Output/MathematicaIO.h"
	"${Output_dir}/Output/MathematicaIO.hxx"
	"${Output_dir}/Output/MathematicaIOExternC.h"
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

include_directories(${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

######################## Additional lines for Wolfram Library link

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/Mathematica")

set (Mathematica_USE_LIBCXX_LIBRARIES ON)

# finds newest Mathematica installation and its default components
find_package(Mathematica)
if (${Mathematica_FOUND})
	message(STATUS "Mathematica Version ${Mathematica_VERSION}")
	# Mathematica_SYSTEM_IDS contains the list of Mathematica platform system IDs that the
	# project is being compiled for. This usually contains just one entry (e.g., "Windows").
	# It may contain multiple entries if we are building a universal binary under Mac OS X
	# (e.g., "MacOSX-x86-64;MacOSX-x86").
	message(STATUS "Mathematica Target System IDs ${Mathematica_SYSTEM_IDS}")
	# Mathematica_HOST_SYSTEM_IDS is the list of Mathematica platform system IDs that can
	# run on the build host. If we are executing the CMake build under a 64-bit version of
	# Windows this would be "Windows-x86-64;Windows". Under a 32-bit version of Windows this
	# would be just "Windows".
	message(STATUS "Mathematica Host System IDs ${Mathematica_HOST_SYSTEM_IDS}")
	# Mathematica_BASE_DIR is the directory for systemwide files to be loaded by Mathematica
	message(STATUS "Mathematica Base Directory ${Mathematica_BASE_DIR}")
	# Mathematica_USERBASE_DIR is the directory for user-specific files to be loaded by Mathematica
	message(STATUS "Mathematica User Base Directory ${Mathematica_USERBASE_DIR}")
endif()
if (${Mathematica_WolframLibrary_FOUND})
	# library found, so let's compile code
	include_directories(${Mathematica_INCLUDE_DIRS})
endif()


# ################## Include the new library link project
option(Executable_AllBase "Executable_AllBase" TRUE)
if(Executable_AllBase)
add_library(MathematicaHFM_AllBase SHARED "MathematicaHFM.cpp" ${Project_Headers})
target_link_libraries(MathematicaHFM_AllBase -Wl, -lstdc++ -lpthread -Wl, ${Mathematica_Libraries})
target_compile_definitions(MathematicaHFM_AllBase PRIVATE AllBaseModels)
endif()

option(Executable_RiemannExtra "Executable_RiemannExtra")
if(Executable_RiemannExtra)
add_library(MathematicaHFM_RiemannExtra SHARED "MathematicaHFM.cpp")
target_link_libraries(MathematicaHFM_RiemannExtra -Wl, -lstdc++ -lpthread -Wl, ${Mathematica_Libraries})
target_compile_definitions(MathematicaHFM_RiemannExtra PRIVATE RiemannExtra)
endif()



# Note from JMM to Erik Bekkers: For OSX compilation, I had to delete all static bindings
#target_link_libraries(MathematicaHFM -static-libgcc -static-libstdc++ -Wl,-Bstatic -lstdc++ -lpthread -Wl,-Bdynamic)



# Headers IDE layout
source_group("DataStructures" FILES ${DataStructures_Headers})
source_group("LinearAlgebra" FILES ${LinearAlgebra_Headers})
source_group("Output" FILES ${Output_Headers})

source_group("Base" FILES ${Base_Headers})
source_group("Base\\Implementation" FILES ${Base_Implementation_Headers})
source_group("Specializations" FILES ${Specializations_Headers})
source_group("Experimental" FILES ${Experimental_Headers})
source_group("ExtraAlgorithms" FILES ${ExtraAlgorithms_Headers})