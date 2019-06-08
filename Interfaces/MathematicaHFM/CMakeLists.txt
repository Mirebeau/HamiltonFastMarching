cmake_minimum_required(VERSION 3.1)
set (CMAKE_CXX_STANDARD 17)

project(MathematicaHFM)

# Include headers (../../JMM_CPPLibs and ../../Headers)
set(DummyBinDir "${CMAKE_CURRENT_BINARY_DIR}/Dummy")
add_subdirectory("../../JMM_CPPLibs" "${DummyBinDir}/JMM_CPPLibs")
add_subdirectory("../../Headers" "${DummyBinDir}/HFMHeaders")

include(${CMAKE_CURRENT_SOURCE_DIR}/../common.cmake)
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

# Single model executables

foreach(model ${ModelNames})
add_library(MathematicaHFM_${model} SHARED "MathematicaHFM.cpp")
target_link_libraries(MathematicaHFM_${model} -Wl, -lstdc++ -lpthread -Wl, ${Mathematica_Libraries})
target_compile_definitions(MathematicaHFM_${model} PRIVATE ModelName=${model})
endforeach(model)


# Optional custom executable

if(CustomExecutable)
add_library(MathematicaHFM_Custom SHARED "MathematicaHFM.cpp" ${Project_Headers})
target_link_libraries(MathematicaHFM_Custom -Wl, -lstdc++ -lpthread -Wl, ${Mathematica_Libraries})
target_compile_definitions(MathematicaHFM_Custom PRIVATE Custom)
endif()


# Note from JMM to Erik Bekkers: For OSX compilation, I had to delete all static bindings
#target_link_libraries(MathematicaHFM -static-libgcc -static-libstdc++ -Wl,-Bstatic -lstdc++ -lpthread -Wl,-Bdynamic)
