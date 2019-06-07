cmake_minimum_required(VERSION 3.1)
set (CMAKE_CXX_STANDARD 17)


project(FileHFM)

# Include headers (../../JMM_CPPLibs and ../../Headers)
set(DummyBinDir "${CMAKE_CURRENT_BINARY_DIR}/Dummy")
add_subdirectory("../../JMM_CPPLibs" "${DummyBinDir}/JMM_CPPLibs")
add_subdirectory("../../Headers" "${DummyBinDir}/HFMHeaders")


include(${CMAKE_CURRENT_SOURCE_DIR}/../common.cmake)
include(ide_headers_layout.cmake)
include_directories(${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

cmake_policy(SET CMP0057 NEW)
if("Seismic3SIMD" IN_LIST ModelNames)
	set(XSIMD_dir "/Users/mirebeau/Dropbox/Programmes/2019/4_Avril/SimdTest/xsimd-master/include" CACHE FILEPATH "xsimd directory")
	include_directories(${XSIMD_dir})
endif()

# Single model executables

foreach(model ${ModelNames})
add_executable(FileHFM_${model} "FileHFM.cpp")
target_compile_definitions(FileHFM_${model} PRIVATE ModelName=${model})
endforeach(model)


# Optional custom executable

if(CustomExecutable)
add_executable(FileHFM_Custom "FileHFM.cpp" ${Project_Headers})
target_compile_definitions(FileHFM_Custom PRIVATE Custom)
endif()

