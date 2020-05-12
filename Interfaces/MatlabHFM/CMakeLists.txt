# ------- THIS DOES NOT BUILD MATLAB HFM ------

# This file is only used for compilation tests with xcode
# Please use the provided Matlab(R) script in the ExampleFiles directory to compile the project using mex

cmake_minimum_required(VERSION 3.1)
set (CMAKE_CXX_STANDARD 17)

project(NotMatlabHFM)

set(DummyBinDir "${CMAKE_CURRENT_BINARY_DIR}/Dummy")
add_subdirectory("../../JMM_CPPLibs" "${DummyBinDir}/JMM_CPPLibs")
add_subdirectory("../../Headers" "${DummyBinDir}/HFMHeaders")

#set a machine dependent path
set(Matlab_dir "/Applications/MATLAB_R2020a.app/extern/include" CACHE FILEPATH "Matlab directory")

include(${CMAKE_CURRENT_SOURCE_DIR}/../common.cmake)
include_directories(${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${Matlab_dir})



# Single model executables

foreach(model ${ModelNames})
add_executable(NotMatlabHFM_${model} "MatlabHFM.cpp")
target_compile_definitions(NotMatlabHFM_${model} PRIVATE ModelName=${model})
endforeach(model)


# Optional custom executable

if(CustomExecutable)
add_executable(NotMatlabHFM_Custom "MatlabHFM.cpp" ${Project_Headers})
target_compile_definitions(NotMatlabHFM_Custom PRIVATE Custom)
endif()

