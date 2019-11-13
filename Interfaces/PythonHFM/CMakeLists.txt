cmake_minimum_required(VERSION 3.3)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_OSX_DEPLOYMENT_TARGET 10.14 CACHE STRING "Minimum OS X deployment version" FORCE)
# Needed for pybind11 2.4.3 under osx. 
#https://stackoverflow.com/questions/34208360/cmake-seems-to-ignore-cmake-osx-deployment-target

project(pybind11_hfm)
FIND_PACKAGE(pybind11)

message("pybind11 version : " ${pybind11_VERSION} ", CMAKE_CXX_COMPILER_ID " ${CMAKE_CXX_COMPILER_ID})
message("CMAKE_OSX_DEPLOYMENT_TARGET " ${CMAKE_OSX_DEPLOYMENT_TARGET})


include(${CMAKE_CURRENT_SOURCE_DIR}/../common.cmake)

INCLUDE_DIRECTORIES(../../Headers ../../JMM_CPPLibs)

# ------- Declare the different builds -------

foreach(model ${ModelNames})
pybind11_add_module(HFM_${model} src/main.cpp)
target_compile_definitions(HFM_${model} PRIVATE ModelName=${model} PythonModuleName=HFM_${model})

# Fixes an issue with pybind11 on clang, see https://github.com/pybind/pybind11/issues/1604
#if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang") 
#	message("Applying fix for pybind11 under AppleClang")
#    target_compile_options(HFM_${model} INTERFACE -fsized-deallocation)
#endif()

endforeach(model)


if(CustomExecutable)
pybind11_add_module(PythonHFM_Custom src/main.cpp ${Project_Headers})
target_compile_definitions(HFM_Custom PRIVATE Custom PythonModuleName=PythonHFM_Custom)

# Fixes an issue with pybind11 on clang, see https://github.com/pybind/pybind11/issues/1604
#if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
#    target_compile_options(HFM_Custom INTERFACE -fsized-deallocation)
#endif()

endif()
