cmake_minimum_required(VERSION 3.1)
set (CMAKE_CXX_STANDARD 17)

project(pybind11_hfm)

FIND_PACKAGE(pybind11)

include(${CMAKE_CURRENT_SOURCE_DIR}/../common.cmake)

INCLUDE_DIRECTORIES(../../Headers ../../JMM_CPPLibs)


# ------- Declare the different builds -------

foreach(model ${ModelNames})
pybind11_add_module(HFM_${model} src/main.cpp)
target_compile_definitions(HFM_${model} PRIVATE ModelName=${model} PythonModuleName=HFM_${model})
endforeach(model)

if(CustomExecutable)
pybind11_add_module(PythonHFM_Custom src/main.cpp ${Project_Headers})
target_compile_definitions(PythonHFM_Custom PRIVATE Custom PythonModuleName=PythonHFM_Custom)
endif()
