cmake_minimum_required(VERSION 3.1)
project(pybind11_hfm)

FIND_PACKAGE(JMM_CPPLibs)
FIND_PACKAGE(pybind11)

include(${CMAKE_CURRENT_SOURCE_DIR}/../common.cmake)

INCLUDE_DIRECTORIES(../../Headers ${JMM_CPPLibs_INCLUDE_DIRS})


# ------- Declare the different builds -------

foreach(model ${ModelNames})
pybind11_add_module(HFM_${model} src/main.cpp)
target_compile_definitions(HFM_${model} PRIVATE ModelName=${model} PythonModuleName=HFM_${model})
endforeach(model)

if(CustomExecutable)
pybind11_add_module(PythonHFM_Custom src/main.cpp ${Project_Headers})
target_compile_definitions(PythonHFM_Custom PRIVATE Custom PythonModuleName=PythonHFM_Custom)
endif()
