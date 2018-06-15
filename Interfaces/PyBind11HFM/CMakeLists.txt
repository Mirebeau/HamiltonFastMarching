cmake_minimum_required(VERSION 2.8.12)
set(PYTHON_EXECUTABLE  "/Users/jean-mariemirebeau/anaconda3/bin/python3.6")

project(pybind11_hfm)
FIND_PACKAGE(pybind11)
include(${CMAKE_CURRENT_SOURCE_DIR}/../common.cmake)

INCLUDE_DIRECTORIES(${PROJECT_DIR})
#pybind11_add_module(pybind11_hfm src/main.cpp)


# ------- Declare the different builds -------

foreach(model ${ModelNames})
pybind11_add_module(PythonHFM_${model} src/main.cpp)
target_compile_definitions(PythonHFM_${model} PRIVATE ModelName=${model} PythonModuleName=PythonHFM_${model})
endforeach(model)

if(CustomExecutable)
pybind11_add_module(PythonHFM_Custom src/main.cpp ${Project_Headers})
target_compile_definitions(PythonHFM_Custom PRIVATE Custom PythonModuleName=PythonHFM_Custom)
endif()
