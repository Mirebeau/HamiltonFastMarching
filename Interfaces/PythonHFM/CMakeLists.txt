cmake_minimum_required(VERSION 2.8.3)
set (CMAKE_CXX_STANDARD 14)

project(PythonHFM)
include(${CMAKE_CURRENT_SOURCE_DIR}/../../common.cmake)

#-------- Link with boost, python, boost-numpy ------

#FILE(COPY hello.py DESTINATION .)

FIND_PACKAGE(PythonInterp)
FIND_PACKAGE(PythonLibs)
if(APPLE AND ${PYTHON_VERSION_MAJOR} EQUAL 3)
    FIND_PACKAGE(Boost COMPONENTS python3)
else()
    FIND_PACKAGE(Boost COMPONENTS python)
endif()

find_library(BoostNumpy NAMES boost_numpy)
INCLUDE_DIRECTORIES(${PROJECT_DIR} ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
LINK_LIBRARIES(${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${BoostNumpy}) # Deprecated but so convenient!

# ------- Declare the different builds -------

foreach(model ${ModelNames})
PYTHON_ADD_MODULE(PythonHFM_${model} PythonHFM.cpp)
target_compile_definitions(PythonHFM_${model} PRIVATE ModelName=${model} PythonModuleName=PythonHFM_${model})
endforeach(model)

if(CustomExecutable)
PYTHON_ADD_MODULE(PythonHFM_Custom PythonHFM.cpp ${Project_Headers})
target_compile_definitions(PythonHFM_Custom PRIVATE Custom PythonModuleName=PythonHFM_Custom)
endif()


# --- Old ----
#option(Executable_RiemannExtra "Executable_RiemannExtra")
#if(Executable_RiemannExtra)
#PYTHON_ADD_MODULE(PythonHFM_RiemannExtra PythonHFM.cpp)
#target_compile_definitions(PythonHFM_RiemannExtra PRIVATE RiemannExtra PythonModuleName=PythonHFM_RiemannExtra)
#endif()
