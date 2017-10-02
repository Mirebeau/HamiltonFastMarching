cmake_minimum_required(VERSION 2.8.3)
set (CMAKE_CXX_STANDARD 14)

project(PythonHFM)

# -------- Header include directories --------

set(JMM_CPPLibs_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../../JMM_CPPLibs-master" CACHE FILEPATH "JMM_CPPLibs directory")
set(LinearAlgebra_dir 		"${JMM_CPPLibs_dir}/LinearAlgebra")
set(Output_dir 				"${JMM_CPPLibs_dir}/Output")
set(DataStructures_dir 		"${JMM_CPPLibs_dir}/DataStructures")

set(DummyBinDir "${CMAKE_CURRENT_BINARY_DIR}/Dummy")

add_subdirectory(${LinearAlgebra_dir} "${DummyBinDir}/LinearAlgebra")
add_subdirectory(${Output_dir} "${DummyBinDir}/Output")
add_subdirectory(${DataStructures_dir} 	"${DummyBinDir}/DataStructures")
add_subdirectory("../../Headers" "${DummyBinDir}/HFMHeaders")

set(Output_Headers 
	${Output_Headers} 
	"${Output_dir}/Output/IO.h"
	"${Output_dir}/Output/IO.hxx"
	"${Output_dir}/Output/BaseIO.h"
	"${Output_dir}/Output/BaseIO.hxx"
	"${Output_dir}/Output/PythonIO.h"
	"${Output_dir}/Output/PythonIO.hxx"
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
option(Executable_AllBase "Executable_AllBase" TRUE)
if(Executable_AllBase)
PYTHON_ADD_MODULE(PythonHFM_AllBase PythonHFM.cpp ${Project_Headers})
target_compile_definitions(PythonHFM_AllBase PRIVATE AllBaseModels PythonModuleName=PythonHFM_AllBase)
endif()

option(Executable_RiemannExtra "Executable_RiemannExtra")
if(Executable_RiemannExtra)
PYTHON_ADD_MODULE(PythonHFM_RiemannExtra PythonHFM.cpp)
target_compile_definitions(PythonHFM_RiemannExtra PRIVATE RiemannExtra PythonModuleName=PythonHFM_RiemannExtra)
endif()

option(Executable_PrescribedCurvature2 "Executable_PrescribedCurvature2")
if(Executable_PrescribedCurvature2)
PYTHON_ADD_MODULE(PythonHFM_PrescribedCurvature2 PythonHFM.cpp)
target_compile_definitions(PythonHFM_PrescribedCurvature2 PRIVATE PrescribedCurvature2 PythonModuleName=PythonHFM_PrescribedCurvature2)
endif()

option(Executable_Isotropic "Executable_Isotropic")
if(Executable_Isotropic)
PYTHON_ADD_MODULE(PythonHFM_Isotropic PythonHFM.cpp)
target_compile_definitions(PythonHFM_Isotropic PRIVATE Isotropic PythonModuleName=PythonHFM_Isotropic)
endif()


# ------- Headers IDE layout -------
source_group("DataStructures" FILES ${DataStructures_Headers})
source_group("LinearAlgebra" FILES ${LinearAlgebra_Headers})
source_group("Output" FILES ${Output_Headers})

source_group("Base" FILES ${Base_Headers})
source_group("Base\\Implementation" FILES ${Base_Implementation_Headers})
source_group("Specializations" FILES ${Specializations_Headers})
source_group("Experimental" FILES ${Experimental_Headers})
source_group("ExtraAlgorithms" FILES ${ExtraAlgorithms_Headers})