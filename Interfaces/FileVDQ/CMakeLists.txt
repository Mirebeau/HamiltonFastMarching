cmake_minimum_required(VERSION 3.1)
set (CMAKE_CXX_STANDARD 17)

project(FileVDQ)

#set(JMM_CPPLibs_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../../../../Github/JMM_CPPLibs" CACHE FILEPATH "JMM_CPPLibs directory")
set(DummyBinDir "${CMAKE_CURRENT_BINARY_DIR}/Dummy")
add_subdirectory("../../JMM_CPPLibs" "${DummyBinDir}/JMM_CPPLibs")
add_subdirectory("Headers")

set(Project_Headers 
 	${Root_Headers}

 	${DataStructures_Headers}
 	${LinearAlgebra_Headers}
 	${Output_Headers}
 	${Macros_Headers}
)
#source_group("Headers" FILES ${Project_Headers})
source_group("DataStructures" FILES ${DataStructures_Headers})
source_group("Macros" FILES ${Macros_Headers})
source_group("LinearAlgebra" FILES ${LinearAlgebra_Headers})
source_group("Output" FILES ${Output_Headers})


include_directories(${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(FileVDQ "FileVDQ.cpp" ${Project_Headers})
