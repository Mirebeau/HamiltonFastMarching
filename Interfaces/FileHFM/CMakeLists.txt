cmake_minimum_required(VERSION 3.1)
set (CMAKE_CXX_STANDARD 14)

project(FileHFM)

include(${CMAKE_CURRENT_SOURCE_DIR}/../common.cmake)
include_directories(${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

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

