set(Root_Headers
	${CMAKE_CURRENT_SOURCE_DIR}/DispatchAndRun.h
	PARENT_SCOPE
	)

add_subdirectory(Base)
set(Base_Headers ${Base_Headers} PARENT_SCOPE)
set(Base_Implementation_Headers ${Base_Implementation_Headers} PARENT_SCOPE)

add_subdirectory(Specializations)
set(Specializations_Headers ${Specializations_Headers} PARENT_SCOPE)
set(Specializations_Implementation_Headers ${Specializations_Implementation_Headers} PARENT_SCOPE)

add_subdirectory(ExtraAlgorithms)
set(ExtraAlgorithms_Headers ${ExtraAlgorithms_Headers} PARENT_SCOPE)
set(ExtraAlgorithms_Implementation_Headers ${ExtraAlgorithms_Implementation_Headers} PARENT_SCOPE)

add_subdirectory(Experimental)
set(Experimental_Headers ${Experimental_Headers} PARENT_SCOPE)
set(Experimental_Implementation_Headers ${Experimental_Implementation_Headers} PARENT_SCOPE)

set(PROJECT_DIR ${PROJECT_DIR} "${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
