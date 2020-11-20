set(ExtraAlgorithms_Headers
	${CMAKE_CURRENT_SOURCE_DIR}/CommonStoppingCriteria.h
	${CMAKE_CURRENT_SOURCE_DIR}/EuclideanPathLength.h
	${CMAKE_CURRENT_SOURCE_DIR}/FirstVariation.h
	${CMAKE_CURRENT_SOURCE_DIR}/TimeDependentFields.h
	${CMAKE_CURRENT_SOURCE_DIR}/VoronoiDiagram.h
	${CMAKE_CURRENT_SOURCE_DIR}/WallObstructionTest.h
#	${CMAKE_CURRENT_SOURCE_DIR}/DynamicFactoring.h
	${CMAKE_CURRENT_SOURCE_DIR}/StaticFactoring.h
	PARENT_SCOPE
	)

add_subdirectory(Implementation)
set(ExtraAlgorithms_Implementation_Headers ${ExtraAlgorithms_Implementation_Headers} PARENT_SCOPE)
