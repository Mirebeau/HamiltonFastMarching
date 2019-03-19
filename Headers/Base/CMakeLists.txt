set(Base_Headers
	${CMAKE_CURRENT_SOURCE_DIR}/BaseGrid.h
	${CMAKE_CURRENT_SOURCE_DIR}/PeriodicGrid.h
	${CMAKE_CURRENT_SOURCE_DIR}/CommonStencil.h
	${CMAKE_CURRENT_SOURCE_DIR}/EulerianStencil.h
	${CMAKE_CURRENT_SOURCE_DIR}/Lagrangian2Stencil.h
	${CMAKE_CURRENT_SOURCE_DIR}/Lagrangian3Stencil.h
	${CMAKE_CURRENT_SOURCE_DIR}/HamiltonFastMarching.h
	${CMAKE_CURRENT_SOURCE_DIR}/HFMInterface.h
	${CMAKE_CURRENT_SOURCE_DIR}/GeodesicODESolver.h
	${CMAKE_CURRENT_SOURCE_DIR}/GeodesicDiscreteSolver.h

	PARENT_SCOPE
	)

add_subdirectory(Implementation)
set(Base_Implementation_Headers ${Base_Implementation_Headers} PARENT_SCOPE)
