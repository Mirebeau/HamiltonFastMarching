set(Base_Implementation_Headers
	${CMAKE_CURRENT_SOURCE_DIR}/EulerianStencil.hxx
	${CMAKE_CURRENT_SOURCE_DIR}/PeriodicGrid.hxx
	${CMAKE_CURRENT_SOURCE_DIR}/EulerianStencil.hxx
	${CMAKE_CURRENT_SOURCE_DIR}/HamiltonFastMarching.hxx
	${CMAKE_CURRENT_SOURCE_DIR}/HFM_ParamDefault.hxx
	${CMAKE_CURRENT_SOURCE_DIR}/HFMInterface.hxx
	${CMAKE_CURRENT_SOURCE_DIR}/GeodesicODESolver.hxx
	${CMAKE_CURRENT_SOURCE_DIR}/GeodesicDiscreteSolver.hxx
	${CMAKE_CURRENT_SOURCE_DIR}/HFM_StencilLag2.hxx
	${CMAKE_CURRENT_SOURCE_DIR}/HFM_StencilLag3.hxx
	${CMAKE_CURRENT_SOURCE_DIR}/HFM_StencilRecompute.hxx
	${CMAKE_CURRENT_SOURCE_DIR}/HFM_StencilShare.hxx

	
	PARENT_SCOPE
)