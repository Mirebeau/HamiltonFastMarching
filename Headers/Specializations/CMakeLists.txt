set(Specializations_Headers
	${CMAKE_CURRENT_SOURCE_DIR}/CommonTraits.h
	${CMAKE_CURRENT_SOURCE_DIR}/Isotropic.h
	${CMAKE_CURRENT_SOURCE_DIR}/Curvature2.h
	${CMAKE_CURRENT_SOURCE_DIR}/Curvature3.h
	${CMAKE_CURRENT_SOURCE_DIR}/Riemannian.h
	${CMAKE_CURRENT_SOURCE_DIR}/QuadLinLag2.h
	${CMAKE_CURRENT_SOURCE_DIR}/QuadLinLag2.h
	PARENT_SCOPE
	)

add_subdirectory(Implementation)
set(Specializations_Implementation_Headers ${Specializations_Implementation_Headers} PARENT_SCOPE)
