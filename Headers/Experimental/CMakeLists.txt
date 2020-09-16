set(Experimental_Headers
	${CMAKE_CURRENT_SOURCE_DIR}/AsymmetricQuadratic.h
	${CMAKE_CURRENT_SOURCE_DIR}/PrescribedCurvature2.h
	${CMAKE_CURRENT_SOURCE_DIR}/RiemannLifted.h
	${CMAKE_CURRENT_SOURCE_DIR}/Quaternionic.h
	${CMAKE_CURRENT_SOURCE_DIR}/ReedsSheppAdaptive2.h
	${CMAKE_CURRENT_SOURCE_DIR}/IsotropicBox.h
	${CMAKE_CURRENT_SOURCE_DIR}/Differentiable.h
	${CMAKE_CURRENT_SOURCE_DIR}/RollingBall.h
	${CMAKE_CURRENT_SOURCE_DIR}/Seismic2.h
	${CMAKE_CURRENT_SOURCE_DIR}/Seismic3.h
	${CMAKE_CURRENT_SOURCE_DIR}/AlignedBillard.h
	${CMAKE_CURRENT_SOURCE_DIR}/TTI.h
	${CMAKE_CURRENT_SOURCE_DIR}/AsymRander.h
	PARENT_SCOPE
	)

add_subdirectory(Implementation)
set(Experimental_Implementation_Headers ${Experimental_Implementation_Headers} PARENT_SCOPE)
