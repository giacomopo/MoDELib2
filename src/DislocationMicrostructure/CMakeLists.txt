target_sources(${PROJECT_NAME} PRIVATE 
	${CMAKE_CURRENT_SOURCE_DIR}/MicrostructureGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/MicrostructureGeneratorBase.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PeriodicDipoleGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PlanarLoopGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PeriodicLoopGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PrismaticLoopGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SphericalInclusionsGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PolyhedronInclusionsGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/StackingFaultTetrahedraGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/FrankLoopsGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/VTKGenerator.cpp
			)


