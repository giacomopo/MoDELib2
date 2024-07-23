target_sources(${PROJECT_NAME} PRIVATE 
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationLoopPatches.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationFieldBase.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/EshelbyInclusionBase.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/SphericalInclusion.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/PolyhedronInclusion.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/StraightDislocationSegment.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/StressStraight.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/VoigtTraits.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/DislocationDynamicsBase.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/ClusterDynamicsParameters.cpp
#		${CMAKE_CURRENT_SOURCE_DIR}/ClusterDynamicsFEM.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/ElasticDeformationFEM.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/MicrostructureContainer.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/MicrostructureBase.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/InclusionMicrostructureBase.cpp
		)
