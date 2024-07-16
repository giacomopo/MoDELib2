
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/MicrostructureGenerator MicrostructureGenerator)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/DDomp DDomp)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/DDqt DDqt)
if(USE_PYBIND11)
	add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/pyMoDELib pyMoDELib)
endif()

