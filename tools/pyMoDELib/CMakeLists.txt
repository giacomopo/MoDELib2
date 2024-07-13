

### pybind11 ###
find_package(Python COMPONENTS Interpreter Development)
if(Python_FOUND)
if(POLICY CMP0057)
  cmake_policy(SET CMP0057 NEW)
endif()
	message(STATUS "Python_VERSION:" ${Python_VERSION})
	find_package(pybind11 CONFIG REQUIRED)
if(pybind11_FOUND)
	message(STATUS "pybind11_VERSION:" ${pybind11_VERSION})
	set(pybind11_LIBRARIES pybind11::embed)
	message(STATUS "pybind11_LIBRARIES: ${pybind11_LIBRARIES}")
	add_definitions(-D_MODEL_PYBIND11_)
else()
	message(WARNING "pybind11 not found, pybind11 is disabled")
endif()
else()
	message(WARNING "Python not found, pybind11 is disabled")
endif()

pybind11_add_module(pyMoDELib SHARED pyMoDELib.cxx)


target_link_libraries(pyMoDELib PUBLIC MoDELib pybind11::module)
