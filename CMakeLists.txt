
# https://mcilloni.ovh/2021/02/09/cxx-cross-clang/
#export LDFLAGS=-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib 
#LDFLAGS="$LDFLAGS -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib/system"


# cmake .. -DCMAKE_CXX_COMPILER=clang++-mp-18

cmake_minimum_required(VERSION 3.1.0)
project(MoDELib VERSION 2.0.0 LANGUAGES CXX)

option(USE_PYBIND11 "enables pybind11" ON)
option(BUILD_TOOLS "builds tools" ON)
set(EIGEN3_INCLUDE_DIRS /opt/local/include/eigen3)
list(APPEND CMAKE_PREFIX_PATH "/Users/Giacomo/Applications/Qt/6.5.2/macos/lib/cmake")

cmake_policy(SET CMP0057 NEW)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -march=native -fopenmp -Ofast -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-parameter")
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

### Boost ### 
find_package(Boost)
if(Boost_FOUND)
	message(STATUS "Boost_VERSION:" ${Boost_VERSION})
	include_directories(${Boost_INCLUDE_DIRS})
else()
	message(WARNING "Boost not found")
endif()

### fftw3 ### find_package is currently broken. Workaround using find_library
find_library(FFTW3_LIBRARIES fftw3)  
if(FFTW3_LIBRARIES)
	cmake_path(GET FFTW3_LIBRARIES PARENT_PATH FFTW3_LIBRARY_DIRS)
	set(FFTW3_INCLUDE_DIRS ${FFTW3_LIBRARY_DIRS}/../include)
	message(STATUS "FFTW3_LIBRARIES:" ${FFTW3_LIBRARIES})
	message(STATUS "FFTW3_LIBRARY_DIRS:" ${FFTW3_LIBRARY_DIRS})
	message(STATUS "FFTW3_INCLUDE_DIRS:" ${FFTW3_INCLUDE_DIRS})
	include_directories(${BOOST_INCLUDE_DIRS})
	include_directories(${FFTW3_INCLUDE_DIRS})
	add_definitions(-D _MODEL_GLIDE_PLANE_NOISE_GENERATOR_)
else()
  message(WARNING "library ${FFTW3_LIBRARIES} not found, noise generator is disabled")
endif()

### pybind11 ###
if(USE_PYBIND11)
	find_package(Python COMPONENTS Interpreter Development)
	if(Python_FOUND)
		message(STATUS "Python_VERSION:" ${Python_VERSION})
		find_package(pybind11 CONFIG REQUIRED)
	else()
		message(WARNING "Python not found, pybind11 is disabled")
	endif()
endif()

if(pybind11_FOUND)
	message(STATUS "pybind11_VERSION:" ${pybind11_VERSION})
	set(pybind11_LIBRARIES pybind11::embed pybind11::headers pybind11::module)
	message(STATUS "pybind11_LIBRARIES: ${pybind11_LIBRARIES}")
	add_definitions(-D_MODEL_PYBIND11_)
	Python_add_library(MoDELib SHARED "")
else()
	message(WARNING "pybind11 is disabled")
	add_library(MoDELib "")
endif()

set(MoDELib_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
include_directories ( ${EIGEN3_INCLUDE_DIRS} )
include_directories(${MoDELib_INCLUDE_DIR}/Utilities/)
include_directories(${MoDELib_INCLUDE_DIR}/IO/)
include_directories(${MoDELib_INCLUDE_DIR}/Math/)
include_directories(${MoDELib_INCLUDE_DIR}/Geometry/)
include_directories(${MoDELib_INCLUDE_DIR}/Mesh/)
include_directories(${MoDELib_INCLUDE_DIR}/Lattices/)
include_directories(${MoDELib_INCLUDE_DIR}/FEM/)
include_directories(${MoDELib_INCLUDE_DIR}/GlidePlanes/)
include_directories(${MoDELib_INCLUDE_DIR}/Quadrature/)
include_directories(${MoDELib_INCLUDE_DIR}/PolycrystallineMaterials/)
include_directories(${MoDELib_INCLUDE_DIR}/DislocationDynamicsIO/)
include_directories(${MoDELib_INCLUDE_DIR}/DislocationMobilities/)
include_directories(${MoDELib_INCLUDE_DIR}/DislocationDynamicsBase/)
include_directories(${MoDELib_INCLUDE_DIR}/DislocationMicrostructure/)
include_directories(${MoDELib_INCLUDE_DIR}/LoopNetwork/)
include_directories(${MoDELib_INCLUDE_DIR}/ParticleInteraction/)
include_directories(${MoDELib_INCLUDE_DIR}/SpaceDecomposition/)
include_directories(${MoDELib_INCLUDE_DIR}/DiscreteCrackMechanics/)
include_directories(${MoDELib_INCLUDE_DIR}/ClusterDynamics/)
include_directories(${MoDELib_INCLUDE_DIR}/ElasticDeformation/)
include_directories(${MoDELib_INCLUDE_DIR}/InclusionMicrostructure/)
include_directories(${MoDELib_INCLUDE_DIR}/DislocationDynamics/)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src/ src)
target_link_libraries(MoDELib PUBLIC ${FFTW3_LIBRARIES} ${pybind11_LIBRARIES})

if(BUILD_TOOLS)
	add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tools tools)
else()
endif()
