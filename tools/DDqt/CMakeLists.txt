### Qt ### #ADD_DEFINITIONS(-DQT_NO_KEYWORDS)
find_package(Qt6 COMPONENTS Widgets OpenGL REQUIRED)
message (STATUS "Qt6_VERSION: ${Qt6_VERSION}")
if (NOT Qt6_FOUND)
  message("Skipping DDqt: ${Qt6_NOT_FOUND_MESSAGE}")
  return()
endif()

### VTK ###
find_package(VTK REQUIRED)
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
if (NOT VTK_FOUND)
  message("Skipping DDvtk: ${VTK_NOT_FOUND_MESSAGE}")
  return()
endif()

add_executable(DDqt DDqt.cxx)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src src)

target_link_libraries(DDqt MoDELib ${FFTW3_LIBRARIES} ${VTK_LIBRARIES} Qt6::Widgets Qt6::OpenGL Qt6::Core)
vtk_module_autoinit(TARGETS DDqt MODULES ${VTK_LIBRARIES} Qt6::Widgets Qt6::OpenGL Qt6::Core)
