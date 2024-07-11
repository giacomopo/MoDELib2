

add_executable(DDqt DDqt.cxx)

target_link_libraries(DDqt MoDELib)
vtk_module_autoinit(TARGETS DDqt MODULES ${VTK_LIBRARIES} )
