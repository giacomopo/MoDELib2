
#if(POLICY CMP0057)
#  cmake_policy(SET CMP0057 NEW)
#endif()
#cmake_policy(SET CMP0057 NEW)

pybind11_add_module(pyMoDELib SHARED pyMoDELib.cxx)
target_link_libraries(pyMoDELib PUBLIC MoDELib pybind11::module)
