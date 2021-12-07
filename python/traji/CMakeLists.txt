# TODO: try to implement pybind11 version

add_cython_target(_lib CXX)
add_library(_lib MODULE ${_lib})
target_link_libraries(_lib traji)
python_extension_module(_lib)
install(TARGETS _lib LIBRARY DESTINATION traji)
