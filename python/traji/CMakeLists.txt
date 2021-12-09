pybind11_add_module(_bindings _bindings.cxx)
target_link_libraries(_bindings PUBLIC traji)
install(TARGETS _bindings LIBRARY DESTINATION traji)
