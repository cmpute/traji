cmake_minimum_required(VERSION 3.5)
project(traji CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(Boost REQUIRED COMPONENTS)
find_package(Eigen3 3.3 REQUIRED)

add_library(traji src/Path.cpp src/Frenet.cpp src/Trajectory.cpp)
# target_include_directories(
#   traji PUBLIC ${Boost_INCLUDE_DIRS} ${Eigen_INCLUDE_DIRS}
# )
set_target_properties(traji PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(traji PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries (traji PUBLIC Boost::boost)
target_link_libraries (traji PUBLIC Eigen3::Eigen)

##### Installation #####

install(
  FILES include/traji.hpp
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include
)
install(TARGETS traji
        DESTINATION lib
        EXPORT trajiTargets)
install(EXPORT trajiTargets
        FILE trajiTargets.cmake
        DESTINATION lib/cmake/traji
)
export(EXPORT trajiTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/trajiTargets.cmake"
)

include(CMakePackageConfigHelpers)
# generate the config file that is includes the exports
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/traji-config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/trajiConfig.cmake"
  INSTALL_DESTINATION "lib/cmake/traji"
  # NO_SET_AND_CHECK_MACRO
  # NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/trajiConfigVersion.cmake"
  VERSION "${traji_VERSION_MAJOR}.${traji_VERSION_MINOR}"
  COMPATIBILITY AnyNewerVersion
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/trajiConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/trajiConfigVersion.cmake
  DESTINATION lib/cmake/traji
)
