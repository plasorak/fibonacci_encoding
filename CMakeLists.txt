project(Compress)
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net XMLIO)
include(${ROOT_USE_FILE})

add_executable(Compress Compress.cxx)
target_link_libraries(Compress ${ROOT_LIBRARIES})
