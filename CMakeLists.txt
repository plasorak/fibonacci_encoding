project(Compress)
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(CMAKE_BUILD_TYPE Debug)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
# set(Gperftools_ROOT_DIR "/Users/plasorak/Dropbox/Sussex_Postdoc/gperftools/")
# set(Gperftools_DIR "/Users/plasorak/Dropbox/Sussex_Postdoc/gperftools/")
string(CONCAT CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/Users/plasorak/Dropbox/Sussex_Postdoc/CompressionStudies/")
#set(Gperftools_ROOT_DIR "/Users/plasorak/Dropbox/Sussex_Postdoc/gperftools/")
find_package(ROOT REQUIRED COMPONENTS RIO Net XMLIO)
include(${ROOT_USE_FILE})
find_package(Gperftools REQUIRED)

add_executable(Compress Compress.cxx)
target_include_directories(Compress PRIVATE ${Gperftools_INCLUDE_DIR})
target_link_libraries(Compress ${ROOT_LIBRARIES} ${Gperftools_PROFILER})
