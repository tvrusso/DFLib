cmake_minimum_required(VERSION 2.6)
project(DFLib)

# set path to additional CMake modules
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

FIND_PACKAGE(Proj)

add_library(DFLib DF_Abstract_Report.cpp DF_Report_Collection.cpp DF_ProjReport_Collection.cpp DF_XY_Point.cpp DF_LatLon_Point.cpp DF_Proj_Point.cpp DF_Proj_Report.cpp Util_Minimization_Methods.cpp gaussian_random.cpp)

include_directories(${DFLib_SOURCE_DIR} ${PROJ_INCLUDE_DIR} )
link_directories(${DFLib_BINARY_DIR} ${PROJ_LIB_DIR})

add_executable(testlsDF_proj testlsDF_proj.cpp)
target_link_libraries(testlsDF_proj DFLib ${PROJ_LIBRARY})

add_executable(SimpleDF SimpleDF.cpp)
target_link_libraries(SimpleDF DFLib ${PROJ_LIBRARY})
