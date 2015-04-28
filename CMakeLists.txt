cmake_minimum_required(VERSION 2.6)
project(DFLib)

# set path to additional CMake modules
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

FIND_PACKAGE(Proj REQUIRED)

add_library(DFLib SHARED DF_Abstract_Report.cpp DF_Report_Collection.cpp DF_ProjReport_Collection.cpp DF_XY_Point.cpp DF_LatLon_Point.cpp DF_Proj_Point.cpp DF_Proj_Report.cpp Util_Minimization_Methods.cpp gaussian_random.cpp)

add_library(DFLibStatic STATIC DF_Abstract_Report.cpp DF_Report_Collection.cpp DF_ProjReport_Collection.cpp DF_XY_Point.cpp DF_LatLon_Point.cpp DF_Proj_Point.cpp DF_Proj_Report.cpp Util_Minimization_Methods.cpp gaussian_random.cpp)

set_target_properties(DFLibStatic PROPERTIES OUTPUT_NAME DFLib)

include_directories(${DFLib_SOURCE_DIR} ${PROJ_INCLUDE_DIR} )
link_directories(${DFLib_BINARY_DIR} ${PROJ_LIB_DIR})

add_executable(testlsDF_proj testlsDF_proj.cpp)
target_link_libraries(testlsDF_proj DFLib ${PROJ_LIBRARY})

add_executable(SimpleDF SimpleDF.cpp)
target_link_libraries(SimpleDF DFLib ${PROJ_LIBRARY})

install(TARGETS SimpleDF testlsDF_proj DFLib DFLibStatic
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

install(FILES  DF_Abstract_Point.hpp DF_Abstract_Report.hpp DF_LatLon_Point.hpp DF_LatLon_Report.hpp DF_ProjReport_Collection.hpp DF_Proj_Point.hpp DF_Proj_Report.hpp DF_Report_Collection.hpp DF_XY_Point.hpp DF_XY_Report.hpp Util_Abstract_Group.hpp Util_Minimization_Methods.hpp Util_Misc.hpp gaussian_random.hpp DFLib_port.h
        DESTINATION include)


