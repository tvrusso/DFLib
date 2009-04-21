lib_LIBRARIES=libDFLib.a
libDFLib_a_SOURCES=DF_Abstract_Report.cpp \
                   DF_Report_Collection.cpp \
                   DF_ProjReport_Collection.cpp \
                   DF_XY_Point.cpp \
                   DF_LatLon_Point.cpp \
                   DF_Proj_Point.cpp \
                   DF_Proj_Report.cpp \
                   Util_Minimization_Methods.cpp \
                   gaussian_random.cpp

bin_PROGRAMS = testlsDFfix testlsDF_ll testlsDF_proj SimpleDF SimpleDF2
testlsDFfix_SOURCES = testlsDFfix.cpp
testlsDFfix_LDADD=-L. -lDFLib
testlsDFfix_DEPENDENCIES=libDFLib.a

testlsDF_ll_SOURCES = testlsDF_ll.cpp
testlsDF_ll_LDADD=-L. -lDFLib
testlsDF_ll_DEPENDENCIES=libDFLib.a

testlsDF_proj_SOURCES = testlsDF_proj.cpp
testlsDF_proj_LDADD=-L. -lDFLib
testlsDF_proj_DEPENDENCIES=libDFLib.a

SimpleDF_SOURCES = SimpleDF.cpp
SimpleDF_LDADD=-L. -lDFLib
SimpleDF_DEPENDENCIES=libDFLib.a

SimpleDF2_SOURCES = SimpleDF2.cpp
SimpleDF2_LDADD=-L. -lDFLib
SimpleDF2_DEPENDENCIES=libDFLib.a