lib_LIBRARIES=libDFLib.a
libDFLib_a_SOURCES=DF_Abstract_Report.cpp \
                   DF_Report_Collection.cpp \
                   DF_XY_Point.cpp \
                   DF_LatLon_Point.cpp \
                   DF_Proj_Point.cpp \
                   Util_Minimization_Methods.cpp \
                   gaussian_random.cpp

bin_PROGRAMS = testlsDFfix testlsDF_ll testlsDF_proj
testlsDFfix_SOURCES = testlsDFfix.cpp
testlsDFfix_LDADD=-L. -lDFLib

testlsDF_ll_SOURCES = testlsDF_ll.cpp
testlsDF_ll_LDADD=-L. -lDFLib

testlsDF_proj_SOURCES = testlsDF_proj.cpp
testlsDF_proj_LDADD=-L. -lDFLib