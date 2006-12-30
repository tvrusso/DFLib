lib_LIBRARIES=libDFLib.a
libDFLib_a_SOURCES=DF_Abstract_Report.cpp \
                   DF_Report_Collection.cpp \
                   DF_XY_Point.cpp \
                   DF_LatLon_Point.cpp \
                   Util_Minimization_Methods.cpp \
                   gaussian_random.cpp

bin_PROGRAMS = testlsDFfix
testlsDFfix_SOURCES = testlsDFfix.cpp
testlsDFfix_LDADD=-lproj -L. -lDFLib