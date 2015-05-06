# DFLib/Python

This is a rudimentary first cut at creating Python bindings for DFLib using
SWIG.

At this time, only the DFLib::ReportCollection class, and the two abstract
interface classes DFLib::Abstract::Point and DFLib::Abstract::Report have
been wrapped.  This enables one to create implementations of the abstract
classes directly in Python, and use them as-is in the report collection to
get DF fixes.

##Quick start:

You must have SWIG, python and a C++ compiler installed.

Build the python module
```
  python setup.py build_ext --inplace [--swig <your swig path>]
```

If your install of SWIG is not called "swig" (i.e. my FreeBSD system,
where it's called "swig3.0") then include the last option and let
setup know how to run swig.

You are now able to import the DFLib python module in your own
scripts, so long as they're in this directory.  I have not yet
progressed to the stage where this process will actually install this
module where it can be found centrally.

There are two test codes in this directory, "pyPointTest.py" and
"pyProjPoint.py".  Both make use of the python module "pyproj" that
you'll need to obtain elsewhere.  That module binds the proj.4
cartographic projection library in Python.

pyPointTest.py implements the equivalent of the DFLib::XY:Point and
DFLib::XY:Report classes in python, and then runs a single DF fix
problem equivalent to a run of the DFLib test program "testlsDFfix"
using the "receivers3" input file and a particular transmitter
location.

pyProjPoint.py implements the equivalent of
the DFLib::LatLon::Point and DFLib::LatLon::Report classes.  These
classes allow the user to specify receiver locations in lat/lon
instead of XY coordinates, and do internal conversions to a usable XY
representation (in this case, Mercator projection).  Its output should
be the same as the output of the pyPointTest.py script.

