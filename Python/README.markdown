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
scripts, so long as your current working directory is this directory.  

If you want to install this centrally instead, do try:
```
  python setup.py build_ext [--swig <your swig path>]
  sudo python setup.py install
```

This works on my system if I make sure the install got the access
permissions right.  Executing "import DFLib" in python when you start
it outside of this directory should give no error, and you should be able
to run the test scripts.

There are two test codes in this directory, "pyPointTest.py" and
"pyProjPoint.py".  Both make use of the python module "pyproj" that
you'll need to obtain elsewhere.  That module binds the proj.4
cartographic projection library in Python.

A third test code, arcpyPoint.py, is similar to pyProjPoint but uses 
the arcpy module from ArcGIS to do the same work.  It is intended as
a demonstration that DFLib can be used within ArcGIS to do DF computations.
It is my intent to add a DFLib-based module to IGT4SAR someday, but 
I caution you not to hold your breath waiting for it.

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

