# DFLib

DFLib is a library of C++ classes for computing the locations of
transmitters from bearings-only measurements.  See
docs/html/index.html for detailed documentation fo the classes.

##Quick start:


Old-style autoconf build:

```
   mkdir build
   cd build
   ../configure
   make
   sudo make install
```

Add configure command line arguments as needed to find PROJ.4 library.
See INSTALL for details.


New style CMake build (cross-platform): See INSTALL_CMake


Two demo programs are built by this process, SimpleDF and
testlsDF_proj.  See their respective entries in docs/html/index.html
for some documentation about them.

DFLib is primarily intended as a set of building blocks for DF fix
programs, not as a finished product itself.  A companion program
[qDF](https://github.com/tvrusso/qDF) is a GUI-based direction finding
program that uses DFLib to compute DF fixes and create KML files that
can be displayed in Google Earth, and can optionally send data for
display or transmission to an external APRS program such as
[Xastir](http://www.xastir.org) that provides server ports.