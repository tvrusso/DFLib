%module(directors="1") DFLib
%feature("director");
%include "std_string.i"
%include "std_vector.i"
namespace std{
  %template(vectord) vector<double>;
 };
%include DFLib_port.h
%include DF_Abstract_Point.i
%include DF_Abstract_Report.i
%include DF_Report_Collection.i
