#!/usr/bin/env python

from distutils.core import setup, Extension

DFLib_module = Extension('_DFLib',
                       sources=['DFLib.i', '../DF_Abstract_Report.cpp','../DF_Report_Collection.cpp', '../Util_Minimization_Methods.cpp'],
                       swig_opts = ['-c++','-I..'],
                       include_dirs = ['..'],
                       )

setup ( name='DFLib',
        version = '0.1',
        author = 'Cuthbert Twillie',
        description = """Toy swig module""",
        ext_modules = [DFLib_module],
        py_modules = ["DFLib"],
        )
