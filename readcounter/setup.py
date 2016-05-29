from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(ext_modules = cythonize(Extension("readcounter",
           sources=["src/readcounter.pyx", "src/cReadCounter.cpp"],
           language="c++",
           extra_compile_args=["-std=c++11", "-O0"],
           extra_link_args=["-std=c++11"]
      ), gdb_debug=True))
