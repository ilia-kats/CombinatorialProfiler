from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize

setup(name='CombinatorialProfiler',
    packages=['combinatorialprofiler'],
    install_requires=['numpy', 'pandas'],
    entry_points={
        'console_scripts': ['CombinatorialProfiler=combinatorialprofiler.CombinatorialProfiler:main'],
        'gui_scripts':['CombinatorialProfilerGUI=combinatorialprofiler.ui.CProfilerGUI:main']
    },
    package_data={
        '': ['*.ui']
    },
    ext_modules = cythonize(Extension("readcounter",
        sources=["combinatorialprofiler/readcounter/readcounter.pyx", "combinatorialprofiler/readcounter/cReadCounter.cpp"],
        language="c++",
        extra_compile_args=["-std=c++14", "-O0"],
        extra_link_args=["-std=c++14"]
    ), gdb_debug=True))
