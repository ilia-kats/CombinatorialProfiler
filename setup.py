from setuptools import setup
from setuptools.extension import Extension
from setuptools.command.develop import develop
from Cython.Build import cythonize

readcounter = Extension("readcounter",
        sources=["combinatorialprofiler/readcounter/readcounter.pyx", "combinatorialprofiler/readcounter/cReadCounter.cpp"],
        language="c++",
        extra_compile_args=["-std=c++14"],
        extra_link_args=["-std=c++14"]
    )

debug = False

class debugmode(develop):
    def __init__(self, dist, **kw):
        super().__init__(dist, **kw)
        debug=True
        readcounter.extra_compile_args.append("-O0")


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
    ext_modules = cythonize(readcounter, gdb_debug=debug),
    cmdclass={
        'develop': debugmode
    }
)
