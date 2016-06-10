from setuptools import setup
from setuptools.extension import Extension
from setuptools.command.develop import develop

readcounter = Extension("readcounter",
        sources=["combinatorialprofiler/readcounter/readcounter.cpp", "combinatorialprofiler/readcounter/cReadCounter.cpp"],
        language="c++",
        extra_compile_args=["-std=c++14"],
        extra_link_args=["-std=c++14"]
    )

class debugmode(develop):
    def __init__(self, dist, **kw):
        super().__init__(dist, **kw)

        from Cython.Build import cythonize
        import os.path
        global readcounter
        readcounter.extra_compile_args.append("-O0")
        path, ext = os.path.splitext(readcounter.sources[0])
        readcounter.sources[0] = path + ".pyx"
        readcounter = cythonize(readcounter)


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
    ext_modules = [readcounter],
    cmdclass={
        'develop': debugmode
    }
)
