import sys
if sys.version < '3.5':
  print('Unsupported Python version: {0:s}.'.format(sys.version))
  print('Supported Python versions are 3.5 or a later 3.x version.')
  sys.exit(1)


from setuptools import setup
from setuptools.extension import Extension
from setuptools.command.develop import develop

from combinatorialprofiler import version

readcounter = Extension("combinatorialprofiler.readcounter",
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
    install_requires=['numpy', 'pandas', 'matplotlib', 'biopython'],
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
    },
    zip_safe = True,
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],

    version=version,
    description='Pipeline for analysis of combinatorial stability profiling data',
    author="Ilia Kats",
    author_email="i.kats@zmbh.uni-heidelberg.de",
    licence="GPLv2"
)
