import sys
if sys.version < '3.4':
  print('Unsupported Python version: {0:s}.'.format(sys.version))
  print('Supported Python versions are 3.4 or a later 3.x version.')
  sys.exit(1)

from distutils.errors import DistutilsPlatformError
from distutils import ccompiler
from setuptools import setup, find_packages
from setuptools.extension import Extension
from setuptools.command.develop import develop

from glob import glob

from combinatorialprofiler import version

readcounter_basepath = "combinatorialprofiler/readcounter/"
readcounter = Extension("combinatorialprofiler.readcounter",
        sources=glob(readcounter_basepath + "*.cpp"),
        language="c++",
        extra_compile_args=["-std=c++14"],
        extra_link_args=["-std=c++14"]
    )
readcounter = [readcounter]

class debugmode(develop):
    def __init__(self, dist, **kw):
        super().__init__(dist, **kw)

        from Cython.Build import cythonize
        import os.path
        global readcounter
        rdcntr = readcounter[0]
        rdcntr.extra_compile_args.append("-O0")
        rdcntr.undef_macros = ['NDEBUG']

        cythonfile = rdcntr.sources.index(readcounter_basepath + "readcounter.cpp")
        path, ext = os.path.splitext(rdcntr.sources[cythonfile])
        rdcntr.sources[cythonfile] = path + ".pyx"
        readcounter[0] = cythonize(readcounter)[0]


dependencies = ['numpy', 'scipy', 'pandas>=0.23', 'matplotlib>=1.5.1', 'scikit-learn>=0.18', 'biopython', 'llist']

def run_setup(deps, with_binary=True):
    entry_points = {
            'gui_scripts':['CombinatorialProfilerGUI=combinatorialprofiler.ui.CProfilerGUI:main']
        }

    if with_binary:
        kw = dict(setup_requires=['pytest_runner'], tests_require=['pytest'], ext_modules=readcounter)
        entry_points['console_scripts'] = ['CombinatorialProfiler=combinatorialprofiler.CombinatorialProfiler:main']
        deps.append('feather-format')
    else:
        deps[deps.index('llist')] = 'pyllist'
        kw = {}
    kw['entry_points'] = entry_points
    setup(name='CombinatorialProfiler',
        packages=find_packages(exclude=['test', 'tests']),
        install_requires=deps,
        package_data={
            '': ['*.ui']
        },
        cmdclass={
            'develop': debugmode
        },
        zip_safe = True,
        version=version,
        description='Pipeline for analysis of combinatorial stability profiling data',
        author="Ilia Kats",
        author_email="i.kats@zmbh.uni-heidelberg.de",
        license="GPLv2",
        dependency_links=['https://github.com/fuzeman/pypy-llist/tarball/master#egg=pyllist-0.1.1'],
        **kw
    )

comp = ccompiler.new_compiler()
binary = True
if comp.compiler_type == 'msvc':
    try:
        comp.initialize()
    except DistutilsPlatformError:
        binary = False

run_setup(dependencies, binary)

