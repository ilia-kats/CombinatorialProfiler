import sys
if sys.version < '3.4':
  print('Unsupported Python version: {0:s}.'.format(sys.version))
  print('Supported Python versions are 3.4 or a later 3.x version.')
  sys.exit(1)


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

class debugmode(develop):
    def __init__(self, dist, **kw):
        super().__init__(dist, **kw)

        from Cython.Build import cythonize
        import os.path
        global readcounter
        readcounter.extra_compile_args.append("-O0")

        print(readcounter.sources)
        cythonfile = readcounter.sources.index(readcounter_basepath + "readcounter.cpp")
        path, ext = os.path.splitext(readcounter.sources[cythonfile])
        readcounter.sources[cythonfile] = path + ".pyx"
        readcounter = cythonize(readcounter)


def run_setup(with_binary=True):
    entry_points = {
            'gui_scripts':['CombinatorialProfilerGUI=combinatorialprofiler.ui.CProfilerGUI:main']
        }

    if with_binary:
        kw = dict(setup_requires=['pytest_runner'], tests_require=['pytest'], ext_modules=[readcounter])
        entry_points['console_scripts'] = ['CombinatorialProfiler=combinatorialprofiler.CombinatorialProfiler:main']
    else:
        kw = {}
    kw['entry_points'] = entry_points
    setup(name='CombinatorialProfiler',
        packages=find_packages(exclude=['test', 'tests']),
        install_requires=['numpy', 'scipy', 'pandas>=0.15', 'matplotlib', 'biopython'],
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
        **kw
    )

try:
    run_setup(True)
except BaseException as e:
    print('*' * 75)
    print(e)
    print('The readcounter extension could not be compiled. Only the GUI will be installed.')
    print('*' * 75)
    run_setup(False)
