import sys
if sys.version < '3.4':
  print('Unsupported Python version: {0:s}.'.format(sys.version))
  print('Supported Python versions are 3.4 or a later 3.x version.')
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


def run_setup(with_binary=True):
    entry_points = {
            'gui_scripts':['CombinatorialProfilerGUI=combinatorialprofiler.ui.CProfilerGUI:main']
        }

    if with_binary:
        kw = dict(setup_requires=['pytest_runner'], tests_require=['pytest'], ext_modules=[readcounter], entry_points=entry_points)
        entry_points['console_scripts'] = ['CombinatorialProfiler=combinatorialprofiler.CombinatorialProfiler:main']
    else:
        kw = {}
    setup(name='CombinatorialProfiler',
        packages=['combinatorialprofiler'],
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
