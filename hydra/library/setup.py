from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

examples_extension = Extension(
    name="cyrand",
    sources=["cyrand.pyx", "cython/rand_index.c"],
    library_dirs=["cython"],
    include_dirs=["cython"])

setup(
    name="cyrand",
    ext_modules=cythonize([examples_extension], 
        gdb_debug=True)
)
