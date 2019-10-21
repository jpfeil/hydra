from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = Extension('rand_index', 
        sources = ['rand_index.pyx'],
        library_dirs=['.'],
        include_dirs=['.'])

setup(name="rand_index", 
      ext_modules = cythonize([ext]))
