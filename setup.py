#!/usr/bin/env python
from setuptools import setup
from pyxarithmetic import __version__

setup(name='pyxarithmetic',
      packages=['pyxarithmetic'],
      version=__version__,
      description='Python package for computing arithmetic with X-ray images',
      author='John ZuHone',
      author_email='jzuhone@gmail.com',
      url='http://github.com/jzuhone/pyxarithmetic',
      install_requires=["six","numpy","astropy","soxs","scipy"],
      classifiers=[
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.5',
          'Topic :: Scientific/Engineering :: Visualization',
      ],
      )
