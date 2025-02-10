# Sourcery install script

import os
import glob
from setuptools import setup
from setuptools.extension import Extension
#from Cython.Distutils import build_ext
#import numpy

setup(name='sourcery',
      version="0.1",
      url=None,
      author='Matt Hilton',
      author_email='matt.hilton@mykolab.com',
      classifiers=[],
      description='Web-based astronomical source list browser and manager.',
      long_description="""Web-based astronomical source list browser and manager.""",
      packages=['sourcery'],
      package_data={'sourcery': ['data/*', 'static/css/*.css', 'templates/*.html']},
      scripts=['bin/sourcery_build_cache', 'bin/sourcery_build_db', 'bin/sourcery_test', 'bin/sourcery_password_hash', 'bin/sourcery_fast_tag', 'bin/sourcery_fetch_skyview'],
      #cmdclass={'build_ext': build_ext},
      #ext_modules=[Extension("sourceryCython", ["sourcery/sourceryCython.pyx"], include_dirs=[numpy.get_include()])]
)
