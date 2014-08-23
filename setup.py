# Sourcery install script

import os
import glob
from distutils.core import setup
from distutils.extension import Extension
import popen2

setup(name='sourcery',
      version="git",
      url=None,
      author='Matt Hilton',
      author_email='matt.hilton@mykolab.com',
      classifiers=[],
      description='Web-based astronomical source list browser and manager.',
      long_description="""Web-based astronomical source list browser and manager.""",
      packages=['sourcery'],
      package_data={'sourcery': ['data/*']},
      scripts=['bin/sourcery_build_cache', 'bin/sourcery_test'],
      #cmdclass={'build_ext': build_ext},
      #ext_modules=[Extension("nemoCython", ["nemo/nemoCython.pyx"])]
)
