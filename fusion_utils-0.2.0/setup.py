#!/usr/bin/env python

from distutils.core import setup

setup(name='fusion_utils',
      version='0.2.0',
      description='Python tools for comparing gene fusion and structural variation results/',
      author='Yuichi Shiraishi',
      author_email='friend1ws@gamil.com',
      url='https://github.com/friend1ws/fusion_utils',
      package_dir = {'': 'lib'},
      packages=['fusion_utils'],
      scripts=['fusion_utils'],
      license='GPL-3'
     )

