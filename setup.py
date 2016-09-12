"""
Create unix package:    python setup.py sdist
"""

import os
import subprocess
import shutil
from setuptools import setup
from os.path import join

subprocess.call('git log --pretty=format:%h -n 1 > msings/data/sha', shell = True)
subprocess.call('git shortlog --format="XXYYXX%h" | grep -c XXYYXX > msings/data/ver', shell = True)

from msings import __version__
from msings.scripts import script

params = {'author': 'Sheena Scroggins',
          'author_email': 'sheena.scroggins@gmail.com',
          'description': script.__doc__.strip(),
          'name': 'msings',
          'packages': ['msings','msings.scripts','msings.subcommands'],
          'package_dir': {'msings': 'msings'},
          'scripts': ['msi'],
          'version': __version__,
          'package_data': {'msings': [join('data',f) for f in ['sha','ver']]},
          'setup_requires':['numpy'],
          'install_requires': [
            'natsort==4.0.4',
            'numpy==1.10.4',
            'wsgiref==0.1.2',
            'xlrd==0.9.4',
            'xlwt==1.0.0',
            'bioscons==0.8.1,'
],
          'dependency_links' : [
              "git+ssh://git@bitbucket.org/uwlabmed/munge.git@master#egg=munge",
          ]
          }

setup(**params)

