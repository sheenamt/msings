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

params = {'author': 'Sheena Todhunter',
          'author_email': 'sheena.todhunter@gmail.com',
          'description': script.__doc__.strip(),
          'name': 'msings',
          'packages': ['msings','msings.scripts','msings.subcommands'],
          'package_dir': {'msings': 'msings'},
          'scripts': ['msi'],
          'version': __version__,
          'package_data': {'msings': [join('data',f) for f in ['sha','ver']]},
          'install_requires': [
              'natsort==3.5.3',
          ],
          }

setup(**params)
