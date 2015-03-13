import glob
from os.path import splitext, split, join
import argparse

def itermodules(subcommands_path, root=__name__):

    commands = [x for x in [splitext(split(p)[1])[0] for p in glob.glob(join(subcommands_path, '*.py'))] if not x.startswith('_')]

    for command in commands:
        yield command, __import__('%s.%s' % (root, command), fromlist=[command])

