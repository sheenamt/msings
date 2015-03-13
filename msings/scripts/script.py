#!/usr/bin/env

"""
Utilities for the msings scripts
"""

import argparse
from argparse import RawDescriptionHelpFormatter
import sys
import os
import logging
from msings import subcommands, __version__ as version

def main(argv):
    action, arguments = parse_arguments(argv)

    loglevel = {
        0: logging.ERROR,
        1: logging.WARNING,
        2: logging.INFO,
        3: logging.DEBUG,
    }.get(arguments.verbosity, logging.DEBUG)

    if arguments.verbosity > 1:
        logformat = '%(levelname)s %(module)s %(lineno)s %(message)s'
    else:
        logformat = '%(message)s'

    # set up logging
    logging.basicConfig(file=sys.stdout, format=logformat, level=loglevel)

    return action(arguments)

def parse_arguments(argv):
    """
    Create the argument parser
    """

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-V', '--version', action='version',
        version = version,
        help = 'Print the version number and exit')

    parser.add_argument('-v', '--verbose',
        action='count', dest='verbosity', default=1,
        help='Increase verbosity of screen output (eg, -v is verbose, '
             '-vv more so)')
    parser.add_argument('-q', '--quiet',
        action='store_const', dest='verbosity', const=0,
        help='Suppress output')

    ##########################
    # Setup all sub-commands #
    ##########################

    subparsers = parser.add_subparsers(dest='subparser_name')

    # Begin help sub-command
    parser_help = subparsers.add_parser(
        'help', help='Detailed help for actions using `help <action>`')
    parser_help.add_argument('action', nargs=1)
    # End help sub-command

    actions = {}

    for name, mod in subcommands.itermodules(os.path.split(subcommands.__file__)[0]):
        # set up subcommand help text. The first line of the dosctring
        # in the module is displayed as the help text in the
        # script-level help message (`script -h`). The entire
        # docstring is displayed in the help message for the
        # individual subcommand ((`script action -h`)).
        subparser = subparsers.add_parser(name,
                                          help = mod.__doc__.lstrip().split('\n', 1)[0],
                                          description = mod.__doc__,
                                          formatter_class = RawDescriptionHelpFormatter)
        mod.build_parser(subparser)
        actions[name] = mod.action

    # Determine we have called ourself (e.g. "help <action>")
    # Set arguments to display help if parameter is set
    #           *or*
    # Set arguments to perform an action with any specified options.
    arguments = parser.parse_args(argv)
    # Determine which action is in play.
    action = arguments.subparser_name

    # Support help <action> by simply having this function call itself and
    # translate the arguments into something that argparse can work with.
    if action == 'help':
        return parse_arguments([str(arguments.action[0]), '-h'])

    return actions[action], arguments
