import sys
import logging
import os
from os import path, mkdir
import unittest

from msings.utils import mkdir

# set up logging for unit tests
verbosity_flag = [x for x in sys.argv if x.startswith('-v')]
verbosity = (verbosity_flag[0] if verbosity_flag else '').count('v')

loglevel = {
    0: logging.WARNING,
    1: logging.INFO,
    2: logging.DEBUG,
}.get(verbosity, logging.DEBUG)

if verbosity > 1:
    logformat = '%(levelname)s %(module)s %(lineno)s %(message)s'
else:
    logformat = '%(message)s'

logging.basicConfig(file=sys.stdout, format=logformat, level=loglevel)
log = logging.getLogger(__name__)

# module data
datadir = 'testfiles'
outputdir = 'test_output'

mkdir(outputdir)

class TestBase(unittest.TestCase):
    """
    Base class for unit tests with methods for defining output
    directories based on method name.
    """

    outputdir = outputdir

    def mkoutdir(self, clobber = True):
        """
        Create outdir as outpudir/module.class.method (destructively
        if clobber is True).
        """
        
        funcname = '.'.join(self.id().split('.')[-3:])
        outdir = path.join(self.outputdir, funcname)
        mkdir(outdir, clobber)
        return outdir
    
class TestCaseSuppressOutput(unittest.TestCase):

    def setUp(self):        
        self.funcname = '_'.join(self.id().split('.')[-2:])
        self.suppress_output = log.getEffectiveLevel() >= logging.INFO
        if self.suppress_output:
            sys.stdout = sys.stderr = open(os.devnull, 'w')
            
    def tearDown(self):
        if self.suppress_output:
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__

