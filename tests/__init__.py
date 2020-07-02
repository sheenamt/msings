import sys
import logging
import os
from os import path, mkdir
import unittest

from msings.utils import mkdir

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

