"""
Genetics NGS MSINGS pipeline
"""
import pprint
import glob
import os
import sys
import csv
from shutil import copyfile
import re

from collections import defaultdict
from itertools import groupby
from os import path
from ConfigParser import SafeConfigParser
import ConfigParser
from SCons.Script import (Decider, Variables, Depends, Alias, Help,
                          Flatten, AllowSubstExceptions, Copy)

Decider('MD5-timestamp')

# declare variables for the environment
vars = Variables(None, ARGUMENTS)
vars.Add(PathVariable('settings', 'configuration file', 'settings.conf'))
vars.Add(PathVariable('data', 'data file', 'data.conf'))
vars.Add(PathVariable('output', 'Path to output directory',
                      'output', PathVariable.PathIsDirCreate))

# Provides access to options prior to instantiation of env object
# below; it's better to access variables through the env object.
varargs = dict({opt.key: opt.default for opt in vars.options}, **vars.args)

# get default locations of input file for pipeline programs from settings.conf
settings = varargs['settings']
if not path.exists(settings):
    print 'Either create "settings.conf" or use "scons settings=settings-example.conf"'
    sys.exit(1)

# get the absolute paths to each input file from the data.conf
data = varargs['data']
if not path.exists(data):
    print 'Either create "data.conf" or use "scons data=data.conf"'
    sys.exit(1)

config = SafeConfigParser(allow_no_value=True)
config.optionxform = str # Make sure keys are case sensitive
config.read(settings)
config.read(data)

# Prepare variables and environment
scons_vars = Variables()
for k, v in config.items('msi_info'):
    scons_vars.Add(k, default=v)

data = {}
for key, bam in config.items('specimen_data'):
    data[key] = bam.split(',')
    assert(len(data[key]) == 1)
    assert(path.exists(bam))

# either we're already operating in an activated virtualenv (which we can activate again without
# causing harm) or use the one defined in settings.conf
venv = os.environ.get('VIRTUAL_ENV') or config.get('msi_info','virtualenv')
if not venv:
    sys.exit('No virtualenv is activated, and none is defined in settings.conf.')
elif not path.exists(venv):
    sys.exit("The virtualenv '{}' does not exist or is not readable".format(venv))
else:
    venv = path.abspath(venv)
    print 'using virtualenv "{}"'.format(venv)
# activate the virtualenv
activate_this = path.join(venv, 'bin', 'activate_this.py')
execfile(activate_this, dict(__file__=activate_this))
os.environ['VIRTUAL_ENV'] = venv

# PATH and Environment (preference given to local venv executables)
PATH = [
    'bin',
    path.join(venv, 'bin'),
    '/usr/local/bin',
    '/usr/bin',
    '/bin'
]

# Prefer the local venv for importing modules
sys.path.insert(0, os.path.join(venv, 'lib', 'python2.7', 'site-packages'))

scons_env = dict(os.environ, PATH=':'.join(PATH))
parent_env = Environment(
    ENV = scons_env,
    variables= scons_vars,
    SHELL='bash',
    time=False)

# Run folder is Date_Machine-run_Project
# 150101_HA0012_OncoPlex92
try:
    seq_run = Dir('.').abspath.split('_')[2]
#Unless this is a validation run, then use the name of the project
except IndexError:
    seq_run = Dir('.').abspath.split('/')[-1]

outputs = []
samples = []  # Sample Directories

for pfx in data.keys():
    env = parent_env.Clone()
    env['specimen'] = pfx
    env['pfxout'] = path.join(path.abspath(env['output']),pfx)
    env['logs'] = path.join(env['pfxout'], 'logs')
    seqs = data[pfx]

    # SNP and INDEL calling through VARSCAN
    mpileup= env.Command(
        target='$pfxout/${specimen}.mpileup',
        source=[seqs,
                '$ref_fasta'],
        action=('samtools mpileup -f ${SOURCES[1]} '
                '-d 100000 '  # max per-BAM depth to avoid excessive memory usage
                '-A '  # count anomalous read pairs
                '-E '  # extended BAQ for higher sensitivity but lower specificity
                '${SOURCES[0]} '
                '> ${TARGET} ')
    )
    Alias('mpileup', mpileup)
#######
##MSI##
#######
    # Call MSI with help of varscan readcounts
    msi_output, log = env.Command(
        target=['$pfxout/${specimen}.msi_output',
                '$logs/readcounts'],
        source=[mpileup,
                '$msi_intervals',
                '$varscan'],
        action=('java -Xmx4g -jar ${SOURCES[2]} readcounts '
                ' ${SOURCES[0]} '
                '--variants-file ${SOURCES[1]} '
                '--min-base-qual 10 '
                '--output-file ${TARGETS[0]} '
                '2> ${TARGETS[1]} ')
    )
    Alias('readcounts', msi_output)
    msi_calls = env.Command(
        target='$pfxout/${specimen}.msi.txt',
        source=[msi_output,
                '$msi_bed'],
        action=('msi analyzer ${SOURCES[0]} ${SOURCES[1]} -o $TARGET')
    )
    Alias('msi_calls', msi_calls)
    msi_analysis = env.Command(
        target='$pfxout/${specimen}.MSI_Analysis.txt',
        source='$msi_baseline',
        action=('msi count_msi_samples '
                '$SOURCE '
                '$pfxout '
                '-m $multiplier '
                '-t $msi_min_threshold $msi_max_threshold '
                '-o $TARGET ')
    )
    env.Depends(msi_analysis, msi_calls)
    Alias('msi_analysis', msi_analysis)

    outputs.append(msi_analysis)
    samples.append(env.subst('$pfxout')) # Add sample directories to sample lis
    
#Run top level scripts
run_env = parent_env.Clone()
run_env['seq_run'] = seq_run
combined_MSI = run_env.Command(
    target='$output/${seq_run}.Combined_MSI.txt',
    source='$msi_baseline',
    action=('msi count_msi_samples '
            '$SOURCE '
            '$output '
            '-m $multiplier '
            '-t $msi_min_threshold $msi_max_threshold '
            '-o $TARGET ')
)
run_env.Depends(combined_MSI, [outputs, samples])
