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

SETTINGS = 'settings.conf'
DATA = 'data.conf'
if not path.exists(SETTINGS):
    sys.exit("Cannot find configuration file '{}'".format(SETTINGS))
if not path.exists(DATA):
    sys.exit("Cannot find configuration file '{}'".format(DATA))

config = SafeConfigParser(allow_no_value=True)
config.optionxform = str # Make sure keys are case sensitive
config.read(SETTINGS)
config.read(DATA)

# Prepare variables and environment
scons_vars = Variables()
for k, v in config.items('msi_info'):
    scons_vars.Add(k, default=v)

data = {}
for key, fastq in config.items('specimen_data'):
    data[key] = fastq.split(',')
    assert(len(data[key]) == 1)
    assert(path.exists(fastq))

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

PERL5LIB = [path.join(venv, 'bin'),
            path.join(venv, 'lib'),
            path.join(venv, 'lib', 'perl5')]

#nproc = int(_params.defaults().get('nproc', 1))
scons_env = dict(os.environ,PATH=':'.join(PATH), PERL5LIB=':'.join(PERL5LIB))
#    THREADS_ALLOC=str(nproc))

from bioscons.slurm import SlurmEnvironment
parent_env = SlurmEnvironment(
    ENV = scons_env,
    variables= scons_vars,
    SHELL='bash',
    time=False)

# Run folder is Project_<Run>
try:
    seq_run = Dir('.').abspath.split('_')[2]
#Unless this is a validation run, then use the name of the project
except IndexError:
    seq_run = Dir('.').abspath.split('/')[-1]

outputs = []
samples = []  # Sample Directories
varscan=os.path.join(venv, 'bin/VarScan.v2.3.7.jar')

for pfx in data.keys():
    env = parent_env.Clone()
    env['specimen'] = pfx
    print "pfx:", pfx
    env['pfxout'] = path.join(path.abspath(env['output']),pfx)
    env['log'] = path.join(env['pfxout'], 'log')

    # Ensure order is R1 followed by R2
    seqs = sorted(data[pfx])
    # create a subdirectory for this specimen
    env['ctrlpath'] = env.subst('$output/$baseline_control')

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
    msi_output, log = env.Command(
        target=['$pfxout/${specimen}.msi_output',
                '$log/readcounts'],
        source=[mpileup,
                '$msi_intervals',
                varscan],
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
    baseline_control = '{}/{}.msi.txt'.format(env['ctrlpath'],env['baseline_control'])
    msi_analysis = env.Command(
        target='$pfxout/${specimen}.MSI_Analysis.txt',
        source=msi_calls,
        action=('msi count_msi_samples '
                '$ctrlpath/${baseline_control}.msi.txt '
                '$SOURCE '
                '-m $multiplier '
                '-t $msi_min_threshold $msi_max_threshold '
                '-o $TARGET ')
    )
    if not pfx == env['baseline_control']:
        env.Depends(msi_analysis, baseline_control)
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
