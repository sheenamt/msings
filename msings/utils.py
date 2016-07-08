from collections import Iterable
import datetime
import re
import os
import shutil
import logging
from collections import namedtuple


from __init__ import __version__

log = logging.getLogger(__name__)
Path = namedtuple('Path', ['dir','fname'])

ASSAYS={'OPX':'OncoPlex',    
        'BRO':'ColoSeq',
        'EPI':'EpiPlex',
        'MRW':'MarrowSeq',
        'IMD':'ImmunoPlex',
        'TESTDATA':'testdata',
        'MSI-PLUS':'msi-plus'}

def walker(dir):
    """Recursively traverse direcory `dir`. For each tuple containing
    (path, dirs, files) return a named tuple with attributes (dir,
    fname) for each file name in `files`.
    """

    for (pth, dirs, files) in os.walk(dir):
        for fname in files:
            yield Path(pth, fname)

def munge_pfx(pfx):
    """
    Get the plate,well, library-version, assay, control 
    and machine-run from the pfx
    """
    output=pfx.split('.')[0].split('_')
    if len(output)==5:
        #New control samples hit this
        keys=['sample_id','well','library-version','control','machine-run']
        pfx_info = dict(zip(keys,output))
        pfx_info['mini-pfx']='{sample_id}_{control}'.format(**pfx_info)
        pfx_info['run']=pfx_info['sample_id'][:-2]
        pfx_info['pfx']='{sample_id}_{well}_{library-version}_{control}_{machine-run}'.format(**pfx_info)
        pfx_info['assay']=[value for key,value in ASSAYS.items() if re.search(key, pfx_info['library-version'])][0]
    elif len(output)==4:
        #New non-control samples hit this
        keys=['sample_id','well','library-version','machine-run']
        pfx_info = dict(zip(keys,output))
        pfx_info['mini-pfx']='{sample_id}'.format(**pfx_info)
        pfx_info['run']=pfx_info['sample_id'][:-2]
        pfx_info['pfx']='{sample_id}_{well}_{library-version}_{machine-run}'.format(**pfx_info)
        pfx_info['assay']=[value for key,value in ASSAYS.items() if re.search(key, pfx_info['library-version'])][0]
    elif len(output)==3:
        #Research non-control samples often hit this
        keys=['sample_id','well','library-version']
        pfx_info = dict(zip(keys,output))
        pfx_info['mini-pfx']='{sample_id}'.format(**pfx_info)
        pfx_info['run']=pfx_info['sample_id'][:-2]
        pfx_info['pfx']='{sample_id}_{well}_{library-version}'.format(**pfx_info)
        pfx_info['assay']=[value for key,value in ASSAYS.items() if re.search(key, pfx_info['library-version'])][0]
    elif len(output)==2:
        #MSI-Plus will hit this
        keys=['sample_id', 'library-version']
        pfx_info = dict(zip(keys,output))
        pfx_info['mini-pfx']='{sample_id}'.format(**pfx_info)
        pfx_info['pfx']='{sample_id}'.format(**pfx_info)
        if re.search('msi',pfx_info['library-version'], re.IGNORECASE):
            pfx_info['assay']=[value for key,value in ASSAYS.items() if re.search(key, pfx_info['library-version'])]
    elif len(output)==1:
        #Only the old LMG/OPX should hit this. 
        pfx=output[0]
        pfx_info={'mini-pfx':pfx,
                  'pfx':pfx,
                  'sample_id':pfx}
        if re.search('LMG', pfx):
            pfx_info['assay']='Coloseq'
        elif re.search('OPX', pfx):
            pfx_info['assay']='OncoPlex'
    else:
        pfx_info={'mini-pfx':pfx,
                  'pfx':pfx,
                  'sample_id':pfx,
                  'assay':'NA'}

    return pfx_info

def multi_split(source, splitlist):
    """
    Function to split a string given a string of multiple split points
    """
    output = []
    atsplit = True
    if source is None:
        return None
    else:
        for char in source:
            if char in splitlist:
                atsplit = True
            else:
                if atsplit:
                    output.append(char)
                    atsplit = False
                else:
                    output[-1] = output[-1] + char
    return output

def mkdir(dirpath, clobber = False):
    """
    Create a (potentially existing) directory without errors. Raise
    OSError if directory can't be created. If clobber is True, remove
    dirpath if it exists.
    """

    if clobber:
        shutil.rmtree(dirpath, ignore_errors = True)

    try:
        os.mkdir(dirpath)
    except OSError, msg:
        pass

    if not os.path.exists(dirpath):
        raise OSError('Failed to create %s' % dirpath)

    return dirpath

    
def munge_path(pth):
    """
    Get date, run, project, machine, assay, prep-type from path
    """
    output=multi_split(pth, '/_')
    #Assuming we want YYMMDD_RUN_PROJECT
    if output[-1]=='output':
        output=output[-4:-1]
    #If old version of data that isn't in a 'output' subfolder
    else:
        output=output[-3:]
    keys=['date','run', 'project']
    pathinfo = dict(zip(keys,output))
    pathinfo['date']=munge_date(pathinfo['date'])
    #Set Machine
    if re.search('HA', pathinfo['run']):
        pathinfo['machine']='hiseq'
    elif re.search('MA', pathinfo['run']):
        pathinfo['machine']='miseq'
    #Set assay
    if re.search('colo', pathinfo['project'].lower()):
        pathinfo['assay']='coloseq'
    elif re.search('onco', pathinfo['project'].lower()):
        pathinfo['assay']='oncoplex'
    elif re.search('epi', pathinfo['project'].lower()):
        pathinfo['assay']='epiplex'
    elif re.search('imm', pathinfo['project'].lower()):
        pathinfo['assay']='immunoplex'
    elif re.search('mrw', pathinfo['project'].lower()):
        pathinfo['assay']='marrowseq'
    #Set prep type
    if re.search('kapa', pathinfo['project'].lower()):
        pathinfo['prep_type']='kapa'
    else:
        pathinfo['prep_type']='sure_select'

    return pathinfo

def munge_date(date):
    """
    Convert new date YYMMDD to YYYY-MM-DD
    """
    #only convert if new date format
    if len(date)<7:
        d = datetime.datetime.strptime(date, '%y%m%d')
        return d.strftime('%Y-%m-%d')
    else:
        return date
